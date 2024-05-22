library(readxl)
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(lemon)

mmP <- 30.97376
mmS <- 32.065
mmFe <- 55.845
mmAl <- 26.981
mmCa <- 40.078
mmMn <- 54.938
mmTi <- 47.867
mmMg <- 24.30506

# First import all raw data. The concentrations of P extraction solutions, and the Dryweight data.
conc_raw <- read_xlsx("input/SequentialExtractions.xlsx")
conc_Bipy_raw <- read_xlsx("input/SequentialExtractionsBipy.xlsx")
wt_raw <- read_xlsx("input/Dryweight.xlsx")
times <- read_xlsx("input/timeseries.xlsx")
kwa_raw <- read_xlsx("input/KWA.xlsx")


# Function to calculate volume based on fraction
get_volume <- function(fraction) {
  case_when(
    fraction == "H2O"   ~ 0.0202,
    fraction %in% c("BD") ~ 0.0415,
    fraction == "NaOH"  ~ 0.035,
    fraction == "Bipy" ~ 0.04,
    TRUE                ~ 0
  )
}
get_mm <- function(parameter) {
  case_when(
    parameter %in% c("SRP","TP", "NRP","NE","P")  ~ mmP,
    parameter %in% c("Fe") ~ mmFe,
    parameter == "S"  ~ mmS,
    parameter == "Mn" ~ mmMn,
    parameter == "Al" ~ mmAl,
    parameter == "Ti" ~ mmTi,
    parameter == "Ca" ~ mmCa,
    parameter == "Mg" ~ mmMg,
    TRUE                ~ 0
  )
}

# Define porosity values for all measurements
porosity <- wt_raw %>%
  group_by(Sample) %>%
  summarise(porosity = mean(Porosity))


#####1st incubation #####
conc <- conc_raw %>%
  left_join(porosity, by = c("Sample")) %>%
  left_join(times, by = c("Sample")) %>%
  mutate(volume = get_volume(Fraction),
         cont = ((Conc - Blank) * volume) / (Weight * (1 - porosity))) %>%
  select(-Blank) %>% 
  transform(Fraction = case_when(
    Fraction == "H2O"   ~ "H2O",
    Fraction == "BD" ~ "BD",
    Fraction == "NaOH" & Parameter == "SRP"  ~ "NaOH SRP",
    Fraction == "NaOH" & Parameter == "TP"  ~ "NaOH TP",
    Fraction == "Bipy" ~ "Bipy",
    TRUE                ~ 0
  ), factor(Fraction, 
                              levels = c("H2O","BD","NaOH")))%>%
  group_by(Time, Fraction, Sample) %>%
  summarise(sd = sd(cont), cont=(mean(cont)), .groups = "drop")

conc$y_pos = NA
conc$y_pos[conc$Fraction == "H2O"] = conc$cont[conc$Fraction == "H2O"]
conc$y_pos[conc$Fraction == "BD"] = conc$cont[conc$Fraction == "H2O"] + 
  conc$cont[conc$Fraction == "BD"]
conc$y_pos[conc$Fraction == "NaOH"] = conc$y_pos[conc$Fraction == "BD"] + 
  conc$cont[conc$Fraction == "NaOH"]
# Combine concentration data with porosity and blanks
# Compute standard deviation for each combination of Time, Fraction, and Sample
# Plot the sequential extractions with error bars

seqextr <- ggplot(conc, aes(x = Time, y = cont, fill = Fraction)) +
  geom_col(position = position_stack(reverse = TRUE), width = 2.0) +
  geom_errorbar(data = conc, aes(ymin = y_pos - sd, ymax = y_pos + sd), width = 1, position = position_dodge(0.1)) +
  scale_fill_manual(values = unlist(colors), labels = legend_labels, name = "Extracted Fraction") +
  labs(y = expression(paste("Extracted P (umol gDW"^{-1}, ")")),
       x = "Time (days)") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16, angle = 0, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.3, "cm")
  )

##### 2nd incubation #####


kwa <- kwa_raw %>% 
  pivot_longer(cols =  Al:Zn, names_to = "Element", values_to = "Conc") %>% 
  mutate(mm = get_mm(Element),
        Total = (0.05*Conc*1000)/(Weight*mm),
         Element = paste("Total_",Element, sep = "")) %>% 
  select(Time,Treatment,Total, Depth, Element) %>% 
  pivot_wider(names_from = Element, values_from = Total)

conc_Bipy <- conc_Bipy_raw %>%
  mutate(volume = get_volume(Fraction),
         cont = ((Conc - Blank) * volume) / (Weight * (1 - Porosity))) %>% 
  select(-Blank) %>% 
  filter(Time != "ts3" | Treatment != "Control" | Nr != 1) %>% 
  filter(Fraction != "Bipy" | Parameter != "SRP") %>% 
  filter(Fraction != "BD" | Parameter != "SRP") %>% 
  mutate(Fraction = case_when(
    Fraction == "H2O"   ~ "H2O",
    Fraction == "BD" ~ "BD",
    Fraction == "NaOH" & Parameter == "SRP"  ~ "NaOHSRP",
    Fraction == "NaOH" & Parameter == "TP"  ~ "NaOHTP",
    Fraction == "Bipy" ~ "Bipy",
    TRUE                ~ as.character(0)
  ))%>% 
  select(Treatment,Depth,Extraction,Time,Nr,cont,Fraction) %>% 
  pivot_wider(names_from = Fraction, values_from = cont)  %>% 
  left_join(kwa, by = c("Treatment","Depth","Time")) %>% 
  mutate(NE = case_when( 
    is.na(Bipy) ~ Total_P-NaOHTP-BD-H2O,
    T ~ Total_P-NaOHTP-Bipy-BD-H2O
                         ),
    NaOHNRP = NaOHTP-NaOHSRP) %>% 
  select(-NaOHTP) %>% 
  pivot_longer(cols= H2O:NaOHNRP, names_to = "Fraction", values_to = "cont", values_drop_na = T)%>% 
group_by( Time,Treatment,Extraction,Fraction,Depth) %>%
  summarise(sd = sd(cont), cont=(mean(cont)))

test <- conc_Bipy_raw %>% filter(Fraction != "H2O", Time != "ts0", Parameter == "TP")

conc_Bipy$y_pos = NA
conc_Bipy$y_pos[conc_Bipy$Fraction == "H2O"] = conc_Bipy$cont[conc_Bipy$Fraction == "H2O"]

conc_Bipy$y_pos[conc_Bipy$Fraction == "Bipy"] = conc_Bipy$cont[conc_Bipy$Fraction == "H2O"] + 
  conc_Bipy$cont[conc_Bipy$Fraction == "Bipy"]

conc_Bipy$y_pos[conc_Bipy$Fraction == "BD"] = conc_Bipy$y_pos[conc_Bipy$Fraction == "Bipy"] + 
  conc_Bipy$cont[conc_Bipy$Fraction == "BD"]
conc_Bipy$y_pos[conc_Bipy$Fraction == "NaOH"] = conc_Bipy$y_pos[conc_Bipy$Fraction == "BD"] + 
  conc_Bipy$cont[conc_Bipy$Fraction == "NaOH"]



# Labels and colors
legend_labels <- c("H2O", "BD", "NaOH SRP", "NaOH NRP")
colors <- list(c("palegreen2", "darkorange2", "cadetblue","coral4"))

legend_labels_Bipy <- c("H2O","Bipy", "BD", "NaOH SRP", "NaOH NRP")
colors_Bipy <- list(c("palegreen2","red", "darkorange2", "cadetblue","coral4","darkgrey"))
fraction_list <- c("H2O","Bipy","BD","NaOHSRP","NaOHNRP", "NE")
legend_labels_Bipy <- c("H2O","Bipy", "BD", "NaOH SRP", "NaOH NRP", "recalcitrant P")

 #### Plots ####




 #Bipy plots

#only Bipy
Bipy_start = filter(conc_Bipy, Fraction =="Bipy", Treatment == "Start") %>% pull(cont)
Bipy_start_line = data.frame(cont = c(Bipy_start[1],Bipy_start[1],Bipy_start[2],Bipy_start[2]),
                             Depth =c(0,2,2,6),
                             Time = c("t=0","t=0","t=0","t=0"),
                             Treatment = factor(c("Start","Start","Start","Start")))

bipyprofile_Aarhus <- ggplot( transform(filter(conc_Bipy, Fraction == "Bipy", Treatment != "Start"),
                  Treatment = factor(Treatment, levels = c("Start","Sulfate", "Control"))), 
        mapping = aes(
          y = cont,
          by = Time,
          color = Treatment,
          x = Depth )) +
  geom_line( aes(
    y = cont,
    x = Depth,
    linetype = Time),
    size = 1)  +
  geom_point(size = 2) +
  geom_path(data = Bipy_start_line,
            aes(linetype = Time),
    size = 1
    )  +
  xlab("Depth (cm)") +
  scale_x_reverse() +
  coord_flip() +
  scale_color_discrete(breaks=c("Start","Control","Sulfate"))+
  scale_linetype_discrete(breaks = c("t=0","ts1","ts2", "ts3"), labels= c("Start","19 days", "48 days", "82 days"))+
  labs(  y = expression(paste("concentration (umol/g dw)")), 
         x = "Depth in cm",
         title = "Bipy profiles" ) +
  theme_bw() +
  theme(axis.title=element_text(size=20),
        plot.title = element_text(size=20, face="bold", hjust = 0.5) ,
        axis.text=element_text(size=14,angle = 0, hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.8, "cm"))


# plot the first timestep
seqextr_Bipy <- ggplot(conc_Bipy, aes(x = Sample, y = cont, fill = Fraction)) +
  geom_col(position = position_stack(reverse = TRUE), width = .7) +
  #geom_errorbar(data = conc_Bipy, aes(ymin = y_pos - sd, ymax = y_pos + sd), width = 1, position = position_dodge(0.1)) +
  scale_fill_manual(values = unlist(colors_Bipy), labels = legend_labels_Bipy, name = "Extracted Fraction") +
  labs(y = expression(paste("Extracted P (umol gDW"^{-1}, ")")),
       x = "Time (days)") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 14, angle = 0, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.3, "cm")
  )

t0plot <- ggplot( transform( filter(conc_Bipy, Time == "ts0", Fraction %in% fraction_list), 
                   Fraction = factor(Fraction, 
                                     levels = fraction_list)), 
        aes( x=Depth, y = cont, fill = Fraction),
        scale_x_discrete(position = 'top')) +
  geom_bar(
    mapping = aes( y = cont, fill = Fraction ),
    stat = "identity", position = position_stack(reverse=TRUE))+
  scale_fill_manual(values = unlist(colors_Bipy), labels = legend_labels_Bipy, name = "Extracted Fraction") +
  coord_flip()+
  scale_x_reverse() + 
  labs(y = expression(paste("Extracted P (",mu,"mol gDW"^-1,")" )),
       x = "Depth (cm)"
  ) +
  theme_bw()+
  facet_grid(.~ Extraction) +
  ylim(0, 250) +
  theme(axis.title=element_text(size=20),plot.title = element_text(hjust = 0.5),
    axis.text=element_text(size=16,angle = 0, hjust = 0.5),
    legend.title = element_text( size = 20),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.3, "cm")) 
show(t0plot)


total_plot <- ggplot( transform( filter(conc_Bipy, Time != "ts0",  Time != "ts2", Extraction != "Bipyridine", Depth != 0 ), 
                                 Extraction = factor(Extraction, levels = c("Normal", "Bipyridine"), labels = c("extraction", "extraction")),
                   Fraction = factor(Fraction, 
                                     levels = c("H2O","Bipy","BD","NaOHSRP","NaOHNRP")),
                   Time = factor(Time, levels = (c("ts0", "ts1","ts2","ts3")), labels = (c("start", "20 days"," more days","98 days"))),
                   Treatment = factor(Treatment, levels = (c("Start", "Control","Sulfate")))), 
        aes( x=Depth, y = cont, fill = Fraction),
        scale_x_discrete(position = 'top')) +
  geom_bar(
    mapping = aes( y = cont, fill = Fraction ),
    stat = "identity", position = position_stack(reverse=TRUE))+
  scale_fill_manual(values = unlist(colors), labels = legend_labels, name = "Extracted Fraction") +
  coord_flip()+
  scale_x_reverse() + 
  labs(y = expression(paste("Extracted P (",mu,"mol gDW"^-1,")" )),
       x = "Depth (cm)"
  ) +
  facet_nested(Treatment ~ Time + Extraction, labeller =label_value) +
  ylim(0, 250) +
  theme_bw()+
  theme(
    axis.title=element_text(size=20),plot.title = element_text(hjust = 0.5),
    axis.text=element_text(size=14,angle = 0, hjust = 0.5),
    legend.title = element_text( size = 20),
    strip.text = element_text(size=10, face="bold", hjust = 0.5) ,
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.3, "cm"),
    legend.position = "top") 

concsum <- conc_Bipy %>% 
  filter(Depth<4.5, Treatment != "Start") %>% 
  group_by(Time,Treatment, Extraction,Fraction) %>% 
  summarise(cont = sum(cont), sd = sqrt(sum(sd^2)))

Psum_Aarhus <- 
  ggplot( transform(filter(concsum, Time != "ts0"),
                      Fraction = factor(Fraction,levels = c("H2O","Bipy","BD","NaOHSRP","NaOHNRP")),
                      Treatment = factor(Treatment, levels = c("Sulfate", "Control"))), 
                              mapping = aes(
                                y = cont,
                                color = Treatment,
                                x = Time)) +
  geom_line(size = 1)  +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = cont - sd, ymax = cont + sd),
    position = position_dodge(width = 0.1),
    width = 0.25
  ) +
  xlab("Time") +
  facet_grid(Fraction ~ Extraction, scales = "free_y")
  labs(  y = expression(paste("concentration (umol/g dw)")), 
         x = "timw",
         title = "Bipy profiles" ) +
  theme_bw() +
  theme(axis.title=element_text(size=20),
        plot.title = element_text(size=20, face="bold", hjust = 0.5) ,
        axis.text=element_text(size=14,angle = 0, hjust = 0.5),
        strip.text = element_text(size=10, face="bold", hjust = 0.5) ,
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.8, "cm"))
  

# Display the plot
print(Psum_Aarhus)
ggsave("figures/fraction_total_Aarhus.png", plot = total_plot, width = 12, height = 7)
ggsave("figures/fraction_t0_Aarhus.png", plot = t0plot, width = 10, height = 5)
