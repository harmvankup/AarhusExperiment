library(readxl)
library(tidyverse)
library(ggplot2)
library(car)
library(FSA)
#library(lemon)




#### values and functions ####
mmP <- 30.97376
mmS <- 32.065
mmFe <- 55.845
mmAl <- 26.981
mmCa <- 40.078
mmMn <- 54.938
mmTi <- 47.867
mmMg <- 24.30506


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

#### First import all raw data. The concentrations of P extraction solutions, and the Dryweight data. ####
conc_seq_raw <- read_xlsx("input/SequentialExtractions.xlsx", sheet = "sediment")
conc_ferrosorb <- read_xlsx("input/SequentialExtractions.xlsx", sheet = "ferrosorb")
conc_Bipy_raw <- read_xlsx("input/SequentialExtractionsBipy.xlsx")
wt_raw <- read_xlsx("input/Dryweight.xlsx")
times <- read_xlsx("input/timeseries.xlsx")
kwa_raw <- read_xlsx("input/KWA.xlsx")


#### calculations ####
# Define porosity values for all measurements
porosity <- wt_raw %>%
  group_by(Sample) %>%
  summarise(Porosity = mean(Porosity))

LOI <- wt_raw %>%   group_by(Timestep,Sample) %>%
  summarise(diff_porosity = abs(diff(Porosity)),
            diff_LOI = abs(diff(LOI)),
            porosity = mean(Porosity),
            LOI = mean(LOI))

meandiff_LOI <- LOI %>%  pull(diff_LOI) %>% mean()
sddiff_LOI <- LOI %>% pull(diff_LOI) %>% sd()
#### ferrosorb calculations ####

ferrosorb<- conc_ferrosorb %>%
  left_join(porosity, by = c("Sample")) %>%
  mutate(volume = get_volume(Fraction),
         Cont =  case_when( 
           Element == "P" ~ ((Conc - Blank) * volume) / (Weight * (1 - porosity)),
           Element == "Fe"~ 1000*(((Conc - Blank) * volume) / (Weight * (1 - porosity)))/mmFe
           )) %>%
  select(-Blank) %>% 
  transform( factor(Fraction, 
            levels = c("H2O","BD","NaOH")))

ferrosorb_calc <- ferrosorb %>% select(-porosity,-volume) %>% 
  pivot_wider(names_from = Element, values_from = c(Conc,Cont)) %>% 
  mutate(FeP = Cont_Fe/Cont_P)

ferrosorb_stat <- ferrosorb_calc %>%  group_by(Fraction) %>%  
  summarise(sd_P = sd(Cont_P), Cont_P=mean(Cont_P), sd_Fe = sd(Cont_Fe), Cont_Fe=mean(Cont_Fe),sd_FeP = sd(FeP), FeP = mean(FeP), .groups = "drop")

#####1st incubation #####
conc_seq_1st <- conc_seq_raw %>%
  left_join(porosity, by = c("Sample")) %>%
  left_join(times, by = c("Sample")) %>%
  mutate(volume = get_volume(Fraction),
         Cont = case_when( Parameter == "Fe" ~ ((Conc -Blank)/mmFe * volume * 1000) / (Weight * (1 - Porosity)),
                           T ~ ((Conc - Blank) * volume) / (Weight * (1 - Porosity)))) %>%
  select(-Blank,-volume,-Conc) %>% 
  pivot_wider(names_from = c(Fraction,Parameter), values_from = Cont)  %>% 
  mutate(
    BD_FeP = BD_Fe/BD_P,
    NaOH_FeP = NaOH_Fe/NaOH_P) %>% 
  pivot_longer(cols= H2O_P:NaOH_FeP, names_to = c("Fraction","Parameter"), values_to = "Cont", values_drop_na = T,names_sep = "_") %>% 
  mutate(
     Fraction = factor(Fraction, levels = c("H2O","BD","NaOH"))) %>% 
  group_by(Time, Fraction, Sample, Parameter) %>% summarise(sd = sd(Cont), Cont=(mean(Cont)), .groups = "drop")

calc_1st <- conc_seq_1st %>% 
  pivot_wider(names_from = c(Fraction,Parameter), values_from = c(Cont,sd)) %>% 
  mutate(NaOHperc = Cont_NaOH_P/(Cont_NaOH_P + Cont_BD_P))

BD_start <- conc_seq_1st %>% filter(Time == 8, Fraction == "BD", Sample == "FeP") %>% pull(Cont)
NaOH_end <- conc_seq_1st %>% filter(Time == 42, Fraction == "NaOH", Sample == "FeP") %>% pull(Cont)
BD_end <- conc_seq_1st %>% filter(Time == 42, Fraction == "BD", Sample == "FeP") %>% pull(Cont)

NaOH_end/BD_start
BD_end/BD_start

conc_porewater_anox <- conc_seq_raw %>%  
  left_join(porosity, by = c("Timestep","Sample")) %>%
  filter(Fraction == "H2O", Parameter == "P") %>% 
  mutate(conc = ((Conc - Blank) * (0.0202+ 0.001*Weight)) / (Weight * porosity))%>% 
  group_by(Sample) %>%  summarise(Conc = mean(conc), sd = sd(conc))

# Combine concentration data with porosity and blanks
# Compute standard deviation for each combination of Time, Fraction, and Sample
# Plot the sequential extractions with error bars

legend_labels_frac_1st <- c("H2O", "BD", "NaOH SRP")

seqextr_1st <- ggplot(conc_seq_1st, aes(x = Time, y = Cont, fill = Fraction)) +
  geom_col(position = position_stack(reverse = TRUE), width = 2.0) +
  #geom_errorbar(data = conc, aes(ymin = y_pos - sd, ymax = y_pos + sd), width = 1, position = position_dodge(0.1)) +
  scale_fill_manual(values = c("palegreen2", "darkorange2", "cadetblue"), labels = legend_labels_frac_1st, name = "Extracted Fraction") +
  labs(y = expression(paste("Extracted P (umol gDW"^{-1}, ")")),
       x = "Time (days)") +
  facet_grid( ~ Sample ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16, angle = 0, hjust = 0.5),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.3, "cm")
  )
show(seqextr_1st)

#####---- 2nd incubation ----#####

Time_Colors <- c("#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442")
Treatment_Colors <- c("#E69F00","#56B4E9")
Treatment_names <- c("Start", "Sulfate treatment", "Control")
Time_labels <- c(
  ts0 = "start" ,
  ts1 =  "19 days",
  ts2 = "62 days",
  ts3 = "96 days"
)
#### KWA ####

kwa <- kwa_raw %>% 
  pivot_longer(cols =  Al:Zn, names_to = "Element", values_to = "Conc") %>% 
  mutate(mm = get_mm(Element),
        Total = (0.05*Conc*1000)/(Weight*mm),
         Element = paste("Cont_Total_",Element, sep = "")) %>% 
  select(Time,Treatment,Total, Depth, Element) %>% 
  pivot_wider(names_from = Element, values_from = Total) %>% 
  mutate(Cont_Total_FeP = Cont_Total_Fe/Cont_Total_P,
         Cont_Total_FeS = Cont_Total_Fe/Cont_Total_S,
         Cont_Total_SFe = Cont_Total_S/Cont_Total_Fe)

totalElements <- c("Total_P","Total_S", "Total_Fe", "Total_Ca","Total_Mg","Total_Mn",  "Total_Al","Total_Ti", "Total_FeP", "Total_FeS","Total_SFe")
Elements <- totalElements %>% str_remove("Total_")


kwa_plotdata <- kwa %>% pivot_longer(Total_Al:Total_SFe, names_to = "Parameter", values_to = "Cont") %>% 
  mutate(Parameter = str_remove(Parameter, "Total_")) %>% 
  filter(Parameter %in% Elements, Time != "ts0") 

plot_kwa <- function(dat,elem,Isumolg = T){
  elem2 <- elem %>% str_remove("Total_")
  plotdat <- dat %>% filter(Parameter %in% elem2) %>% 
    transform( Parameter = factor(Parameter, levels = elem2),
               Treatment = factor(Treatment, levels = c("Start","Sulfate","Control"), labels = Treatment_names),
               Time = factor(Time, levels = c("ts0","ts1","ts2","ts3"),labels = Time_labels ))
  plot <-  ggplot( plotdat, 
                   mapping =  aes( x=Depth, y = Cont, color = Treatment, shape = Treatment)) +
    geom_line()+
    geom_point(size= 2) +
    scale_color_manual(values = Treatment_Colors, name = "treatment")+
    scale_shape_manual(values = c(19,17,15), name = "treatment") +
    coord_flip()+
    scale_x_reverse(limits = c(6,0)) + 
    #scale_y_Continuous(breaks = equal_breaks(3,0.13,0)) +
    labs(y = if(Isumolg == T) expression(paste(" concentration (", mu,"mol gDW"^-1,") " ))
             else "Molar ratio",
         x = "Depth (cm)"
    ) +
    facet_grid(Time ~ Parameter, scales = "free_x") +
    theme_classic()+
    theme(
      axis.title=element_text(size=10),plot.title = element_text(hjust = 0.5),
      axis.text=element_text(size=10,angle = 0, hjust = 0.5),
      panel.background = element_rect(colour = "grey20",linewidth = 1),
      strip.text = element_text(size=13, face="bold", hjust = 0.5) ,
      strip.background = element_blank(),
      legend.position = "bottom",
      legend.title = element_text( size = 10),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, "cm"),
      legend.background = element_rect(colour='grey')) 
  return(plot)
}


kwa_plot1 <- plot_kwa(c("Total_P","Total_S", "Total_Ca","Total_Mg","Total_Mn", "Total_Fe"))

ratio_plot <- plot_kwa(c("Total_FeP", "Total_FeS", "Total_SFe"))
ggsave(kwa_plot1, file="figures/elements.png",width = 10, height = 5)
show(ratio_plot)
ggsave(ratio_plot, file="figures/ratio.png",width = 10, height = 5)
#### Parallel Sequential extractions. ####

conc_extractions <- conc_Bipy_raw %>%
  mutate(volume = get_volume(Fraction),
         Cont = case_when( Parameter == "Fe" ~ ((Conc - Blank)/mmFe * volume * 1000) / (Weight * (1 - Porosity)),
           T ~ ((Conc - Blank) * volume) / (Weight * (1 - Porosity)))) %>% 
  select(-Blank) %>% 
  filter(Time != "ts3" | Treatment != "Control" | Depth != 5.4, # remove outliers
         Time != "ts3" | Treatment != "Control" | Core != 3)%>% #filter() %>% 
  filter(Fraction != "BD" | Parameter != "SRP",
         Fraction != "Bipy" | Parameter != "SRP") %>% 
  select(Treatment,Depth,Core,Nr,Extraction,Time,Cont,Fraction, Parameter) %>% 
  pivot_wider(names_from = c(Fraction,Parameter), values_from = Cont)  %>% 
  mutate(
    BD_FeP = BD_Fe/BD_TP,
    NaOH_FeP = NaOH_Fe/NaOH_TP,
    Bipy_FeP = Bipy_Fe/Bipy_TP) %>% 
  pivot_longer(cols= H2O_SRP:Bipy_FeP, names_to = c("Fraction","Parameter"), values_to = "Cont", values_drop_na = T,names_sep = "_")


#### some statistics ####
conc_test_bipy <- conc_extractions %>% pivot_wider(names_from = c(Fraction,Extraction), values_from = Cont)%>% 
  group_by( Time,Treatment) %>% summarize(n =n(), test = t.test(Bipy_Bipyridine, BipySRP_Bipyridine, paired = TRUE)$p.value) # significantly different so we leave out the BipySRP




conc__norm_test <- conc_extractions %>% 
  filter(Extraction == "Normal",
         Treatment != "Start", 
         Parameter %in% c("SRP","TP"), 
         Time != "ts2",
         Depth != 0) %>% 
  group_by(Fraction,Parameter, Time, Treatment, Depth) %>% 
  do({
    data_group <- .  # This stores the current group data
    n <- nrow(data_group)
    shapiro_test <- if (n >= 3) {
      shapiro_test_result <- shapiro.test(data_group$Cont)
      shapiro_test_result$p.value
    } 
    else NA  
    data.frame(
      n = n, 
      p_value = shapiro_test,
      sd = sd(data_group$Cont), 
      Cont = mean(data_group$Cont)
    )
  }) %>% 
  ungroup()

leveneTest(Cont ~ Time, data = conc_extractions)


conc_test_group <- conc_extractions  %>% 
  filter(Treatment != "Start", Extraction == "Normal", Time != "ts2", Depth != 0) %>% 
  mutate(group = paste(Time,Treatment)) 

conc_test_data <- conc_test_group %>% 
  group_by( Parameter,Fraction, Depth)%>% 
  do({
    # Perform Kruskal-Wallis test for each group
    kruskal_result <- kruskal.test(Cont ~ group, data = .)
    anova_result <- rev(summary(aov(Cont ~ group, data = .))[[1]])[[1]][1]
    signif <- case_when(
      kruskal_result$p.value<0.001 ~ "***",
      kruskal_result$p.value<0.01 ~ "**",
      kruskal_result$p.value<0.05 ~ "*",
              T ~ "")
    # Create a data frame with the results to return
    data.frame(
      n = nrow(.),
      p_value = kruskal_result$p.value,
      statistic = kruskal_result$statistic,
      df = kruskal_result$parameter,
      significance = signif ,
      anova = anova_result
      
    )
  }) %>% 
  ungroup()  # Ungroup at the end if needed

conc_test_wilcox_time <-
  conc_extractions  %>% 
  filter(Treatment != "Start", Extraction == "Normal", Time != "ts2", Depth != 0, Parameter %in% c("SRP","TP"), Depth<5) %>% 
  select(Parameter,Fraction, Depth, Treatment, Cont, Time) %>% 
  group_by( Parameter,Fraction, Depth, Treatment)%>% 
  do({
    data_group <- . 
    # Perform  test for each group
    ts1 <- data_group %>% filter(Time == "ts1") %>% pull(Cont)
    ts3 <- data_group %>% filter(Time == "ts3") %>% pull(Cont)
    wilcox_result <- wilcox.test(ts1, ts3, exact=F)
    
    signif <- case_when(
      wilcox_result$p.value<0.001 ~ "***",
      wilcox_result$p.value<0.01 ~ "**",
      wilcox_result$p.value<0.05 ~ "*",
      T ~ "")
    # Create a data frame with the results to return
    data.frame(
      n = nrow(.),
      p_value = wilcox_result$p.value,
      statistic = wilcox_result$statistic,
      #df = wilcox_result$parameter,
      significance = signif,
      ts1 =length(ts1),
      ts3 =length(ts3)
    )
  }) 

conc_test_wilcox_treatment <-
  conc_extractions  %>% 
  filter(Treatment != "Start", 
         Extraction == "Normal", 
         Time != "ts2", Depth != 0, 
         Parameter %in% c("SRP","TP"), 
         Depth<5) %>% 
  select(Parameter,Fraction, Depth, Treatment, Cont, Time) %>% 
  group_by( Parameter,Fraction, Depth, Time) %>%
  do({
    data_group <- . 
    # Perform  test for each group
    t1 <- data_group %>% filter(Treatment == "Sulfate") %>% pull(Cont)
    t3 <- data_group %>% filter(Treatment == "Control") %>% pull(Cont)
    wilcox_result <- wilcox.test( t3,t1, exact=F, paired =F)
    stat <- wilcox_result$statistic %>% as.numeric()
    signif <- case_when(
      wilcox_result$p.value<0.001 ~ "***",
      wilcox_result$p.value<0.01 ~ "**",
      wilcox_result$p.value<0.05 ~ "*",
      T ~ "")
    # Create a data frame with the results to return
    data.frame(
      n = nrow(.),
      p_value = wilcox_result$p.value,
      statistic = stat,
      #df = wilcox_result$parameter,
      significance = signif,
      t1 =length(t1),
      t3 =length(t3)
    )
  }) 

conc_test_wilcox_treatment[[7]]

conc_extractions  %>% 
  filter(Treatment == "Sulfate", 
         Extraction == "Normal", 
         Time == "ts3",  
         Parameter %in% c("SRP"), 
         Fraction == "NaOH",
         Depth == 0.6) %>% pull(Cont) 

conc_extractions  %>% 
  filter(Treatment == "Control", 
         Extraction == "Normal", 
         Time == "ts3",  
         Parameter %in% c("SRP"), 
         Fraction == "NaOH",
         Depth == 0.6) %>% pull(Cont)

wilcox_result <- wilcox.test(conc_extractions  %>% 
                               filter(Treatment == "Control", 
                                      Extraction == "Normal", 
                                      Time == "ts3",  
                                      Parameter %in% c("SRP"), 
                                      Fraction == "NaOH",
                                      Depth == 0.6) %>% pull(Cont),conc_extractions  %>% 
                               filter(Treatment == "Sulfate", 
                                      Extraction == "Normal", 
                                      Time == "ts3",  
                                      Parameter %in% c("SRP"), 
                                      Fraction == "NaOH",
                                      Depth == 0.6) %>% pull(Cont),exact=F , paired =F)$statistic %>% as.numeric()

data <- conc_test_group %>% filter(Parameter == "TP", Fraction == "BD", Depth == 1.8)
anova_result <- aov(Cont ~ group, data = data)
summary(anova_result)
TukeyHSD(anova_result)
# Summary of ANOVA
test <- rev(summary(anova_result)[[1]])[[1]][1]

dunn_results <- function(p,f,d){
data <- conc_test_group %>% filter(Parameter == p, Fraction == f, Depth == d)
dunnTest(Cont ~ group, data = data, method = "bh")
}  

pw_results <- function(p,f,d){
  data <- conc_test_group %>% filter(Parameter == p, Fraction == f, Depth == d)
  pairwise.wilcox.test(data$Cont, data$group)
} 



pw_results("SRP","NaOH",0.6)
dunn_results("SRP","NaOH",0.6)
dunn_results("SRP","NaOH",1.8)
dunn_results("SRP","NaOH",3.0)

dunn_results("TP","BD",1.8)

dunn_results("TP","NaOH",0.6)
dunn_results("TP","NaOH",1.8)
dunn_results("TP","NaOH",3.0)

# Print the results


#### calculations ####
conc_data <- conc_extractions %>% 
group_by( Time,Treatment,Extraction,Fraction,Parameter,Depth) %>%
summarise(sd = sd(Cont), n = n(), Cont=(mean(Cont))) %>% 
pivot_wider(names_from = c(Fraction,Parameter), values_from = c(Cont, sd,n))  %>% 
left_join(kwa) %>% 
  mutate(
    Cont_NE_P = case_when( 
      is.na(Cont_Bipy_TP) ~ Cont_Total_P - Cont_NaOH_TP - Cont_BD_TP-Cont_H2O_SRP,
      T ~ Cont_Total_P- Cont_NaOH_TP - Cont_Bipy_TP - Cont_BD_TP-Cont_H2O_SRP), 
    Cont_NaOH_NRP = Cont_NaOH_TP-Cont_NaOH_SRP,
    sd_NE_P = case_when( 
      is.na(Cont_Bipy_TP) ~ sqrt( 0#sd_Total_P
      + sd_NaOH_TP^2 + sd_BD_TP^2 + sd_H2O_SRP^2),
      T ~ sqrt( 0#sd_P
      + sd_NaOH_TP^2 + sd_Bipy_TP^2 + sd_BD_TP^2 + sd_H2O_SRP^2 )),
    sd_NaOH_NRP = sqrt(sd_NaOH_TP^2+sd_NaOH_SRP^2)) %>% 
  pivot_longer(cols= Cont_BD_Fe:sd_NaOH_NRP, names_to = c("aspect","Fraction","Parameter"), values_to = "value", values_drop_na = T,names_sep = "_")%>% 
  pivot_wider(names_from = aspect, values_from = value) #select(-NaOHTP) %>% 

conc_porewater <- conc_Bipy_raw %>% filter(Fraction == "H2O", Parameter == "SRP") %>% 
  mutate(conc = ((Conc - Blank) * (0.0202+ 0.001*Weight)) / (Weight * Porosity)) %>% 
  group_by(Treatment, Depth, Time) %>%  summarise(Conc = mean(conc), sd = sd(conc))




 #### Plots ####
# Labels and colors
legend_labels <- c("H2O", "BD", "NaOH SRP", "NaOH NRP", "recalcitrant P")
colors <- list(c("palegreen2", "darkorange2", "cadetblue","coral4","darkgrey"))

legend_labels_Bipy <- c("H2O","Bipy", "BD", "NaOH SRP", "NaOH NRP")
colors_Bipy <- list(c("palegreen2","red", "darkorange2", "cadetblue","coral4","darkgrey"))

legend_labels_Bipy <- c("H2O","Bipy", "BD", "NaOH SRP", "NaOH NRP", "recalcitrant P")

plotdata <- transform(filter(conc_data, Extraction == "Normal" )) %>% 
  pivot_wider(names_from = Fraction, values_from = c(Cont, sd, Total))

#### Bipy plots ####

#only Bipy
Bipy_start = filter(conc_extractions, Fraction =="Bipy", Treatment == "Start") %>% pull(Cont)
Bipy_start_line = data.frame(Cont = c(Bipy_start[1],Bipy_start[1],Bipy_start[2],Bipy_start[2]),
                             Depth =c(0,2,2,6),
                             Time = c("t=0","t=0","t=0","t=0"),
                             Treatment = factor(c("Start","Start","Start","Start")))

bipyprofile_Aarhus <- ggplot( transform(filter(conc_data, Extraction == "Normal" ) 
                                        %>% pivot_wider(names_from = Fraction, values_from = c(Cont, sd, Total)),
                  Treatment = factor(Treatment, levels = c("Start","Sulfate", "Control"))), 
        mapping = aes(
          y = Cont_NaOHSRP,
          by = Depth,
          color = Treatment,
          x = Total_P )) +
  #geom_line( size = 1)  +
  geom_point(aes(shape = Time), size = 2) +
  #geom_path(data = Bipy_start_line,aes(linetype = Time),size = 1)  +
  #geom_boxplot(aes(y = Cont,x = Depth)) +
  xlab("Depth (cm)") +
  scale_x_reverse() +
  coord_flip() +
  #facet_grid( Time ~Treatment ) +
  #scale_color_discrete(breaks=c("Start","Control","Sulfate"))+
  #scale_linetype_discrete(breaks = c("t=0","ts1","ts2", "ts3"), labels= c("Start","19 days", "48 days", "82 days"))+
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
show(bipyprofile_Aarhus)
ggsave("figures/allfractions_Aarhus.png", plot = bipyprofile_Aarhus, width = 20, height = 5)

#### plot the ratios  ####
plot_data_ratios <- conc_data %>% filter(Extraction == "Normal",
                                             Time != "ts0",  Time != "ts2", 
                                             Parameter %in% c("FeP","FeS","SFe") ) %>% 
  mutate(Parameter = paste0(Fraction,Parameter))

plot_ratios <- plot_kwa(plot_data_ratios,c("NaOHFeP","BDFeP","TotalFeP","TotalFeS","TotalSFe"),F)

show(plot_ratios)


t0plot <- ggplot( transform( filter(conc_data, Time == "ts0", Fraction %in% fraction_list), 
                   Fraction = factor(Fraction, 
                                     levels = fraction_list)), 
        aes( x=Depth, y = Cont, fill = Fraction),
        scale_x_discrete(position = 'top')) +
  geom_bar(
    mapping = aes( y = Cont, fill = Fraction ),
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

#### Extractions, No Bipy #####

fraction_list <- c("H2OSRP","BipyTP","BDTP","NaOHSRP","NaOHNRP", "NEP")
seq_plot_data <-  
  conc_data %>% 
    mutate(Fraction = paste0(Fraction,Parameter)) %>% 
    filter( Time != "ts0",  Time != "ts2", 
                                   Extraction != "Bipyridine", 
                                   Depth != 0 , Fraction %in% fraction_list ) %>%  
    select(Time, Treatment, Extraction, Depth, Fraction, Cont, sd, n) %>% 
    transform( 
           Extraction = factor(Extraction, levels = c("Normal", "Bipyridine"), labels = c("Normal", "Bipyridine")),
           Fraction = factor(Fraction, 
                             levels = fraction_list),
           Time = factor(Time, levels = (c("ts0", "ts1","ts2","ts3")), labels = (c("start", "19 days"," more days","96 days"))),
           Treatment = factor(Treatment, levels = (c("Start", "Control","Sulfate"))))

seq_plot_data$y_pos = NA
seq_plot_data$y_pos[seq_plot_data$Fraction == "H2OSRP"] = seq_plot_data$Cont[seq_plot_data$Fraction == "H2OSRP"]

seq_plot_data$y_pos[seq_plot_data$Fraction == "BDTP"] = seq_plot_data$y_pos[seq_plot_data$Fraction == "H2OSRP"] + 
  seq_plot_data$Cont[seq_plot_data$Fraction == "BDTP"]
seq_plot_data$y_pos[seq_plot_data$Fraction == "NaOHSRP"] = seq_plot_data$y_pos[seq_plot_data$Fraction == "BDTP"] + 
  seq_plot_data$Cont[seq_plot_data$Fraction == "NaOHSRP"]
seq_plot_data$y_pos[seq_plot_data$Fraction == "NaOHNRP"] = seq_plot_data$y_pos[seq_plot_data$Fraction == "NaOHSRP"] + 
  seq_plot_data$Cont[seq_plot_data$Fraction == "NaOHNRP"]
seq_plot_data$y_pos[seq_plot_data$Fraction == "NEP"] = seq_plot_data$y_pos[seq_plot_data$Fraction == "NaOHNRP"] + 
  seq_plot_data$Cont[seq_plot_data$Fraction == "NEP"]

total_plot <- ggplot( seq_plot_data, 
        aes( x=Depth, y = Cont, fill = Fraction),
        scale_x_discrete(position = 'top')) +
  geom_bar(
    mapping = aes( y = Cont, fill = Fraction ),
    stat = "identity", position = position_stack(reverse=TRUE))+
  scale_fill_manual(values = unlist(colors), labels = legend_labels, name = "Extracted Fraction") +
  geom_errorbar(data = seq_plot_data, aes(ymin = y_pos - sd, ymax = y_pos + sd), width = 0.5, position = position_dodge(0.2))+
  coord_flip()+
  scale_x_reverse() + 
  labs(y = expression(paste("Extracted P (",mu,"mol gDW"^-1,")" )),
       x = "Depth (cm)"
  ) +
  facet_grid(Treatment ~ Time, labeller =label_value) +
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

concsum <- conc_extractions %>% 
  filter(Depth<4.5, Treatment != "Start",Fraction != "BipySRP") %>% 
  group_by(Extraction, Core, Time, Treatment,Fraction,Depth) %>%
  summarise( Cont=(mean(Cont)), sd = 0, .groups = "drop") %>% 
  pivot_wider(names_from = Fraction, values_from = c(Cont,sd))  %>% 
  left_join(kwa)  %>% 
  mutate(Cont_NE = case_when( 
    is.na(Cont_Bipy) ~ Cont_P-Cont_NaOHTP-Cont_BD-Cont_H2O,
    T ~ Cont_P- Cont_NaOHTP - Cont_Bipy - Cont_BD-Cont_H2O
  ), Cont_NaOHNRP = Cont_NaOHTP-Cont_NaOHSRP) %>% 
  pivot_longer(cols= Cont_BD:Cont_NaOHNRP, names_to = c("Parameter","Fraction"), values_to = "value", values_drop_na = T,names_sep = "_")%>% 
  pivot_wider(names_from = Parameter, values_from = value) %>%  
  filter( Fraction %in% c("H2O","Bipy","BD","NaOHSRP","NaOHNRP", "NE","P")) %>% 
group_by( Core, Time,Treatment,Extraction, Fraction) %>%
  summarise(Cont = sum(Cont))

concsum %>% filter(Extraction != "Bipyridine") %>% view()

Psum_Aarhus <- 
  ggplot( transform(filter(concsum, Extraction != "Bipyridine") ,
                      Fraction = factor(Fraction,levels = c("H2O","Bipy","BD","NaOHSRP","NaOHNRP", "NE","P")),
                    Core = factor(Core,levels = c(1,2,3,4,5)),
                      Treatment = factor(Treatment, levels = c("Sulfate", "Control"))), 
                              mapping = aes(
                                y = Cont,
                                fill = treatment,
                                x = Core)) +
  geom_col(mapping = aes(
    y = Cont,
    fill = Treatment), position = position_dodge())+
  #geom_errorbar(aes(ymin = Cont - sd, ymax = Cont + sd),position = position_dodge(width = 0.1),width = 0.25) +
  xlab("Time") +
  facet_wrap( ~ Fraction , ncol = 2,scales = "free_y" ) +
  labs(  y = expression(paste("concentration (umol/g dw)")), 
         x = "time",
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
print(total_plot)
ggsave("figures/fraction_total_Aarhus.png", plot = total_plot, width = 12, height = 7)
ggsave("figures/fraction_t0_Aarhus.png", plot = t0plot, width = 10, height = 5)
