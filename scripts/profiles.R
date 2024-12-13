
library(readxl)
library(tidyverse)
library(ggplot2)
library(scales)
library(cowplot)
library(ggh4x)
#### First import all raw data. The concentrations of P extraction solutions, and the Dryweight data. ####
profiles_sulfide_raw <- read_xlsx("input/Profiles.xlsx", sheet = "Sulfide")

profiles_oxygen_raw <- read_xlsx("input/Profiles.xlsx", sheet = "Oxygen")

correction <- read_xlsx("input/Profiles.xlsx", sheet = "Corrections")

porosity <- read_xlsx("input/Profiles.xlsx", sheet = "Porosity")

Treatment_Colors <- c("#E69F00","#56B4E9")

#### Calculations ####
profiles_oxygen <- profiles_oxygen_raw %>% 
  filter(Depth < -3000) %>% 
  group_by(Core, Timestep,Nr) %>% 
  summarise( Conc = mean(Signal),sd = sd(Signal), n = length(Signal)) %>% 
  right_join(profiles_oxygen_raw) %>% 
  remove_missing() %>% 
  filter(Core != "O1" | Timestep != "tp1" ) %>% 
  mutate(Treatment = substr(Core,1,1),
         DO = (Signal/Conc)*100,
         slope = (lead(DO) - DO) / (lead(Depth) - Depth),
         por = case_when( (slope/lag(slope)) <= 0 ~ 1,
         T ~  (((slope/lag(slope))+3)^(-1))*4 
         ),
         rate = case_when( (slope/lag(slope)) <= 0 ~ 0,
         T ~ ((slope-lag(slope))/(Depth-lag(Depth))))
         ) 

profiles_oxygen_corrected <- profiles_oxygen %>% 
  left_join(correction) %>% 
  mutate( Depthcor = Depth-Corr)

#### Statistics, Oxygen ####
depth_test <- profiles_oxygen_corrected %>% 
 ddply(.(Core, Treatment, Timestep, Nr), 
  summarise,
  depth = first(pull(subset(tibble(DO,Depthcor), DO < 1), Depthcor))) %>% 
  ddply(.(Treatment, Timestep), 
        summarise,
        sd = sd(depth, na.rm = T), depth = mean(depth, na.rm = T)) 

O <- profiles_oxygen_corrected %>% 
  ddply(.(Core, Treatment, Timestep, Nr), 
        summarise,
        depth = first(pull(subset(tibble(DO,Depthcor), DO < 1), Depthcor))) %>% 
  filter(Timestep == "ts3", Treatment == "O") %>% pull(depth)

S <- profiles_oxygen_corrected %>% 
  ddply(.(Core, Treatment, Timestep, Nr), 
        summarise,
        depth = first(pull(subset(tibble(DO,Depthcor), DO < 1), Depthcor))) %>% 
  filter(Timestep == "ts3", Treatment == "S") %>% pull(depth)
  

 wilcox.test(O, S, exact=F)

profiles_oxygen_test <- ggplot(  profiles_oxygen %>% 
                                   filter(Timestep == "ts3", 
                                          Core == "O3" , 
                                          Nr == 1
                                                            ) %>% 
                                   pivot_longer(cols = DO:rate, names_to =  "type", values_to = "value"), 
                                mapping = aes(
                                  y = value,
                                  x = Depth,
                                  label = Depth)) +
  geom_text(hjust=0, vjust=0)  +
  geom_point(size = 2) +
  scale_x_reverse() +
  coord_flip() +
  facet_grid(.~ type, scales = "free") +
  labs(  y = expression(paste("concentration ",mu,"M")), 
         title = "dissolved oxygen porewater profiles" ) +
  theme_bw() +
  theme( text = element_text(size = 10))

print(profiles_oxygen_test)

#### Figures Oxygen ####
profiles_oxygen_figdata <- profiles_oxygen_corrected %>% 
  filter(Core != "Ocontrol", Core!= "Scontrol") %>% 
  group_by(Treatment, Depthcor, Timestep) %>% 
  summarise( 
  sd = sd(DO), DO = mean(DO)) %>% filter(Timestep == "ts3")%>% 
  transform(Treatment = factor(Treatment, levels = c("S", "O"), labels = c("sulfate treatment","control")),
            Depthcor = Depthcor/10000)

profiles_oxygen_fig <- ggplot( profiles_oxygen_figdata , 
                     mapping = aes(
                       y = DO,
                       x = Depthcor,
                       color = Treatment)) +
  #geom_line()  +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = DO - sd, ymax = DO + sd),position = position_dodge(width = 0.01), width = 0.05)+
  scale_color_manual(values = Treatment_Colors, name = "treatment")+
  scale_x_reverse() +
  coord_flip() +
  labs(  y = expression(paste("concentration %")), 
         x = expression(paste("Depth in cm")),
         title = "dissolved oxygen porewater profiles" ) +
  #facet_grid(Treatment ~ Timestep, scales = "free_x")+
  theme_bw() +
  theme( text = element_text(size = 10))
print(profiles_oxygen_fig)

ggsave(profiles_oxygen_fig, file="figures/profiles_oxygen.png",width = 20, height = 20)
ggsave(profiles_oxygen_fig + theme_void(), file="figures/profiles_oxygen_simple.png", width = 10, height = 20)

#### Calculations sulfide ####
profiles_sulfide_raw %>% filter(time== "ts3", Core %in% c("s3","s5"), Depth <= 2000) %>% view()
 profiles_conc <- profiles_sulfide_raw %>% 
   mutate(Treatment = substr(Core,1,1),
          #Conc = case_when(Conc <0 ~ 0, T ~ Conc)
          ) %>% 
   remove_missing() %>% 
   mutate(totsulf = Conc + ((10^(-6.98))/(10^(-pH)))*Conc ) %>% 
   filter(Core!= "ocontrol", Core != "scontrol") %>% 
   group_by(Depth, Treatment, time) %>% 
   summarise(Conc = mean(totsulf),sd = sd(totsulf))
   
#### Figures sulfide ####
 legend_names <- c(
   ts0 = "start" ,
   ts1 =  "19 days",
   ts2 = "62 days",
   ts3 = "96 days"
 )
 
figure_data <- profiles_conc %>% 
  filter(Depth >= 0) %>% 
  transform(Treatment = factor(Treatment, levels = c("s", "o"), labels = c("sulfate treatment","control")),
            Depth = Depth/10000)

profiles <- ggplot(  figure_data, 
                     mapping = aes(
                       y = Conc,
                       x = Depth,
                       color = Treatment)) +
  #geom_line()  +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = Conc - sd, ymax = Conc + sd),
    position = position_dodge(width = 0.01),
    width = 0.05)+
  xlab("Depth in cm)") +
  scale_x_reverse(limits = c(4.5, -0.5)) +
  scale_color_manual(values = Treatment_Colors, name = "treatment")+
  coord_flip() +
  labs(  y = expression(paste("concentration ",mu,"M")), 
         x = "Depth in cm",
         title = "dissolved sulfide porewater profiles" ) +
  facet_grid(.~ time, labeller = labeller(time = legend_names))+
  theme_bw() +
  theme(      text = element_text(size = 10))
print(profiles)

ggsave(profiles, file="figures/profiles.png",width = 10, height = 5)
show(kwaplots)

plot_sulf <- profiles + theme(legend.position = "none")
plot_ox <- profiles_oxygen_fig + theme(legend.position = "none")
plot_dgt <- Fe_profiles + theme(legend.position = "none")

legend <- get_legend(
  # create some space to the left of the legend
  profiles + theme(legend.box.margin = margin(0, 0, 0, 12))
)

top_row <- plot_grid(
  plot_sulf, plot_dgt,
  labels = c("A", "B"),
  align = "h",
  axis = "bt",
  rel_widths = c(2, 1),
  nrow = 1
)
bottom_row <- plot_grid(
  plot_ox,legend, 
  labels = c("C", ""),
  rel_widths = c(2, 1),
  nrow = 1
)

gridplot <- plot_grid(top_row, bottom_row, labels = c(""), ncol = 1)
ggsave("figures/all_profiles.tiff", plot = gridplot, device = "tiff",
       width = 174, height = 200,  units = "mm", limitsize = FALSE, dpi = 600)

#### calculations flux and reaction rates ####

depth_interpolation <-function(time, core, nr){

  concentration_data <-  profiles_oxygen_corrected %>% filter(Timestep == time, Core == core, Nr == nr)
  porewater_data <- porosity %>% filter(Core == core, Nr == 1)
fine_depth <- concentration_data$Depthcor/10000

# Interpolate porewater values to the finer grid
porewater_interpolated <- ifelse(fine_depth < 0, 1, 
                                 approx(
  x = porewater_data$Depth,        # Original depths (coarse grid)
  y = porewater_data$Porosity,    # Porewater values
  xout = fine_depth,               # Interpolated depths (fine grid)
  method = "linear",                # Linear interpolation,'
  rule = 2
)$y
)

# Combine into a single data frame
aligned_data <- data.frame(
  Depth = fine_depth,
  Porosity = porewater_interpolated,
  Biodiffusity = rep(0, length(fine_depth)),
  Irrigation = rep(0, length(fine_depth)),
  Concentration = concentration_data$DO
  
)

# View the combined data frame
print(aligned_data)
}

test <- depth_interpolation("tp1","O2",1)
