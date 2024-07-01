
library(readxl)
library(tidyverse)
library(ggplot2)

# First import all raw data. The concentrations of P extraction solutions, and the Dryweight data.
conc_raw <- read_xlsx("input/Profiles.xlsx")

conc_raw %>% filter(time== "ts3", Core %in% c("s3","s5"), Depth <= 2000) %>% view()
 conc <- conc_raw %>% 
   mutate(Treatment = substr(Core,1,1),
          #Conc = case_when(Conc <0 ~ 0, T ~ Conc)
          ) %>% 
   remove_missing() %>% 
   mutate(totsulf = Conc + ((10^(-6.98))/(10^(-pH)))*Conc ) %>% 
   filter(Core!= "ocontrol", Core != "scontrol") %>% 
   group_by(Depth, Treatment, time) %>% 
   summarise(Conc = mean(totsulf),sd = sd(totsulf))
   

 legend_names <- c(
   ts0 = "start" ,
   ts1 =  "20 days",
   ts2 = "62 days",
   ts3 = "98 days"
 )
 
figure_data <- conc %>% 
  filter(Depth >= -0) %>% 
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
  scale_x_reverse() +
  coord_flip() +
  labs(  y = expression(paste("concentration ",mu,"M")), 
         x = "Depth in cm",
         title = "dissolved sulfide porewater profiles" ) +
  facet_grid(.~ time, labeller = labeller(time = legend_names))+
  theme_bw() +
  theme(      aspect.ratio = 2,
              axis.title=element_text(size=10),
              plot.title = element_text(size=20, face="bold", hjust = 0.5) ,
              axis.text=element_text(size=10,angle = 0, hjust = 0.5),
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10, face = "bold"),
              legend.key.size = unit(1.5, "cm"))
print(profiles)

ggsave(profiles, file="figures/profiles.png",width = 10, height = 5)
show(kwaplots)
