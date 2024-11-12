library(readxl)
library(tidyverse)
library(ggplot2)

xrd_raw <- read_xlsx("input/XRD.xlsx")
#depths <- read_xlsx("input_data/depths.xlsx")

# Deepvalues <- xrd_raw %>% filter(Sample == "DEEP") %>%   pivot_longer(cols = 2:4, names_to = "peak", values_to = "Deep") %>%  select(!Sample)
#AAR_t0_treated_1mmDS_3mmAS_2.5Fefilt_2.5secsol_1h_10rpm_10-90_0.014_0.6s-step_1h_874.brml #1

xrd <- xrd_raw %>% 
  separate(Scan, into = c("Lake", "Time",  "Treatment","Depth","DS","AS","Fefilt","SS","meastime","rpm","range","stepsize","steplength","col","nr"), sep = "_") %>% 
  mutate(Depth = case_when(
    Depth == "0" ~ 0,
    Depth == "0-1.2" ~ 0.6,
    Depth == "1.2-2.4" ~ 1.8,
    Depth == "2.4-3.6" ~ 3.0,
    Depth == "3.6-4.8" ~ 4.2,
    Depth == "4.8-6" ~ 5.4
  ),
  Treatment = str_remove(Treatment, "[15]"))
  

xrd_stat <- xrd %>% 
  filter(Time == "t0",Treatment == "treated") %>% 
  mutate(Start = NetArea) %>% 
  select(Name,Start) %>% 
  right_join(xrd) %>% 
  mutate(Relative = NetArea/Start,
         Loss = (Start-NetArea)/Start*100)

xrd_stat %>% 
  filter(Time != "t0") %>% 
  select(Time,Depth,Treatment,Name,NetArea,Relative,Loss) %>% view()


  color_box <- tibble(xmin = c(2,0),
                      xmax = c(8,2.5),
                      ymin = c(-Inf,-Inf),
                      ymax = c(Inf,Inf),
                      fill = factor(c("D" ,"S"), levels = c("D","S")),
                      Treatment = factor(c("A","A"),levels = c("A","B","C")))

xrdplot <- ggplot(xrd, 
                  mapping = aes(
                    y = NetArea,
                    x = Depth,
                    color = Treatment
                  )) +
  geom_line() +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin = meannormarea - sd, ymax = meannormarea + sd), width = .5, position = position_dodge(0.05)) +
  xlab("Depth (cm)") +
  scale_x_reverse() +
  scale_y_continuous(
    name = expression(paste(" rel. peak area"))
  ) +
  coord_flip() +
  facet_grid( Time ~Name, scales = "free")+
  labs(
    x = "Depth in cm",
    title = "",
    shape = "Measurement"  # Adding legend title
  ) +
  #scale_shape_manual(name = "Legend",labels = c("XRD relative peak area", "Mössbauer spectroscopy"),values = c(1, 2, 3)) + 
  theme_bw(base_size = 8) +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 8, angle = 0, hjust = 0.5))
show(xrdplot)

ggsave( "figures/XRDplot.tiff", plot = xrdplot, device = "tiff",
width = 84, height = 84,  units = "mm", limitsize = FALSE, dpi = 600)
