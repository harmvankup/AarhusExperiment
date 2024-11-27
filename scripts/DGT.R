library(haven)
library(tidyverse)
library(ggplot2)

t0 <-  ISOdate(2023,7,05)
tDGT <-  ISOdate(2023,7,20)
t1 <- ISOdate(2023,7,24)
t2 <- ISOdate(2023,9,5)
t3 <- ISOdate(2023,10,9)
t1-tDGT 
t2-tDGT 
t3 - tDGT
t3 - t0
Fe_raw <- read_xlsx("input/DGT_Fe.xlsx")

DGT_1 <-   read.delim(
  "input/microplates/20231130_P_DGT_1.txt", sep = ",", header = FALSE, colClasses = "character"
) 
DGT_2<-   read.delim(
  "input/microplates/20231130_P_DGT_2.txt", sep = ",", header = FALSE, colClasses = "character"
) 
DGT_3<-   read.delim(
  "input/microplates/20231130_P_DGT_3.txt", sep = ",", header = FALSE, colClasses = "character"
) 
DGT_11<-   read.delim(
  "input/microplates/20231201_P_DGT_1.txt", sep = ",", header = FALSE, colClasses = "character"
) 
DGT_22<-   read.delim(
  "input/microplates/20231201_P_DGT_2.txt", sep = ",", header = FALSE, colClasses = "character"
) 
DGT_33<-   read.delim(
  "input/microplates/20231201_P_DGT_3.txt", sep = ",", header = FALSE, colClasses = "character"
) 

mp_data_transform <-  function(dat){
df <- dat  %>% as.data.frame() %>% filter(!str_detect(V1,"SP"))
n <- c(1:12)
nam <- as.character(n)
l <-  list()
  for(x in 1:12){
    x1 <- names(df)[(x*2)-1]
    x2 <- names(df)[x*2]
    li <- unite(df, x, x1:x2, sep = '.', remove = TRUE) %>% select(x) %>% mutate(x = as.numeric(x))
    l[[x]] <- li
}
t <-as.data.frame(l) 
names(t) <- c(1:12)
t <- cbind(col = c("A","B","C","D","E","F","G","H"),t)
return(t)


}

extract_values <- function(dat, dil, DE, GH){
  cali <- dat %>% filter(col == "A" | col == "B") %>% 
    select(1:9) %>% pivot_longer(2:9,names_to = "nr") %>% 
    group_by(nr) %>%   summarise( value = mean(value)) %>% 
    cbind("uM" = c(0,0.5,2.5,5,10,15,20,25))
    
  linear <- lm( cali$value ~ cali$uM )
  intercept <- linear$coefficients[1]
  slope <- linear$coefficients[2]
  
  S <- dat %>% select(1:11) %>% 
    filter(col == "D" | col == "E") %>% 
    pivot_longer(2:11,names_to = "nr") %>% 
    mutate(nmcm2 = case_when(col == "D" & nr %in% c(1:5) ~ ((0.8+1.8*0.094)*10^(0)*(value-intercept)*dil/(0.9*slope*1.8)),
                          T ~  ((0.8+1.8*0.094*0.5)*10^(0)*(value-intercept)*dil/(0.9*0.5*slope*1.8))))%>% pull(nmcm2) 
  O <- dat %>% select(1:11) %>% 
    filter(col == "G" | col == "H") %>% 
    pivot_longer(2:11,names_to = "nr") %>% 
    mutate(nmcm2 = case_when(col == "G" & nr %in% c(1:5) ~ ((0.8+1.8*0.094)*10^(0)*(value-intercept)/(0.9*slope*1.8))*dil,
                        T ~  ((0.8+1.8*0.094*0.5)*10^(0)*(value-intercept)/(0.9*0.5*slope*1.8)*dil))) %>% pull(nmcm2) 
  
  
  t <-cbind( S[1:15],  O[1:15]) %>% as.data.frame() 
  names(t) <- c(DE, GH)
  test <- t
    #
  return(t)
}

P_1 <- mp_data_transform(DGT_1) %>% extract_values(2,"S1","O1")
#P_2 <- mp_data_transform(DGT_2) %>% extract_values("S2","O2")
#P_3 <- mp_data_transform(DGT_3) %>% extract_values("S1","O1")
P_11 <- mp_data_transform(DGT_11) %>% extract_values(20,"S2","O2")
P_22 <- mp_data_transform(DGT_22) %>% extract_values(20,"S3","S4")
P_33 <- mp_data_transform(DGT_33) %>% extract_values(20,"O3","O4")

DGT_P_table <- cbind(P_1,P_11,P_22,P_33, depth = c(-5,-4,-3,-2,-1,-0.5,0,.5,1,1.5,2,2.5,3,3.5,4))

DGT_P_data <- DGT_P_table %>% 
  pivot_longer(S1:O4, names_to = "core", values_to = "nmol_p_cm2") %>% 
  filter(!(core == "O3" & depth == 2)) %>% 
  mutate(treatment = substr(core,1,1),
         time = case_when( substr(core,2,2) == "1" ~ 4,
                           substr(core,2,2) == "2" ~ 47,
                           substr(core,2,2) %in% c("3","4") ~ 81),
         uM = (nmol_p_cm2*0.094)/(time*24*60*60*5.26*10^(-6)) ,
         nmol_p_cm2_s = nmol_p_cm2/(time*24*60*60),
         mol_p_cm2_y = 10^(-9)*nmol_p_cm2/(time/365)
         ) %>% 
  transform(treatment = factor(treatment, 
                               levels = c("S", "O"), 
                               labels = c("sulfate treatment","control")), 
            depth = depth) %>% 
    pivot_longer(cols = c("nmol_p_cm2","uM","nmol_p_cm2_s","mol_p_cm2_y"), names_to = "unit")


profiles <- ggplot(  DGT_P_data, 
                     mapping = aes(
                       y = value,
                       x = depth,
                       color = treatment, by = core)) +
  geom_line()  +
  geom_point(size = 2) +
  xlab("Depth in cm)") +
  scale_color_manual(values = Treatment_Colors)+
  scale_x_reverse() +
  coord_flip() +
  labs(  y = expression(paste("concentration ",n,"mol/cm"^2)), 
         x = "Depth in cm",
         title = " porewater profiles" ) +
  facet_grid(time ~ unit, scales = "free")+
  theme_bw() +
  theme(      aspect.ratio = 1,
              axis.title=element_text(size=10),
              plot.title = element_text(size=20, face="bold", hjust = 0.5) ,
              axis.text=element_text(size=10,angle = 0, hjust = 0.5),
              legend.title = element_text(size = 0),
              legend.position = "bottom",
              legend.text = element_text(size = 10, face = "bold"),
              legend.key.size = unit(1.5, "cm"))
print(profiles)
ggsave("figures/profile_DGT_P_Aarhus.png", plot = profiles, width = 10, height = 10)


DGT_Fe_table <- Fe_raw %>% mutate(SurfConc = case_when(ID %in% c(1:5) ~ (0.8+1.8*0.094)*10^(3)*FeFF/(0.9*1*mmFe*1.8),
                                     T ~   (0.8+1.8*0.094*0.5)*10^(3)*FeFF/(0.9*0.5*mmFe*1.8)),
                                  PwConc = (SurfConc*0.094)/(81*24*60*60*5.32*10^(-6)) ) %>% 
  select(ID, Core, SurfConc, PwConc) %>% 
  pivot_wider(names_from = Core, values_from = c(PwConc,SurfConc)) %>% 
  cbind(depth = c(-5,-4,-3,-2,-1,-0.5,0,.5,1,1.5,2,2.5,3,3.5,4))


Treatment_names <- c("Sulfate treatment", "Control")

DGT_Fe_data <- DGT_Fe_table %>% 
  pivot_longer(PwConc_O4:SurfConc_S4, names_to = c("Concentration","core"), values_to = "value", values_drop_na = T,names_sep = "_") %>% 
  mutate(treatment = substr(core,1,1),
         time = case_when( substr(core,2,2) == "1" ~ 4,
                           substr(core,2,2) == "2" ~ 47,
                           substr(core,2,2) %in% c("3","4") ~ 81),
          depth = case_when(core == "S3" ~ depth + 0.5,
                            core == "S4" ~ depth + 0.5,
                            core == "O3" ~ depth + 0.5,
                            core == "O4"~ depth + 0.5)) %>% 
  transform(treatment = factor(treatment, levels = c("S","O"), labels = Treatment_names))

DGT_Fe_error <-  DGT_Fe_data %>% group_by(treatment,depth) %>% reframe(diff = diff(value)) %>% 
group_by(treatment) %>% summarize( meandiff = mean(diff),n = length(diff), RMS = sqrt(sum(diff^2)/(2*n)))

DGT_Fe_mean <-  DGT_Fe_data %>% group_by(treatment,depth,Concentration) %>% reframe(conc = mean(value)) 


Treatment_Colors <- c("#E69F00","#56B4E9")
Fe_profiles <- ggplot(  DGT_Fe_mean, 
                     mapping = aes(
                       y = conc,
                       x = depth,
                       color = treatment)) +
  geom_line()  +
  geom_point(size = 2) +
  xlab("Depth in cm)") +
  scale_color_manual(values = Treatment_Colors)+
  scale_x_reverse() +
  coord_flip() +
  labs(  y = expression(paste("concentration in nmol/cm2")), 
         x = "Depth in cm",
         title = " DGT Fe profile" ) +
  facet_grid(.~ Concentration, scales = "free")+
  theme_bw() +
  theme(      aspect.ratio = 2,
              axis.title=element_text(size=10),
              plot.title = element_text(size=20, face="bold", hjust = 0.5) ,
              axis.text=element_text(size=10,angle = 0, hjust = 0.5),
              legend.title = element_text(size = 0),
              legend.position = "bottom",
              legend.text = element_text(size = 10, face = "bold"),
              legend.key.size = unit(1.5, "cm"))
print(Fe_profiles)

calc <- function(cap,mag,mm){
  cap*10^(mag)/(mm*3.8)
  }

calc(3000,0,mmP)


histogram <- ggplot(  DGT_P_data, 
                     mapping = aes(
                       x = value,
                       fill = treatment, by = core)) +
geom_histogram()+
  facet_grid(.~ time, scales = "free")

show(histogram)

M <- 300*0.8 # in nmol
D <- 5.27*10^-6 # in cm2 s-1
t <- 47*(60*60*24)# in sec
A <- 1.8*0.5 # in cm2
Vgel <- A*0.1
co_Mvst <-  1*(D*A*(60*60*24))/(0.094)
co_Cvst <- (M*0.094)/(D*A*(60*60*24))
co_CvsM <- (0.094)/(D*1.8*0.5*48)

M/(A*t)
7/(A*4)
50*(D*A*(60*60*24)*4)/(0.094)

(7*0.094)/(D*1.8*0.5*4*(60*60*24))
(M*0.094)/(D*A*t)

(D*1.8*0.5*t*0.05)/(0.094)
(M*0.094)/(D*1.8*0.5*0.05)

curve(co_Mvst*x, from = 0, to = 5, col = "blue", lwd = 2, main = "Plot of y = x^2")
#read_xpt(  "input/microplates/20231130_P_DGT_1.xpt")

DGT_1 %>% subset( c(2,4,6,8,10))

