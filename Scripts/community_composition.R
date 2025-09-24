#NMDS and Adonis test for Mar Lodge moths
#Patrick Cook and Katie Powell
#Last modified 23/01/2025

#Intro####

#Clear workspace
rm(list=ls())

#Load packages
library(tidyverse)
library(vegan) #for running NMDS
library(reshape2) #for reshaping dataframes
library(ggplot2) #fore creating graphs
library(ggforce) #for drawing ellipses
library(ggrepel) #for moving labels in plots
library(ggpubr)
library(sf)
library(rnrfa)
library(ape)
library(glmmTMB)
library(performance)
library(DHARMa)
library(lsmeans)
library(lme4)
library(MuMIn)
library(multcomp)
## Get raw data
moth <- read.csv("Data/Moth_data.csv", header=TRUE)

#Run an NMDS ordination#########################################################
#Take long version of moth dataset and re-arrange into wide format based on abundance

moth_wide <- dcast(moth, Plot + Treatment ~ Taxon, value.var="Quantity", fun.aggregate = sum)
moth_wide <- moth_wide %>%
  mutate(Treatment = recode(Treatment, "moorland" = 'Moorland', 
                            regeneration = 'Early Successional Woodland', woodland = 'Mature Woodland'))

head(moth_wide)
moth_wide<- distinct(moth_wide)
moth_wide[is.na(moth_wide)]<- 0

#run the NMDS
set.seed(1)
NMDS_moth<-metaMDS(moth_wide[,3:ncol(moth_wide)],
                   autotransform = TRUE,
                   distance = "bray", k=3)# when k=2 stress value was above 0.2 so have raised to k=3

str(NMDS_moth)# check the structure of the NMDS output

NMDS_moth$stress #print the stress score with k=3
#Stress score= 0.19 which is below 0.2 so ok to proceed

stressplot(object=NMDS_moth)#create a stress plot of observed dissimilarity and ordination distance

#create a data frame of site NMDS scores
data.scores <- as.data.frame(vegan:::scores.metaMDS(NMDS_moth, "sites")) #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  #create a column of site names, from the rownames of data.scores
data.scores$treatment <- moth_wide$Treatment  #add the habitat zone information
data.scores$treatment <- factor(data.scores$treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

head(data.scores)  #look at the data
moth_wide$Treatment
#create a data frame of species NMDS scores
species.scores <- as.data.frame(vegan:::scores.metaMDS(NMDS_moth, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

#create a two panel plot of site and species scores from the NMDS
p1<-ggplot() +
  geom_point(data=data.scores, mapping= aes(x=NMDS1, y=NMDS2, colour= treatment)) + #create points of site locations
  geom_mark_ellipse(data=data.scores, mapping= aes(x=NMDS1, y=NMDS2, colour= treatment)) + #create an ellipse around sites
  scale_colour_manual(values=c("Moorland" = "#CC0FF9", "Early Successional Woodland" = "#4DFD5E", "Mature Woodland"= "#20A430")) + #set colours for each treatment
  geom_abline(intercept=0, slope=0, linetype="dashed", colour="black")+ #create a horizontal line through 0
  geom_vline(aes(xintercept=0), linetype="dashed", colour="black")+ #create a vertical lime through 0
  labs(x = "NMDS1", y = "NMDS2") + #change x and y labels
  xlim(-1.25, 1.25) + #change x limit of graph
  ylim(-1.25,1.25) + #change y limit of graph
  theme_bw() +
  theme(panel.background = element_blank(), #blank background for the graph
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position="inside",
        legend.position.inside = c(.75,.9),
        legend.title=element_blank(),
        legend.text = element_text(size=14),
        legend.key.size = unit(0.75, 'cm'),
        axis.text=element_text(size=18), 
        axis.title = element_text(size=18)) 
p2<-ggplot() +
  geom_segment(data=species.scores, mapping= aes(x=0, y=0, xend= NMDS1, yend= NMDS2),
               arrow = arrow(length=unit(0.01, "npc")),
               colour= "black") + #create arrow segments for species scores stating at the origin of the graph
  geom_text_repel(data=species.scores, mapping= aes(label=species, x=NMDS1, y=NMDS2))+ #create labels for species and disperse labels to minimise overlap
  geom_abline(intercept=0, slope=0, linetype="dashed", colour="black")+ #create a horizontal line through 0
  geom_vline(aes(xintercept=0), linetype="dashed", colour="black")+ #create a vertical lime through 0
  annotation_custom(grobTree(textGrob("Stress Score=0.19", x=0.72,  y=0.02, hjust=0,
                                      gp=gpar(col="black", fontsize=12))))+
  labs(x = "NMDS1", y = "NMDS2") + #change x and y labels
  xlim(-1.25, 1.25) + #change x limit of graph
  ylim(-1.25,1.25) + #change y limit of graph
  theme_bw() +
  theme(panel.background = element_blank(), #blank background for the graph
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=18), 
        axis.title = element_text(size=18)) 

p3 <- ggarrange(p1, p2, labels=c("(a)", "(b)"), font.label=list(size=18))

ggsave(p3, file="Graphs/NMDS.png", height=8, width=15)



################## Run adonis function to determine differences between habitat zones #############

### First create distance matrix, using Bray Curtis dissimilarities
head(moth_wide)

#######Compare moorland and regeneration zone
temp<- moth_wide%>%
  filter(Treatment=="Moorland"|Treatment=="Early Successional Woodland")
dist_mat<- as.matrix(vegdist(temp[3:ncol(temp)], method="bray"))
ad_test1<- adonis2(dist_mat~Treatment, data= temp, permutations=10000, method="bray")#perform an adonis test

ad_test1
ad_test1$`Pr(>F)`[1]

#p<0.001, moorland and regen are different

temp

#######Compare moorland and woodland zone
temp<- moth_wide%>%
  filter(Treatment=="Moorland"|Treatment=="Mature Woodland")
dist_mat<- as.matrix(vegdist(temp[3:ncol(temp)], method="bray"))
ad_test2<- adonis2(dist_mat~Treatment, data= temp, permutations=10000, method="bray")#perform an adonis test

ad_test2
ad_test2$`Pr(>F)`[1]

#p<0.001, moorland and woodland zones are different




################# compare woodland and regeneration zones
temp<- moth_wide%>%
  filter(Treatment=="Early Successional Woodland"|Treatment=="Mature Woodland")
head(temp)
dist_mat<- as.matrix(vegdist(temp[3:ncol(temp)], method="bray"))
ad_test3<- adonis2(dist_mat~Treatment, data= temp, permutations=10000, method="bray")#perform an adonis test
ad_test3
ad_test3$`Pr(>F)`[1]

#p<0.001, regeneration and woodland zones are different
ad_test1
ad_test2
ad_test3



################## Accounting for spatial autocorrelation?

#################### CREATE DISTANCE MATRICES ##############

##### Calculate geographic distances between pairs of sites
head(moth)
## convert gridref to easting and northing and add to moth dataframe
coords<- osg_parse(moth$Gridref, coord_system="BNG") 
moth$easting<- coords$easting
moth$northing<-coords$northing
head(moth)
plots<-moth%>% distinct(Plot, Gridref, easting, northing)%>%arrange(Plot)
## convert moth to spatial object
plots_sf<- st_as_sf(plots, coords=c("easting", "northing"), crs=27700)
head(plots_sf)
## calculate distances between plots and turn into a matrix
geo_dist<- st_distance(plots_sf)/1000 ## (convert to km by /1000)
geo_dist<- matrix(geo_dist, nrow=nrow(plots), ncol=nrow(plots))
row.names(geo_dist)<- plots$Plot
colnames(geo_dist)<- plots$Plot  


##### Calculate compositional differences between each pair of sites
comp_dist<- as.matrix(vegdist(moth_wide[3:ncol(moth_wide)], method="bray"))
comp_dist

##### Code an artificial matrix for differences in treatment
unique(moth_wide$Treatment)
## Moorland = 1, Early Successional Woodland = 2, Mature Woodland = 3
moth_wide<- moth_wide %>% 
  mutate(Treatment_code=recode(Treatment, 
                               "Moorland" = 1, 
                               "Early Successional Woodland" = 2, 
                               "Mature Woodland" = 3 ))%>%
  select(Plot, Treatment, Treatment_code, 3:ncol(.))
head(moth_wide)
treat_dist<- as.matrix(vegdist(moth_wide[3], method="euclidean"))
treat_dist


##############  calculate Moran's I value - looking for SAC between each NMDS axis and geographic distance

## Calculate inverse distance weights from the geo_dist matrix
geo_dist_inv<- 1/(geo_dist) ## This turns distance into positive weights - i.e., higher values = closer together
geo_dist[1:6,1:6]
geo_dist_inv[1:6,1:6]
diag(geo_dist_inv)<- 0 ## replace Infinity values!
head(data.scores$NMDS1)
head(data.scores)
geo_dist

## Calculate Moran's I values for:
## 1. NMDS1
Moran.I(data.scores$NMDS1, geo_dist_inv, scaled=TRUE, alternative="two.sided") ## cannot reject null hypothesis of spatial-autocorrelation along NMDS1
## Moran's I = 0.48 (expected = -0.02); sd=0.07, p<0.001

## 2. NMDS2
Moran.I(data.scores$NMDS2, geo_dist_inv, scaled=TRUE, alternative="two.sided") ## Cannot reject null hypothesis of spatial-autocorrelation along NMDS2
## Moran's I = 0.57 (expected = -0.02); sd=0.07, p<0.001
##

################# Mantel tests: Correlation between matrices ##################################

##Compare geo_dist matrix with a matrix of compositional dissimilarities
set.seed(1)
mantel(comp_dist, geo_dist) ##mantel R= 0.256, p<0.001
#Significant positive correlation between geographic distance and compositional differences

## Test correlation between compositional differences and treatment differences
mantel(comp_dist, treat_dist) ## Mantel R = 0.349, p<0.001 # Significant and positive - and stronger than geo


## Test correlation between compositional differences and treatment differences
mantel(treat_dist, geo_dist) ## Mantel R = 0.412, p<0.001 # Significant and positive - and strong


## Do a partial mantel to test correlation between composition and treatment, whilst accounting for geographic distance
mantel.partial(comp_dist, treat_dist, geo_dist)  ## Mantel R= 0.276, - so is weakened from 0.349 by the effect of geography
                                                 ## p<0.001    

############### Run some simple models

## Distance models
mod1<-lm(comp_dist~treat_dist+geo_dist, data=plot_data)
summary(mod1)
## Significant effect of both treatment distance and geographic distance on community dissimilarity, with a greater effect of treatment 


## Modelling NMDS axes as a function of treatment
hist(data.scores$NMDS1)
hist(data.scores$NMDS2)

mod2<- glmmTMB(NMDS1~treatment + (1|site), family="gaussian", data=data.scores)
summary(mod2)  

testDispersion(mod2) # looks good
simulationOutput <- simulateResiduals(fittedModel = mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(mod2) # all <3
r.squaredGLMM(mod2) # 35%

# Pairwise differences in management effects
glht1 <- glht(mod2, mcp(treatment="Tukey"))
summary(glht1) 
# no significant difference between regen and moorland
# significant difference between woodland and moorland
# significant difference between woodland and regen

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Broadleaf richness"







############### xy plots of distance values
head(comp_dist)
plot_data<- geo_dist %>% 
  as.data.frame() %>%
  mutate(plot1=rownames(.))%>% 
  gather("plot2", "geo_dist", 1:(ncol(.)-1))
plot_data<- comp_dist %>% as.data.frame() %>%
  mutate(plot1=rownames(.))%>% gather("plot2", "comp_dist", 1:(ncol(.)-1))%>%
  left_join(plot_data, by=c("plot1", "plot2"))
plot_data<- treat_dist %>% as.data.frame() %>%
  mutate(plot1=rownames(.))%>% gather("plot2", "treat_dist", 1:(ncol(.)-1))%>%
  left_join(plot_data, by=c("plot1", "plot2"))%>%
  filter(treat_dist!=0&comp_dist!=0&geo_dist!=0)
head(plot_data)
ggplot(plot_data)+
  geom_point(aes(x=geo_dist, y=comp_dist, col=as.factor(treat_dist)))+
  theme_bw()+
  labs(x="Geographic Distance (km)", y="Compositional Dissimilarity (Bray-Curtis)")+
  scale_colour_discrete(name="Step difference in treatment")

