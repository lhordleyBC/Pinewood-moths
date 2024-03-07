#NMDS and Adonis test for Mar Lodge moths
#Patrick Cook
#Last modified 07/03/2024


#Intro####

#Clear workspace
rm(list=ls())

#set working directory
setwd("C:/Users/pcook/OneDrive - Butterfly Conservation/Homeworking/BC/Marvellous Moths of Mar Lodge/Stats")

#Load packages
library(tidyverse)
library(dplyr)
library(vegan) #for running NMDS
library(reshape2) #for reshaping dataframes
library(ggplot2) #fore creating graphs
library(gridExtra) #for creating grids
library(ggforce) #for drawing ellipses
library(ggrepel) #for moving labels in plots
library(grid)


#read in layers
moth <- read.csv("data.csv") #moth data
plots <- read.csv("plots.csv") #plot information

#Run an NMDS ordination#########################################################
#Take long version of moth dataset and re-arrange to be a matrix based on abundance
moth_wide <- dcast(moth, plot + treatment ~ scientific_name, value.var="abundance", fun.aggregate = sum)
sort(rowSums(moth_wide[,3:ncol(moth_wide)]>0)) #sort columns so they are not duplicated for each plot

#run the NMDS
NMDS_moth<-metaMDS(moth_wide[,3:ncol(moth_wide)],
                   distance = "bray", k=3)# when k=2 stress value was above 0.2 so have raised to k=3

str(NMDS_moth)# check the structure of the NMDS output

NMDS_moth$stress #print the stress score with k=3
#Stress score= 0.19 which is below 0.2 so ok to proceed

stressplot(object=NMDS_moth)#create a stress plot of observed dissimilarity and ordination distance


#create a data frame of site NMDS scores
data.scores <- as.data.frame(vegan:::scores.metaMDS(NMDS_moth, "sites")) #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  #create a column of site names, from the rownames of data.scores
data.scores$treatment <- moth_wide$treatment  #add the habitat zone information
head(data.scores)  #look at the data

#create a data frame of species NMDS scores
species.scores <- as.data.frame(vegan:::scores.metaMDS(NMDS_moth, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

#create a two panel plot of site and species scores from the NMDS
p1<-ggplot() +
  geom_point(data=data.scores, mapping= aes(x=NMDS1, y=NMDS2, colour= treatment)) + #create points of site locations
  #geom_segment(data=species.scores, mapping= aes(x=0, y=0, xend= NMDS1, yend= NMDS2),
  #arrow = arrow(length=unit(0.01, "npc")),
  #colour= "black") + #create arrow segments for species scores stating at the origin of the graph
  #geom_text_repel(data=pollinator.species.scores, mapping= aes(label=species, x=NMDS1, y=NMDS2))+ #create labels for species and disperse labels to minimise overlap
  geom_mark_ellipse(data=data.scores, mapping= aes(x=NMDS1, y=NMDS2, colour= treatment)) + #create an ellipse around sites
  scale_colour_manual(values=c("moorland" = "#CC0FF9", "regeneration" = "#4DFD5E", "woodland"= "#20A430")) + #set colours for each treatment
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
        legend.position = "none") 

p2<-ggplot() +
  #geom_point(data=pollinator.data.scores, mapping= aes(x=NMDS1, y=NMDS2, colour= treatment)) + #create points of site locations
  geom_segment(data=species.scores, mapping= aes(x=0, y=0, xend= NMDS1, yend= NMDS2),
               arrow = arrow(length=unit(0.01, "npc")),
               colour= "black") + #create arrow segments for species scores stating at the origin of the graph
  geom_text_repel(data=species.scores, mapping= aes(label=species, x=NMDS1, y=NMDS2))+ #create labels for species and disperse labels to minimise overlap
  #geom_mark_ellipse(data=data.scores, mapping= aes(x=NMDS1, y=NMDS2, colour= treatment)) + 
  scale_colour_manual(values=c("moorland" = "#CC0FF9", "regeneration" = "#4DFD5E", "woodland"= "#20A430")) + #set colours for each treatment
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
        legend.position = "none") 

grid.arrange(p1, p2, ncol=2, nrow= 1)


#Run adonis function to determine differences between habitat zones#############

#Comapre moorland and regeneration zone
matrix1 <- moth_wide%>%
  filter(treatment=="moorland"|treatment=="regeneration")

matrix2<- matrix1[,3:ncol(matrix1)]
ad_test1<- adonis2(matrix2~treatment, data= matrix1, permutations=10000)#perform an adonis test
ad_test1
ad_test1$`Pr(>F)`[1]

#p<0.001, moorland and regen are different

#Comapre moorland and woodland zone
matrix3 <- moth_wide%>%
  filter(treatment=="moorland"|treatment=="woodland")

matrix4<- matrix3[,3:ncol(matrix3)]
ad_test2<- adonis2(matrix4~treatment, data= matrix3, permutations=10000)#perform an adonis test
ad_test2
ad_test2$`Pr(>F)`[1]

#p<0.001, moorland and woodland zones are different

#comapre woodland and regeneration zones
matrix5 <- moth_wide%>%
  filter(treatment=="woodland"|treatment=="regeneration")

matrix6<- matrix5[,3:ncol(matrix1)]
ad_test3<- adonis2(matrix6~treatment, data= matrix5, permutations=10000)#perform an adonis test
ad_test3
ad_test3$`Pr(>F)`[1]

#p<0.001, regeneration and woodland zones are different
