#Clear workspace################################################################
rm(list=ls())

#Load packages####

library(tidyverse)
library(dplyr) #for data wrangling
library(ggplot2) #for maming graphs
library(vegan) #for community analysis
library(reshape2) #for reshaping datasets
library(grid) #for creating a grid plot
library(gridExtra) 
library(iNEXT)

#ste wd####
#load cleaned species data####
df_moth<- read.delim("Data/moth_rar.txt")

#Carry out rarefaction
df_moth_rar <- iNEXT(df_moth, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95,
                        nboot=50) #run rarefaction analysis, with q=0 as species richness, and abundance


#create a graph to plot the moth rarefaction curve
p1<- ggiNEXT(df_moth_rar, type=2, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE) +
  theme_classic() +
  scale_colour_manual(values=c("#CC0FF9", "#4DFD5E", "#20A430")) + #change the colour of the line
  scale_fill_manual(values=c("#CC0FF9", "#4DFD5E", "#20A430")) + #change the colour of the bands
  theme(legend.position = "right") + #no legend
  theme(axis.text.x=element_text(size=12), #changes size of the element text
        axis.text.y=element_text(size=12), #changes size of x axis label text
        axis.title=element_text(size=12)) #changes size of the y axis label text
  # labs(x="Number of Individuals", y=" Number of Moth Species")

p1
