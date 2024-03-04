##########################
## User: Lisbeth Hordley
## Date: November 2023
## Info: Pinewood moth analysis - structure and plant cover models

rm(list = ls())

library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(lme4)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(ggeffects)
library(multcomp)
library(broom)
library(ggpubr)
library(lsmeans)
library(lubridate)
library(performance)
library(vegan)
library(factoextra)
options(scipen=999)

# read in data
moths_final <- read.csv("Data/Moths_data_final.csv", header=TRUE)
habitat <- read.csv("Data/Habitat_data.csv", header=TRUE)

colSums(is.na(habitat)) # for some reason there are 5 rows of all NAs so remove these
habitat <- na.omit(habitat)
str(habitat)
habitat <- habitat %>% 
  mutate_at(c(30:35), as.numeric)

## Structural explanatory variables:
# Conifer cover (conifer_canopy)
# Broadleaf cover (broadleaf_canopy)
# Regeneration cover (tree_regeneration)
# Ground layer cover (ground_layer)
# Complexity score (complexity_score)
# Coefficient of structural variance (coefficient of variance of ground layer height columns)

# first calculate the mean and SD across the ground layer height columns before calculating SD and CV
habitat$ground_layer_height_sd <- apply(habitat[, c(30:35)],1,sd)
habitat$ground_layer_height_mean <- rowMeans(habitat[,c(30:35)])
habitat$ground_layer_height_cv <- habitat$ground_layer_height_sd / habitat$ground_layer_height_mean
# CV = variability of ground layer height at each plot 

## Plant cover explanatory variables:
# Calluna vulgaris
# Erica cinerea
# Vaccinium myrtilus
# Vaccinium vitis idaea
# Grass
# Forb 
# Pinus sylvestris
# Betula spp

# put both datasets together
moths_final <- merge(moths_final, habitat, by.x="Plot", by.y="plot")

## First check for correlations between covariates and structural features
round(cor(moths_final[,c("Min_temp","Wind_speed", "conifer_canopy", "broadleaf_canopy", 
                         "tree_regeneration", "ground_layer", "complexity_score", "ground_layer_height_cv")]),3) 
# all correlations are <0.7 - keep all variables in models

################### STRUCTURAL FEATURE AND MOTH MODELS ################### 

## Abundance 

ggplot(moths_final, aes(conifer_canopy, richness))+
  geom_point()+ 
  geom_smooth()+
  theme_classic()

ggplot(moths_final, aes(broadleaf_canopy, richness))+
  geom_point()+ 
  # geom_smooth()+
  theme_classic()

ggplot(moths_final, aes(ground_layer, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(complexity_score, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(ground_layer_height_cv, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()


##########################################################################################

abund_struc <- glmmTMB(tot_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + scale(tree_regeneration) +
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(ground_layer_height_cv) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]

AIC(abund_struc)

# minimum temperature (positive)
# wind speed (negative)
# conifer canopy (positive)
# broadleaf canopy (positive)
# ground layer (positive)
# complexity score (positive)

## check model assumptions
testDispersion(abund_struc) ## some underdispersion.. not when using genpois
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) ## quantile deviations
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## complexity score needs removed from this model

abund_broadleaf <- ggpredict(abund_struc, terms="broadleaf_canopy") 
abund_conifer <- ggpredict(abund_struc, terms="conifer_canopy") 
abund_ground <- ggpredict(abund_struc, terms="ground_layer")
abund_complex <- ggpredict(abund_struc, terms="complexity_score") 

abund_broadleaf$predictor <- "Broadleaf canopy cover"
abund_conifer$predictor <- "Conifer canopy cover"
abund_ground$predictor <- "Ground layer cover"
abund_complex$predictor <- "Vertical complexity score"

abund_broadleaf$response <- "Total abundance"
abund_conifer$response <- "Total abundance"
abund_ground$response <- "Total abundance"
abund_complex$response <- "Total abundance"

abund_structure <- rbind(abund_broadleaf, abund_conifer, abund_ground, abund_complex)

## Woodland abundance
abund_struc <- glmmTMB(woodland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + 
                         scale(tree_regeneration) + (1|visit) + (1|Sub_site), data=moths_final, 
                         family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]

AIC(abund_struc)

# wind speed (negative)
# conifer canopy (positive)
# broadleaf canopy (positive)
# ground layer (positive)

## check model assumptions
testDispersion(abund_struc) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) 
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## all <3

wood_abund_broadleaf <- ggpredict(abund_struc, terms="broadleaf_canopy") 
wood_abund_conifer <- ggpredict(abund_struc, terms="conifer_canopy") 
wood_abund_ground <- ggpredict(abund_struc, terms="ground_layer")

wood_abund_broadleaf$predictor <- "Broadleaf canopy cover"
wood_abund_conifer$predictor <- "Conifer canopy cover"
wood_abund_ground$predictor <- "Ground layer cover"

wood_abund_broadleaf$response <- "Woodland total abundance"
wood_abund_conifer$response <- "Woodland total abundance"
wood_abund_ground$response <- "Woodland total abundance"

abund_structure <- rbind(abund_structure, wood_abund_broadleaf, wood_abund_conifer, wood_abund_ground)


## Moorland abundance
abund_struc <- glmmTMB(moorland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(tree_regeneration) + scale(ground_layer_height_cv) + (1|visit) + (1|Sub_site), data=moths_final, 
                         family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]


AIC(abund_struc)

# wind speed (negative)
# conifer canopy (positive)
# broadleaf canopy (positive)
# ground layer (positive)
# complexity score (positive)

## check model assumptions
testDispersion(abund_struc) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) 
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## complexity score needs removed

moor_abund_broadleaf <- ggpredict(abund_struc, terms="broadleaf_canopy") 
moor_abund_conifer <- ggpredict(abund_struc, terms="conifer_canopy") 
moor_abund_ground <- ggpredict(abund_struc, terms="ground_layer")
moor_abund_complex <- ggpredict(abund_struc, terms="complexity_score")

moor_abund_broadleaf$predictor <- "Broadleaf canopy cover"
moor_abund_conifer$predictor <- "Conifer canopy cover"
moor_abund_ground$predictor <- "Ground layer cover"
moor_abund_complex$predictor <- "Vertical complexity score"

moor_abund_broadleaf$response <- "Moorland total abundance"
moor_abund_conifer$response <- "Moorland total abundance"
moor_abund_ground$response <- "Moorland total abundance"
moor_abund_complex$response <- "Moorland total abundance"

abund_structure <- rbind(abund_structure, moor_abund_broadleaf, moor_abund_conifer, 
                         moor_abund_ground, moor_abund_complex)


## Grassland abundance
abund_struc <- glmmTMB(grassland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) +
                         scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]


AIC(abund_struc)

# minimum temperature (positive)
# wind speed (negative)
# ground layer (positive)

## check model assumptions
testDispersion(abund_struc) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) 
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## all low correlations

ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Conifer abundance
abund_struc <- glmmTMB(conifer_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]


AIC(abund_struc)

# minimum temperature (positive)
# wind speed (negative)
# complexity score (positive)

## check model assumptions
testDispersion(abund_struc) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) 
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## all low correlations

ggpredict(abund_struc, terms="complexity_score") %>% plot(add.data=TRUE)


## Broadleaf abundance
abund_struc <- glmmTMB(broadleaf_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]


AIC(abund_struc)

# conifer canopy (positive)
# broadleaf canopy (negative)
# ground layer (positive)

## check model assumptions
testDispersion(abund_struc) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) 
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## all low correlations

ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Shrub abundance
abund_struc <- glmmTMB(shrub_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) +
                         scale(tree_regeneration) + scale(ground_layer_height_cv) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

summary(abund_struc)

abund_struc_sum <- as.data.frame(summary(abund_struc)$coefficients$cond)
abund_struc_sum$parameters <- row.names(abund_struc_sum)
row.names(abund_struc_sum) <- 1:nrow(abund_struc_sum)
abund_struc_sum$`z value` <- NULL
colnames(abund_struc_sum)[3] <- "p_value"
abund_struc_sum <- as.data.frame(abund_struc_sum)
abund_struc_sum$significance <- case_when(
  abund_struc_sum$p_value>=0.05 ~ "ns",
  abund_struc_sum$p_value<0.05 & abund_struc_sum$p_value>=0.01 ~ "*",
  abund_struc_sum$p_value<0.01 & abund_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_struc_sum <- abund_struc_sum[,c(4,1,2,3,5)]

AIC(abund_struc)

# wind speed (negative)
# conifer canopy (positive)
# broadleaf canopy (negative)
# ground layer (positive)

## check model assumptions
testDispersion(abund_struc) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_struc, plot = F)
plot(simulationOutput) 
# no assumptions violated using genpois

## check for multicolinearity
check_collinearity(abund_struc) ## complexity score needs removed

ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)



###########################################################################################################


## Richness

rich_struc <- glmer(richness ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# minimum temperature (positive)
# wind speed (negative)
# conifer canopy (positive)
# ground layer (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## quantile deviations
testZeroInflation(simulationOutput)
## check for multicolinearity
check_collinearity(rich_struc) ## all low correlations

rich_conifer <- ggpredict(rich_struc, terms="conifer_canopy")
rich_ground <- ggpredict(rich_struc, terms="ground_layer")

rich_conifer$predictor <- "Conifer canopy cover"
rich_ground$predictor <- "Ground layer cover"

rich_conifer$response <- "Total species richness"
rich_ground$response <- "Total species richness"

rich_structure <- rbind(rich_conifer, rich_ground)


## Woodland richness

rich_struc <- glmer(woodland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# minimum temperature (positive)
# wind speed (negative)
# conifer canopy (positive)
# broadleaf canopy (positive)
# ground layer (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## quantile deviations
testZeroInflation(simulationOutput)

## check for multicolinearity
check_collinearity(rich_struc) ## all low correlations

wood_rich_broadleaf <- ggpredict(rich_struc, terms="broadleaf_canopy")
wood_rich_conifer <- ggpredict(rich_struc, terms="conifer_canopy")
wood_rich_ground <- ggpredict(rich_struc, terms="ground_layer")

wood_rich_broadleaf$predictor <- "Broadleaf canopy cover"
wood_rich_conifer$predictor <- "Conifer canopy cover"
wood_rich_ground$predictor <- "Ground layer cover"

wood_rich_broadleaf$response <- "Woodland species richness"
wood_rich_conifer$response <- "Woodland species richness"
wood_rich_ground$response <- "Woodland species richness"

rich_structure <- rbind(rich_structure, wood_rich_broadleaf, wood_rich_conifer, wood_rich_ground)

## Moorland richness
rich_struc <- glmer(moorland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) +
                      scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# minimum temperature (positive)
# wind speed (negative)
# ground layer (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]

## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## quantile deviations
testZeroInflation(simulationOutput)

## check for multicolinearity
check_collinearity(rich_struc) ## all low correlations

moor_rich_ground <- ggpredict(rich_struc, terms="ground_layer")
moor_rich_ground$predictor <- "Ground layer cover"
moor_rich_ground$response <- "Moorland species richness"

rich_structure <- rbind(rich_structure, moor_rich_ground)

## Grassland richness
rich_struc <- glmer(grassland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# minimum temperature (positive)
# wind speed (negative)
# ground layer (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## quantile deviations
testZeroInflation(simulationOutput)

## check for multicolinearity
check_collinearity(rich_struc) ## all low correlations

ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Conifer richness
rich_struc <- glmer(conifer_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(tree_regeneration) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# minimum temperature (positive)
# wind speed (negative)
# complexity score (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## quantile deviations
testZeroInflation(simulationOutput)

## check for multicolinearity
check_collinearity(rich_struc) ## complexity score needs removed

ggpredict(rich_struc, terms="complexity_score") %>% plot(add.data=TRUE)

## Broadleaf richness
rich_struc <- glmer(broadleaf_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# minimum temperature (positive)
# conifer cover (positive)
# broadleaf cover (positive)
# ground layer (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## quantile deviations
testZeroInflation(simulationOutput)

## check for multicolinearity
check_collinearity(rich_struc) ## all low correlations

ggpredict(rich_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Shrub richness
rich_struc <- glmer(shrub_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(tree_regeneration) + scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_struc)
# wind speed (negative)
# lunar cycle (waxing crescent) (negative)
# conifer cover (positive)
# ground layer (positive)

rich_struc_sum <- as.data.frame(coef(summary(rich_struc)))
rich_struc_sum$parameters <- row.names(rich_struc_sum)
row.names(rich_struc_sum) <- 1:nrow(rich_struc_sum)
rich_struc_sum$`z value` <- NULL
colnames(rich_struc_sum)[3] <- "p_value"
rich_struc_sum <- as.data.frame(rich_struc_sum)
rich_struc_sum$significance <- case_when(
  rich_struc_sum$p_value>=0.05 ~ "ns",
  rich_struc_sum$p_value<0.05 & rich_struc_sum$p_value>=0.01 ~ "*",
  rich_struc_sum$p_value<0.01 & rich_struc_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_struc_sum <- rich_struc_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_struc) 
simulationOutput <- simulateResiduals(fittedModel = rich_struc, plot = F)
plot(simulationOutput) ## no assumptions violated
testZeroInflation(simulationOutput)

## check for multicolinearity
check_collinearity(rich_struc) ## all low correlations

ggpredict(rich_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)



############# Supplementary material plot

structure_predictions <- rbind(abund_structure, rich_structure)
lookup <- c("Total abundance"="Total \nabundance", "Total species richness"="Total \nspecies \nrichness", "Woodland total abundance"="Woodland \ntotal \nabundance",
            "Woodland species richness"="Woodland \nspecies \nrichness", "Moorland total abundance"="Moorland \ntotal \nabundance",
            "Moorland species richness"="Moorland \nspecies \nrichness")
structure_predictions$response <- as.character(lookup[structure_predictions$response])
structure_predictions$response  = factor(structure_predictions$response, levels=c("Total \nabundance", 
                                                                                  "Total \nspecies \nrichness", "Woodland \ntotal \nabundance",
                                                                                  "Woodland \nspecies \nrichness", "Moorland \ntotal \nabundance",
                                                                                  "Moorland \nspecies \nrichness"))


abund_broad_raw <- moths_final[,c("tot_abund", "broadleaf_canopy")]
abund_broad_raw$response <- "Total \nabundance"
abund_broad_raw$predictor <- "Broadleaf canopy cover"
colnames(abund_broad_raw)[1:2] <- c("y", "x")

abund_conifer_raw <- moths_final[,c("tot_abund", "conifer_canopy")]
abund_conifer_raw$response <- "Total \nabundance"
abund_conifer_raw$predictor <- "Conifer canopy cover"
colnames(abund_conifer_raw)[1:2] <- c("y", "x")

abund_ground_raw <- moths_final[,c("tot_abund", "ground_layer")]
abund_ground_raw$response <- "Total \nabundance"
abund_ground_raw$predictor <- "Ground layer cover"
colnames(abund_ground_raw)[1:2] <- c("y", "x")

abund_complex_raw <- moths_final[,c("tot_abund", "complexity_score")]
abund_complex_raw$response <- "Total \nabundance"
abund_complex_raw$predictor <- "Vertical complexity score"
colnames(abund_complex_raw)[1:2] <- c("y", "x")


wood_abund_broad_raw <- moths_final[,c("woodland_abund", "broadleaf_canopy")]
wood_abund_broad_raw$response <- "Woodland \ntotal \nabundance"
wood_abund_broad_raw$predictor <- "Broadleaf canopy cover"
colnames(wood_abund_broad_raw)[1:2] <- c("y", "x")

wood_abund_conifer_raw <- moths_final[,c("woodland_abund", "conifer_canopy")]
wood_abund_conifer_raw$response <- "Woodland \ntotal \nabundance"
wood_abund_conifer_raw$predictor <- "Conifer canopy cover"
colnames(wood_abund_conifer_raw)[1:2] <- c("y", "x")

wood_abund_ground_raw <- moths_final[,c("woodland_abund", "ground_layer")]
wood_abund_ground_raw$response <- "Woodland \ntotal \nabundance"
wood_abund_ground_raw$predictor <- "Ground layer cover"
colnames(wood_abund_ground_raw)[1:2] <- c("y", "x")


moor_abund_broad_raw <- moths_final[,c("moorland_abund", "broadleaf_canopy")]
moor_abund_broad_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_broad_raw$predictor <- "Broadleaf canopy cover"
colnames(moor_abund_broad_raw)[1:2] <- c("y", "x")

moor_abund_conifer_raw <- moths_final[,c("moorland_abund", "conifer_canopy")]
moor_abund_conifer_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_conifer_raw$predictor <- "Conifer canopy cover"
colnames(moor_abund_conifer_raw)[1:2] <- c("y", "x")

moor_abund_ground_raw <- moths_final[,c("moorland_abund", "ground_layer")]
moor_abund_ground_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_ground_raw$predictor <- "Ground layer cover"
colnames(moor_abund_ground_raw)[1:2] <- c("y", "x")

moor_abund_complex_raw <- moths_final[,c("moorland_abund", "complexity_score")]
moor_abund_complex_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_complex_raw$predictor <- "Vertical complexity score"
colnames(moor_abund_complex_raw)[1:2] <- c("y", "x")

rich_conifer_raw <- moths_final[,c("richness", "conifer_canopy")]
rich_conifer_raw$response <- "Total \nspecies \nrichness"
rich_conifer_raw$predictor <- "Conifer canopy cover"
colnames(rich_conifer_raw)[1:2] <- c("y", "x")

rich_ground_raw <- moths_final[,c("richness", "ground_layer")]
rich_ground_raw$response <- "Total \nspecies \nrichness"
rich_ground_raw$predictor <- "Ground layer cover"
colnames(rich_ground_raw)[1:2] <- c("y", "x")


wood_rich_broad_raw <- moths_final[,c("woodland_rich", "broadleaf_canopy")]
wood_rich_broad_raw$response <- "Woodland \nspecies \nrichness"
wood_rich_broad_raw$predictor <- "Broadleaf canopy cover"
colnames(wood_rich_broad_raw)[1:2] <- c("y", "x")

wood_rich_conifer_raw <- moths_final[,c("woodland_rich", "conifer_canopy")]
wood_rich_conifer_raw$response <- "Woodland \nspecies \nrichness"
wood_rich_conifer_raw$predictor <- "Conifer canopy cover"
colnames(wood_rich_conifer_raw)[1:2] <- c("y", "x")

wood_rich_ground_raw <- moths_final[,c("woodland_rich", "ground_layer")]
wood_rich_ground_raw$response <- "Woodland \nspecies \nrichness"
wood_rich_ground_raw$predictor <- "Ground layer cover"
colnames(wood_rich_ground_raw)[1:2] <- c("y", "x")

moor_rich_ground_raw <- moths_final[,c("moorland_rich", "ground_layer")]
moor_rich_ground_raw$response <- "Moorland \nspecies \nrichness"
moor_rich_ground_raw$predictor <- "Ground layer cover"
colnames(moor_rich_ground_raw)[1:2] <- c("y", "x")

raw_dat <- rbind(abund_broad_raw, abund_conifer_raw, abund_ground_raw, abund_complex_raw,
                 wood_abund_broad_raw, wood_abund_conifer_raw, wood_abund_ground_raw,
                 moor_abund_broad_raw, moor_abund_conifer_raw, moor_abund_ground_raw, moor_abund_complex_raw,
                 rich_conifer_raw, rich_ground_raw, wood_rich_broad_raw, wood_rich_conifer_raw,
                 wood_rich_ground_raw, moor_rich_ground_raw)
raw_dat$response  = factor(raw_dat$response, levels=c("Total \nabundance", 
                                                      "Total \nspecies \nrichness", "Woodland \ntotal \nabundance",
                                                      "Woodland \nspecies \nrichness", "Moorland \ntotal \nabundance",
                                                      "Moorland \nspecies \nrichness"))

plot <- ggplot()+
  geom_point(data=raw_dat, aes(x=x, y=y), colour="lightgrey")+
  geom_line(data=structure_predictions, aes(x=x, y=predicted))+
  geom_ribbon(data=structure_predictions, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2)+
  facet_grid(response ~ predictor, scales="free")+
  labs(x="", y="")+
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0))

a <- c("aaaa
        #aa#
        aaa#
        aaa#
        aaa#
        ##a#")

plot <- remove_facets(plot, a)
plot
ggsave(plot, file="Graphs/Structure_supplementary_plot.png", height=5, width=7)




##########################################################################################################
##########################################################################################################


################### PLANT COVER MODELS ################### 

## Plant cover explanatory variables:
# Calluna vulgaris
# Erica cinerea
# Vaccinium myrtilus
# Vaccinium vitis idaea
# Grass
# Forb 
# Pinus sylvestris
# Betula spp

## correlations between covariates and plant cover variables
round(cor(moths_final[,c("Min_temp","Wind_speed", "calluna_vulgaris", "erica_cinerea",
                         "vaccinium_myrtilus", "vaccinium_vitis.idaea",
                         "grass", "forb", "pinus_sylvestris", "betula_spp")]),3) 
# vaccinium_myrtilus and vaccinium_vitis.idaea r=0.61
cor.test(moths_final$vaccinium_vitis.idaea, moths_final$vaccinium_myrtilus) # r=0.61, p<0.001

cor.test(moths_final$richness, moths_final$vaccinium_vitis.idaea) # r=0.2 p=0.038
cor.test(moths_final$richness, moths_final$vaccinium_myrtilus) # r=0.42, p<0.001
# remove vaccinium_vitis.idaea from richness models
cor.test(moths_final$tot_abund, moths_final$vaccinium_vitis.idaea) # r=0.37, p<0.001
cor.test(moths_final$tot_abund, moths_final$vaccinium_myrtilus) # r=0.51, p<0.001
# remove vaccinium_vitis.idaea from abundance models

# high VIF values in models for calluna and vaccinium - one needs to be removed
# correlation between them is -0.58
cor.test(moths_final$tot_abund, moths_final$calluna_vulgaris) # r=-0.26, p<0.001
cor.test(moths_final$tot_abund, moths_final$vaccinium_myrtilus) # r=0.51, p<0.001
# remove calluna vulgaris from abundance models - not ideal
cor.test(moths_final$richness, moths_final$calluna_vulgaris) # r=-0.19, p=0.05
cor.test(moths_final$richness, moths_final$vaccinium_myrtilus) # r=0.42, p<0.001
# remove calluna vulgaris from richness models too

ggplot(moths_final, aes(x=calluna_vulgaris, y=vaccinium_myrtilus))+
  geom_point() # negatively correlated

## Abundance

ggplot(moths_final, aes(calluna_vulgaris, tot_abund))+
  geom_point()+ 
  geom_smooth()+
  theme_classic()

ggplot(moths_final, aes(erica_cinerea, tot_abund))+
  geom_point()+ 
  geom_smooth()+
  theme_classic()

ggplot(moths_final, aes(vaccinium_myrtilus, tot_abund))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(grass, tot_abund))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(forb, tot_abund))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()
# two big outliers

ggplot(moths_final, aes(pinus_sylvestris, tot_abund))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(betula_spp, tot_abund))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()


abund_plant <- glmmTMB(tot_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")
summary(abund_plant)
# min temp (positive)
# wind speed (negative)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]

## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3

abund_betula <- ggpredict(abund_plant, terms="betula_spp") 
abund_pine <- ggpredict(abund_plant, terms="pinus_sylvestris") 
abund_vaccinium <- ggpredict(abund_plant, terms="vaccinium_myrtilus")

abund_betula$predictor <- "Betula spp"
abund_pine$predictor <- "Pinus sylvestris"
abund_vaccinium$predictor <- "Vaccinium myrtilus"

abund_betula$response <- "Total \nabundance"
abund_pine$response <- "Total \nabundance"
abund_vaccinium$response <- "Total \nabundance"

abund_plants <- rbind(abund_betula, abund_pine, abund_vaccinium)

## Woodland abundance
abund_plant <- glmmTMB(woodland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")
summary(abund_plant)
# wind speed (negative)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3

wood_abund_betula <- ggpredict(abund_plant, terms="betula_spp") 
wood_abund_pine <- ggpredict(abund_plant, terms="pinus_sylvestris") 
wood_abund_vaccinium <- ggpredict(abund_plant, terms="vaccinium_myrtilus")

wood_abund_betula$predictor <- "Betula spp"
wood_abund_pine$predictor <- "Pinus sylvestris"
wood_abund_vaccinium$predictor <- "Vaccinium myrtilus"

wood_abund_betula$response <- "Woodland \ntotal \nabundance"
wood_abund_pine$response <- "Woodland \ntotal \nabundance"
wood_abund_vaccinium$response <- "Woodland \ntotal \nabundance"

abund_plants <- rbind(abund_plants, wood_abund_betula, wood_abund_pine, wood_abund_vaccinium)


## Moorland abundance
abund_plant <- glmmTMB(moorland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")
summary(abund_plant)
# wind speed (negative)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3


moor_abund_betula <- ggpredict(abund_plant, terms="betula_spp") 
moor_abund_pine <- ggpredict(abund_plant, terms="pinus_sylvestris") 
moor_abund_vaccinium <- ggpredict(abund_plant, terms="vaccinium_myrtilus")

moor_abund_betula$predictor <- "Betula spp"
moor_abund_pine$predictor <- "Pinus sylvestris"
moor_abund_vaccinium$predictor <- "Vaccinium myrtilus"

moor_abund_betula$response <- "Moorland \ntotal \nabundance"
moor_abund_pine$response <- "Moorland \ntotal \nabundance"
moor_abund_vaccinium$response <- "Moorland \ntotal \nabundance"

abund_plants <- rbind(abund_plants, moor_abund_betula, moor_abund_pine, moor_abund_vaccinium)


## Grassland abundance
abund_plant <- glmmTMB(grassland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")
summary(abund_plant)
# minimum temperature (positive)
# wind speed (negative)
# forb (positive)

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3

ggpredict(abund_plant, terms="forb") %>% plot(add.data=TRUE) # potentially being driven by outliers

## Conifer abundance
abund_plant <- glmmTMB(conifer_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, ziformula=~1, family="genpois", na.action = "na.fail")
# model assumptions not great with poisson
# but convergence errors with genpois.. 
# fixed with zero inflated genpois 

summary(abund_plant)
# nothing significant

## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]


## Broadleaf abundance
abund_plant <- glmmTMB(broadleaf_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")
summary(abund_plant)

# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3

## Shrub abundance
abund_plant <- glmmTMB(shrub_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                         scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                         (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")
summary(abund_plant)

# wind speed (negative)
# lunar cycle (waxing crescent) (negative)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

abund_plant_sum <- as.data.frame(summary(abund_plant)$coefficients$cond)
abund_plant_sum$parameters <- row.names(abund_plant_sum)
row.names(abund_plant_sum) <- 1:nrow(abund_plant_sum)
abund_plant_sum$`z value` <- NULL
colnames(abund_plant_sum)[3] <- "p_value"
abund_plant_sum <- as.data.frame(abund_plant_sum)
abund_plant_sum$significance <- case_when(
  abund_plant_sum$p_value>=0.05 ~ "ns",
  abund_plant_sum$p_value<0.05 & abund_plant_sum$p_value>=0.01 ~ "*",
  abund_plant_sum$p_value<0.01 & abund_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
abund_plant_sum <- abund_plant_sum[,c(4,1,2,3,5)]

## check model assumptions
testDispersion(abund_plant) 
simulationOutput <- simulateResiduals(fittedModel = abund_plant, plot = F)
plot(simulationOutput) ## no assumptions violated with genpois
## check for multicolinearity
check_collinearity(abund_plant) ## all <3


# Richness 

ggplot(moths_final, aes(calluna_vulgaris, richness))+
  geom_point()+ 
  geom_smooth()+
  theme_classic()

ggplot(moths_final, aes(erica_cinerea, richness))+
  geom_point()+ 
  geom_smooth()+
  theme_classic()
# possibly 4 outliers

ggplot(moths_final, aes(vaccinium_myrtilus, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(grass, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(forb, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()
# two big outliers

ggplot(moths_final, aes(pinus_sylvestris, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()

ggplot(moths_final, aes(betula_spp, richness))+
  geom_point()+ 
  stat_smooth()+
  theme_classic()


rich_plant <- glmer(richness ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# wind speed (negative)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)

rich_betula <- ggpredict(rich_plant, terms="betula_spp") 
rich_pine <- ggpredict(rich_plant, terms="pinus_sylvestris") 
rich_vaccinium <- ggpredict(rich_plant, terms="vaccinium_myrtilus")

rich_betula$predictor <- "Betula spp"
rich_pine$predictor <- "Pinus sylvestris"
rich_vaccinium$predictor <- "Vaccinium myrtilus"

rich_betula$response <- "Total \nspecies \nrichness"
rich_pine$response <- "Total \nspecies \nrichness"
rich_vaccinium$response <- "Total \nspecies \nrichness"

rich_plants <- rbind(rich_betula, rich_pine, rich_vaccinium)


## Woodland richness
rich_plant <- glmer(woodland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# wind speed (negative)
# vaccinium myrtilus (positive)
# forb (positive)
# pinus sylvestris (positive)
# betula (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)


wood_rich_betula <- ggpredict(rich_plant, terms="betula_spp") 
wood_rich_forbs <- ggpredict(rich_plant, terms="forb") 
wood_rich_vaccinium <- ggpredict(rich_plant, terms="vaccinium_myrtilus")

wood_rich_betula$predictor <- "Betula spp"
wood_rich_forbs$predictor <- "Forbs"
wood_rich_vaccinium$predictor <- "Vaccinium myrtilus"

wood_rich_betula$response <- "Woodland \nspecies \nrichness"
wood_rich_forbs$response <- "Woodland \nspecies \nrichness"
wood_rich_vaccinium$response <- "Woodland \nspecies \nrichness"

rich_plants <- rbind(rich_plants, wood_rich_betula, wood_rich_forbs, wood_rich_vaccinium)

## Moorland richness
rich_plant <- glmer(moorland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# wind speed (negative)
# pinus sylvestris (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)

moor_rich_pine <- ggpredict(rich_plant, terms="pinus_sylvestris") 
moor_rich_pine$predictor <- "Pinus sylvestris"
moor_rich_pine$response <- "Moorland \nspecies \nrichness"
rich_plants <- rbind(rich_plants, moor_rich_pine)


## Grassland richness
rich_plant <- glmer(grassland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# wind speed (negative)
# grass (negative)
# forb (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)

ggpredict(rich_plant, terms="grass") %>% plot(add.data=TRUE)

## Conifer richness
rich_plant <- glmer(conifer_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# wind speed (negative)
# pinus sylvestris (positive)
# betula spp (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)

## Broadleaf richness
rich_plant <- glmer(broadleaf_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)

## Shrub richness
rich_plant <- glmer(shrub_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(erica_cinerea) + scale(vaccinium_myrtilus) + 
                      scale(grass) + scale(forb) + scale(pinus_sylvestris) + scale(betula_spp) + (1|visit) + 
                      (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
summary(rich_plant)
# min temp (positive)
# wind speed (negative)
# vaccinium myrtilus (positive)
# pinus sylvestris (positive)
# betula spp (positive)

rich_plant_sum <- as.data.frame(coef(summary(rich_plant)))
rich_plant_sum$parameters <- row.names(rich_plant_sum)
row.names(rich_plant_sum) <- 1:nrow(rich_plant_sum)
rich_plant_sum$`z value` <- NULL
colnames(rich_plant_sum)[3] <- "p_value"
rich_plant_sum <- as.data.frame(rich_plant_sum)
rich_plant_sum$significance <- case_when(
  rich_plant_sum$p_value>=0.05 ~ "ns",
  rich_plant_sum$p_value<0.05 & rich_plant_sum$p_value>=0.01 ~ "*",
  rich_plant_sum$p_value<0.01 & rich_plant_sum$p_value>=0.001 ~ "**",
  TRUE ~ "***"
)
rich_plant_sum <- rich_plant_sum[,c(4,1,2,3,5)]


## check model assumptions
testDispersion(rich_plant) 
simulationOutput <- simulateResiduals(fittedModel = rich_plant, plot = F)
plot(simulationOutput) ## quantile deviations
## check for multicolinearity
car::vif(rich_plant) ## all good (under 3)




############# Supplementary material plots

plant_predictions <- rbind(abund_plants, rich_plants)
plant_predictions$response  = factor(plant_predictions$response, levels=c("Total \nabundance", 
                                                                          "Total \nspecies \nrichness", "Woodland \ntotal \nabundance",
                                                                          "Woodland \nspecies \nrichness", "Moorland \ntotal \nabundance",
                                                                          "Moorland \nspecies \nrichness"))


abund_betula_raw <- moths_final[,c("tot_abund", "betula_spp")]
abund_betula_raw$response <- "Total \nabundance"
abund_betula_raw$predictor <- "Betula spp"
colnames(abund_betula_raw)[1:2] <- c("y", "x")

abund_pine_raw <- moths_final[,c("tot_abund", "pinus_sylvestris")]
abund_pine_raw$response <- "Total \nabundance"
abund_pine_raw$predictor <- "Pinus sylvestris"
colnames(abund_pine_raw)[1:2] <- c("y", "x")

abund_vaccinium_raw <- moths_final[,c("tot_abund", "vaccinium_myrtilus")]
abund_vaccinium_raw$response <- "Total \nabundance"
abund_vaccinium_raw$predictor <- "Vaccinium myrtilus"
colnames(abund_vaccinium_raw)[1:2] <- c("y", "x")

wood_abund_betula_raw <- moths_final[,c("woodland_abund", "betula_spp")]
wood_abund_betula_raw$response <- "Woodland \ntotal \nabundance"
wood_abund_betula_raw$predictor <- "Betula spp"
colnames(wood_abund_betula_raw)[1:2] <- c("y", "x")

wood_abund_pine_raw <- moths_final[,c("woodland_abund", "pinus_sylvestris")]
wood_abund_pine_raw$response <- "Woodland \ntotal \nabundance"
wood_abund_pine_raw$predictor <- "Pinus sylvestris"
colnames(wood_abund_pine_raw)[1:2] <- c("y", "x")

wood_abund_vaccinium_raw <- moths_final[,c("woodland_abund", "vaccinium_myrtilus")]
wood_abund_vaccinium_raw$response <- "Woodland \ntotal \nabundance"
wood_abund_vaccinium_raw$predictor <- "Vaccinium myrtilus"
colnames(wood_abund_vaccinium_raw)[1:2] <- c("y", "x")

moor_abund_betula_raw <- moths_final[,c("moorland_abund", "betula_spp")]
moor_abund_betula_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_betula_raw$predictor <- "Betula spp"
colnames(moor_abund_betula_raw)[1:2] <- c("y", "x")

moor_abund_pine_raw <- moths_final[,c("moorland_abund", "pinus_sylvestris")]
moor_abund_pine_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_pine_raw$predictor <- "Pinus sylvestris"
colnames(moor_abund_pine_raw)[1:2] <- c("y", "x")

moor_abund_vaccinium_raw <- moths_final[,c("moorland_abund", "vaccinium_myrtilus")]
moor_abund_vaccinium_raw$response <- "Moorland \ntotal \nabundance"
moor_abund_vaccinium_raw$predictor <- "Vaccinium myrtilus"
colnames(moor_abund_vaccinium_raw)[1:2] <- c("y", "x")

rich_betula_raw <- moths_final[,c("richness", "betula_spp")]
rich_betula_raw$response <- "Total \nspecies \nrichness"
rich_betula_raw$predictor <- "Betula spp"
colnames(rich_betula_raw)[1:2] <- c("y", "x")

rich_pine_raw <- moths_final[,c("richness", "pinus_sylvestris")]
rich_pine_raw$response <- "Total \nspecies \nrichness"
rich_pine_raw$predictor <- "Pinus sylvestris"
colnames(rich_pine_raw)[1:2] <- c("y", "x")

rich_vaccinium_raw <- moths_final[,c("richness", "vaccinium_myrtilus")]
rich_vaccinium_raw$response <- "Total \nspecies \nrichness"
rich_vaccinium_raw$predictor <- "Vaccinium myrtilus"
colnames(rich_vaccinium_raw)[1:2] <- c("y", "x")

wood_rich_betula_raw <- moths_final[,c("woodland_rich", "betula_spp")]
wood_rich_betula_raw$response <- "Woodland \nspecies \nrichness"
wood_rich_betula_raw$predictor <- "Betula spp"
colnames(wood_rich_betula_raw)[1:2] <- c("y", "x")

wood_rich_forb_raw <- moths_final[,c("woodland_rich", "forb")]
wood_rich_forb_raw$response <- "Woodland \nspecies \nrichness"
wood_rich_forb_raw$predictor <- "Forbs"
colnames(wood_rich_forb_raw)[1:2] <- c("y", "x")

wood_rich_vaccinium_raw <- moths_final[,c("woodland_rich", "vaccinium_myrtilus")]
wood_rich_vaccinium_raw$response <- "Woodland \nspecies \nrichness"
wood_rich_vaccinium_raw$predictor <- "Vaccinium myrtilus"
colnames(wood_rich_vaccinium_raw)[1:2] <- c("y", "x")

moor_rich_pine_raw <- moths_final[,c("moorland_rich", "pinus_sylvestris")]
moor_rich_pine_raw$response <- "Moorland \nspecies \nrichness"
moor_rich_pine_raw$predictor <- "Pinus sylvestris"
colnames(moor_rich_pine_raw)[1:2] <- c("y", "x")


raw_dat <- rbind(abund_betula_raw, abund_pine_raw, abund_vaccinium_raw, wood_abund_betula_raw, 
                 wood_abund_pine_raw, wood_abund_vaccinium_raw, moor_abund_betula_raw, 
                 moor_abund_pine_raw, moor_abund_vaccinium_raw, rich_betula_raw, 
                 rich_pine_raw, rich_vaccinium_raw, wood_rich_betula_raw, wood_rich_forb_raw,
                 wood_rich_vaccinium_raw, moor_rich_pine_raw)
raw_dat$response  = factor(raw_dat$response, levels=c("Total \nabundance",
                                                      "Total \nspecies \nrichness", "Woodland \ntotal \nabundance",
                                                      "Woodland \nspecies \nrichness", "Moorland \ntotal \nabundance",
                                                      "Moorland \nspecies \nrichness"))

plot2 <- ggplot()+
  geom_point(data=raw_dat, aes(x=x, y=y), colour="lightgrey")+
  geom_line(data=plant_predictions, aes(x=x, y=predicted))+
  geom_ribbon(data=plant_predictions, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2)+
  facet_grid(response ~ predictor, scales="free")+
  labs(x="", y="")+
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0))
plot2
a <- c("a#aa
        a#aa
        a#aa
        aa#a
        a#aa
        ##a#")

plot2 <- remove_facets(plot2, a)
plot2
ggsave(plot2, file="Graphs/Plant_cover_supplementary_plot.png", height=5, width=7)

