
##########################
## User: Lisbeth Hordley
## Date: November 2023
## Info: Pinewood moth analysis - differences between treatments

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
library(cowplot)
library(flextable)
library(stringr)
options(scipen=999)

moths <- read.csv("Data/Moth_data.csv", header=TRUE)

# remove NAs in Quantity - abundance not recorded for that species 
# (these are all aggregate spp - so they aren't included in richness either so can be removed completely)
moths <- moths %>% drop_na(Quantity)
moths <- moths %>% filter(!is.na(Quantity))

moth_traits <- read.csv("Data/moth_traits.csv", header=TRUE)

## Keep the following traits:
# X1..Woodland
# X3_moorland
# moorland_specialist
# X4_grassland
# coniferous_trees
# broadleaf_trees
# shrubs_dwarf_shrubs

moth_traits <- moth_traits[,c("scientific_name", "X1..Woodland", "X3_moorland", "moorland_specialist", "X4_grassland", 
                              "coniferous_trees", "broadleaf_trees", "shrubs_dwarf_shrubs")]
moth_traits[moth_traits == 3] <- NA
moth_traits[is.na(moth_traits)] <- 0
## merge in with moth data
moths <- merge(moths, moth_traits, by.x="Taxon", by.y="scientific_name", all.x=TRUE)

#**************** kp addition: rename columns *****************
names(moths)
moths<- moths %>% rename(woodland=17, moorland=18, grassland=20)

# calculate species richness and total abundance for each site and visit and for each guild
# note: richness calculated using common name to ensure 'moth spp' are not included in this
# but they are included in total abundance


moths_final <- moths %>% 
  filter(!is.na(Quantity))%>%
  group_by(Plot, Date, Treatment, Treatment_code, Sub_site, Cloud_cover, Min_temp, 
           Wind_speed, Wind_direction, Lunar_cycle, Precipitation) %>%
  summarise(### RICHNESS
            richness = n_distinct(Common_name),
            woodland_rich = n_distinct(Common_name[woodland==1], na.rm=TRUE),
            moorland_rich = n_distinct(Common_name[moorland==1], na.rm=TRUE),
            moorland_spec_rich = n_distinct(Common_name[moorland_specialist==1], na.rm=TRUE),
            grassland_rich = n_distinct(Common_name[grassland==1], na.rm=TRUE),
            conifer_rich = n_distinct(Common_name[coniferous_trees==1], na.rm=TRUE),
            broadleaf_rich = n_distinct(Common_name[broadleaf_trees==1], na.rm=TRUE),
            shrub_rich = n_distinct(Common_name[shrubs_dwarf_shrubs==1], na.rm=TRUE),
            ### ABUNDANCE
            tot_abund=sum(Quantity),
            woodland_abund = sum(Quantity[woodland==1], na.rm=TRUE),
            moorland_abund = sum(Quantity[moorland==1], na.rm=TRUE),
            moorland_spec_abund = sum(Quantity[moorland_specialist==1], na.rm=TRUE),
            grassland_abund = sum(Quantity[grassland==1], na.rm=TRUE),
            conifer_abund = sum(Quantity[coniferous_trees==1], na.rm=TRUE),
            broadleaf_abund = sum(Quantity[broadleaf_trees==1], na.rm=TRUE),
            shrub_abund = sum(Quantity[shrubs_dwarf_shrubs==1], na.rm=TRUE),
            #**************** kp addition: calculate shannon weiner diversity index *****************
            ### ALPHA DIVERSITY
            tot_div=diversity(Quantity, index="shannon"),
            woodland_div = diversity(Quantity[woodland%in%c(1)], index="shannon"),
            moorland_div = diversity(Quantity[moorland%in%c(1)], index="shannon"),
            moorland_spec_div = diversity(Quantity[moorland_specialist%in%c(1)], index="shannon"),
            grassland_div = diversity(Quantity[grassland%in%c(1)], index="shannon"),
            conifer_div = diversity(Quantity[coniferous_trees%in%c(1)], index="shannon"),
            broadleaf_div = diversity(Quantity[broadleaf_trees%in%c(1)], index="shannon"),
            shrub_div = diversity(Quantity[shrubs_dwarf_shrubs%in%c(1)], index="shannon")
            )

head(moths_final)%>%as.data.frame()

moths_final <- data.frame(moths_final)
colSums(is.na(moths_final)) # no NAs
str(moths_final)
moths_final$tot_abund <- as.numeric(moths_final$tot_abund)
moths_final$woodland_abund <- as.numeric(moths_final$woodland_abund)
moths_final$grassland_abund <- as.numeric(moths_final$grassland_abund)
moths_final$moorland_abund <- as.numeric(moths_final$moorland_abund)
moths_final$conifer_abund <- as.numeric(moths_final$conifer_abund)
moths_final$broadleaf_abund <- as.numeric(moths_final$broadleaf_abund)
moths_final$shrub_abund <- as.numeric(moths_final$shrub_abund)
moths_final$moorland_spec_abund <- as.numeric(moths_final$moorland_spec_abund)

## create visit variable - early (1) and late (2)
# trap 8 and 39 were empty on the first visit - data we have need to be recorded as second visit
moths_final$Date <- dmy(moths_final$Date)
moths_final$year <- year(moths_final$Date)
moths_final$month <- month(moths_final$Date)
moths_final$day <- day(moths_final$Date)

x <- moths_final %>% group_by(Plot) %>% distinct(month)

moths_final <- moths_final %>% group_by(Plot) %>%
  mutate(visit=ifelse(Date==min(Date), 1, 2))
# just need to fix traps 8 and 39 to be visit 2, not 1 (as traps were empty on visit 1)
moths_final$visit[moths_final$Plot == 8] <- 2
moths_final$visit[moths_final$Plot == 39] <- 2

#**************** kp amendment: change "Regeneration" to "Early Successional Woodland"  *****************
moths_final <- moths_final %>%
  mutate(Treatment = recode(Treatment, "moorland" = 'Moorland', 
                            regeneration = 'Early Successional Woodland', woodland = 'Mature Woodland'))


moths_final$Treatment

## save file
write.csv(moths_final, file="Data/Moths_data_final.csv", row.names=FALSE)


###############################################################################################################

moths_final <- read.csv("Data/Moths_data_final.csv", header=TRUE)

round(cor(moths_final[,c("richness","Cloud_cover","Min_temp","Wind_speed")]),3) # min temp and cloud cover highly correlated (r=0.8)
round(cor(moths_final[,c("tot_abund","Cloud_cover","Min_temp","Wind_speed")]),3) # min temp and cloud cover highly correlated (r=0.8)
cor.test(moths_final$richness, moths_final$Min_temp) # r=0.35, p<0.001
cor.test(moths_final$richness, moths_final$Cloud_cover) # r=0.34, p<0.001
cor.test(moths_final$tot_abund, moths_final$Min_temp) # r=0.23, p=0.018
cor.test(moths_final$tot_abund, moths_final$Cloud_cover) # r=0.2, p=0.04
#**************** kp addition: add cor tests for diversity  *****************
round(cor(moths_final[,c("tot_div","Cloud_cover","Min_temp","Wind_speed")]),3)
cor.test(moths_final$tot_div, moths_final$Min_temp) # r=0.38, p<0.001
cor.test(moths_final$tot_div, moths_final$Cloud_cover) # r=0.38, p<0.001
# keep minimum temperature and remove cloud cover

# look at relationships with categorical covariates

# wind direction
ggplot(moths_final, aes(Wind_direction, richness))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
ggplot(moths_final, aes(Wind_direction, tot_abund))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
#**************** kp addition:
ggplot(moths_final, aes(Wind_direction, tot_div))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()

# lunar cycle
ggplot(moths_final, aes(Lunar_cycle, richness))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
ggplot(moths_final, aes(Lunar_cycle, tot_abund))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
#**************** kp addition:
ggplot(moths_final, aes(Lunar_cycle, tot_div))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()

# precipitation
ggplot(moths_final, aes(Precipitation, richness))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
ggplot(moths_final, aes(Precipitation, tot_abund))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
#**************** kp addition:
ggplot(moths_final, aes(Precipitation, tot_div))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()

# Include all covariates apart from cloud cover 

ggplot(moths_final, aes(Treatment, tot_abund))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
ggplot(moths_final, aes(Treatment, richness))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
#**************** kp addition:
ggplot(moths_final, aes(Treatment, tot_div))+
  geom_boxplot()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
###########################################################
## First look at significant differences in treatment 
richness_mod <- glmer(richness ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                        (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(richness_mod)

testDispersion(richness_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(richness_mod) # all <3
r.squaredGLMM(richness_mod) # 40%

# Pairwise differences in management effects
glht1 <- glht(richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Total richness"
treatment_final <- rich_treatment

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = richness)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Total species richness")+
  ylim(0,40)+
  xlab("Treatment") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
rich_treatment_p
# no difference in species richness between moorland and regen 
# but woodland has significantly more species than moorland and regen

## woodland richness

wood_richness_mod <- glmer(woodland_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                             (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(wood_richness_mod)

testDispersion(wood_richness_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = wood_richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(wood_richness_mod) # all <3
r.squaredGLMM(wood_richness_mod) # 46%

# Pairwise differences in management effects
glht1 <- glht(wood_richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Woodland richness"
treatment_final <- rbind(treatment_final, rich_treatment)


## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

wood_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = woodland_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Woodland species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
wood_rich_treatment_p
# same as overall richness

## moorland richness

moor_richness_mod <- glmer(moorland_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                             (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(moor_richness_mod)

testDispersion(moor_richness_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = moor_richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(moor_richness_mod) # all <3
r.squaredGLMM(moor_richness_mod) # 37%

# Pairwise differences in management effects
glht1 <- glht(moor_richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# significant difference between woodland and moorland 
# no significant difference between woodland and regen
# no signficant difference between regen and moorland (only just sig.)

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Moorland richness"
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

moor_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = moorland_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
moor_rich_treatment_p
# no sig difference between regen and woodland now 

## moorland specialist richness

moor_richness_mod2 <- glmmTMB(moorland_spec_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                                (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(moor_richness_mod2)

testDispersion(moor_richness_mod2) # looks good
simulationOutput <- simulateResiduals(fittedModel = moor_richness_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(moor_richness_mod2) # all <3
r.squaredGLMM(moor_richness_mod2) # 28%

# Pairwise differences in management effects
glht1 <- lsmeans(moor_richness_mod2, c("Treatment")) 
summary(glht1) 
cld(glht1)

# no significant difference between woodland and moorland
# no significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
glht1 <- lsmeans(moor_richness_mod2, c("Treatment")) %>% pairs
rich_treatment <- summary(glht1)
rich_treatment$response <- "Moorland specialist richness"
rich_treatment$term <- "Treatment"
rich_treatment$null.value <- 0
names(rich_treatment)
colnames(rich_treatment)[3] <- "std.error"
#rich_treatment$df <- NULL
#rich_treatment$t.ratio <- NULL
rich_treatment$statistic <- NA
names(rich_treatment)
colnames(rich_treatment)[6] <- "adj.p.value"
names(treatment_final)
rich_treatment <- rich_treatment[colnames(treatment_final)]
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "a")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

moor_rich_treatment_p2 <- ggplot(mdata1, aes(x=Treatment, y = moorland_spec_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=max(ddply(moths_final, "Treatment", summarise, fivenum(moorland_spec_rich)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland specialist species richness")+
  xlab("Treatment") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
moor_rich_treatment_p2
# no sig difference between all treatments
ggsave(moor_rich_treatment_p, file="Graphs/Moorland_specialists_treatment.png")

## grassland richness

grass_richness_mod <- glmer(grassland_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                              (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(grass_richness_mod)

testDispersion(grass_richness_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = grass_richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(grass_richness_mod) # all <3
r.squaredGLMM(grass_richness_mod) # 35.5%

# Pairwise differences in management effects
glht1 <- glht(grass_richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# no significant differences

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Grassland richness"
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

grass_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = grassland_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Grassland species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
grass_rich_treatment_p
# no sig difference between all treatments

## conifer richness

conifer_richness_mod <- glmer(conifer_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                                (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(conifer_richness_mod)

testDispersion(conifer_richness_mod) # possible underdispersion..
simulationOutput <- simulateResiduals(fittedModel = conifer_richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(conifer_richness_mod) # all <3
r.squaredGLMM(conifer_richness_mod) # 31%

# Pairwise differences in management effects
glht1 <- glht(conifer_richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# significant difference between regen and moorland
# significant difference between woodland and moorland
# no significant difference between woodland and regen

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Conifer richness"
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

conifer_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = conifer_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Conifer species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
conifer_rich_treatment_p
# higher richness in regen and woodland compared to moorland

## broadleaf richness

broadleaf_richness_mod <- glmer(broadleaf_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                                  (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(broadleaf_richness_mod)

testDispersion(broadleaf_richness_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = broadleaf_richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(broadleaf_richness_mod) # all <3
r.squaredGLMM(broadleaf_richness_mod) # 35%

# Pairwise differences in management effects
glht1 <- glht(broadleaf_richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# no significant difference between regen and moorland
# significant difference between woodland and moorland
# significant difference between woodland and regen

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Broadleaf richness"
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

broadleaf_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = broadleaf_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Broadleaf species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
broadleaf_rich_treatment_p
# higher richness in woodland compared to regen and moorland

## shrub richness

shrub_richness_mod <- glmer(shrub_rich ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                              (1|visit) + (1|Sub_site), data = moths_final, family = poisson(), na.action = "na.fail")
summary(shrub_richness_mod)

testDispersion(shrub_richness_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = shrub_richness_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(shrub_richness_mod) # all <3
r.squaredGLMM(shrub_richness_mod) # 34%

# Pairwise differences in management effects
glht1 <- glht(shrub_richness_mod, mcp(Treatment="Tukey"))
summary(glht1) 
# no significant difference between regen and moorland
# significant difference between woodland and moorland
# significant difference between woodland and regen

# produce summary
CI <- summary(glht1)
rich_treatment <- data.frame(tidy(CI))
rich_treatment$response <- "Shrub richness"
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

shrub_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = shrub_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Dwarf shrub species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
shrub_rich_treatment_p
# higher richness in woodland compared to regen and moorland

### Put all richness plots together 
richness_plots <- ggarrange(rich_treatment_p, wood_rich_treatment_p, moor_rich_treatment_p, grass_rich_treatment_p, 
                            conifer_rich_treatment_p, broadleaf_rich_treatment_p, shrub_rich_treatment_p, labels=c("(a)",
                                                                                                                   "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"),  hjust=0.05)
richness_plots
ggsave(richness_plots, file="Graphs/Moth_richness_treatment.png", height=10, width=10)
ggsave(richness_plots, file="Graphs/Moth_richness_treatment.pdf", height=10, width=10)

treatment_final$term <- NULL
treatment_final$null.value <- NULL
treatment_final$statistic <- NULL
treatment_final <- treatment_final[,c(5,1:4)]
write.csv(treatment_final, file="Pinewood_moths_treatment_richness_results.csv", row.names=FALSE)

#################################################################################################

abund_mod <- glmmTMB(tot_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                       (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(abund_mod)

testDispersion(abund_mod) # looks ok - maybe slightly underdispersed? Much better with genpois
simulationOutput <- simulateResiduals(fittedModel = abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(abund_mod) # all <3
r.squaredGLMM(abund_mod) # 39%

# Pairwise differences in management effects
glht1 <- lsmeans(abund_mod, c("Treatment")) 
summary(glht1) 
cld(glht1)

# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "b", "c")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = tot_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Total abundance")+
  ylim(0,125)+
  xlab("Treatment") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
abund_treatment_p
# no difference in abundance between moorland and regen 
# but woodland has significantly higher abundance than moorland and regen

glht1 <- lsmeans(abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Total abundance"
treatment_final2 <- abund_treatment

## woodland abundance 

wood_abund_mod <- glmmTMB(woodland_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                            (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(wood_abund_mod)

testDispersion(wood_abund_mod) # looks ok - maybe slightly underdispersed? 
simulationOutput <- simulateResiduals(fittedModel = wood_abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(wood_abund_mod) # all <3
r.squaredGLMM(wood_abund_mod) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(wood_abund_mod, c("Treatment"))
summary(glht1) 
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "b", "c")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

wood_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = woodland_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Woodland total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
wood_abund_treatment_p

glht1 <- lsmeans(wood_abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Woodland abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)

## Moorland abundance

moor_abund_mod <- glmmTMB(moorland_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                            (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(moor_abund_mod)

testDispersion(moor_abund_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = moor_abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(moor_abund_mod) # all <3
r.squaredGLMM(moor_abund_mod) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(moor_abund_mod, c("Treatment"))
summary(glht1) 
cld(glht1)
# significant difference between woodland and moorland
# significant difference between woodland and regen
# significant difference between regen and moorland

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "b", "c")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

moor_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = moorland_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
moor_abund_treatment_p

glht1 <- lsmeans(moor_abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Moorland abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)


## Moorland specialist abundance

moor_abund_mod2 <- glmmTMB(moorland_spec_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                             (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(moor_abund_mod2)

testDispersion(moor_abund_mod2) # looks good
simulationOutput <- simulateResiduals(fittedModel = moor_abund_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(moor_abund_mod2) # all <3
r.squaredGLMM(moor_abund_mod2) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(moor_abund_mod2, c("Treatment"))
summary(glht1) 
cld(glht1)
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "a")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

moor_abund_treatment_p2 <- ggplot(mdata1, aes(x=Treatment, y = moorland_spec_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=max(ddply(moths_final, "Treatment", summarise, fivenum(moorland_spec_abund)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland specialist total abundance")+
  xlab("Treatment") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
moor_abund_treatment_p2
ggsave(moor_abund_treatment_p2, file="Graphs/Moorland_spec_abund_treatment.png")

glht1 <- lsmeans(moor_abund_mod2, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Moorland specialist abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)


## Grassland abundance

grass_abund_mod <- glmmTMB(grassland_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                             (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(grass_abund_mod)

testDispersion(grass_abund_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = grass_abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(grass_abund_mod) # all <3
r.squaredGLMM(grass_abund_mod) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(grass_abund_mod, c("Treatment"))
summary(glht1) 
cld(glht1)
# no significant differences

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "a")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

grass_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = grassland_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Grassland total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
grass_abund_treatment_p

glht1 <- lsmeans(grass_abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Grassland abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)

## Conifer abundance
conifer_abund_mod <- glmmTMB(conifer_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                               (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(conifer_abund_mod)

testDispersion(conifer_abund_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = conifer_abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(conifer_abund_mod) # all <3
r.squaredGLMM(conifer_abund_mod) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(conifer_abund_mod, c("Treatment"))
summary(glht1) 
cld(glht1)
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "b")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

conifer_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = conifer_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Conifer total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
conifer_abund_treatment_p

glht1 <- lsmeans(conifer_abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Conifer abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)

## Broadleaf abundance
broadleaf_abund_mod <- glmmTMB(broadleaf_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                                 (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(broadleaf_abund_mod)

testDispersion(broadleaf_abund_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = broadleaf_abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(broadleaf_abund_mod) # all <3
r.squaredGLMM(broadleaf_abund_mod) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(broadleaf_abund_mod, c("Treatment"))

summary(glht1) 
cld(glht1)
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "b")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

broadleaf_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = broadleaf_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Broadleaf total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
broadleaf_abund_treatment_p

glht1 <- lsmeans(broadleaf_abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Broadleaf abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)

## Shrub abundance

shrub_abund_mod <- glmmTMB(shrub_abund ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                             (1|visit) + (1|Sub_site), data = moths_final, family = "genpois", na.action = "na.fail")
summary(shrub_abund_mod)

testDispersion(shrub_abund_mod) # looks good
simulationOutput <- simulateResiduals(fittedModel = shrub_abund_mod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(shrub_abund_mod) # all <3
r.squaredGLMM(shrub_abund_mod) # 34%

# Pairwise differences in management effects
glht1 <- lsmeans(shrub_abund_mod, c("Treatment"))

summary(glht1) 
cld(glht1)
# significant difference between woodland and moorland
# significant difference between woodland and regen
# no signficant difference between regen and moorland

# produce summary
CI <- summary(glht1)
abund_treatment <- data.frame(tidy(CI))

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "b")
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

shrub_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = shrub_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Dwarf shrub total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
shrub_abund_treatment_p

glht1 <- lsmeans(shrub_abund_mod, c("Treatment")) %>% pairs
abund_treatment <- summary(glht1)
abund_treatment$response <- "Shrub abundance"
treatment_final2 <- rbind(treatment_final2, abund_treatment)

## put all abundance plots together

abundance_plots <- ggarrange(abund_treatment_p, wood_abund_treatment_p, moor_abund_treatment_p, grass_abund_treatment_p, 
                             conifer_abund_treatment_p, broadleaf_abund_treatment_p, shrub_abund_treatment_p, 
                             labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"))
abundance_plots
ggsave(abundance_plots, file="Graphs/Moth_abundance_treatment.png", height=10, width=10)
ggsave(abundance_plots, file="Graphs/Moth_abundance_treatment.pdf", height=10, width=10)

treatment_final2$df <- NULL
treatment_final2$z.ratio <- NULL
treatment_final2$t.ratio <- NULL
treatment_final2 <- treatment_final2[,c(5,1:4)]
write.csv(treatment_final2, file="Pinewood_moths_treatment_abundance_results.csv", row.names=FALSE)

####################################################################################
#**************** kp addition: models and summaries for diversity  *****************
#*
################# Total Diversity

dmod <- glmer(tot_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                        (1|visit) + (1|Sub_site), data = moths_final, family=gaussian, na.action = "na.fail")

summary(dmod)

testDispersion(dmod) # looks good
simulationOutput <- simulateResiduals(fittedModel = dmod, plot = F)
plot(simulationOutput) 

summary(dmod)
car::vif(dmod) # all <3
r.squaredGLMM(dmod) # 

# Pairwise differences in management effects
glht1 <- glht(dmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Alpha Diversity (Shannon Index)"
treatment_final <- d_treatment
treatment_final
## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
# Change the order of levels
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))


d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = tot_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment,y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Total diversity (Shannon)")+
 # ylim(0,40)+
  xlab("Treatment") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
d_treatment_p+   scale_x_discrete(labels = function(x) str_wrap(x, width = 10))


############### Woodland species diversity model

hist(moths_final$woodland_div)
wooddmod <- glmer(woodland_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                             (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")
summary(wooddmod)

testDispersion(wooddmod) # looks good
simulationOutput <- simulateResiduals(fittedModel = wooddmod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(wooddmod) # all <3
r.squaredGLMM(wooddmod) # 46%

# Pairwise differences in management effects
glht1 <- glht(wooddmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Woodland Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)


## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

wood_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = woodland_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Woodland diversity (Shannon)")+
  xlab("Treatment") +
 # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
wood_d_treatment_p


############### moorland species diversity model
## 
hist(moths_final$moorland_div)
moordmod <- glmer(moorland_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                    (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")
summary(moordmod)

testDispersion(moordmod) # looks good
simulationOutput <- simulateResiduals(fittedModel = moordmod, plot = F)
plot(simulationOutput) ## no assumptions violated

car::vif(moordmod) # all <3
r.squaredGLMM(moordmod) # 

# Pairwise differences in management effects
glht1 <- glht(moordmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Moorland Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)


## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

moor_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = moorland_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland diversity (Shannon)")+
  xlab("Treatment") +
  # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
moor_d_treatment_p


############### Moorland specialist species diversity model
## woodland richness
hist(moths_final$moorland_spec_div)
hist(log(moths_final$moorland_spec_div))
moorsdmod <- glmer(moorland_spec_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                    (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")
summary(moorsdmod)

testDispersion(moorsdmod) # looks good
simulationOutput <- simulateResiduals(fittedModel = moorsdmod, plot = F)
plot(simulationOutput) ##


summary(moorsdmod)


car::vif(moorsdmod) # all <3
r.squaredGLMM(moorsdmod) #

# Pairwise differences in management effects
glht1 <- glht(moorsdmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Moorland Specialist Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)


## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

moorspec_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = moorland_spec_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=2, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland specialist diversity (Shannon)")+
  xlab("Treatment") +
  # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
moorspec_d_treatment_p
ggsave(moor_rich_treatment_p, file="Graphs/Moorland_specialists_treatment_diversity.png")


############### grassland species diversity model
hist(moths_final$grassland_div)
hist(log(moths_final$grassland_div))
grassdmod <- glmer(grassland_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                     (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")
summary(grassdmod)

testDispersion(grassdmod) # looks good
simulationOutput <- simulateResiduals(fittedModel = grassdmod, plot = F)
plot(simulationOutput) 
moths_final$Treatment<- as.factor(moths_final$Treatment)
### Some assumptions violated - try a zero inflated model
grassdmod <- glmmTMB(grassland_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                     (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail", ziformula=~1)

testDispersion(grassdmod) # looks good
simulationOutput <- simulateResiduals(fittedModel = grassdmod, plot = F)
plot(simulationOutput) 

r.squaredGLMM(grassdmod) #

# Pairwise differences in management effects
glht1 <- glht(grassdmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Grassland Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)
treatment_final

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

grass_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = grassland_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Grassland diversity (Shannon)")+
  xlab("Treatment") +
  # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
grass_d_treatment_p




############### conifer species diversity model
hist(moths_final$conifer_div)
coniferdmod <- glmer(conifer_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                     (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")
summary(coniferdmod)

testDispersion(coniferdmod) 
simulationOutput <- simulateResiduals(fittedModel = coniferdmod, plot = F)
plot(simulationOutput) ##


car::vif(coniferdmod) # all <3
r.squaredGLMM(coniferdmod) #

# Pairwise differences in management effects
glht1 <- glht(coniferdmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Conifer Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)

treatment_final
## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

conifer_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = conifer_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Conifer diversity (Shannon)")+
  xlab("Treatment") +
  # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
conifer_d_treatment_p


############### shrub species diversity model
hist(moths_final$shrub_div)
#hist(log(moths_final$grassland_div))
shrubdmod <- glmer(shrub_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                       (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")

summary(shrubdmod)

testDispersion(shrubdmod) 
simulationOutput <- simulateResiduals(fittedModel = shrubdmod, plot = F)
plot(simulationOutput) ##


car::vif(shrubdmod) # all <3
r.squaredGLMM(shrubdmod) #

# Pairwise differences in management effects
glht1 <- glht(shrubdmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Shrub Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)

treatment_final
## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

shrub_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = shrub_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Dwarf shrub diversity (Shannon)")+
  xlab("Treatment") +
  # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
shrub_d_treatment_p




############### broadleaf species diversity model
hist(moths_final$broadleaf_div)
#hist(log(moths_final$grassland_div))
broaddmod <- glmer(broadleaf_div ~ Treatment + scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                     (1|visit) + (1|Sub_site), data = moths_final, family = gaussian, na.action = "na.fail")

summary(broaddmod)

testDispersion(broaddmod) 
simulationOutput <- simulateResiduals(fittedModel = broaddmod, plot = F)
plot(simulationOutput) ##


car::vif(broaddmod) # all <3
r.squaredGLMM(broaddmod) #

# Pairwise differences in management effects
glht1 <- glht(broaddmod, mcp(Treatment="Tukey"))
summary(glht1) 


# produce summary
CI <- summary(glht1)
d_treatment <- data.frame(tidy(CI))
d_treatment$response <- "Broadleaf Alpha Diversity (Shannon Index)"
treatment_final <- rbind(treatment_final, d_treatment)


## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)
mdata1$Treatment <- factor(mdata1$Treatment, levels = c("Moorland", "Early Successional Woodland", "Mature Woodland"))

broad_d_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = broadleaf_div)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=3.5, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Broadleaf diversity (Shannon)")+
  xlab("Treatment") +
  # ylim(0,40)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+   theme_classic()
broad_d_treatment_p

## put all diversity plots together

d_plots <- ggarrange(d_treatment_p, wood_d_treatment_p, moor_d_treatment_p, grass_d_treatment_p, 
                             conifer_d_treatment_p, broad_d_treatment_p, shrub_d_treatment_p, 
                             labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"))
d_plots
ggsave(d_plots, file="Graphs/Moth_diversity_treatment.png", height=10, width=10)
ggsave(d_plots, file="Graphs/Moth_diversity_treatment.pdf", height=10, width=10)

#treatment_final$df <- NULL
#treatment_final$z.ratio <- NULL
#treatment_final$t.ratio <- NULL
#treatment_final <- treatment_final[,c(5,1:4)]
treatment_final
write.csv(treatment_final, file="Pinewood_moths_treatment_diversity_results.csv", row.names=FALSE)


d_plots1 <- ggarrange(moor_rich_treatment_p2, moor_abund_treatment_p2, moorspec_d_treatment_p, nrow=1)#, 
                    # labels=c("(a)", "(b)", "(c)",))
d_plots1
ggsave(d_plots1, file="Graphs/Moth_moorland_spec_treatment.png", height=4, width=10)

##*********  KP addition: Produce table from treatment comparison results **************
l<-list.files(pattern="Pinewood")
abund<-read.csv(paste(l[1]))
div<- read.csv(paste(l[2]))
rich<-read.csv(paste(l[3]))
abund
div
rich
names(abund)
names(div)
names(rich)
div<- div %>% rename(SE=std.error, p.value=adj.p.value)
rich<- rich %>% rename(SE=std.error, p.value=adj.p.value)
div<- div[colnames(abund)]
rich<-rich[colnames(abund)]
head(abund)
head(div)
head(rich)
table<- rbind(rich, abund, div)
table
table<-table %>% mutate_if(is.numeric, function(x) round(x, digits=3))
table$p.value[table$p.value<=0.001] <- "<0.001"
table$p.value[table$p.value<=0.01&table$p.value!="<0.001"] <- "<0.01"
table$p.value[table$p.value<=0.05&!table$p.value%in%c("<0.001", "<0.01")] <- "<0.01"
table<- table%>%flextable() 
table
save_as_docx(table, 
             path="TableS4.docx",
             align = "left")
