##########################
## User: Lisbeth Hordley
## Date: November 2023
## Info: Pinewood moth analysis

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

moths <- read.csv("Data/Moth_data.csv", header=TRUE)

# remove NAs in Quantity - abundance not recorded for that species 
# (these are all aggregate spp - so they aren't included in richness either so can be removed completely)
moths <- moths %>% drop_na(Quantity)

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


# calculate species richness and total abundance for each site and visit and for each guild
# note: richness calculated using common name to ensure 'moth spp' are not included in this
# but they are included in total abundance
moths_final <- moths %>% 
  group_by(Plot, Date, Treatment, Treatment_code, Sub_site, Cloud_cover, Min_temp, 
           Wind_speed, Wind_direction, Lunar_cycle, Precipitation) %>%
  summarise(richness = n_distinct(Common_name),
            woodland_rich = n_distinct(Common_name[X1..Woodland==1], na.rm=TRUE),
            moorland_rich = n_distinct(Common_name[X3_moorland==1], na.rm=TRUE),
            moorland_spec_rich = n_distinct(Common_name[moorland_specialist==1], na.rm=TRUE),
            grassland_rich = n_distinct(Common_name[X4_grassland==1], na.rm=TRUE),
            conifer_rich = n_distinct(Common_name[coniferous_trees==1], na.rm=TRUE),
            broadleaf_rich = n_distinct(Common_name[broadleaf_trees==1], na.rm=TRUE),
            shrub_rich = n_distinct(Common_name[shrubs_dwarf_shrubs==1], na.rm=TRUE),
            tot_abund=sum(Quantity),
            woodland_abund = sum(Quantity[X1..Woodland==1], na.rm=TRUE),
            moorland_abund = sum(Quantity[X3_moorland==1], na.rm=TRUE),
            moorland_spec_abund = sum(Quantity[moorland_specialist==1], na.rm=TRUE),
            grassland_abund = sum(Quantity[X4_grassland==1], na.rm=TRUE),
            conifer_abund = sum(Quantity[coniferous_trees==1], na.rm=TRUE),
            broadleaf_abund = sum(Quantity[broadleaf_trees==1], na.rm=TRUE),
            shrub_abund = sum(Quantity[shrubs_dwarf_shrubs==1], na.rm=TRUE))

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

moths_final <- moths_final %>%
  mutate(Treatment = recode(Treatment, "moorland" = 'Moorland', 
                            regeneration = 'Regeneration', woodland = 'Woodland'))

round(cor(moths_final[,c("richness","Cloud_cover","Min_temp","Wind_speed")]),3) # min temp and cloud cover highly correlated (r=0.8)
round(cor(moths_final[,c("tot_abund","Cloud_cover","Min_temp","Wind_speed")]),3) # min temp and cloud cover highly correlated (r=0.8)
cor.test(moths_final$richness, moths_final$Min_temp) # r=0.35, p<0.001
cor.test(moths_final$richness, moths_final$Cloud_cover) # r=0.34, p<0.001
cor.test(moths_final$tot_abund, moths_final$Min_temp) # r=0.23, p=0.018
cor.test(moths_final$tot_abund, moths_final$Cloud_cover) # r=0.2, p=0.04
# keep minimum temperature and remove cloud cover

# look at relationships with categorical covariates

# wind direction
ggplot(moths_final, aes(Wind_direction, richness))+
  geom_boxplot()+ 
  theme_classic()
ggplot(moths_final, aes(Wind_direction, tot_abund))+
  geom_boxplot()+ 
  theme_classic()

# lunar cycle
ggplot(moths_final, aes(Lunar_cycle, richness))+
  geom_boxplot()+ 
  theme_classic()
ggplot(moths_final, aes(Lunar_cycle, tot_abund))+
  geom_boxplot()+ 
  theme_classic()

# precipitation
ggplot(moths_final, aes(Precipitation, richness))+
  geom_boxplot()+ 
  theme_classic()
ggplot(moths_final, aes(Precipitation, tot_abund))+
  geom_boxplot()+ 
  theme_classic()

# Include all covariates apart from cloud cover 

ggplot(moths_final, aes(Treatment, tot_abund))+
  geom_boxplot()+ 
  theme_classic()

ggplot(moths_final, aes(Treatment, richness))+
  geom_boxplot()+ 
  theme_classic()

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

rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = richness)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Total species richness")+
  ylim(0,40)+
  xlab("Treatment") +
  theme_classic()
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

wood_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = woodland_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Woodland species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  theme_classic()
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

moor_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = moorland_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  theme_classic()
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
colnames(rich_treatment)[3] <- "std.error"
rich_treatment$df <- NULL
rich_treatment$t.ratio <- NULL
rich_treatment$statistic <- NA
colnames(rich_treatment)[4] <- "adj.p.value"
rich_treatment <- rich_treatment[,c(6,1,7,2,3,8,4,5)]
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = cld(glht1)$Treatment,
                      L = cld(glht1)$.group)
letters$L <- c("a", "a", "a")
mdata1 <- merge(moths_final, letters)

moor_rich_treatment_p2 <- ggplot(mdata1, aes(x=Treatment, y = moorland_spec_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=max(ddply(moths_final, "Treatment", summarise, fivenum(moorland_spec_rich)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland specialist species richness")+
  xlab("Treatment") +
  theme_classic()
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
rich_treatment$response <- "Gassland richness"
treatment_final <- rbind(treatment_final, rich_treatment)

## plot result
letters <- data.frame(Treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(moths_final, letters)

grass_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = grassland_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Grassland species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  theme_classic()
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

conifer_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = conifer_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Conifer species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  theme_classic()
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

broadleaf_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = broadleaf_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Broadleaf species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  theme_classic()
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

shrub_rich_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = shrub_rich)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=40, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Dwarf shrub species richness")+
  xlab("Treatment") +
  ylim(0,40)+
  theme_classic()
shrub_rich_treatment_p
# higher richness in woodland compared to regen and moorland

### Put all richness plots together 
richness_plots <- ggarrange(rich_treatment_p, wood_rich_treatment_p, moor_rich_treatment_p, grass_rich_treatment_p, 
                            conifer_rich_treatment_p, broadleaf_rich_treatment_p, shrub_rich_treatment_p, labels=c("(a)",
                            "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"),  hjust=0.05)
richness_plots
ggsave(richness_plots, file="Graphs/Moth_richness_treatment.png", height=10, width=10)

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

abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = tot_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Total abundance")+
  ylim(0,125)+
  xlab("Treatment") +
  theme_classic()
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

wood_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = woodland_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Woodland total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  theme_classic()
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

moor_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = moorland_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  theme_classic()
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

moor_abund_treatment_p2 <- ggplot(mdata1, aes(x=Treatment, y = moorland_spec_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=max(ddply(moths_final, "Treatment", summarise, fivenum(moorland_spec_abund)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Moorland specialist total abundance")+
  xlab("Treatment") +
  theme_classic()
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

grass_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = grassland_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Grassland total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  theme_classic()
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

conifer_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = conifer_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Conifer total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  theme_classic()
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

broadleaf_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = broadleaf_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Broadleaf total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  theme_classic()
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

shrub_abund_treatment_p <- ggplot(mdata1, aes(x=Treatment, y = shrub_abund)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0), colour="darkgrey") +
  geom_text(aes(Treatment, y=125, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Dwarf shrub total abundance")+
  xlab("Treatment") +
  ylim(0,125)+
  theme_classic()
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

treatment_final2$df <- NULL
treatment_final2$z.ratio <- NULL
treatment_final2$t.ratio <- NULL
treatment_final2 <- treatment_final2[,c(5,1:4)]
write.csv(treatment_final2, file="Pinewood_moths_treatment_abundance_results.csv", row.names=FALSE)



######################################################################################################################
######################################################################################################################



# read in habitat data
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
# Mean vegetation height

# first calculate the mean and SD across the ground layer height columns before calculating SD and CV
habitat$ground_layer_height_sd <- apply(habitat[, c(30:35)],1,sd)
habitat$ground_layer_height_mean <- rowMeans(habitat[,c(30:35)])
habitat$ground_layer_height_cv <- habitat$ground_layer_height_sd / habitat$ground_layer_height_mean
# CV = variability of ground layer height at each plot 

## Plant cover explanatory variables:
# Calluna vulgaris
# Erica cinerea
# Erica tetralix ?
# Vaccinium myrtilus
# Vaccinium vitis idaea
# Empetrum nigrum X
# Grass
# Forb ?
# Pinus sylvestris
# Betula spp
# Sorbus aucuparia X
# Salix spp X

## all variables are % cover to the nearest 5% ## 

# think it might be worth dredging this model? To determine which plant cover is most important 
# should do the same for the structural features too 
# or maybe backwards selection based on p-value and AIC
# need to decide whether to force the model to retain covariates or not 

# put both datasets together

moths_final <- merge(moths_final, habitat, by.x="Plot", by.y="plot")

## First check for correlations between covariates and structural features
round(cor(moths_final[,c("Min_temp","Wind_speed", "conifer_canopy", "broadleaf_canopy", "ground_layer_height_mean", 
                         "tree_regeneration", "ground_layer", "complexity_score", "ground_layer_height_cv")]),3) 
# tree regen and complexity score (0.66) - higher % cover of tree regen is also more complex
cor.test(moths_final$tree_regeneration, moths_final$complexity_score) # r=0.66, p<0.001

cor.test(moths_final$richness, moths_final$tree_regeneration) # r=0.027, p=0.8
cor.test(moths_final$richness, moths_final$complexity_score) # r=0.28, p=0.003
# remove tree regen from richness models
cor.test(moths_final$tot_abund, moths_final$tree_regeneration) # r=-0.08, p=0.4
cor.test(moths_final$tot_abund, moths_final$complexity_score) # r=0.25, p=0.009
# remove tree regen from abundance models

#### PCA of structural cover variables ####
pca_struc <- prcomp(moths_final[,c("conifer_canopy", "broadleaf_canopy", "tree_regeneration", "ground_layer",
                                   "complexity_score", "ground_layer_height_cv")])
summary(pca_struc) ## PC1 = 24.4%, PC2 = 22.8%, PC3 = 21.5%
struc_eigen <- get_eigenvalue(pca_struc)
struc_eigen ## first 4 dimensions have >1 eigenvalue - not ideal. Need first 3 at minimum
pca_struc$rotation
# PC1: higher tree regen, lower conifer canopy and lower ground layer cover
# PC2: higher conifer canopy and lower ground layer cover
# PC3: higher tree regen
# PC4: higher broadleaf canopy

struc_pca <- cbind(pca_struc$x[,1:2], moths_final[,3]) %>% as.data.frame()
struc_pca$PC1 <- as.numeric(struc_pca$PC1)
struc_pca$PC2 <- as.numeric(struc_pca$PC2)
colnames(struc_pca)[3] <- "Treatment"
p1 <- ggplot(struc_pca, aes(PC1, PC2, colour = Treatment)) +
  geom_point(size = 3, aes(shape = Treatment)) +
  scale_shape_manual(values=c(0,1,2,3,4))+
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), linewidth = 1)+
  # scale_colour_manual(values=c("#FF0000", "#0000FF", "darkgrey", "#00FF00", "#FF9900"))+
  # labs(x="Tree height", y="Mean vegetation height vs \npresence of False Brome \nand Calamagrostis")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
p1
# PC values are duplicated - due to lack of variance in structural data
# e.g. moorland plots may have the exact same values, like 80% for ground layer
# and 0 for canopy cover variables, therefore exact same PC values
# not sure this is ideal.. check plant data next

#### PCA of plant cover variables ####
pca_plant <- prcomp(moths_final[,c("calluna_vulgaris", "erica_cinerea", "erica_tetralix", 
                                   "vaccinium_myrtilus", "vaccinium_vitis.idaea",
                                   "grass", "forb", "pinus_sylvestris", "betula_spp")])
summary(pca_plant) ## PC1 = 31.2%, PC2 = 25.1%, PC3 = 22.1%
plant_eigen <- get_eigenvalue(pca_plant)
plant_eigen ## all dimensions have >1 eigenvalue - not ideal. Need first 4 at absolute minimum
pca_plant$rotation
# PC1: higher calluna vulgaris, lower cavvinium myrtilus and lower pine
# PC2: lower grass, higher vaccinium myrtilus and vitis idaea
# PC3: lower pine and culluna vulgaris, higher vaccinium vitis idaea and grass
# PC4: higher birch and vaccinium vitis idaea and lower myrtilus

plant_pca <- cbind(pca_plant$x[,1:2], moths_final[,3]) %>% as.data.frame()
plant_pca$PC1 <- as.numeric(plant_pca$PC1)
plant_pca$PC2 <- as.numeric(plant_pca$PC2)
colnames(plant_pca)[3] <- "Treatment"
p2 <- ggplot(plant_pca, aes(PC1, PC2, colour = Treatment)) +
  geom_point(size = 3, aes(shape = Treatment)) +
  scale_shape_manual(values=c(0,1,2,3,4))+
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), linewidth = 1)+
  # scale_colour_manual(values=c("#FF0000", "#0000FF", "darkgrey", "#00FF00", "#FF9900"))+
  # labs(x="Tree height", y="Mean vegetation height vs \npresence of False Brome \nand Calamagrostis")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
p2
# still quite a lot of duplicated values in the dataframe.. 
# plot does look better though
# concerned about needing at least 4 dimensions to capture variation though
# probably due to lack of variance in plant cover values


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

abund_struc <- glmmTMB(tot_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                       scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                       scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                       scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, 
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
check_collinearity(abund_struc) ## all good (under 3)

a <- ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE) # a lot of zeros in broadleaf canopy
b <- ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
c <- ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)
d <- ggpredict(abund_struc, terms="complexity_score") %>% plot(add.data=TRUE)

struc <- ggarrange(a,b,c,d)
struc
ggsave(struc, file="Graphs/Abundance_structural_features.png", height=6, width=8)


## Woodland abundance
abund_struc <- glmmTMB(woodland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

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
check_collinearity(abund_struc) ## lunar cycle is 3.09, but the rest are <3

ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE) # a lot of zeros in broadleaf canopy
ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)


## Moorland abundance
abund_struc <- glmmTMB(moorland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

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
check_collinearity(abund_struc) ## all <3

ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE) # a lot of zeros in broadleaf canopy
ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="complexity_score") %>% plot(add.data=TRUE)


## Grassland abundance
abund_struc <- glmmTMB(grassland_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

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
check_collinearity(abund_struc) ## lunar cycle 3.26

ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Conifer abundance
abund_struc <- glmmTMB(conifer_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

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
check_collinearity(abund_struc) ## lunar cycle 3.26

ggpredict(abund_struc, terms="complexity_score") %>% plot(add.data=TRUE)


## Broadleaf abundance
abund_struc <- glmmTMB(broadleaf_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

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
check_collinearity(abund_struc) ## all <3

ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Shrub abundance
abund_struc <- glmmTMB(shrub_abund ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                         scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                         scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                         scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="genpois", na.action = "na.fail")

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
check_collinearity(abund_struc) ## all <3

ggpredict(abund_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(abund_struc, terms="ground_layer") %>% plot(add.data=TRUE)



###########################################################################################################


## Richness

rich_struc <- glmer(richness ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                       scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                       scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

# try dredging model
rich_struc_d <- dredge(rich_struc)
rich_struc_topmods <- subset(rich_struc_d, delta <6)
rich_struc_avgmod <- model.avg(rich_struc_topmods, fit=TRUE) 
summary(rich_struc_avgmod)

a <- ggpredict(rich_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
b <- ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

struc2 <- ggarrange(a,b)
ggsave(struc2, file="Graphs/Richness_structural_features.png", height=5, width=8)


## Woodland richness

rich_struc <- glmer(woodland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

ggpredict(rich_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Moorland richness
rich_struc <- glmer(moorland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Grassland richness
rich_struc <- glmer(grassland_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Conifer richness
rich_struc <- glmer(conifer_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

ggpredict(rich_struc, terms="complexity_score") %>% plot(add.data=TRUE)

## Broadleaf richness
rich_struc <- glmer(broadleaf_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

ggpredict(rich_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="broadleaf_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)

## Shrub richness
rich_struc <- glmer(shrub_rich ~ scale(Min_temp) + scale(Wind_speed) + Lunar_cycle + 
                      scale(conifer_canopy) + scale(broadleaf_canopy) + scale(ground_layer) + 
                      scale(complexity_score) + scale(ground_layer_height_cv) + (1|visit) + 
                      scale(ground_layer_height_mean) + (1|Sub_site), data=moths_final, family="poisson", na.action = "na.fail")
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
car::vif(rich_struc) ## all good (under 3)

ggpredict(rich_struc, terms="conifer_canopy") %>% plot(add.data=TRUE)
ggpredict(rich_struc, terms="ground_layer") %>% plot(add.data=TRUE)


##########################################################################################################
##########################################################################################################


################### PLANT COVER MODELS ################### 

# Calluna vulgaris
# Erica cinerea
# Erica tetralix ?
# Vaccinium myrtilus
# Vaccinium vitis idaea X
# Empetrum nigrum X
# Grass
# Forb ?
# Pinus sylvestris
# Betula spp
# Sorbus aucuparia X
# Salix spp X

## correlations between covariates and plant cover variables
round(cor(moths_final[,c("Min_temp","Wind_speed", "calluna_vulgaris", "erica_cinerea",
                         "erica_tetralix", "vaccinium_myrtilus", "vaccinium_vitis.idaea",
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

a <- ggpredict(abund_plant, terms="vaccinium_myrtilus") %>% plot(add.data=TRUE)
b <- ggpredict(abund_plant, terms="pinus_sylvestris") %>% plot(add.data=TRUE)
c <- ggpredict(abund_plant, terms="betula_spp") %>% plot(add.data=TRUE)
abund_plants <- ggarrange(a,b,c)
ggsave(abund_plants, file="Graphs/Abundance_plant_cover.png", height=6, width=8)

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

c <- ggpredict(rich_plant, terms="vaccinium_myrtilus") %>% plot(add.data=TRUE)
d <- ggpredict(rich_plant, terms="pinus_sylvestris") %>% plot(add.data=TRUE)
e <- ggpredict(rich_plant, terms="betula_spp") %>% plot(add.data=TRUE)
rich_plants <- ggarrange(c,d,e, ncol=2, nrow=2)
ggsave(rich_plants, file="Graphs/Richness_plant_cover.png", height=8, width=8)


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
