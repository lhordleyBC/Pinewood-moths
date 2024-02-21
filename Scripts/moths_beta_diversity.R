library(betapart)
rm(list = ls())

moths <- read.csv("Data/Moth_data.csv", header=TRUE)

# remove NAs in Quantity - abundance not recorded for that species 
# (these are all aggregate spp - so they aren't included in richness either so can be removed completely)
moths <- moths %>% drop_na(Quantity)

### Moorland ###
moths_moorland <- moths[moths$Treatment=="moorland",]
moths_moorland <- moths_moorland[,c(4,6)]

library(reshape2)
moths_moorland_pa <- dcast(moths_moorland, Plot~Common_name, length)
moths_moorland_pa[moths_moorland_pa>0] <-1

moths_moorland_pa <- moths_moorland_pa[,!names(moths_moorland_pa) %in% c("Plot", "Var.2")]

moorland_pa_core <- betapart.core(moths_moorland_pa)
moorland_pa_multi <- beta.multi(moorland_pa_core)
glimpse(moorland_pa_multi)
# total dissimilarity = 87.9%
# mainly comprised of turnover (72.9%)
# only small amount of nestedness (15%)

### Regen ###
moths_regen <- moths[moths$Treatment=="regeneration",]
moths_regen <- moths_regen[,c(4,6)]

library(reshape2)
moths_regen_pa <- dcast(moths_regen, Plot~Common_name, length)
moths_regen_pa[moths_regen_pa>0] <-1

moths_regen_pa <- moths_regen_pa[,!names(moths_regen_pa) %in% c("Plot", "Var.2")]

regen_pa_core <- betapart.core(moths_regen_pa)
regen_pa_multi <- beta.multi(regen_pa_core)
glimpse(regen_pa_multi)
# total dissimilarity = 88%
# mainly comprised of turnover (82.4%)
# only small amount of nestedness (5.6%)

### Woodland ###
moths_wood <- moths[moths$Treatment=="woodland",]
moths_wood <- moths_wood[,c(4,6)]

library(reshape2)
moths_wood_pa <- dcast(moths_wood, Plot~Common_name, length)
moths_wood_pa[moths_wood_pa>0] <-1

moths_wood_pa <- moths_wood_pa[,!names(moths_wood_pa) %in% c("Plot", "Var.2")]

wood_pa_core <- betapart.core(moths_wood_pa)
wood_pa_multi <- beta.multi(wood_pa_core)
glimpse(wood_pa_multi)
# total dissimilarity = 85.8%
# mainly comprised of turnover (74.3%)
# only small amount of nestedness (11.5%)

##### Species composition changes within treatments are mainly comprised of turnover
##### rather than nestedness 

# sampling across equal sites

moor_pa_samp <- beta.sample(moorland_pa_core)
regen_pa_samp <- beta.sample(regen_pa_core)
wood_pa_samp <- beta.sample(wood_pa_core)

dist_moor <- moor_pa_samp$sampled.values
dist_regen <- regen_pa_samp$sampled.values
dist_wood <- wood_pa_samp$sampled.values

dist_moor$Treatment <- "Moorland"
dist_regen$Treatment <- "Regeneration"
dist_wood$Treatment <- "Woodland"

dist_all <- rbind(dist_moor, dist_regen, dist_wood)
dist_all2 <- melt(setDT(dist_all), id.vars = "Treatment", variable.name = "Beta_div")
lookup <- c("beta.SNE"="Nestedness", "beta.SIM"="Turnover", "beta.SOR"="Total dissimilarity")
dist_all2$Beta_div <- as.character(lookup[dist_all2$Beta_div])
dist_all2$Beta_div  = factor(dist_all2$Beta_div, levels=c("Nestedness", 
                                                          "Turnover", "Total dissimilarity"))

div_bar <- ggplot(dist_all2, aes(x = Treatment, y = value, fill = Beta_div)) +
  geom_col()+
  labs(y="Dissimilarity value")+
  theme(legend.title=element_blank())
div_bar
ggsave(div_bar, file="Graphs/Beta_diversity_treatment_barchart.png", height=4, width=6)
# SOR = total dissimilarity
# SIM = turnover
# SNE = nestedness 


moor_pa_samp <- beta.sample(moorland_pa_core, sites=10, sample=100)
regen_pa_samp <- beta.sample(regen_pa_core, sites=10, sample=100)
wood_pa_samp <- beta.sample(wood_pa_core, sites=10, sample=100)

dist_moor <- moor_pa_samp$sampled.values
dist_regen <- regen_pa_samp$sampled.values
dist_wood <- wood_pa_samp$sampled.values

plot(density(dist_moor$beta.SOR), xlim= c(0, 0.9), ylim= c(0, 40), col="blue", xlab = 'Beta Diversity', main='', lwd=2)
lines(density(dist_moor$beta.SNE), col="blue", lty=2, lwd=2)
lines(density(dist_moor$beta.SIM), col="blue", lty=3, lwd=2)

lines(density(dist_regen$beta.SOR), col='red', lwd=2)
lines(density(dist_regen$beta.SNE), col='red', lty=2, lwd=2)
lines(density(dist_regen$beta.SIM), col='red', lty=3, lwd=2)

lines(density(dist_wood$beta.SOR), col='green', lwd=2)
lines(density(dist_wood$beta.SNE), col='green', lty=2, lwd=2)
lines(density(dist_wood$beta.SIM), col='green', lty=3, lwd=2)

# blue is moorland
# red is regen
# green is woodland

# solid line = total dissimilarity
# dashed line = nestedness
# dotted line = turnover

## regen has higher turnover than both woodland and moorland
## regen also has lower nestedness than woodland and moorland


### look at pairwise samples
pair_moorland <- beta.pair(moths_moorland_pa)
pair_moorland$beta.sim
pair_moorland$beta.sne
pair_moorland$beta.sor

plot(hclust(pair_moorland$beta.sim, method="average"), hang=-1, main='', sub='', xlab='')
# turnover

## All this compares sites within each treatment, but I assume we want to 
## compare between each treatment, e.g. how similar are sites across the 3 
## treatments? Maybe just take out treatment and treat as one big group
## to do the pairwise samples? 



