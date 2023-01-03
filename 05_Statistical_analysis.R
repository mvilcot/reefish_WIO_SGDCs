## -------------------------------------------------------------------------- ##
## Codes for the paper: Spatial genetic differentiation correlates with species 
## assemblage turnover across tropical reef fish lineages - Vilcot et al. 2023
## Script name: 05_Statistical_analysis
## Date: 02/01/2023
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## -------------------------------------------------------------------------- ##

#### -------- Setup ------------------------------------------------------------
source("00_Setup.R")

## Read diversity indices
Data_full <- read.csv("Results/Table_diversity_full.csv")
list_sp <- as.character(Data_full$Gen_s) # list of the 20 studied species
list_species <- as.character(Data_full$Genus_species) # list of the 20 studied species
list_pop <- c("MF",	"MV",	"MY",	"SC")

## Time-calibrated phylogeny of fishes and resctriction to the 20 studied species
rownames(Data_full) <- Data_full$Genus_species
TreeFull <- read.tree("Data/actinopt_full.trees")[[1]]
TreeSpecies <- keep.tip(TreeFull, list_species)
rownames(Data_full) <- Data_full$Gen_s



#### -------- Lineage-based Gamma SGDC -----------------------------------------
## Linear model
reg1 <- summary(lm(Data_full$Ht ~ Data_full$SR_Gamma))
reg1

## Mixed model
mod <- lmer(Ht ~ SR_Gamma + (1 | Family), data=Data_full) 
summary(mod) 
R2_reg1 <- r2beta (model = mod, partial = FALSE, method = "sgv") # mixed model R2
R2_reg1

# PLot
gg1 <- ggplot(Data_full, aes(y = Ht, x = SR_Gamma)) +
  geom_point(aes(color = Family), size = 2) +
  xlab(TeX("$\\gamma$ species diversity"))+
  ylab(TeX("$\\gamma$ genetic diversity"))+
  geom_text_repel(aes(label = Gen_s), size=2.5) +
  annotate('text', x=160, y=0.330, hjust = 1, vjust = 1, size=3.5,
           label=paste0("LM adjusted R² = ", round(reg1$adj.r.squared, 3), "\nP = ", round(coef(reg1)[2,4], 3))) +
  theme_Publication() +
  theme(legend.position="none") +
  labs(tag = "(c)")



#### -------- Lineage-based Alpha SGDC -----------------------------------------

## Linear model
reg2 <- summary(lm(Data_full$Hs ~ Data_full$SR_Alpha))
reg2

## Mixed model
mod <- lmer(Hs ~ SR_Alpha + (1 | Family), data=Data_full) 
summary(mod) 
R2_reg2 <- r2beta (model = mod, partial = FALSE, method = "sgv") 
R2_reg2

# PLot
gg2 <- ggplot(Data_full, aes(y = Hs, x = SR_Alpha)) +
  geom_point(aes(color = Family), size = 2) +
  xlab(TeX("$\\bar{\\alpha}$ species diversity"))+
  ylab(TeX("$\\bar{\\alpha}$ genetic diversity"))+
  geom_text_repel(aes(label = Gen_s), size=2.5) +
  annotate('text', x=124, y=0.33, hjust = 1, vjust=1, size=3.5,
           label=paste0("LM adjusted R² = ", round(reg2$adj.r.squared, 3), "\nP = ", round(coef(reg2)[2,4], 3))) +
  theme_Publication() +
  theme(legend.position="none") +
  labs(tag = "(b)")



#### -------- Lineage-based Beta SGDC ------------------------------------------

## Linear regression
reg3 <- summary(lm(log(Data_full$GstPPHedrick) ~ Data_full$JTU_multi))
reg3

## Mixed model
mod <- lmer(log(GstPPHedrick) ~ JTU_multi + (1 | Family), data=Data_full) 
summary(mod)
R2_reg3 <- r2beta (model = mod, partial = FALSE, method = "sgv")
R2_reg3

# Plot
gg3 <- ggplot(Data_full, aes(y = log(GstPPHedrick), x = JTU_multi)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(aes(color = Family), size = 2) +
  theme_Publication() +
  xlab(TeX("$\\beta_{jtu}$ species diversity"))+
  ylab(TeX("log ($\\beta$ genetic diversity)"))+
  geom_text_repel(aes(label = Gen_s), size=2.5) +
  annotate('text', x=0.05, y=-1.2, hjust = 0, vjust = 1, size=3.5,
           label=paste0("LM adjusted R² = ", round(reg3$adj.r.squared, 3), ", P = ", round(coef(reg3)[2,4], 3))) +
  coord_fixed(ratio = 0.06) +
  labs(tag = "(a)")


 
## Plot all together (Figure 2)
pdf("Figures/SGDC_Gamma_Alpha_Beta_multispecies.pdf", width = 9, height = 9)
plot_grid(gg3, plot_grid(gg2, gg1), ncol = 1) 
dev.off()





## Same with overall Jaccard beta dissimilarity (Figure S2)
reg4 <- summary(lm(log(Data_full$GstPPHedrick) ~ Data_full$JAC_multi))
reg4

mod <- lmer(log(GstPPHedrick) ~ JAC_multi + (1 | Family), data=Data_full) 
summary(mod)
R2_reg4 <- r2beta (model = mod, partial = FALSE, method = "sgv") 
R2_reg4

ggplot(Data_full, aes(y = log(GstPPHedrick), x = JAC_multi)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(aes(color = Family), size=2) +
  geom_text_repel(aes(label = Gen_s), size=2.5) +
  theme_Publication() +
  xlab(TeX("$\\beta_{jac}$ species diversity"))+
  ylab(TeX("log($\\beta$ genetic diversity)"))+
  annotate('text', x=0.28, y=-1.5, hjust=0, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg4$adj.r.squared, 3), ", p = ", round(coef(reg4)[2,4], 3)))

ggsave(filename = "Figures/SGDC_Beta_multispecies_JAC.pdf", width = 8, height = 5)





## Same without Parapercis hexophtalma (Figure S3)
Data_full2 <- Data_full[!(Data_full$Gen_s == "Par_h"),]
reg5 <- summary(lm(log(Data_full2$GstPPHedrick) ~ Data_full2$JTU_multi))
reg5

mod <- lmer(log(GstPPHedrick) ~ JTU_multi + (1 | Family), data=Data_full2) 
summary(mod) 
R2_reg5 <- r2beta(model = mod, partial = FALSE, method = "sgv") # te donne le R2 du modèle mixte 
R2_reg5

ggplot(Data_full2, aes(y = log(GstPPHedrick), x = JTU_multi)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(aes(color = Family), size=2) +
  geom_text_repel(aes(label = Gen_s), size=2.5) +
  theme_Publication() +
  xlab(TeX("$\\beta_{jtu}$ species diversity"))+
  ylab(TeX("log($\\beta$ genetic diversity)"))+
  annotate('text', x=0.05, y=-2.6, hjust=0, vjust=1,
           label=paste0("LM adjusted R² = ", sprintf("%.3f", reg5$adj.r.squared), ", P = ", sprintf("%.3f", coef(reg5)[2,4])))
ggsave(filename = "Figures/SGDC_Beta_multispecies_without_Parapapercis_hexophtalma.pdf", width = 8, height = 5)








#### -------- Beta SGDC by species ---------------------------------------------

## Read pairwise beta values
listJaccard <- readRDS("Intermediate/Jaccard_DistMat_sf.RDS")
listGst <- readRDS("Intermediate/GD_beta_DistMat.RDS") 

## Create empty data frame
plot_list_beta <- list()
BetaSGDC <- data.frame()

## Loop on species
for (i in 1:nrow(Data_full)){
  species <- as.character(Data_full[i, "Gen_s"])
  family <- as.character(Data_full[i, "Family"])
  distGst <- listGst[[species]]$GstPP_mean
  distJaccard <- listJaccard[[family]]$beta.jtu
  set.seed(999)
  
  ## Mantel test
  Mantel <- vegan::mantel(xdis = distGst, ydis = distJaccard)
  BetaSGDC[species, "Mantel_p"] <- Mantel$signif
  BetaSGDC[species, "Mantel_r"] <- Mantel$statistic
  
  ## Procruste test
  Protest <- protest(distGst, distJaccard)
  BetaSGDC[species,"Procruste_p"] <- Protest$signif
  BetaSGDC[species,"Procruste_r"] <- sqrt(1-Protest$ss)
  
  ## Keep only lower triangle of the matrix
  distGst <- as.matrix(distGst)
  distJaccard <- as.matrix(distJaccard) 
  distGst[upper.tri(distGst, diag = T)] <- NA
  distJaccard[upper.tri(distJaccard, diag = T)] <- NA
  
  ## Melt ant merge both species and genetic beta diversities
  merged <- merge(melt(distGst, value.name = "GstPP"),
                  melt(distJaccard, value.name = "JTU"))
  merged <- merged[complete.cases(merged), ] #remove NA
  merged$sites <- paste0(merged$Var1, "-", merged$Var2)
  merged$group <- ifelse(grepl("MV", merged$sites), "inter", "intra")

  ## Create SGDC plot for each species
  plot_list_beta[[i]] <-
    ggplot(merged, aes(JTU, GstPP)) +
    geom_point(aes(color=group), show.legend = FALSE, size=3) +
    scale_color_manual(values=c("#6BD800", "darkorange")) +
    theme_Publication() +
    xlab(TeX("$\\beta_{jtu}$ species diversity")) +
    ylab(TeX("$\\beta$ genetic diversity")) +
    geom_text_repel(aes(label = sites), size=3) +
    annotate('text', x=min(merged$JTU), y=max(merged$GstPP)+0.3*(max(merged$GstPP) - min(merged$GstPP)), 
             hjust = 0, vjust = 1,
             # label=paste0("Procruste R² = ", round(sqrt(1-Protest$ss), 4), "\np = ", round(Protest$signif, 4))) +
             label=paste0("r Mantel = ", round(Mantel$statistic, 4), "\np = ", round(Mantel$signif, 4))) +
    ggtitle(gsub("_", " ", as.character(Data_full[i, "Genus_species"])))
}

## Put all plots into one file
g <- arrangeGrob(grobs=plot_list_beta, ncol=2)
ggsave("Figures/SGDC_Beta_by_species.pdf", g, width = 10, height = 46, limitsize = F)

## Merge diversity dataset and SGDCs values
BetaSGDC$Gen_s <- rownames(BetaSGDC)
Data_full <- merge(Data_full, BetaSGDC) #to use later
BetaSGDC <- merge(Data_full[, c( "Gen_s", "Genus_species", "Family")], BetaSGDC, by = "Gen_s")
write.csv(BetaSGDC, "Results/SGDC_Beta_by_species.csv", row.names = F)


### ------ Test H0: observed p-value is not different from the null distributions of p-values ------ ###
pvalmean <- mean(BetaSGDC$Mantel_p)
null.distrib <- c()
for (i in 1:99999){
  null.distrib[i] <- mean(runif(nrow(BetaSGDC)))
}
pvaltest <- length(null.distrib[which(null.distrib <= pvalmean)])/99999 # 2-sided test p-value 

ggplot() +
  geom_histogram(aes(null.distrib), bins = 40, color="black") +
  geom_vline(xintercept = pvalmean, linetype = "dashed") +
  theme_Publication() +
  xlab("null distribution of p-value") +
  ylab("count") +
  annotate('text', x=0.79, y=9500, hjust=1, vjust=1,
           label=paste0("mean observed p-value = ", round(pvalmean, 4), "\np-value (test) < 0.001"))





#### -------- Beta SGDC ~ dispersal traits -------------------------------------

Data_Ucrit <- read.csv("Data/Ucrit_values_Fisher2005.csv")
Data_Ucrit <- Data_Ucrit[!is.na(Data_Ucrit$Ucrit_Family),] 
Data_Ucrit <- merge(Data_Ucrit, Data_full[,c("Gen_s", "Mantel_r", "Procruste_r")], by="Gen_s")



### ------ Compare Procruste with Mantel test ------ ###
reg <- summary(lm(BetaSGDC$Mantel_r ~ BetaSGDC$Procruste_r))
reg
ggplot(BetaSGDC, aes(y = Mantel_r, x = Procruste_r)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(size=2) +
  theme_Publication() +
  geom_text_repel(aes(label = Gen_s), size=4) +
  annotate('text', x=0.6, y=1, hjust=0, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg$adj.r.squared, 4), ", p = ", round(coef(reg)[2,4], 4)))

reg <- summary(lm(BetaSGDC$Mantel_p ~ BetaSGDC$Procruste_p))
reg
ggplot(BetaSGDC, aes(y = Mantel_p, x = Procruste_p)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(size=2) +
  theme_Publication() +
  geom_text_repel(aes(label = Gen_s), size=4) +
  annotate('text', x=0, y=1.1, hjust=0, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg$adj.r.squared, 4), ", p = ", round(coef(reg)[2,4], 4)))

cor(BetaSGDC$Mantel_r, BetaSGDC$Procruste_r, method = "spearman")




### ------ Mantel - PLD ------ ### 
## Linear model
reg6 <- summary(lm(Data_full$Mantel_r ~ log(Data_full$PLD_days)))
reg6

## PGLS
rownames(Data_full) <- Data_full$Genus_species
regGLS6 <- summary(gls(Mantel_r ~ log(PLD_days), correlation = corBrownian(phy = TreeSpecies, form =~Genus_species),
                       data = Data_full, method = "ML"))

## Plot
gg1 <- ggplot(Data_full, aes(log(PLD_days), Mantel_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3) +
  theme_Publication() +
  theme(axis.title=element_text(size=18))+
  xlab("log(PLD)")+
  ylab("Mantel r")+
  annotate('text', x=4.38, y=1.3, hjust=1, vjust=1, size=4.5,
           label=paste0("LM adjusted R² = ", round(reg6$adj.r.squared, 3), 
                        ", P = ", round(coef(reg6)[2,4], 4), 
                        "\nPGLS P = ", round(coef(regGLS6)[2,4], 3)))+
  labs(tag = "(a)")


### ------ Mantel - Ucrit ------ ### 
## Linear model
reg7 <- summary(lm(Data_Ucrit$Mantel_r ~ Data_Ucrit$Ucrit_Family))
reg7

## PGLS
rownames(Data_Ucrit) <- Data_Ucrit$Genus_species
TreeUcrit <- keep.tip(TreeSpecies, as.character(Data_Ucrit$Genus_species))
regGLS7 <- summary(gls(Mantel_r ~ Ucrit_Family, correlation = corBrownian(phy = TreeUcrit, form =~Genus_species),
                       data = Data_Ucrit, method = "ML"))
regGLS7

## Plot
gg2 <- ggplot(Data_Ucrit, aes(Ucrit_Family, Mantel_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3, max.overlaps = 20) +
  theme_Publication() +
  theme(axis.title=element_text(size=18))+
  xlab("U-crit")+
  ylab("Mantel r")+
  annotate('text', x=78, y=1.5, hjust=1, vjust=1, size=4.5,
           label=paste0("LM adjusted R² = ", round(reg7$adj.r.squared, 3), 
                        ", P = ", round(coef(reg7)[2,4], 5), 
                        "\nPGLS P = ", round(coef(regGLS7)[2,4], 3)))+
  labs(tag = "(b)")


gg1 / gg2 
ggsave("Figures/Dispersal_rMantel_comparison_PLD-Ucrit_bigger.pdf", width = 8, height = 12)





### ------ Procruste - PLD ------ ### 
## Linear model
reg8 <- summary(lm(Data_full$Procruste_r ~ log(Data_full$PLD_days)))
reg8

## PGLS 
rownames(Data_full) <- Data_full$Genus_species
regGLS8 <- summary(gls(Procruste_r ~ log(PLD_days), correlation = corBrownian(phy = TreeSpecies, form =~Genus_species),
                       data = Data_full, method = "ML"))
regGLS8

## Plot 
gg3 <- ggplot(Data_full, aes(log(PLD_days), Procruste_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3) +
  theme_Publication() +
  theme(axis.title=element_text(size=18))+
  xlab("log(PLD)")+
  ylab("Procruste r")+
  annotate('text', x=4.38, y=1.05, hjust=1, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg8$adj.r.squared, 3), 
                        ", P = ", round(coef(reg8)[2,4], 3), 
                        "\nPGLS P = ", round(coef(regGLS8)[2,4], 3))) +
  labs(tag = "(a)")



### ------ Procruste - Ucrit ------ ### 
## Linear model
reg9 <- summary(lm(Data_Ucrit$Procruste_r ~ Data_Ucrit$Ucrit_Family))
reg9

## PGLS 
rownames(Data_Ucrit) <- Data_Ucrit$Genus_species
TreeUcrit <- keep.tip(TreeSpecies, as.character(Data_Ucrit$Genus_species))
regGLS9 <- summary(gls(Procruste_r ~ Ucrit_Family, correlation = corBrownian(phy = TreeUcrit, form =~Genus_species),
                       data = Data_Ucrit, method = "ML"))
regGLS9

## Plot 
gg4 <- ggplot(Data_Ucrit, aes(Ucrit_Family, Procruste_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3, max.overlaps = 15) +
  theme_Publication() +
  theme(axis.title=element_text(size=18))+
  xlab("U-crit")+
  ylab("Procruste r")+
  annotate('text', x=78, y=1.1, hjust=1, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg9$adj.r.squared, 3), 
                        ", P = ", round(coef(reg9)[2,4], 4), 
                        "\nPGLS P = ", round(coef(regGLS9)[2,4], 3)))+
  labs(tag = "(b)")


gg3 / gg4
ggsave("Figures/Dispersal_rProcruste_comparison_PLD-Ucrit_bigger.pdf", width = 8, height = 12)




#### -------- Least cost distance ----------------------------------------------
# Geographic distance is computed as least-cost path distances, i.e. the length of the shortest path within water 
# between two sites.

coord_pop <- read.csv("Data/Sites_coordinates.csv", row.names = 1) # Population mean coordinates
Bathy <- getNOAA.bathy(lon1 = min(coord_pop[,"Longitude"])-5, lon2 = max(coord_pop[,"Longitude"])+5,
                       lat1 = min(coord_pop[,"Latitude"])-5, lat2 = max(coord_pop[,"Latitude"])+5,
                       resolution = 1)
saveRDS(Bathy, "Intermediate/Bathymetry.RDS")
Bathy <- readRDS("Intermediate/Bathymetry.RDS")
blues <- colorRampPalette(c("purple", "blue", "cadetblue1", "cadetblue2", "white"))

pdf(file = "Figures/Bathymetry_plot.pdf", height=7, width=9.5)
plot.bathy(Bathy, step=1000, deepest.isobath = -14000, shallowest.isobath = 0, image = TRUE, bpal = blues(20))
scaleBathy(Bathy, deg = 5, x = "topleft", inset = 5)
points(coord_pop[,"Longitude"], coord_pop[,"Latitude"],
       pch = 21, col = "black", bg = "red", cex = 1.3)
dev.off()

Trans <- trans.mat(Bathy, min.depth=-5, max.depth=NULL)
dist_geo <- lc.dist(Trans, coord_pop, res=c("dist","path"))
write.csv(as.matrix(dist_geo), file = "Intermediate/Least_cost_distance.csv")




#### -------- MRM - Isolation by barrier ---------------------------------------

### ------ Get pairwise distances into a melt table ------ ###
combin <- combn(list_pop, 2)
Combinaisons <- paste0(combin[1,], ".", combin[2,])
dist_geo_allSpecies <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), Combinaisons))
dist_geo <- as.matrix(read.csv(file = "Intermediate/Least_cost_distance.csv", row.names = 1))

for (species in list_sp){
  dist_geo_allSpecies[species,] <- c(dist_geo[combin[1,1], combin[2,1]], 
                                     dist_geo[combin[1,2], combin[2,2]], 
                                     dist_geo[combin[1,3], combin[2,3]], 
                                     dist_geo[combin[1,4], combin[2,4]], 
                                     dist_geo[combin[1,5], combin[2,5]], 
                                     dist_geo[combin[1,6], combin[2,6]])
}
dist_geo_allSpecies$Gen_s <- rownames(dist_geo_allSpecies)
dist_geo_melt <- melt(dist_geo_allSpecies, variable.name = "sites", value.name = "dist_geo")


### ------ Melt pairwise beta diversity indices ------ ###
pairwise_genet <- Data_full[, c("Gen_s", "Genus_species", "Family", "GstPP_MF.MV", "GstPP_MF.MY", "GstPP_MF.SC", "GstPP_MV.MY", "GstPP_MV.SC", "GstPP_MY.SC")]
pairwise_taxo <- Data_full[, c("Gen_s", "Genus_species", "Family", "JTU_MF.MV", "JTU_MF.MY", "JTU_MF.SC", "JTU_MV.MY", "JTU_MV.SC", "JTU_MY.SC")]

colnames(pairwise_genet) <- gsub("GstPP_", "", colnames(pairwise_genet))
colnames(pairwise_taxo) <- gsub("JTU_", "", colnames(pairwise_taxo))
melt_genet <- melt(pairwise_genet, variable.name = "sites", value.name = "GstPP")
melt_taxo <- melt(pairwise_taxo, variable.name = "sites", value.name = "JTU")
Data_melt <- merge(melt_genet, melt_taxo)
Data_melt <- merge(Data_melt, dist_geo_melt)

## Group in-between the barrier or not
Data_melt$group <- ifelse(grepl("MV", Data_melt$sites), "inter", "intra")




### ------ MRM ------ ###
Data_scaled <- Data_melt
Data_scaled$group <- ifelse(grepl("MV", Data_scaled$sites), 1, 0)
Data_scaled$logGstPP <- scale(log(Data_scaled$GstPP)) #GstPP log here
Data_scaled$JTU <- scale(Data_scaled$JTU)
Data_scaled$dist_geo <- scale(Data_scaled$dist_geo)

MRM(formula = logGstPP ~ dist_geo + group, data = Data_scaled, nperm = 99999)
MRM(formula = JTU ~ dist_geo + group, data = Data_scaled, nperm = 99999)


