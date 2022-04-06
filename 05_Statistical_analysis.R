## ---------------------------------------------------------------------------------------------------------------- ##
## Script name: 05_Statistical_analysis
## Codes for the paper: Species-genetic diversity correlations in tropical reef fishes
## Date: 28 July 2020
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## Aim: Assess the relationships between:
##        - alpha and gamma genetic diveristy
##        - genetic and species diversity for each diversity component (alpha, beta, gamma).
## ---------------------------------------------------------------------------------------------------------------- ##

#### -------- Setup ---------------------------------------------------------------------------------------------
source("00_Setup.R")

# Diversity measures
Data_full <- read.csv("Results/Table_diversity_full.csv")
list_sp <- as.character(Data_full$Gen_s) # list of the 20 studied species
list_species <- as.character(Data_full$Genus_species) # list of the 20 studied species
list_pop <- c("MF",	"MV",	"MY",	"SC")

# Time-calibrated phylogeny of fishes and resctriction to the 20 studied species
rownames(Data_full) <- Data_full$Genus_species
TreeFull <- read.tree("Data/actinopt_full.trees")[[1]]
TreeSpecies <- keep.tip(TreeFull, list_species)
rownames(Data_full) <- Data_full$Gen_s

# Donati <- readRDS("../Donati_2021/Biorxiv_Donati_dryad/Data/F_dataset_2.RDS")
# Test <- merge(Data_full, Donati[,c("Abd_cum", "Abd_relative", "Mean_Adb")], by = 0)


#### -------- Gamma SGDC multi-species --------------------------------------------------------------------------------------------
## Linear model
reg1 <- summary(lm(Data_full$Ht ~ Data_full$SR_Gamma))
reg1

## Mixed model
mod <- lmer(Ht ~ SR_Gamma + (1 | Family), data=Data_full) # facteur random 
summary(mod) # pente et p.value. Tu verras que ça change rien au t de student du modèle linéaire simple avec lm()
R2_reg1 <- r2beta (model = mod, partial = FALSE, method = "sgv") # te donne le R2 du modèle mixte 
R2_reg1

# PLot gamma-SGDC munti-species
gg1 <-
  ggplot(Data_full, aes(y = Ht, x = SR_Gamma)) +
  geom_point(aes(color = Family), size = 2) +
  theme_Publication() +
  xlab(TeX("$\\gamma$ species richness$"))+
  ylab(TeX("$\\gamma$ genetic diversity$"))+
  annotate('text', x=160, y=0.330, hjust = 1, vjust = 1,
           label=paste0("LM R² = ", round(reg1$adj.r.squared, 4), "\np-value = ", round(coef(reg1)[2,4], 4))) +
  theme(legend.position="none") +
  labs(tag = "c")
# ggsave(filename = "Results/Figures/SGDC_Gamma_multispecies.pdf", width = 8, height = 5)
# ggsave(filename = "Results/Figures/SGDC_Gamma_multispecies.png", width = 8, height = 5)



#### -------- Alpha SGDC multi-species -------------------------------------------------------------------------------------------
# Assessing the relationships between genetic and species diversity by applying a linear model
# between each component of GD (alpha, beta, gamma) and its respective SD. Because of non-independence 
# between observations in regressions where species diversity was replicated for species belonging to that family, 
# significance of linear models is confirmed by linear mixed model adding family as a random effect. 

## Linear model
reg2 <- summary(lm(Data_full$Hs ~ Data_full$SR_Alpha))
reg2

## Mixed model
mod <- lmer(Hs ~ SR_Alpha + (1 | Family), data=Data_full) # facteur random 
summary(mod) # pente et p.value. Tu verras que ça change rien au t de student du modèle linéaire simple avec lm()
R2_reg2 <- r2beta (model = mod, partial = FALSE, method = "sgv") # te donne le R2 du modèle mixte 
R2_reg2

# PLot alpha-SGDC munti-species
gg2 <-
  ggplot(Data_full, aes(y = Hs, x = SR_Alpha)) +
  # geom_smooth(method='lm', formula=y~x) +
  geom_point(aes(color = Family), size = 2) +
  theme_Publication() +
  xlab(TeX("$\\bar{\\alpha}$ species richness$"))+
  ylab(TeX("$\\bar{\\alpha}$ genetic diversity$"))+
  annotate('text', x=124, y=0.33, hjust = 1, vjust=1, size=4,
           label=paste0("LM R² = ", round(reg2$adj.r.squared, 4), "\np-value = ", round(coef(reg2)[2,4], 4))) +
  theme(legend.position="none") +
  labs(tag = "b")
# ggsave(filename = "Results/Figures/SGDC_Alpha_multispecies.pdf", width = 8, height = 5)
# ggsave(filename = "Results/Figures/SGDC_Alpha_multispecies.png", width = 8, height = 5)




# #### -------- Alpha SGDC by species ---------------------------------------------------------------------------------- #
# rownames(Data_full) <- Data_full$Gen_s
# plot_list_alpha <- list()
# AlphaSGDC <- data.frame()
# list_alpha_Hs <- c()
#   
# for (i in 1:nrow(Data_full)){
#   species <- list_sp[i]
#   
#   alphagenet <- t(Data_full[species, c("Hs_MF", "Hs_MY", "Hs_MV", "Hs_SC")])
#   rownames(alphagenet) <- gsub("Hs_", "", rownames(alphagenet))
#   alphataxo <- t(Data_full[species, c("SR_MF", "SR_MY", "SR_MV", "SR_SC")])
#   rownames(alphataxo) <- gsub("SR_", "", rownames(alphataxo))
#   merged <- data.frame(alphagenet, alphataxo)
#   colnames(merged) <- c("Hs", "SR")
#   # list_alpha_Hs <- append(list_alpha_Hs, max(alphagenet, na.rm=T) - min(alphagenet, na.rm = T))
#   
#   reg <- summary(lm(merged$Hs ~ merged$SR))
#   AlphaSGDC[species, "R2"] <- reg$adj.r.squared
#   AlphaSGDC[species, "pvalue"] <- coef(reg)[2,4]
# 
#   plot_list_alpha[[i]] <-
#   ggplot(merged, aes(SR, Hs)) +
#     # geom_smooth(method='lm', formula=y~x, color = "grey10") +
#     geom_point() +
#     theme_Publication() +
#     xlab(TeX("$\\alpha$ species richness$")) +
#     ylab(TeX("$\\alpha$ genetic diversity$")) +
#     # ylim(c(min(merged$Hs, na.rm = T), max(merged$Hs, na.rm = T)+0.0025))+
#     geom_text_repel(aes(label = rownames(merged)), size=4) +
#     annotate('text', x=min(merged$SR), y=max(merged$Hs, na.rm = T)+0.0025, hjust = 0, vjust = 1, size=4,
#              label=paste0("LM R² = ", round(reg$adj.r.squared, 4), "\np-value = ", round(coef(reg)[2,4], 4))) +
#     ggtitle(gsub("_", " ", as.character(Data_full[i, "Genus_species"])))
# }
# 
# g <- arrangeGrob(grobs=plot_list_alpha, ncol=2)
# ggsave("Results/Figures/SGDC_Alpha_by_species.png", g, width = 10, height = 46, limitsize = F)
# ggsave("Results/Figures/SGDC_Alpha_by_species.pdf", g, width = 10, height = 46, limitsize = F)
# 
# # Save LM alpha by species results
# AlphaSGDC$Gen_s <- rownames(AlphaSGDC)
# AlphaSGDC <- merge(Data_full[, c( "Gen_s", "Genus_species", "Family")], AlphaSGDC, by = "Gen_s")
# write.csv(AlphaSGDC, "Results/SGDC_Alpha_by_species.csv", row.names = F)
# 
# 
# # Hypothesis H0: observed p-value is not different from the null distributions of p-valuespvalmean <- mean(mtest.table$pval)
# pvalmean <- mean(AlphaSGDC$pvalue)
# null.distrib <- c()
# for (i in 1:99999){
#   null.distrib[i] <- mean(runif(nrow(AlphaSGDC)))
# }
# 
# hist(null.distrib, breaks = 20, main = " Test of the overall (over the species) residual correlation", xlim = c(0,1)) # null distribution of p-values
# abline(v = pvalmean, col = "red") # mean value of observed r
# length(null.distrib[which(null.distrib <= pvalmean)])/99999 # test p-value 





#### -------- Beta SGDC multi-species--------------------------------------------------------------------------------------

## Linear regression
reg3 <- summary(lm(log(Data_full$GstPPHedrick) ~ Data_full$JTU_multi))
reg3

## Mixed model
mod <- lmer(log(GstPPHedrick) ~ JTU_multi + (1 | Family), data=Data_full) # facteur random 
summary(mod) # pente et p.value. Tu verras que ça change rien au t de student du modèle linéaire simple avec lm()
R2_reg3 <- r2beta (model = mod, partial = FALSE, method = "sgv") # te donne le R2 du modèle mixte 
R2_reg3


## Plot (figure 2)
# getPalette = colorRampPalette(brewer.pal(12, "Dark2"))
# paletteblind <- safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                             # "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

gg3 <- ggplot(Data_full, aes(y = log(GstPPHedrick), x = JTU_multi)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(aes(color = Family), size = 2) +
  theme_Publication() +
  xlab(TeX("$\\beta_{multi}$ species turnover$"))+
  ylab(TeX("$\\beta_{multi}$ genetic differentiation$"))+
  # geom_text_repel(aes(label = Gen_s), size=4) +
  annotate('text', x=0.05, y=-1.2, hjust = 0, vjust = 1, 
           label=paste0("LM R² = ", round(reg3$adj.r.squared, 4), "; p-value = ", round(coef(reg3)[2,4], 4),
                        "\nGLMM pseudo-R² = ", round(R2_reg3$Rsq, 4))) +
  # theme(legend.position="none") +
  # scale_color_brewer(palette = "Paired") +
  coord_fixed(ratio = 0.06) +
  labs(tag = "a")
  # scale_color_manual(values = manualcolors)
  # scale_color_manual(values = getPalette(12))
# ggsave(filename = "Results/Figures/SGDC_Beta_multispecies.pdf", width = 8, height = 5)
# ggsave(filename = "Results/Figures/SGDC_Beta_multispecies.png", width = 8, height = 5)

# Method 1 plot
library(cowplot)
pdf("Results/Figures/SGDC_Gamma_Alpha_Beta_multispecies_test2.pdf", width = 9, height = 9)
png("Results/Figures/SGDC_Gamma_Alpha_Beta_multispecies_test2.png", width = 9, height = 9, units = 'in', res = 100)
plot_grid(gg3, plot_grid(gg2, gg1), ncol = 1) 
dev.off()

# Method 2 plot
gg3 / (gg2 + gg1)
ggsave("Results/Figures/SGDC_Gamma_Alpha_Beta_multispecies_test.png", width = 9, height = 9)
ggsave("Results/Figures/SGDC_Gamma_Alpha_Beta_multispecies_test.pdf", width = 9, height = 9)




## Same without Parapercis hexophtalma (Figure S..)
Data_full2 <- Data_full[!(Data_full$Gen_s == "Par_h"),]
reg4 <- summary(lm(log(Data_full2$GstPPHedrick) ~ Data_full2$JTU_multi))
reg4

mod <- lmer(log(GstPPHedrick) ~ JTU_multi + (1 | Family), data=Data_full2) # facteur random 
summary(mod) # pente et p.value. Tu verras que ça change rien au t de student du modèle linéaire simple avec lm()
R2_reg4 <- r2beta (model = mod, partial = FALSE, method = "sgv") # te donne le R2 du modèle mixte 
R2_reg4

ggplot(Data_full2, aes(y = log(GstPPHedrick), x = JTU_multi)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(aes(color = Family), size=2) +
  theme_Publication() +
  xlab(TeX("$\\beta_{multi}$ species turnover$"))+
  ylab(TeX("$\\beta_{multi}$ genetic differentiation$"))+
  # geom_text(aes(label = Gen_s), size=4) +
  annotate('text', x=0.05, y=-2.6, hjust=0, vjust=1,
           label=paste0("LM R² = ", round(reg4$adj.r.squared, 4), "; p-value = ", round(coef(reg4)[2,4], 4),
                        "\nGLMM pseudo-R² = ", round(R2_reg4$Rsq, 4)))
ggsave(filename = "Results/Figures/SGDC_Beta_multispecies_without_Parapapercis_hexophtalma.pdf", width = 8, height = 5)
ggsave(filename = "Results/Figures/SGDC_Beta_multispecies_without_Parapapercis_hexophtalma.png", width = 8, height = 5)





## Same with full Jaccard beta dissimilarity
reg5 <- summary(lm(log(Data_full$GstPPHedrick) ~ Data_full$JAC_multi))
reg5

mod <- lmer(log(GstPPHedrick) ~ JTU_multi + (1 | Family), data=Data_full) # facteur random 
summary(mod) # pente et p.value. Tu verras que ça change rien au t de student du modèle linéaire simple avec lm()
R2_reg5 <- r2beta (model = mod, partial = FALSE, method = "sgv") # te donne le R2 du modèle mixte 
R2_reg5

ggplot(Data_full, aes(y = log(GstPPHedrick), x = JAC_multi)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(aes(color = Family), size=2) +
  theme_Publication() +
  xlab(TeX("$\\beta_{multi}$ species dissimilarity"))+
  ylab(TeX("$\\beta_{multi}$ genetic differentiation$"))+
  # geom_text(aes(label = Gen_s), size=4) +
  annotate('text', x=0.28, y=-1.5, hjust=0, vjust=1,
           label=paste0("LM R² = ", round(reg5$adj.r.squared, 4), "; p-value = ", round(coef(reg5)[2,4], 4),
                        "\nGLMM pseudo-R² = ", round(R2_reg5$Rsq, 4)))
ggsave(filename = "Results/Figures/SGDC_Beta_multispecies_JAC.pdf", width = 8, height = 5)
ggsave(filename = "Results/Figures/SGDC_Beta_multispecies_JAC.png", width = 8, height = 5)






#### -------- Beta SGDC by species -----------------------------------------------------------------------------------
# Computation for each species of the correlation coefficient (r) associated with Procrustes superimposition 
# between genetic distance matrix (pairwise b-GD) and the corresponding 
# species distance matrix (pairwise b-SD) between the four sites.

listJaccard <- readRDS("Intermediate/Jaccard_DistMat.RDS")
listGst <- readRDS("Intermediate/GD_beta_DistMat.RDS") 

plot_list_beta <- list()
BetaSGDC <- data.frame()

for (i in 1:nrow(Data_full)){
  species <- as.character(Data_full[i, "Gen_s"])
  family <- as.character(Data_full[i, "Family"])
  distGst <- listGst[[species]]$GstPP_mean
  distJaccard <- listJaccard[[family]]$beta.jtu
  set.seed(999)
  
  Mantel <- vegan::mantel(xdis = distGst, ydis = distJaccard)
  BetaSGDC[species, "Mantel_p"] <- Mantel$signif
  BetaSGDC[species, "Mantel_r"] <- Mantel$statistic
  
  Protest <- protest(distGst, distJaccard)
  BetaSGDC[species,"Procruste_p"] <- Protest$signif
  BetaSGDC[species,"Procruste_r"] <- sqrt(1-Protest$ss)
  
  # Keep only lower triangle
  distGst <- as.matrix(distGst)
  distJaccard <- as.matrix(distJaccard)
  distGst[upper.tri(distGst, diag = T)] <- NA
  distJaccard[upper.tri(distJaccard, diag = T)] <- NA
  
  merged <- merge(melt(distGst, value.name = "GstPP"),
                  melt(distJaccard, value.name = "JTU"))
  merged <- merged[complete.cases(merged), ] #remove NA
  merged$sites <- paste0(merged$Var1, "-", merged$Var2)
  merged$group <- ifelse(grepl("MV", merged$sites), "inter", "intra")

  
  plot_list_beta[[i]] <-
    ggplot(merged, aes(JTU, GstPP)) +
    geom_point(aes(color=group), show.legend = FALSE, size=3) +
    scale_color_manual(values=c("#6BD800", "darkorange")) +
    theme_Publication() +
    xlab(TeX("$\\beta_{pair}$ species turnover$")) +
    ylab(TeX("$\\beta_{pair}$ genetic differentiation$")) +
    geom_text_repel(aes(label = sites), size=3) +
    annotate('text', x=min(merged$JTU), y=max(merged$GstPP)+0.3*(max(merged$GstPP) - min(merged$GstPP)), 
             hjust = 0, vjust = 1,
             # label=paste0("Procruste R² = ", round(sqrt(1-Protest$ss), 4), "\np-value = ", round(Protest$signif, 4))) +
             label=paste0("r Mantel = ", round(Mantel$statistic, 4), "\np-value = ", round(Mantel$signif, 4))) +
    ggtitle(gsub("_", " ", as.character(Data_full[i, "Genus_species"])))
}

# g <- arrangeGrob(grobs=plot_list_beta, ncol=2)
# ggsave("Results/Figures/SGDC_Beta_by_species.png", g, width = 10, height = 46, limitsize = F)
# ggsave("Results/Figures/SGDC_Beta_by_species.pdf", g, width = 10, height = 46, limitsize = F)

BetaSGDC$Gen_s <- rownames(BetaSGDC)
Data_full <- merge(Data_full, BetaSGDC) #to use later
BetaSGDC <- merge(Data_full[, c( "Gen_s", "Genus_species", "Family")], BetaSGDC, by = "Gen_s")
write.csv(BetaSGDC, "Results/SGDC_Beta_by_species.csv", row.names = F)


# Hypothesis H0: observed p-value is not different from the null distributions of p-valuespvalmean <- mean(mtest.table$pval)
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
  annotate('text', x=0.79, y=9500, hjust=1, vjust=1,
           label=paste0("mean observed p-value = ", round(pvalmean, 4), "\np-value (test) = ", round(pvaltest, 4)))
ggsave("Results/Figures/SGDC_Beta_by_species_hist_p_values_random.png", width = 8, height = 5)
ggsave("Results/Figures/SGDC_Beta_by_species_hist_p_values_random.pdf", width = 8, height = 5)
  

## Test Procruste vs. Mantel ##
reg <- summary(lm(BetaSGDC$Mantel_r ~ BetaSGDC$Procruste_r))
reg
ggplot(BetaSGDC, aes(y = Mantel_r, x = Procruste_r)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(size=2) +
  theme_Publication() +
  geom_text_repel(aes(label = Gen_s), size=4) +
  annotate('text', x=0.6, y=1, hjust=0, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg$adj.r.squared, 4), "; p-value = ", round(coef(reg)[2,4], 4)))+
  ggsave("Results/Figures/Test_Procruste_Mantel_betaSGDC_rcoeff.png", width = 8, height = 5)

reg <- summary(lm(BetaSGDC$Mantel_p ~ BetaSGDC$Procruste_p))
reg
ggplot(BetaSGDC, aes(y = Mantel_p, x = Procruste_p)) +
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(size=2) +
  theme_Publication() +
  geom_text_repel(aes(label = Gen_s), size=4) +
  annotate('text', x=0, y=1.1, hjust=0, vjust=1,
           label=paste0("LM adjusted R² = ", round(reg$adj.r.squared, 4), "; p-value = ", round(coef(reg)[2,4], 4)))+
  ggsave("Results/Figures/Test_Procruste_Mantel_betaSGDC_pvalue.png", width = 8, height = 5)

cor(BetaSGDC$Mantel_r, BetaSGDC$Procruste_r,
    method = "spearman"
)


#### -------- Beta SGDC ~ dispersal traits -----------------------------------------------------------------------
# To examine the effects of species traits on SGDCs strength, a linear model and a PGLS are implemented 
# between the Procrustes r of each species and dispersal traits

Data_Ucrit <- read.csv("Data/Ucrit_values.csv")
Data_Ucrit <- read.csv("Data/Ucrit_values_Maurine.csv")
Data_Ucrit <- merge(Data_full, Data_Ucrit, by = "Genus_species")
Data_Ucrit <- Data_Ucrit[!is.na(Data_Ucrit$Ucrit_Family),]

Data_Ucrit <- read.csv("Data/tropical_ucrit_data/output_mean_ucrit_family.csv")




## PLD - Mantel
reg6 <- summary(lm(Data_full$Mantel_r ~ log(Data_full$PLD_days)))
reg6

rownames(Data_full) <- Data_full$Genus_species
regGLS6 <- summary(gls(Mantel_r ~ log(PLD_days), correlation = corBrownian(phy = TreeSpecies),
                       data = Data_full, method = "ML"))


gg1 <- ggplot(Data_full, aes(log(PLD_days), Mantel_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3) +
  theme_Publication() +
  xlab("log(PLD)")+
  ylab("Mantel r")+
  annotate('text', x=4.38, y=1.3, hjust=1, vjust=1,
           label=paste0("LM R² = ", round(reg6$adj.r.squared, 4), 
                        "; p-value = ", round(coef(reg6)[2,4], 4), 
                        "\nPGLS p-value = ", round(coef(regGLS6)[2,4], 4)))+
  labs(tag = "a")
# ggsave(filename = "Results/Figures/Dispersal_rMantel_logPLD.png", width = 8, height = 6)
# ggsave(filename = "Results/Figures/Dispersal_rMantel_logPLD.pdf", width = 8, height = 6)



## Ucrit - Mantel
reg7 <- summary(lm(Data_Ucrit$Mantel_r ~ Data_Ucrit$Ucrit_Family))
reg7

rownames(Data_Ucrit) <- Data_Ucrit$Genus_species
TreeUcrit <- keep.tip(TreeSpecies, as.character(Data_Ucrit$Genus_species))
regGLS7 <- summary(gls(Mantel_r ~ Ucrit_Family, correlation = corBrownian(phy = TreeUcrit),
                       data = Data_Ucrit, method = "ML"))
regGLS7

gg2 <- ggplot(Data_Ucrit, aes(Ucrit_Family, Mantel_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3) +
  theme_Publication() +
  xlab("U-crit")+
  ylab("Mantel r")+
  annotate('text', x=78, y=1.5, hjust=1, vjust=1,
           label=paste0("LM R² = ", round(reg7$adj.r.squared, 4), 
                        "; p-value = ", round(coef(reg7)[2,4], 5), 
                        "\nPGLS p-value = ", round(coef(regGLS7)[2,4], 4)))+
  labs(tag = "b")
# ggsave(filename = "Results/Figures/Dispersal_rMantel_Ucrit.png", width = 8, height = 6)
# ggsave(filename = "Results/Figures/Dispersal_rMantel_Ucrit.pdf", width = 8, height = 6)


gg1 / gg2 
ggsave("Results/Figures/Dispersal_rMantel_comparison_PLD-Ucrit.png", width = 8, height = 12)
ggsave("Results/Figures/Dispersal_rMantel_comparison_PLD-Ucrit.pdf", width = 8, height = 12)





## PLD - Procruste
reg8 <- summary(lm(Data_full$Procruste_r ~ log(Data_full$PLD_days)))
reg8

rownames(Data_full) <- Data_full$Genus_species
regGLS8 <- summary(gls(Procruste_r ~ log(PLD_days), correlation = corBrownian(phy = TreeSpecies),
                       data = Data_full, method = "ML"))
regGLS8

gg3 <- ggplot(Data_full, aes(log(PLD_days), Procruste_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3) +
  theme_Publication() +
  xlab("log(PLD)")+
  ylab("Procruste r")+
  annotate('text', x=4.38, y=1.05, hjust=1, vjust=1,
           label=paste0("LM R² = ", round(reg8$adj.r.squared, 4), 
                        "; p-value = ", round(coef(reg8)[2,4], 4), 
                        "\nPGLS p-value = ", round(coef(regGLS8)[2,4], 4))) +
  labs(tag = "a")
# ggsave(filename = "Results/Figures/Dispersal_rProcruste_logPLD.png", width = 8, height = 6)
# ggsave(filename = "Results/Figures/Dispersal_rProcruste_logPLD.pdf", width = 8, height = 6)



## Ucrit - Procruste
reg9 <- summary(lm(Data_Ucrit$Procruste_r ~ Data_Ucrit$Ucrit_Family))
reg9

rownames(Data_Ucrit) <- Data_Ucrit$Genus_species
TreeUcrit <- keep.tip(TreeSpecies, as.character(Data_Ucrit$Genus_species))
regGLS9 <- summary(gls(Procruste_r ~ Ucrit_Family, correlation = corBrownian(phy = TreeUcrit),
                       data = Data_Ucrit, method = "ML"))
regGLS9

gg4 <- ggplot(Data_Ucrit, aes(Ucrit_Family, Procruste_r)) +
  geom_smooth(method='lm', formula=y~x, colour = "grey35") + 
  geom_point(size=2, color="grey10") +
  geom_text_repel(aes(label = Gen_s), size=3) +
  theme_Publication() +
  xlab("U-crit")+
  ylab("Procruste r")+
  annotate('text', x=78, y=1.1, hjust=1, vjust=1,
           label=paste0("LM R² = ", round(reg9$adj.r.squared, 4), 
                        "; p-value = ", round(coef(reg9)[2,4], 5), 
                        "\nPGLS p-value = ", round(coef(regGLS9)[2,4], 4)))+
  labs(tag = "b")
# ggsave(filename = "Results/Figures/Dispersal_rMantel_Ucrit.png", width = 8, height = 6)
# ggsave(filename = "Results/Figures/Dispersal_rMantel_Ucrit.pdf", width = 8, height = 6)


gg3 / gg4
ggsave("Results/Figures/Dispersal_rProcruste_comparison_PLD-Ucrit.png", width = 8, height = 12)
ggsave("Results/Figures/Dispersal_rProcruste_comparison_PLD-Ucrit.pdf", width = 8, height = 12)




#### -------- Least cost distance -----------------------------------------------------------------------------------
# Geographic distance is computed as least-cost path distances, i.e. the length of the shortest path within water 
# between two sites.

coord_pop <- read.csv("Data/Sites_coordinates.csv", row.names = 1) # Population mean coordinates
Bathy <- getNOAA.bathy(lon1 = min(coord_pop[,"Longitude"])-5, lon2 = max(coord_pop[,"Longitude"])+5,
                       lat1 = min(coord_pop[,"Latitude"])-5, lat2 = max(coord_pop[,"Latitude"])+5,
                       resolution = 1)
saveRDS(Bathy, "Intermediate/Bathymetry.RDS")
Bathy <- readRDS("Intermediate/Bathymetry.RDS")
blues <- colorRampPalette(c("purple", "blue", "cadetblue1", "cadetblue2", "white"))

pdf(file = "Results/Figures/Bathymetry_plot.pdf", height=7, width=9.5)
plot.bathy(Bathy, step=1000, deepest.isobath = -14000, shallowest.isobath = 0, image = TRUE, bpal = blues(20))
scaleBathy(Bathy, deg = 5, x = "topleft", inset = 5)
points(coord_pop[,"Longitude"], coord_pop[,"Latitude"],
       pch = 21, col = "black", bg = "red", cex = 1.3)
dev.off()

Trans <- trans.mat(Bathy, min.depth=-5, max.depth=NULL)
dist_geo <- lc.dist(Trans, coord_pop, res=c("dist","path"))
write.csv(as.matrix(dist_geo), file = "Intermediate/Least_cost_distance.csv")

#### -------- Isolation by distance (IBD) ----------------------------------------------------------------------------------
# Isolation by distance (IBD) is tested in a multispecies context, by performing a linear model between 
# pairwise genetic distance and marine geographic distance. As data are spatially non-independent, to test for significance of the correlation, a randomization test is
# performed with 99,999 permutations as explained above.

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


## Melt pairwise data table
pairwise_genet <- Data_full[, c("Gen_s", "Genus_species", "Family", "GstPP_MF.MV", "GstPP_MF.MY", "GstPP_MF.SC", "GstPP_MV.MY", "GstPP_MV.SC", "GstPP_MY.SC")]
pairwise_taxo <- Data_full[, c("Gen_s", "Genus_species", "Family", "JTU_MF.MV", "JTU_MF.MY", "JTU_MF.SC", "JTU_MV.MY", "JTU_MV.SC", "JTU_MY.SC")]

colnames(pairwise_genet) <- gsub("GstPP_", "", colnames(pairwise_genet))
colnames(pairwise_taxo) <- gsub("JTU_", "", colnames(pairwise_taxo))
melt_genet <- melt(pairwise_genet, variable.name = "sites", value.name = "GstPP")
melt_taxo <- melt(pairwise_taxo, variable.name = "sites", value.name = "JTU")
Data_melt <- merge(melt_genet, melt_taxo)
Data_melt <- merge(Data_melt, dist_geo_melt)
Data_melt$group <- ifelse(grepl("MV", Data_melt$sites), "inter", "intra")


## Linear regression
reglin <- summary(lm(data = Data_melt, formula = log(GstPP) ~ dist_geo)) # log
reglin

reglin <- summary(lm(data = Data_melt[Data_melt$group == "intra",], formula = log(GstPP) ~ dist_geo)) # log
reglin


## Randomization test
X <- Data_melt$dist_geo
Y <- log(Data_melt$GstPP)
reg_real <- summary(lm(Y ~ X))
R2real <- reg_real$adj.r.squared
Rreal <- cor(X,Y)
COEFreal <- coef(reg_real)[2,1]

COEFrandom <- list()
for (i in 1:9999){
  Xsample <- sample(X) # X randomization
  reg_random <- summary(lm(Y ~ Xsample))
  COEFrandom[i] <- coef(reg_random)[2,1] # slope
}
COEFrandom <- unlist(COEFrandom)
pvalue <- sum( abs(COEFrandom) >= abs(COEFreal))/length(COEFrandom)
pvalue


ggplot(as.data.frame(COEFrandom), aes(COEFrandom)) + # Histogram
  geom_histogram(color="black", bins = 50) +
  geom_vline(aes(xintercept=COEFreal, color = "real_coef"), linetype="dashed", size=1) +
  ggtitle("Muti-species IBD", subtitle = paste0("p-value = ", round(pvalue, 10)))
# ggsave("Results/Figures/IBD_randomization_test.pdf", width = 8, height = 6)


## IBD plot (Figure 3)
ggplot(Data_melt, aes(x = dist_geo, y = log(GstPP))) +
  geom_smooth(method='lm', formula=y~x, color="grey35") +
  geom_point(aes(color=group)) +
  scale_color_manual(values=c("#6BD800", "darkorange")) +
  theme_Publication() +
  xlab("geographic distance")+
  ylab(TeX("$\\beta_{pair}$ genetic differentiation$"))+
  annotate('text', x=800, y=-0.5, hjust = 0, vjust=1,
           label=paste0("LM R² = ", round(reglin$adj.r.squared, 4),
                        "\nrandomization p-value = ", round(pvalue, 4)))
ggsave(filename = "Results/Figures/IBD_multi_species.png", width = 8, height = 6)
ggsave(filename = "Results/Figures/IBD_multi_species.pdf", width = 8, height = 6)


## T test
t.test(log(Data_melt[Data_melt$group == "intra",]$GstPP), 
       log(Data_melt[Data_melt$group == "inter",]$GstPP))



#### -------- Isolation by barrier (IBB) - MRM ----------------------------------------------------------------------------------
# Isolation by barrier (IBB) is tested as we assume that Maldives can be potentially isolated 
# from the other sampling sites. The difference of pairwise beta genetic diversity between Maldives and the 
# three Southern-Western sites (Mafia Island, Mayotte, and, Seychelles), compared to pairwise beta genetic 
# diversity values among Southern-Western sites is examined by a t.test.

Data_scaled <- Data_melt
Data_scaled$group <- ifelse(grepl("MV", Data_scaled$sites), 1, 0)
Data_scaled$GstPP <- scale(log(Data_scaled$GstPP))
Data_scaled$JTU <- scale(Data_scaled$JTU)
Data_scaled$dist_geo <- scale(Data_scaled$dist_geo)

MRM(formula = GstPP ~ dist_geo + group, data = Data_scaled, nperm = 99999)
MRM(formula = JTU ~ dist_geo + group, data = Data_scaled, nperm = 99999)

MRM_list <- list()
for (species in list_sp){
  MRM_list[[species]]["GstPP"] <- MRM(formula = GstPP ~ dist_geo + group, data = Data_scaled[Data_scaled$Gen_s == species, ], nperm = 99999)
  MRM_list[[species]]["JTU"] <- MRM(formula = JTU ~ dist_geo + group, data = Data_scaled[Data_scaled$Gen_s == species, ], nperm = 99999)
}

MRM_table <- melt(MRM_list)
MRM_table <- MRM_table[MRM_table$Var1 != "Int",]
MRM_table <- MRM_table[!(MRM_table$Var2 %in% c("GstPP", "JTU")),]
MRM_table <- MRM_table[,c("L1", "L2", "Var1", "value")]
colnames(MRM_table) <- c("Gen_s", "metric", "variable", "pval")


# #### -------- Tables ----------------------------------------------------------------------------------
# 
# Data_full <- read.csv("Results/")
# 
# Table1 <- Data_full[,c("Gen_s", "Genus_species", "Family", "Body_size_cm", "PLD_days", "Hs", "Ht", "GstPPHedrick", "SR_Alpha", "SR_Gamma", "JTU_multi")]
# colnames(Table1) <- c("Code", "Species", "Family", "Body size", "PLD", "alpha (Hs)", "gamma (Ht)", "beta Gst''", "alpha SR", "gamma SR", "beta (JTU)")
# write.csv(Table1, "Results/Table1.csv", row.names = F)
# 
# 
# TableS3 <- Data_full[,c("Gen_s", "Genus_species", "Family", 
#                         "Hs", "Hs_MF", "Hs_MV", "Hs_MY", "Hs_SC", 
#                         "SR_Alpha", "SR_MF", "SR_MV", "SR_MY", "SR_SC")]
# TableS3 <- TableS3[order(TableS3$Family),]
# write.csv(TableS3, "Results/TableS3.csv", row.names = F)
# 
# TableS4 <- Data_full[,c("Gen_s", "Genus_species", "Family", 
#                         "GstPPHedrick", "GstPP_MF.MV", "GstPP_MF.MY", "GstPP_MF.SC", "GstPP_MV.MY", "GstPP_MV.SC", "GstPP_MY.SC")]
# TableS4 <- TableS4[order(TableS4$Family),]
# write.csv(TableS4, "Results/TableS4.csv", row.names = F)
# 
# TableS5 <- Data_full[,c("Family", 
#                         "JAC_multi", "JAC_MF.MV", "JAC_MF.MY", "JAC_MF.SC", "JAC_MV.MY", "JAC_MV.SC", "JAC_MY.SC",       
#                         "JTU_multi", "JTU_MF.MV", "JTU_MF.MY", "JTU_MF.SC", "JTU_MV.MY", "JTU_MV.SC", "JTU_MY.SC")]
# TableS5 <- TableS5[order(TableS5$Family),]
# TableS5 <- TableS5[!duplicated(TableS5$Family),]
# write.csv(TableS5, "Results/TableS5.csv", row.names = F)
# 
# 
# colnames(Data_full)
# 




