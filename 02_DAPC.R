## ---------------------------------------------------------------------------------------------------------------- ##
## Script name: 02_DAPC
## Codes for the paper: Species-genetic diversity correlations in tropical reef fishes
## Date: 28 July 2020
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## Aim: To investigate the genetic structure at individual level, we applied 
##      a Discriminant Analysis of Principal Components (DAPC), a multivariate statistical approach
##      for which variance in the sample is partitioned into a betweengroup and a within-group components, 
##      in an effort to maximize discrimination between groups. DAPC was applied for each species 
##      to one of the randomized dataset using the dapc function implemented in the R package “adegenet”.
##      In sight of multi-species comparison, four clusters were defined as the four sampling sites 
##      and the same number of PCA axis were retained in the discriminant analysis. As Caranx melampygus
##      is only represented by 15 individuals, we choose to retain five PCA axis to avoid overfitting.
## ---------------------------------------------------------------------------------------------------------------- ##

#### -------- Setup ----------------------------------------------------------------------------------------------
source("00_Setup.R")

Database <- read.csv("Data/Dataset_species_traits.csv")
list_sp <- as.character(Database$Gen_s) # list of the 20 studied species
myCol <- c("#f8766d","#76ee00","#00bfc4","#febe01") #Best


#### -------- Loop by species ----------------------------------------------------------------------------------------
n_pca=5
par(oma=c(0,0,0,0))
# DAPC_var <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("PC1", "PC2", "PC3"))
for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  print(paste0("Species n°", sp, "/20: ", species)) 

  Genlight <- readRDS(paste0("Data/SNPs_ReeFISH/", species, "/radiator_genlight_", species, ".RDS"))  #Choose the first sampling
  Genlight$pop <- as.factor(as.character(Genlight$pop)) #so as population levels are in alphabetic order
  # Genlight <- Genlight[,sample(nLoc(Genlight), 4479, replace=F)]
  
  #### -------- DAPC -------------------------------------------------------------------------------------------------
  set.seed(999) # Setting a seed for a consistent result
  dapc <- dapc(Genlight, pop = Genlight$pop, n.da = 6, n.pca = n_pca)
  
  saveRDS(dapc, paste0("Intermediate/02_DAPC//DAPC_result_", n_pca, "PC_", species, ".RDS"))
  print.dapc(dapc) # decomment if needed
  summary.dapc(dapc)
  predict.dapc(dapc)

  #### -------- Plots ------------------------------------------------------------------------------------------------
  pdf(file = paste0("Results/Figures/DAPC/DAPC_plot_", n_pca, "PC_", species, ".pdf"), width = 6, height = 8)

  layout(matrix(c(1,2), nrow = 2, ncol = 1), heights=c(2.4,0.8))
  par(oma=c(1,1,2,1))
  scatter(dapc, scree.pca = F, legend = T, clabel = F, 
                posi.leg = "bottomleft", col=myCol,)
  par(mar=c(1,4,1,0))
  compoplot(dapc, lab=NULL, col=myCol, border = NA, legend = F, show.lab=F, cex.lab = 0.8)
  mtext(Database$Genus_species[sp], outer=TRUE, cex=1.5, line=0.5)
  
  dev.off()
  
  # DAPC_var[species,] <- c(dapc$eig/sum(dapc$eig))
}
# write.csv(DAPC_var, "Results/01_DAPC_percentage_axis_variance.csv", row.names = T)


