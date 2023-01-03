## -------------------------------------------------------------------------- ##
## Codes for the paper: Spatial genetic differentiation correlates with species 
## assemblage turnover across tropical reef fish lineages - Vilcot et al. 2023
## Script name: 02_DAPC
## Date: 02/01/2023
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## -------------------------------------------------------------------------- ##

#### -------- Setup ------------------------------------------------------------
source("00_Setup.R")

## list the 20 species studied
Database <- read.csv("Data/Dataset_species_traits.csv")
list_sp <- as.character(Database$Gen_s) 

## Setup colors
myCol <- c("#f8766d","#76ee00","#00bfc4","#febe01")
par(oma=c(0,0,0,0))



#### -------- DAPC by species --------------------------------------------------
## Set the number of PCA axis to retain
n_pca=5 

## Loop on each species
for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  print(paste0("Species n°", sp, "/20: ", species)) 

  ## Read SNPs data
  Genlight <- readRDS(paste0("Data/SNPs_reefish/", species, "/radiator_genlight_", species, ".RDS"))
  levels(Genlight$pop) <- c("MV","SC","MF","MY") # relevel populations

  ## DAPC computation
  set.seed(999) # Setting a seed for a consistent result
  dapc <- dapc(Genlight, pop = Genlight$pop, n.da = 6, n.pca = n_pca)
  
  ## Save DAPC result
  saveRDS(dapc, paste0("Intermediate/02_DAPC/DAPC_result_", n_pca, "PC_", species, ".RDS"))
  dapc <- readRDS(paste0("Intermediate/02_DAPC/DAPC_result_", n_pca, "PC_", species, ".RDS"))
  print.dapc(dapc) # decomment if needed
  summary.dapc(dapc)
  predict.dapc(dapc)

  ## Plot
  pdf(file = paste0("Figures/DAPC/DAPC_plot_", n_pca, "PC_", species, ".pdf"), width = 6, height = 8)
  layout(matrix(c(1,2), nrow = 2, ncol = 1), heights=c(2.4,0.8))
  par(oma=c(1,1,2,1))
  scatter(dapc, scree.pca = F, legend = T, clabel = F, 
                posi.leg = "bottomleft", col=myCol,)
  par(mar=c(1,4,1,0))
  compoplot(dapc, lab=NULL, col=myCol, border = NA, legend = F, show.lab=F, cex.lab = 0.8)
  mtext(Database$Genus_species[sp], outer=TRUE, cex=1.5, line=0.5)
  dev.off()

}


