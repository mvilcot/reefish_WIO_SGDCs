## -------------------------------------------------------------------------- ##
## Codes for the paper: Spatial genetic differentiation correlates with species 
## assemblage turnover across tropical reef fish lineages - Vilcot et al. 2023
## Script name: 02_Genetic_diversity
## Date: 02/01/2023
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## -------------------------------------------------------------------------- ##

#### -------- Setup ------------------------------------------------------------
source("00_Setup.R")

## List of species 
list_sp <- list.files("Data/SNPs_reefish/")

## List of sampling sites
list_pop <- c("MF",	"MV",	"MY",	"SC")


#### -------- Compute genetic diversity for the 999 resampled SNPs -------------

## Loop on species
for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  
  ## Read resampled SNPs data
  cat("Data loading for species:", species, ", species ", sp, "/", length(list_sp), "\n")
  SNP_resampled <- readRDS(paste0("Intermediate/01_SNP_resampled/Data_999resampling_4479SNP_10ind_", species, ".rds"))
  
  ## Basic stats for the 999 resamplings
  cat("Basic stats computation for species:", species, "\n")
  GenetDiv999 <- mclapply(1:999, mc.cores=20, function(i){ # Parallel computation on 20 cores
    cat("Basic stats computation:",i,"\n")
    BS <- basic.stats(SNP_resampled[[i]][[3]], diploid=TRUE, digits=4) # [[3]] for Hierfstat
    list(c(BS$overall),BS[-7])
  })
  saveRDS(GenetDiv999, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_999_",species,".rds"))
  
  ## Mean and sd over the 999 resamplings
  GenetDivMean <- apply(do.call(rbind,lapply(GenetDiv999, function(e){e[[1]]})),2,mean,na.rm=T)
  GenetDivSd <- apply(do.call(rbind,lapply(GenetDiv999, function(e){e[[1]]})),2,sd,na.rm=T)
  saveRDS(GenetDivMean, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_mean_", species, ".rds"))
  saveRDS(GenetDivSd, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_sd_", species, ".rds"))
  
  ## Hedrick Gst for the 999 resamplings
  cat("Hedrick G''st computation for species:", species, "\n")
  GenetDivHedrick999 <- mclapply(1:999, mc.cores=20, function(i){ # Parallel computation on 20 cores
    Mygenind <- gl2gi(SNP_resampled[[i]][[1]]) # [[2]] for Genlight
    Gst_Hedrick(Mygenind)
  })
  saveRDS(GenetDivHedrick999, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_GstPPHedrick_999_", species, ".rds"))
  # rm(SNP_resampled)
}



#### -------- Population Alpha diversity (= by sampling site) ------------------

Alpha_list <- list()

for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  
  ## Read Basic stats computations over 999 resamplings
  cat("Alpha population values:", species, ", species ", sp, "/", length(list_sp),"\n")
  GenetDiv999 <- readRDS(paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_999_", species, ".rds"))
  
  ## mean an sd of Hs for population
  Hs_999 <- do.call(rbind, lapply(1:999, function(x) {
      cat(x, "\n")
      apply(GenetDiv999[[x]][[2]]$Hs, 2, mean, na.rm = T)
    }))
  Alpha_list[[species]]$Hs_mean <- apply(Hs_999, 2, mean)
  Alpha_list[[species]]$Hs_sd <- apply(Hs_999, 2, sd)
}

saveRDS(Alpha_list, file="Intermediate/GD_alpha_population.RDS")




#### -------- Pairwise Beta diversity ------------------------------------------

Beta_list <- list()

for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  
  ## Read resampled SNPs data
  cat("Beta pairwise values:", species,", species ", sp, "/",length(list_sp),"\n")
  SNP_resampled <- readRDS(paste0("Intermediate/01_SNP_resampled/Data_999resampling_4479SNP_10ind_", species, ".rds"))
  
  ## Pairwise Gst'' 999 computations
  GstPP_999 <- mclapply(1:999, mc.cores=20, function(i){
    Mygenind <- gl2gi(SNP_resampled[[i]][[1]]) # [[1]] for genlight
    distGstPP <- as.matrix(pairwise_Gst_Hedrick(Mygenind)) # transform dist object to matrix
    colnames(distGstPP) <- rownames(distGstPP) <- levels(SNP_resampled[[i]][[3]][,"pop"]) # Name rows and columns as pops
    distGstPP
  })
  saveRDS(GstPP_999, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_DistMatGstPP999_", species, ".rds"))
  
  ## Pairwise Gst'' mean & sd
  npop <- nrow(GstPP_999[[1]])
  GstPP_mean <- apply(array(unlist(GstPP_999), c(npop, npop, length(GstPP_999))), c(1,2), mean)
  GstPP_sd <- apply(array(unlist(GstPP_999), c(npop, npop, length(GstPP_999))), c(1,2), sd)
  colnames(GstPP_mean) <- colnames(GstPP_sd) <- colnames(GstPP_999[[1]])
  rownames(GstPP_mean) <- rownames(GstPP_sd) <- rownames(GstPP_999[[1]])
  GstPP_mean <- GstPP_mean[order(rownames(GstPP_mean)), order(colnames(GstPP_mean))]
  GstPP_sd <- GstPP_sd[order(rownames(GstPP_sd)), order(colnames(GstPP_sd))]
  Beta_list[[species]]$GstPP_mean <- as.dist(GstPP_mean)
  Beta_list[[species]]$GstPP_sd <- as.dist(GstPP_sd)
  
  
  ## Pairwise Fst 99 computations
  Fst_999 <- mclapply(1:99, mc.cores=20, function(i){
    distFst <- as.matrix(genet.dist(SNP_resampled[[i]][[3]], diploid=TRUE, method="Nei87"))
    colnames(distFst) <- rownames(distFst) <- levels(SNP_resampled[[i]][[3]][,"pop"])
    distFst
  })
  saveRDS(Fst_999, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_DistMatFst999_", species, ".rds"))
  
  ## Pairwise Fst mean & sd
  Fst_mean <- apply(array(unlist(Fst_999), c(npop, npop, length(Fst_999))), c(1,2), mean)
  Fst_sd <- apply(array(unlist(Fst_999), c(npop, npop, length(Fst_999))), c(1,2), sd)
  colnames(Fst_mean) <- colnames(Fst_sd) <- colnames(Fst_999[[1]])
  rownames(Fst_mean) <- rownames(Fst_sd) <- rownames(Fst_999[[1]])
  Fst_mean <- Fst_mean[order(rownames(Fst_mean)), order(colnames(Fst_mean))]
  Fst_sd <- Fst_sd[order(rownames(Fst_sd)), order(colnames(Fst_sd))]
  Beta_list[[species]]$Fst_mean <- as.dist(Fst_mean)
  Beta_list[[species]]$Fst_sd <- as.dist(Fst_sd)
  
  rm(GstPP_99, Fst_999, GstPP_mean, GstPP_sd, Fst_mean, Fst_sd, SNP_resampled)
  
  gc()
}

saveRDS(Beta_list, "Intermediate/GD_beta_DistMat.RDS")




#### -------- Merge all genetic diversity indices into one table ---------------

### ----- Multi-site values ----- ###
## Read basic stats values
F2load <- list.files("Intermediate/03_Genetic_diversity/", pattern = "Genetic_diversity_mean", full.names = T)
Genet_multi <- data.frame(do.call(rbind,lapply(F2load, function(i){readRDS(i)})))

## Rename columns from species name
rownames(Genet_multi) <- gsub(F2load, pattern = "Intermediate/03_Genetic_diversity/Genetic_diversity_mean_|.rds", replacement = "")
Genet_multi$Gen_s <- rownames(Genet_multi)

## Add traits indices
Traits <- read.csv("Data/Dataset_species_traits.csv")
Genet_multi <- merge(Traits, Genet_multi, by = "Gen_s", all = T)
rownames(Genet_multi) <- Genet_multi$Gen_s

## Add Hedrick Gst values
for (species in list_sp){
  print(species)
  GenetDivHedrick999 <- readRDS(paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_GstPPHedrick_999_",species,".rds"))
  Genet_multi[species, "GstPPHedrick"] <- apply(do.call(rbind,lapply(GenetDivHedrick999, function(e){e$global})),2,mean,na.rm=T)
}



### ----- Population Alpha values ----- ###
Alpha_list <- readRDS(file="Intermediate/GD_alpha_population.RDS")
for (metric in c("Hs", "Js")){
  Genet_population <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(Genet_population) <- list_pop
  for (species in list_sp){
    pops <- names(Alpha_list[[species]][[paste0(metric, "_mean")]])
    for (pop in pops){
      Genet_population[species, pop] <- Alpha_list[[species]][[paste0(metric, "_mean")]][[pop]]
    }
  }
  colnames(Genet_population) <- paste0(metric, "_", colnames(Genet_population))
  Genet_population$Gen_s <- rownames(Genet_population)
  
  Genet_multi <- merge(Genet_multi, Genet_population, by = "Gen_s")
}


### ----- Pairwise Beta values ----- ###
Beta_list <- readRDS("Intermediate/GD_beta_DistMat.RDS")
for (metric in c("Fst", "GstPP")){
  Genet_pairwise <- data.frame(t(combn(list_pop, 2)))
  for (species in list_sp){
    Beta_mean <- as.matrix(Beta_list[[species]][[paste0(metric,"_mean")]])
    
    melt <- t(combn(colnames(Beta_mean), 2))
    melt <- data.frame(melt, species=Beta_mean[melt])
    colnames(melt)[3] <- species
    Genet_pairwise <- merge(Genet_pairwise, melt, all = T)
  }
  
  Genet_pairwise$sites <- paste0(Genet_pairwise$X1, '.', Genet_pairwise$X2)
  Genet_pairwise <- data.frame(Genet_pairwise, row.names = Genet_pairwise$sites)
  Genet_pairwise <- data.frame(t(Genet_pairwise[, 3:(ncol(Genet_pairwise)-1)]))
  colnames(Genet_pairwise) <- paste0(metric, "_", colnames(Genet_pairwise))
  Genet_pairwise$Gen_s <- rownames(Genet_pairwise)
  
  Diversity_genet <- merge(Genet_multi, Genet_pairwise, by = "Gen_s")
}

## Write whole genetic data set
write.csv(Diversity_genet, "Intermediate/Table_Genetic_Diversity_full.csv", row.names=F)


