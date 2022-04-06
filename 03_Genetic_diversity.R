## ---------------------------------------------------------------------------------------------------------------- ##
## Script name: 02_Genetic_diversity
## Codes for the paper: Species-genetic diversity correlations in tropical reef fishes
## Date: 28 July 2020
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## Aim: Measure alpha, beta, and gamma genetic diversity for each of the 99 randomizations
##      and then compute the mean and standard deviation (sd) across the 99 randomizations.
##        - global alpha genetic diversity: average of the mean Hs obtained between the four sites
##            (i.e. the mean genetic alpha diversity of the Western Indian Ocean by species)
##        - global gamma genetic diversity: global expected heterozygosity (Hs)
##        - global beta genetic diversity: global G”ST, a corrected version of Hedrick’s G’ST
##        - pairwise beta genetic diversity: pairwise G”ST
## ---------------------------------------------------------------------------------------------------------------- ##

#### -------- Setup ---------------------------------------------------------------------------------------------
source("00_Setup.R")

list_sp <- list.files("Data/SNPs_ReeFISH/")
list_pop <- c("MF",	"MV",	"MY",	"SC")


#### -------- Compute genetic diversity on the 999 resampling ---------------------------------------------------------------------------------------------

for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  cat("##### Data loading for species:", species, ", species ", sp, "/", length(list_sp), "#####","\n")
  SNP_resampled <- readRDS(paste0("Intermediate/01_SNP_resampled/Data_999resampling_4479SNP_10ind_", species, ".rds"))
  
  cat("##### Basic stats computation and giaggotti framework application for species:", species, "#####","\n")
  GenetDiv999 <- mclapply(1:999, mc.cores=20, function(i){
    
    cat("Basic stats computation:",i,"\n")
    BS <- basic.stats(SNP_resampled[[i]][[3]],diploid=TRUE,digits=4)
    Bs_o <- BS$overall
    
    Gst <- (Bs_o["Ht"]-Bs_o["Hs"])/Bs_o["Ht"];k=4 # ATTENTION ICI on devrait avoir k=3 pour certaines espèces...
    Gst_Max <- (k-1)*(1-Bs_o["Hs"])/ (k-1+Bs_o["Hs"])
    
    GstP <- Gst/Gst_Max
    GstPP <- k*(Bs_o["Ht"]-Bs_o["Hs"])/ ((k*Bs_o["Ht"]-Bs_o["Hs"])*(1-Bs_o["Hs"]))
    
    cat("Jost multiplicatif framework application:",i,"\n")
    Hst <- (Bs_o["Ht"]-Bs_o["Hs"])/(1-Bs_o["Hs"])
    Jost <- c(Jt=1/(1-Bs_o["Ht"]),Js=1/(1-Bs_o["Hs"]),Jst=1/(1-Hst))
    
    list(c(BS$overall,Jost,Gst=Gst,GstP=GstP,GstPP=GstPP),BS[-7])
    
  })
  cat("##### Results compeltion for :", species, "#####","\n")
  saveRDS(GenetDiv999, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_999_",species,".rds"))
  
  GenetDivMean <- apply(do.call(rbind,lapply(GenetDiv999, function(e){e[[1]]})),2,mean,na.rm=T)
  GenetDivSd <- apply(do.call(rbind,lapply(GenetDiv999, function(e){e[[1]]})),2,sd,na.rm=T)
  saveRDS(GenetDivMean, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_mean_", species, ".rds"))
  saveRDS(GenetDivSd, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_sd_", species, ".rds"))
  
  GenetDivHedrick999 <- mclapply(1:2, mc.cores=1, function(i){
    Mygenind <- gl2gi(SNP_resampled[[i]][[1]])
    Gst_Hedrick(Mygenind)
  })
  saveRDS(GenetDivHedrick999, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_GstPPHedrick_999_", species, ".rds"))
  # rm(SNP_resampled)
}



#### -------- Alpha computation by site ---------------------------------------------------------------------------------------------
Alpha_list <- list()
for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  
  ## PREVIOUS METHOD
  cat("##### Alpha population values:", species, ", species ", sp, "/", length(list_sp), "#####","\n")
  GenetDiv999 <- readRDS(paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_999_", species, ".rds"))
  Hs_999 <- do.call(rbind,lapply(1:999,function(x){cat(x,"\n");apply(GenetDiv999[[x]][[2]]$Hs,2,mean,na.rm=T)}))
  Alpha_list[[species]]$Hs_mean <- apply(Hs_999, 2, mean)
  Alpha_list[[species]]$Hs_sd <- apply(Hs_999, 2, sd)

  Js_999 <- 1/(1-Hs_999)
  Alpha_list[[species]]$Js_mean <- apply(Js_999, 2, mean)
  Alpha_list[[species]]$Js_sd <- apply(Js_999, 2, sd)
  
  # ## New method
  # Alpha_list[[species]]$Hs <- Hs(SNP_resampled[[1]][[2]])
  
} # end of i

saveRDS(Alpha_list, file="Intermediate/GD_alpha_population.RDS")
# saveRDS(Alpha_list, file="Intermediate/Population_Alpha_adegenet.RDS")
# Alpha_list <- readRDS(file="Intermediate/GD_alpha_population.RDS")



#### -------- Beta calculation pairwise ---------------------------------------------------------------------------------------------
Beta_list <- list()
# Beta_list <- readRDS("Intermediate/DistMat_GstPP.RDS")
for (sp in 1:length(list_sp)){
  species <- list_sp[sp]
  cat("##### Beta pairwise values:", species,", species ", sp, "/",length(list_sp),"#####","\n")
  SNP_resampled <- readRDS(paste0("Intermediate/01_SNP_resampled/Data_999resampling_4479SNP_10ind_", species, ".rds"))
  
  # GstPP 99 computations
  GstPP_99 <- mclapply(1:99, mc.cores=20, function(i){
    Mygenind <- gl2gi(SNP_resampled[[i]][[1]])
    distGstPP <- as.matrix(pairwise_Gst_Hedrick(Mygenind))
    colnames(distGstPP) <- rownames(distGstPP) <- levels(SNP_resampled[[i]][[3]][,"pop"])
    distFst
  })
  saveRDS(GstPP_99, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_DistMatGstPP99_", species, ".rds"))
  
  # GstPP mean & sd
  npop <- nrow(GstPP_99[[1]])
  GstPP_mean <- apply(array(unlist(GstPP_99), c(npop, npop, length(GstPP_99))), c(1,2), mean)
  GstPP_sd <- apply(array(unlist(GstPP_99), c(npop, npop, length(GstPP_99))), c(1,2), sd)
  colnames(GstPP_mean) <- colnames(GstPP_sd) <- colnames(GstPP_99[[1]])
  rownames(GstPP_mean) <- rownames(GstPP_sd) <- rownames(GstPP_99[[1]])
  GstPP_mean <- GstPP_mean[order(rownames(GstPP_mean)), order(colnames(GstPP_mean))]
  GstPP_sd <- GstPP_sd[order(rownames(GstPP_sd)), order(colnames(GstPP_sd))]
  Beta_list[[species]]$GstPP_mean <- as.dist(GstPP_mean)
  Beta_list[[species]]$GstPP_sd <- as.dist(GstPP_sd)
  
  # Fst 99 computations
  Fst_99 <- mclapply(1:99, mc.cores=20, function(i){
    distFst <- as.matrix(genet.dist(SNP_resampled[[i]][[3]], diploid=TRUE, method="Nei87"))
    colnames(distFst) <- rownames(distFst) <- levels(SNP_resampled[[i]][[3]][,"pop"])
    distFst
  })
  saveRDS(Fst_99, paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_DistMatFst99_", species, ".rds"))
  
  # Fst mean & sd
  Fst_mean <- apply(array(unlist(Fst_99), c(npop, npop, length(Fst_99))), c(1,2), mean)
  Fst_sd <- apply(array(unlist(Fst_99), c(npop, npop, length(Fst_99))), c(1,2), sd)
  colnames(Fst_mean) <- colnames(Fst_sd) <- colnames(Fst_99[[1]])
  rownames(Fst_mean) <- rownames(Fst_sd) <- rownames(Fst_99[[1]])
  Fst_mean <- Fst_mean[order(rownames(Fst_mean)), order(colnames(Fst_mean))]
  Fst_sd <- Fst_sd[order(rownames(Fst_sd)), order(colnames(Fst_sd))]
  Beta_list[[species]]$Fst_mean <- as.dist(Fst_mean)
  Beta_list[[species]]$Fst_sd <- as.dist(Fst_sd)
  
  rm(GstPP_99, Fst_99, GstPP_mean, GstPP_sd, Fst_mean, Fst_sd, SNP_resampled)
  gc()
}

saveRDS(Beta_list, "Intermediate/GD_beta_DistMat.RDS")




#### -------- Genetic diversity dataset compilation --------------------------------------------------------

## Multi-site/regional
F2load <- list.files("Intermediate/03_Genetic_diversity/", pattern = "Genetic_diversity_mean", full.names = T)
Genet_global <- data.frame(do.call(rbind,lapply(F2load, function(i){readRDS(i)})))
rownames(Genet_global) <- gsub(F2load, pattern = "Intermediate/03_Genetic_diversity/Genetic_diversity_mean_|.rds", replacement = "")
Genet_global$Gen_s <- rownames(Genet_global)

Traits <- read.csv("Data/Dataset_species_traits.csv")
Genet_global <- merge(Traits, Genet_global, by = "Gen_s", all = T)

rownames(Genet_global) <- Genet_global$Gen_s
colnames(Genet_global) <- gsub(colnames(Genet_global), pattern = ".Ht|.Hs", replacement = "")
Genet_global["Car_m",c("Dst","Dstp","Fst","Fstp","Dest","Gst","GstP", "GstPP")] <- 0
Genet_global["Car_m","Hs"] <- Genet_global["Car_m","Ht"]
Genet_global["Car_m","Jst"] <- 1

for (species in list_sp){
  print(species)
  GenetDivHedrick999 <- readRDS(paste0("Intermediate/03_Genetic_diversity/Genetic_diversity_GstPPHedrick_999_",species,".rds"))
  Genet_global[species, "GstPPHedrick"] <- apply(do.call(rbind,lapply(GenetDivHedrick999, function(e){e$global})),2,mean,na.rm=T)
}



# Population alpha
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
  
  Genet_global <- merge(Genet_global, Genet_population, by = "Gen_s")
}


## Pairwise beta
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
  
  ## Merge into one data frame
  Genet_global <- merge(Genet_global, Genet_pairwise, by = "Gen_s")
}

## Write whole dataset
write.csv(Genet_global, "Intermediate/Table_Genetic_Diversity_full.csv", row.names=F)


