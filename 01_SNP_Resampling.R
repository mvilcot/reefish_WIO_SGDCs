## -------------------------------------------------------------------------- ##
## Codes for the paper: Spatial genetic differentiation correlates with species 
## assemblage turnover across tropical reef fish lineages - Vilcot et al. 2023
## Script name: 01_SNPs_resampling
## Date: 02/01/2023
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## -------------------------------------------------------------------------- ##

#### -------- Setup ------------------------------------------------------------
source("00_Setup.R")

list_sp <- list.files("Data/SNPs_reefish/")


#### -------- Resampling of individuals and SNPs -------------------------------
for (species in list_sp){
  ## Read SNP dataset
  cat("##### Data resampling species:", species ,"#####","\n")
  BaseGenlight <- readRDS(paste0("Data/SNPs_reefish/", species, "/radiator_genlight_", species, ".RDS"))
  
  ## Individuals resampling, 999 times
  SNP_resampled <- mclapply(1:999, mc.cores=20, function(i){
    cat("i=", i, "\n")
    new_sample_list <- list()
    for (pop in levels(BaseGenlight$pop)) {
      # if more than 10 individuals in a population, sample 10 individuals
      if(length(which(BaseGenlight@pop == pop)) > 10) { 
        sample_new <- sample(grep(pop, BaseGenlight$ind.names), 10, replace=F) }
      
      # if less than 10 individuals in a population, just keep all individuals
      if(length(which(BaseGenlight@pop == pop)) <= 10) {
        sample_new <- grep(pop, BaseGenlight$ind.names) }
      new_sample_list <- append(new_sample_list, sample_new)
    }
    names_to_keep <- BaseGenlight$ind.names[unlist(new_sample_list)]
    Genlight <- gl.keep.ind(BaseGenlight, ind.list = names_to_keep)
    
    ## SNPs resampling
    Genlight <- Genlight[,sample(nLoc(Genlight), 4479, replace=F)]
    
    ## Convert to other formats
    Genpop <- genind2genpop(gl2gi(Genlight,v=0),quiet=T) 
    Hierfstat <- genind2hierfstat(dat=gl2gi(Genlight),pop=NULL) 
    list(Genlight, Genpop, Hierfstat) # list of different formats
  })
  
  ## Save 999 resampled dataset for each species
  saveRDS(SNP_resampled, 
          paste0("Intermediate/01_SNP_resampled/Data_999resampling_4479SNP_10ind_", species, ".rds")) 
  rm(SNP_resampled)
}

