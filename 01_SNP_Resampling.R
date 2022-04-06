## ---------------------------------------------------------------------------------------------------------------- ##
## Script name: 01_SNPs_resampling
## Codes for the paper: Species-genetic diversity correlations in tropical reef fishes
## Date: 28 July 2020
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## Aim: To account for differences in site sampling success across species, 
##      the number of individuals used for genetic diversity measurements is randomly sampled 
##      to a maximum number of 10 individuals per sampling site 99 times (median and tradeoff
##      value of the overall sampling).
##      As genomic diversity is also directly related to the number of SNPs, for each of the
##      99 samplings, the filtered SNPs is randomly down-sampled to the lowest common number 
##      of SNPs found across all species (4470 SNPs).
##      This 99x data is then stored in a genlight list for each species.
## ---------------------------------------------------------------------------------------------------------------- ##

#### -------- Setup ---------------------------------------------------------------------------------------------
source("00_Setup.R")

list_sp <- list.files("Data/SNPs_ReeFISH/")


#### -------- SNPs and individuals resampling ------------------------------------------------------------------------
for (species in list_sp){
  BaseGenlight <- readRDS(paste0("Data/SNPs_ReeFISH/", species, "/radiator_genlight_", species, ".RDS"))
  
  cat("##### Data resampling and giaggotti framework application for species:", species ,"#####","\n")
  SNP_resampled <- mclapply(1:999, mc.cores=20, function(i){
    cat("i=", i, "\n")
    new_sample_list <- list()
    for (pop in levels(BaseGenlight$pop)) {
      # if more than 10 individuals in a population, sample 10 individuals names
      if(length(which(BaseGenlight@pop == pop)) > 10) { 
        sample_new <- sample(grep(pop, BaseGenlight$ind.names), 10, replace=F) }
      # if less than 10 individuals in a population, just collect all individual names
      if(length(which(BaseGenlight@pop == pop)) <= 10) {
        sample_new <- grep(pop, BaseGenlight$ind.names) }
      new_sample_list <- append(new_sample_list, sample_new)
    }
    names_to_keep <- BaseGenlight$ind.names[unlist(new_sample_list)]
    Genlight <- gl.keep.ind(BaseGenlight, ind.list = names_to_keep)
    
    Genlight <- Genlight[,sample(nLoc(Genlight), 4479, replace=F)]
    Genpop <- genind2genpop(gl2gi(Genlight,v=0),quiet=T) # conversion
    Hierfstat <- genind2hierfstat(dat=gl2gi(Genlight),pop=NULL)
    list(Genlight, Genpop, Hierfstat)
  })
  saveRDS(SNP_resampled, paste0("Intermediate/01_SNP_resampled/Data_999resampling_4479SNP_10ind_", species, ".rds")) 
  rm(SNP_resampled)
}

