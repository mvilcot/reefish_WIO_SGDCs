## -------------------------------------------------------------------------- ##
## Codes for the paper: Spatial genetic differentiation correlates with species 
## assemblage turnover across tropical reef fish lineages - Vilcot et al. 2023
## Script name: 04_Taxonomic_diversity
## Date: 02/01/2023
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## -------------------------------------------------------------------------- ##

#### -------- Setup ------------------------------------------------------------
source("00_Setup.R")

## Read data set 
Database <- read.csv("Data/Dataset_species_traits.csv")
list_species <- as.character(Database$Genus_species) # list of the 20 studied species
list_families <- sort(unique(Database$Family))
list_pop <- c("MF",	"MV",	"MY",	"SC")




#### -------- Buffer around each site shapefile --------------------------------
## Mafia Island
MF_reef <- st_read("Data/ShapeFiles/Mafia_reef.shp")
MF_reef_sp <- as_Spatial(MF_reef$geom) #Transform to spatial object
MF_reef_spPlanar <- spTransform(MF_reef_sp, CRS( "+init=epsg:4326" ))
BufferSp_MF <- gBuffer(MF_reef_spPlanar, width = 1)
Buffer_MF <- st_as_sf(BufferSp_MF) #Transform back to sf
saveRDS(Buffer_MF, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MF.RDS")

## Maldives
MV_inter <- st_read("Data/ShapeFiles/Maldive_inter.shp")
MV_interJSON <- geojson_json(MV_inter, geometry = "polygon", group = "group")
MV_interSP <- geojson_sp(MV_interJSON)
MV_interSP <- crop(MV_interSP, extent(72.55, 73.78, 3.33, 4.78)) # crop MV shape file to sampling locations
MV_inter <- st_as_sf(MV_interSP) #Transform back to sf the crop file
MV_interSIMP <- ms_simplify(MV_interSP, keep = 0.001) #simplify MV
MV_reef_spPlanar <- spTransform(MV_interSIMP, CRS( "+init=epsg:4326" ))
BufferSp_MV <- gBuffer(MV_reef_spPlanar, width = 1)
Buffer_MV <- st_as_sf(BufferSp_MV) #Transform back to sf
saveRDS(Buffer_MV, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MV.RDS")

## Mayotte
MY_reef <- st_read("Data/ShapeFiles/Mayotte_inter.shp")
MY_reef_sp <- as_Spatial(MY_reef$geom) #Transform to spatial object
MY_reef_spPlanar <- spTransform(MY_reef_sp, CRS( "+init=epsg:4326" ))
BufferSp_MY <- gBuffer(MY_reef_spPlanar, width = 1)
Buffer_MY <- st_as_sf(BufferSp_MY) #Transform back to sf
saveRDS(Buffer_MY, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MY.RDS")

## Seychelles
SC_reef <- st_read("Data/ShapeFiles/Seychelle_reef.shp")
SC_reef_sp <- as_Spatial(SC_reef$geom) #Transform to spatial object
SC_reef_spPlanar <- spTransform(SC_reef_sp, CRS( "+init=epsg:4326" ))
BufferSp_SC <- gBuffer(SC_reef_spPlanar, width = 1)
Buffer_SC <- st_as_sf(BufferSp_SC) #Transform back to sf
saveRDS(Buffer_SC, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_SC.RDS")



#### -------- Species Presence by sampling site --------------------------------
# ## Read sampling site buffers
# Buffer_MF <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MF.RDS")
# Buffer_MV <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MV.RDS")
# Buffer_MY <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MY.RDS")
# Buffer_SC <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_SC.RDS")

## use sf
sf::sf_use_s2(FALSE) 

## Create list of number of species in each family
nb_species_glob <- list()

## loop on the familes
for (i in 1:length(list_families)){
  family <- list_families[i]
  print(paste0("Family n°", i, "/12: ", family)) 
  
  ## Read global species presence data
  PAmatrix_global <- readRDS(paste0("Data/Presence_matrix/Mat_", family, ".RDS"))
  nb_species_glob[family] <- length(PAmatrix_global) - 2
  PAmatrix_globalsf <- st_as_sf(PAmatrix_global, coords = c("Longitude", "Latitude"), crs=4326, remove=FALSE)
  
  ### Restriction of presence data to each sampling site buffer
  ## Mafia Island
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_MF) 
  test_inter2 <- sapply(test_inter, length) 
  inter_MFbuff <- PAmatrix_globalsf[test_inter2>0,] 
  inter_MFbuff_df <- inter_MFbuff
  st_geometry(inter_MFbuff_df) <- NULL #convert sf to dataframe
  
  ## Maldives
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_MV) 
  test_inter2 <- sapply(test_inter, length)
  inter_MVbuff <- PAmatrix_globalsf[test_inter2>0,] 
  inter_MVbuff_df <- inter_MVbuff
  st_geometry(inter_MVbuff_df) <- NULL #convert sf to dataframe
  
  ## Mayotte
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_MY) 
  test_inter2 <- sapply(test_inter, length)
  inter_MYbuff <- PAmatrix_globalsf[test_inter2>0,] 
  inter_MYbuff_df <- inter_MYbuff
  st_geometry(inter_MYbuff_df) <- NULL #convert sf to dataframe
  
  ## Seychelles
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_SC) 
  test_inter2 <- sapply(test_inter, length) 
  inter_SCbuff <- PAmatrix_globalsf[test_inter2>0,] 
  inter_SCbuff_df <- inter_SCbuff
  st_geometry(inter_SCbuff_df) <- NULL #convert sf to dataframe
  
  ## Create the presence matrix by sampling site
  PAmatrix_FAMbuff <- as.data.frame(setNames(replicate(length(PAmatrix_global), numeric(0), simplify = F), colnames(PAmatrix_global)))
  PAmatrix_FAMbuff["MF", ] <- apply(inter_MFbuff_df, 2, max)
  PAmatrix_FAMbuff["MV", ] <- apply(inter_MVbuff_df, 2, max)
  PAmatrix_FAMbuff["MY", ] <- apply(inter_MYbuff_df, 2, max)
  PAmatrix_FAMbuff["SC", ] <- apply(inter_SCbuff_df, 2, max)
  
  ## Save 
  saveRDS(PAmatrix_FAMbuff[,-c(1,2)], paste0("Intermediate/04_PAmatrix_by_site/Species_PAmatrix_", family, "_sf.RDS"))
}



#### -------- Compute diversity indices ----------------------------------------

## Create pairwise data
combin <- combn(list_pop, 2)
Combinaisons <- paste(combin[1,], combin[2,], sep=".")

## Create empty data frames
Taxo_multi <- data.frame() # Multi-site values
Taxo_pairwise_jac <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), Combinaisons)) # Pairwise beta Jaccard 
Taxo_pairwise_jtu <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), Combinaisons)) # Pairwise beta Jaccard turnover
Taxo_population <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F), c("MF", "MV", "MY", "SC"))) # Population Alpha

## Create empty list for Jaccard output
Jaccard_multi <- list()
Jaccard_pair <- list()

## Loop on families
for (family in list_families){
  
  ## Read Presence data by site
  print(family)
  PAmatrix_FAM <- readRDS(paste0("Intermediate/04_PAmatrix_by_site/Species_PAmatrix_", family, "_sf.RDS"))

  
  ### ------ Taxonomic diversity computations ------ ### 
  ## Gamma diversity (= species richness on the 4 sites)
  PAmatrix_FAMsum <- data.frame(t(apply(PAmatrix_FAM, 2, max)))
  Taxo_multi[family, "SR_Gamma"] <- rowSums(PAmatrix_FAMsum)

  ## Alpha diversity (= species richness on each site)
  alphaDiv <- rowSums(PAmatrix_FAM)
  Taxo_population[family,] <- alphaDiv
  Taxo_multi[family, "SR_Alpha"] <- mean(alphaDiv)

  ## Beta diversity (multiple-site and pairwise)
  Jaccard_multi[[family]] <- beta.multi(PAmatrix_FAM, index.family="jaccard")
  Taxo_multi[family,"JAC_multi"] <- Jaccard_multi[[family]]$beta.JAC # Jaccard
  Taxo_multi[family,"JTU_multi"] <- Jaccard_multi[[family]]$beta.JTU # Jaccard turnover component

  Jaccard_pair[[family]] <- beta.pair(PAmatrix_FAM, index.family="jaccard")
  Taxo_multi[family,"JAC_pair"] <- mean(Jaccard_pair[[family]]$beta.jac) # Jaccard 
  Taxo_multi[family,"JTU_pair"] <- mean(Jaccard_pair[[family]]$beta.jtu) # Jaccard turnover component


  
  #### ------ Pairwise beta diversity table ------ ### 
  DistJAC <- as.matrix(Jaccard_pair[[family]][["beta.jac"]])
  Taxo_pairwise_jac[family,] <- c(DistJAC[combin[1,1], combin[2,1]],
                                       DistJAC[combin[1,2], combin[2,2]],
                                       DistJAC[combin[1,3], combin[2,3]],
                                       DistJAC[combin[1,4], combin[2,4]],
                                       DistJAC[combin[1,5], combin[2,5]],
                                       DistJAC[combin[1,6], combin[2,6]])
  
  DistJTU <- as.matrix(Jaccard_pair[[family]][["beta.jtu"]])
  Taxo_pairwise_jtu[family,] <- c(DistJTU[combin[1,1], combin[2,1]],
                                       DistJTU[combin[1,2], combin[2,2]],
                                       DistJTU[combin[1,3], combin[2,3]],
                                       DistJTU[combin[1,4], combin[2,4]],
                                       DistJTU[combin[1,5], combin[2,5]],
                                       DistJTU[combin[1,6], combin[2,6]])
  }

saveRDS(Jaccard_pair, "Intermediate/Jaccard_DistMat_sf.RDS")




#### -------- Merge all species diversity indices into one table ---------------
## Rename columns for Species Richness, JAC and JTU tables
colnames(Taxo_population) <- paste0("SR_", colnames(Taxo_population))
colnames(Taxo_pairwise_jac) <- paste0("JAC_", colnames(Taxo_pairwise_jac))
colnames(Taxo_pairwise_jtu) <- paste0("JTU_", colnames(Taxo_pairwise_jtu))

## Rename rows with family names before merging tables
Taxo_multi$Family <- rownames(Taxo_multi)
Taxo_population$Family <- rownames(Taxo_population)
Taxo_pairwise_jac$Family <- rownames(Taxo_pairwise_jac)
Taxo_pairwise_jtu$Family <- rownames(Taxo_pairwise_jtu)

## Merge all species diversity tables
Diversity_taxo <- merge(Taxo_multi, Taxo_population, by = "Family", sort = FALSE)
Diversity_taxo <- merge(Diversity_taxo, Taxo_pairwise_jac, by = "Family", sort = FALSE)
Diversity_taxo <- merge(Diversity_taxo, Taxo_pairwise_jtu, by = "Family", sort = FALSE)

## Save
write.csv(Diversity_taxo, "Intermediate/Table_Taxonomic_Diversity_full.csv", row.names = FALSE)




#### -------- Combine genetic and species tables -------------------------------
## Read tables
Diversity_taxo <- read.csv(file = "Intermediate/Table_Taxonomic_Diversity_full.csv") # Species diversity 
Diversity_genet <- read.csv(file = "Intermediate/Table_Genetic_Diversity_full.csv") # Genetic diversity

## Merge tables and order by species name
Diversity <- merge(Diversity_genet, Diversity_taxo)
Diversity <- Diversity[order(Diversity$Gen_s), ]

## Reorder columns
Diversity <- Diversity[,c(2,3,1,4:ncol(Diversity))]

## Save
write.csv(Diversity, "Results/Table_diversity_full.csv", row.names = F)







