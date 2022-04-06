## ---------------------------------------------------------------------------------------------------------------- ##
## Script name: 04_Taxonomic_diversity
## Codes for the paper: Species-genetic diversity correlations in tropical reef fishes
## Date: 28 July 2020
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## Aim: From a worldwide 2,446 fish species occurrence data of the 12 families, data is restricted to a 1° distance 
##      from the reefs of the four sampling areas (gbuffer). It allows to produce a species presence/absence matrix 
##      by samplng area.
##      From this dataset is computed:
##        - population alpha species diversity: the number of species present on each site
##        - global alpha species diversity: average of population alpha species diversity
##        - global gamma species diversity: the number of species present on the 4 sites
##        - pairwise beta species diversity: Jaccard total dissimilarity index (Bjac) 
##                                           & the turnover component of Jaccard dissimilarity index (Bjtu)
##        - global (multisite) beta species diversity: Bjac & Bjtu 
## ---------------------------------------------------------------------------------------------------------------- ##

#### -------- Setup ---------------------------------------------------------------------------------------------
source("00_Setup.R")

Database <- read.csv("Data/Dataset_species_traits.csv")
list_species <- as.character(Database$Genus_species) # list of the 20 studied species
list_families <- c(sort(as.character(unique(Database$Family))), "12Families")
list_pop <- c("MF",	"MV",	"MY",	"SC")
myCol <- c("#f8766d","#76ee00","#00bfc4","#febe01")


#### -------- Buffer shapefile around 4 islands ---------------------------------------------------------------------
MF_reef <- st_read("Data/ShapeFiles/Mafia_reef.shp")
MF_reef_sp <- as_Spatial(MF_reef$geom) #Transform to spatial object
MF_reef_spPlanar <- spTransform(MF_reef_sp, CRS( "+init=epsg:4326" ))
BufferSp_MF <- gBuffer(MF_reef_spPlanar, width = 1)
Buffer_MF <- st_as_sf(BufferSp_MF) #Transform back to sf
saveRDS(Buffer_MF, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MF.RDS")

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

MY_reef <- st_read("Data/ShapeFiles/Mayotte_inter.shp")
MY_reef_sp <- as_Spatial(MY_reef$geom) #Transform to spatial object
MY_reef_spPlanar <- spTransform(MY_reef_sp, CRS( "+init=epsg:4326" ))
BufferSp_MY <- gBuffer(MY_reef_spPlanar, width = 1)
Buffer_MY <- st_as_sf(BufferSp_MY) #Transform back to sf
saveRDS(Buffer_MY, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MY.RDS")

SC_reef <- st_read("Data/ShapeFiles/Seychelle_reef.shp")
SC_reef_sp <- as_Spatial(SC_reef$geom) #Transform to spatial object
SC_reef_spPlanar <- spTransform(SC_reef_sp, CRS( "+init=epsg:4326" ))
BufferSp_SC <- gBuffer(SC_reef_spPlanar, width = 1)
Buffer_SC <- st_as_sf(BufferSp_SC) #Transform back to sf
saveRDS(Buffer_SC, "Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_SC.RDS")



#### -------- Species Presence by island -------------------------------------------------------------
Buffer_MF <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MF.RDS")
Buffer_MV <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MV.RDS")
Buffer_MY <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_MY.RDS")
Buffer_SC <- readRDS("Intermediate/04_ShapeFiles_buffered/Buffer_shapefile_SC.RDS")

nb_species_glob <- list()
for (i in 1:length(list_families)){
  family <- list_families[i]
  print(paste0("Family n°", i, "/12: ", family)) 
  PAmatrix_global <- readRDS(paste0("Data/Presence_matrix/Mat_", family, ".RDS"))
  nb_species_glob[family] <- length(PAmatrix_global) - 2
  PAmatrix_globalsf <- st_as_sf(PAmatrix_global, coords = c("Longitude", "Latitude"), crs=4326, remove=FALSE)
  
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_MF) # test si intersection (renvoie une liste)
  test_inter2 <- sapply(test_inter, length) # test_inter2 = nb d'element intersecte (normalemt 0 ou 1)
  inter_MFbuff <- PAmatrix_globalsf[test_inter2>0,] # filtrer sf_data : eliminer test_inter2=0 
  inter_MFbuff_df <- inter_MFbuff
  st_geometry(inter_MFbuff_df) <- NULL #convert sf to dataframe
  
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_MV) # test si intersection (renvoie une liste)
  test_inter2 <- sapply(test_inter, length) # test_inter2 = nb d'element intersecte (normalemt 0 ou 1)
  inter_MVbuff <- PAmatrix_globalsf[test_inter2>0,] # filtrer sf_data : eliminer test_inter2=0 
  inter_MVbuff_df <- inter_MVbuff
  st_geometry(inter_MVbuff_df) <- NULL #convert sf to dataframe
  
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_MY) # test si intersection (renvoie une liste)
  test_inter2 <- sapply(test_inter, length) # test_inter2 = nb d'element intersecte (normalemt 0 ou 1)
  inter_MYbuff <- PAmatrix_globalsf[test_inter2>0,] # filtrer sf_data : eliminer test_inter2=0 
  inter_MYbuff_df <- inter_MYbuff
  st_geometry(inter_MYbuff_df) <- NULL #convert sf to dataframe
  
  test_inter <- st_intersects(PAmatrix_globalsf, Buffer_SC) # test si intersection (renvoie une liste)
  test_inter2 <- sapply(test_inter, length) # test_inter2 = nb d'element intersecte (normalemt 0 ou 1)
  inter_SCbuff <- PAmatrix_globalsf[test_inter2>0,] # filtrer sf_data : eliminer test_inter2=0 
  inter_SCbuff_df <- inter_SCbuff
  st_geometry(inter_SCbuff_df) <- NULL #convert sf to dataframe
  
  PAmatrix_FAMbuff <- as.data.frame(setNames(replicate(length(PAmatrix_global), numeric(0), simplify = F), colnames(PAmatrix_global)))
  PAmatrix_FAMbuff["MF", ] <- apply(inter_MFbuff_df, 2, max)
  PAmatrix_FAMbuff["MV", ] <- apply(inter_MVbuff_df, 2, max)
  PAmatrix_FAMbuff["MY", ] <- apply(inter_MYbuff_df, 2, max)
  PAmatrix_FAMbuff["SC", ] <- apply(inter_SCbuff_df, 2, max)
  
  saveRDS(PAmatrix_FAMbuff[,-c(1,2)], paste0("Intermediate/04_PAmatrix_by_site/Species_PAmatrix_", family, ".RDS"))
}


#### -------- Venn Diagram & table nb species ------------------------------------------------------------------------
### !!!!!!!!!!!xxxxxxxxxxxxxxxxxx ICI, redondant avec calcul alpha et gamma div !!! ##############
plot_list <- list()
nb_species <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(nb_species) <- c("MF", "MV", "MY", "SC", "regional")
for (i in c(1:length(list_families))){
  family <- list_families[i]
  print(paste0("Family n°", i, "/12: ", family)) 
  PAmatrix_FAMbuff <- readRDS(paste0("Intermediate/04_PAmatrix_by_site/Species_PAmatrix_", family, ".RDS"))
  PAmatrix <- as.data.frame(t(PAmatrix_FAMbuff))
  PAlist <- list(MF = row.names(PAmatrix[PAmatrix$MF == 1,]),
                  MV = row.names(PAmatrix[PAmatrix$MV == 1,]),
                  MY = row.names(PAmatrix[PAmatrix$MY == 1,]),
                  SC = row.names(PAmatrix[PAmatrix$SC == 1,]))
  test <- unique(c(PAlist$MF, PAlist$MV, PAlist$MY, PAlist$SC))
  PAlist$regional <- test
  nb_species[family,] <- lengths(PAlist)
  
  ## ggvenn
  # plot_list[[family]] <-
  #   ggvenn(PAlist, fill_color = myCol, stroke_size = 0.3, set_name_size = 6) +
  #   ggtitle(family) +
  #   annotate('text', x=-2, y=2, hjust=0, vjust=0, size=4,
  #            label=paste("regional SR =", as.character(length(PAlist$regional))))
  
  ## or ggvenndiagram
  plot_list[[family]] <-
  ggVennDiagram(PAlist[1:4], label_alpha = 0, set_color = myCol) +
    scale_fill_gradient(low="white", high = "orange") +
    ggtitle(family) +
    annotate('text', x=0.1, y=0.9, hjust=0, vjust=0, size=4,
             label=paste("regional SR =", as.character(length(PAlist$regional))))
}

g <- arrangeGrob(grobs=plot_list, ncol=2)
ggsave("Results/Figures/SR_Venn_diagram.pdf", g, width = 10, height = 46, limitsize = F)

# nb_species$global <- unlist(nb_species_glob)
write.csv(nb_species[,c(5,1:4)], "Results/TableS2_nb_species_by_site.csv")


#### -------- Loop on family ----------------------------------------------------------------------------------------
combin <- combn(list_pop, 2)
Combinaisons <- paste(combin[1,], combin[2,], sep=".")
Taxo_global <- data.frame()
Taxo_pairwise_jac <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), Combinaisons))
Taxo_pairwise_jtu <- as.data.frame(setNames(replicate(6,numeric(0), simplify = F), Combinaisons))
Taxo_population <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F), c("MF", "MV", "MY", "SC")))

Jaccard_multi <- list()
Jaccard_pair <- list()

for (family in list_families){
  PAmatrix_FAM <- readRDS(paste0("Intermediate/04_PAmatrix_by_site/Species_PAmatrix_", family, ".RDS"))

  #### -------- Taxonomic diversity computations --------------------------------------------------------------------
  ## Gamma diversity
  PAmatrix_FAMsum <- data.frame(t(apply(PAmatrix_FAM, 2, max)))
  Taxo_global[family, "SR_Gamma"] <- rowSums(PAmatrix_FAMsum)
  
  ## Alpha diversity
  alphaDiv <- rowSums(PAmatrix_FAM)
  Taxo_population[family,] <- alphaDiv
  Taxo_global[family, "SR_Alpha"] <- mean(alphaDiv)
  
  ## Beta diversity
  Jaccard_multi[[family]] <- beta.multi(PAmatrix_FAM, index.family="jaccard")
  Taxo_global[family,"JAC_multi"] <- Jaccard_multi[[family]]$beta.JAC
  Taxo_global[family,"JTU_multi"] <- Jaccard_multi[[family]]$beta.JTU
  
  Jaccard_pair[[family]] <- beta.pair(PAmatrix_FAM, index.family="jaccard")
  Taxo_global[family,"JAC_pair"] <- mean(Jaccard_pair[[family]]$beta.jac)
  Taxo_global[family,"JTU_pair"] <- mean(Jaccard_pair[[family]]$beta.jtu)

  
  #### -------- Pairwise beta diversity table -----------------------------------------------------------------------
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

saveRDS(Jaccard_pair, "Intermediate/Jaccard_DistMat.RDS")




#### -------- One table with all taxnonomic diversity values --------------------------------------------------------
colnames(Taxo_population) <- paste0("SR_", colnames(Taxo_population))
colnames(Taxo_pairwise_jac) <- paste0("JAC_", colnames(Taxo_pairwise_jac))
colnames(Taxo_pairwise_jtu) <- paste0("JTU_", colnames(Taxo_pairwise_jtu))
Taxo_global$Family <- rownames(Taxo_global)
Taxo_population$Family <- rownames(Taxo_population)
Taxo_pairwise_jac$Family <- rownames(Taxo_pairwise_jac)
Taxo_pairwise_jtu$Family <- rownames(Taxo_pairwise_jtu)

Data_taxo <- merge(Taxo_global, Taxo_population, by = "Family", sort = FALSE)
Data_taxo <- merge(Data_taxo, Taxo_pairwise_jac, by = "Family", sort = FALSE)
Data_taxo <- merge(Data_taxo, Taxo_pairwise_jtu, by = "Family", sort = FALSE)

write.csv(Data_taxo, file = "Intermediate/Table_Taxonomic_Diversity_full.csv", row.names = FALSE)



#### -------- Combined Genetic and Taxonomic tables ------------------------------------------------------------------
Data_taxo <- read.csv(file = "Intermediate/Table_Taxonomic_Diversity_full.csv")
Data_genet <- read.csv(file = "Intermediate/Table_Genetic_Diversity_full.csv")
Data_full <- merge(Data_genet, Data_taxo[-nrow(Data_taxo),])
Data_full <- Data_full[order(Data_full$Gen_s),]
Data_full <- Data_full[,c(2,3,1,4:ncol(Data_full))]

Data_taxoCommunity <- data.frame(Data_taxo[Data_taxo$Family=="12Families",c("SR_Gamma", "SR_Alpha", "JAC_multi", "JTU_multi",
                                                                              "SR_MF", "SR_MV", "SR_MY", "SR_SC")])
colnames(Data_taxoCommunity) <- paste0("12Families_", colnames(Data_taxoCommunity))
Data_full <- cbind(Data_full, as.data.frame(lapply(Data_taxoCommunity, rep, 20)))

write.csv(Data_full, "Results/Table_diversity_full.csv", row.names = F)
