## -------------------------------------------------------------------------- ##
## Codes for the paper: Spatial genetic differentiation correlates with species 
## assemblage turnover across tropical reef fish lineages - Vilcot et al. 2023
## Script name: 00_Setup
## Date: 02/01/2023
## Author: Maurine Vilcot
## Email: maurine.vilcot@ens-lyon.fr
## -------------------------------------------------------------------------- ##

#### -------- Create arborescence ----------------------------------------------

dir.create("Figures")
dir.create("Figures/DAPC")
dir.create("Intermediate")
dir.create("Intermediate/01_SNP_resampled")
dir.create("Intermediate/02_DAPC")
dir.create("Intermediate/03_Genetic_diversity")
dir.create("Intermediate/04_PAmatrix_by_site")
dir.create("Intermediate/04_ShapeFiles_buffered")
dir.create("Results")



#### -------- Load libraries ---------------------------------------------------

## Genetic diversity
library(adegenet)
library(hierfstat)
library(dartR)
library(mmod)

## Species diversity
library(ape)
library(betapart)


## R code
library(parallel) # mclapply
library(reshape2) #melt

## GIS
library(geojsonio) #geojson_json
library(raster) #crop #extent
library(rgeos) #gbuffer
library(rmapshaper) #ms_simplify
library(sf) #st_as_sf & st_read
library(marmap) #lc.dist
library(rgdal) #spTransform

## Statistics
library(vegan) #procrustes #mantel
library(ecodist) #MRM
library(nlme) #gls
library(r2glmm) #r2beta

## Plot
library(ggplot2)
library(lme4) #lmer modele mixte
library(cowplot)
library(ggrepel)
library(gridExtra)
library(patchwork)


# library(RColorBrewer)
# library(latex2exp)
# library(lmerTest) # pour avoir les p.values
# library(picante) #pd
# library(geiger) #name.check



#### -------- Plot function ----------------------------------------------------

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90, vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="darkgrey"),
            axis.line.y = element_line(colour="darkgrey"),
            axis.ticks = element_line(colour="darkgrey"),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.key.size= unit(0.5, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


