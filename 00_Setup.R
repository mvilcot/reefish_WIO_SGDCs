
#### -------- Load libraries ---------------------------------------------------------------------------------------------

library(ade4) #mantel.rtest, #dudi.pca
library(vegan) #procrustes
library(factoextra) #fviz_eig
library(nlme) #gls
library(lme4) #lmer modele mixte
library(r2glmm) #r2beta
library(lmerTest) # pour avoir les p.values
library(marmap) #lc.dist

library(ape) #read.tree
library(betapart) #beta.multi beta.pair phylo.beta.multi
library(dartR) # To keep some individuls (gl.keep.ind) and drop some loci (gl.drop.loc) from a genlight
library(hierfstat)
library(adegenet)
library(mmod) #Gst_Hedrick & pairwise_Gst_Hedrick
library(ecodist) #MRM

library(sf) #st_as_sf & st_read
library(rgdal) #spTransform
library(rgeos) #gbuffer
library(geojsonio) #geojson_json
library(raster) #crop #extent
library(rmapshaper) #ms_simplify
library(picante) #pd
library(geiger) #name.check

library(ggplot2)
library(patchwork)
library(ggVennDiagram)
library(ggrepel)
library(gridExtra)
library(RColorBrewer)
library(parallel)
library(reshape2) #melt
library(latex2exp)
# library(tidyr)


#### -------- Customed function ---------------------------------------------------------------------------------------------

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
            # legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


