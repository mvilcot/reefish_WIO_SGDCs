# REEFISH SGDCs
Codes for the paper:  
__Spatial genetic differentiation correlates with species assemblage turnover across tropical reef fish lineages__  
Maurine Vilcot, maurine.vilcot@gmail.com


# Data set
The paper explores species-genetic diversity correlations (SGDCs) across 20 reef fish species in four sampling sites of the Western Indian Ocean:
Maldives (MV), Mayotte (MY), Mafia Island (MF) and Seychelles (SC).

### Genetic SNPs data set
Filtered SNPs data for 20 species, from 4479 to 38931 SNPs. From Donati et al. (2021).

### Species presence data set
Global species presence data. From Albouy et al. (2019).


# 01 SNP Resampling
999 down-sampling of each of the 20 species SNPs data to:  
- 10 individuals per sampling site
- 4479 SNPs (i.e. the lowest common number of SNPs found across all species)


# 02 DAPC
Application of a Discriminant Analysis of Principal Components (DAPC) to each of the 20 species SNPs data, with sampling site as a prior information.


# 03 Genetic diversity
Computation of genetic diversity 
- β-GD: Hedrick G''st 
- α-SD: local expected heterozygosity (Hs) 
- γ-SD: regional expected heterozygosity (Ht) 


# 04 Species diversity
Computation of species diversity  
- β-SD: Species dissimilarity between sites, with Jaccard’s dissimilarity index  
- α-SD: Local species richness, i.e. at each site  
- γ-SD: Regional species richness, i.e. from the four sites  


# 05 Statistical analysis
- α, β, and γ lineage-based SGDCs: linear models  
- β-SGDCs for each of the 20 species separatedly: Mantel correlation test and Procruste analysis  
- Relation between β-SGDCs and dispersal (PLD or U-crit): linear models and PGLS  
- MRM relating β-GD or β-SD to (i) geographic distance between sites and a (ii) binary distance metric reflecting whether the studied sites are separated or not from the Maldives by the Monsoon Drift  


