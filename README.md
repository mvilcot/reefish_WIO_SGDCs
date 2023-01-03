# Reef fish SGDCs
Codes for the paper: Spatial genetic differentiation correlates with species assemblage turnover across tropical reef fish lineages  
Maurine Vilcot, maurine.vilcot@ens-lyon.fr


# Data set
The paper explores species-genetic diversity correlations (SGDCs) across 20 reef fish species in four sampling sites of the Western Indian Ocean:
Maldives (MV), Mayotte (MY), Mafia Island (MF) and Seychelles (SC).

### Genetic SNPs data set
20 species, from 4479 to ... SNPs (Donati et al., 2021)

### Species presence data set
from (Albouy et al., 2019)


# 01 SNP Resampling
"As genomic diversity is directly related to the number of SNPs (Moragues et al., 2010), we randomly down-sampled 999 times the SNPs to the lowest common number of SNPs found across all species (i.e., n = 4479 SNPs in Pseudanthias squamipinnis; Table S2). To best account for differences in site sampling success across species, the number of individuals used for genetic diversity measurements was also randomly sampled 999 times to a maximum number of 10 individuals per sampling site (median value of the overall sampling). "


# 02 DAPC
Application of a Discriminant Analysis of Principal Components (DAPC) to each of the 20 species SNPs data, with sampling site as a prior information.


# 03 Genetic diversity
Computation of genetic diversity 
- β-GD: Hedrick G''st 
- α-SD: local expected heterozygosity (Hs) 
- γ-SD: regional expected heterozygosity (Ht) 


# 04 Species diversity
Computation of species diversity 
- β-SD: Species dissimilarity between sites 
	βjac-SD using Jaccard’s dissimilarity index 
	βjtu-SD  using turnover component of Jaccard’s dissimilarity index 
- α-SD: Local species richness, i.e. at each site 
- γ-SD: Regional species richness, i.e. from the four sites 


# 05 Statistical analysis
a. Application of α, β, and γ lineage-based SGDCs: linear models  
b. β-SGDCs for each of the 20 species separatedly: Mantel correlation test and Procruste analysis  
c. Influence of dispersal on β-SGDCs: linear models, and PGLS  
	PLD ~ β-SGDCs  
	U-crit ~ β-SGDCs  
c. MRM relating β-GD or β-SD to (i) geographic distance between sites and a (ii) binary distance metric reflecting whether the studied sites are separated or not from the Maldives by the Monsoon Drift  


