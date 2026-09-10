# soybean_microbiome_16S_analysis
R scripts for field and greenhouse soybean 16S rRNA amplicon sequencing and phenotypic analyses investigating soybean microbiome responses to iron deficiency chlorosis (IDC).

# Script organization
Scripts are numbered according to the analysis workflow:

•	“01_DADA2.R” - DADA2 processing of 16S rRNA amplicon sequences from bulk soil, soybean rhizosphere, and endosphere samples collected from field sites (Prosper, Casselton, Colfax, and Leonard; 2021-2023), and from greenhouse rhizosphere samples.

•	“02_preparing_phyloseq_field_trials.R” - preparation of combined field phyloseq object (Prosper, Casselton, Colfax and Leonard; 2021-2023).

•	“03_bacterial_community_analysis_across_all_sites_years.R” - alpha diversity, beta diversity and statistical analyses of soybean bacterial communities across the four field sites over the three years.

•	“04_IDC_scores_ND17009GT_across_all_sites_years.R” - assessing IDC scores of ND17009GT, IDC-sensitive soybean genotype across the sites and years.

•	“05_IDC_scores_genotypes_across_all_sites_years.R” - assessing IDC scores of IDC-sensitive and tolerant soybean genotypes across the sites and years.

•	“06_bacterial_community_analysis_across_genotypes_Leonard.R” - alpha diversity, beta diversity and statistical analyses of soybean microbiome across IDC-sensitive and tolerant genotypes at Leonard over the years.

•	“07_IDC_scores_across_genotypes_Leonard.R” - assessing IDC scores of IDC-sensitive and tolerant soybean genotypes at Leonard over the years.

•	“08_bacterial_community_analysis_Leonard_Casselton_2022.R” - alpha diversity, beta diversity, compositional and statistical analyses of soybean microbiome under different fertilizer treatments at Leonard and Casselton in 2022.

•	“09_soybean_phenotypic_analysis_Leonard_Casselton_2022.R” - assessing IDC scores, chlorophyll content and shoot dry weight of IDC-sensitive and tolerant soybean genotypes under different fertilizer treatments at Leonard and Casselton in 2022.

•	“10_bacterial_community_analysis_Leonard_Casselton_2023.R” - alpha diversity, beta diversity, compositional and statistical analyses of soybean microbiome under different fertilizer treatments at Leonard and Casselton in 2023.

•	“11_soybean_phenotypic_analysis_Leonard_Casselton_2023.R” - assessing IDC scores, chlorophyll content and shoot dry weight of IDC-sensitive and tolerant soybean genotypes under different fertilizer treatments at Leonard and Casselton in 2023.

•	“12_preparing_phyloseq_greenhouse.R” - Preparation of the greenhouse phyloseq object.

•	“13_rhizosphere_bacterial_community_analysis_greenhouse.R” - alpha diversity, beta diversity, compositional and statistical analyses of soybean microbiome under different fertilizer treatments from the greenhouse experiment

•	“14_soybean_phenotypic_analysis_greenhouse.R” - assessing chlorophyll content and shoot dry weight of IDC-sensitive and tolerant soybean genotypes under different soil iron deficient and fertilizer treatments in the greenhouse experiment.

# Taxonomic classification
Bacterial 16S rRNA amplicon sequences were taxonomically classified using the SILVA v132 reference database.

# Data availability
Raw sequencing data and large data files are not included in this repository. Sequencing data accession information will be provided with the associated publication.
