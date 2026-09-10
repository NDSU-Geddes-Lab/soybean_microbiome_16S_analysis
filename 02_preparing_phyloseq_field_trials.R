# Preparation of the field 16S rRNA phyloseq object
# Merges DADA2 output with metadata, assigns ASV IDs, and creates the rarefied dataset

# Load required packages
library(phyloseq)
library(readr)
library(Biostrings)

# Load phyloseq object from DADA2 processing
ps <- read_rds("ps.noeuk.rds")

# Load metadata
metadata <- import_qiime_sample_data("Metadata_noeuk.txt")

# Merge phyloseq object with metadata
physeq <- merge_phyloseq(ps, sample_data(metadata))
sample_data(physeq)

# Save phyloseq object with metadata
saveRDS(physeq, file.path("PS.all_field_data(2021_to_2023)_withoutDNAstringset.rds"))

# Add DNA sequences and assign ASV IDs
dna <- Biostrings::DNAStringSet(taxa_names(physeq))
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna)
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))
physeq


# Save final phyloseq object
saveRDS(physeq, file.path("Condo7_PS.all_field_data(2021_to_2023)"))

# Rarefy to 8,000 reads per sample
set.seed(777)
ps.rare <- rarefy_even_depth(physeq, sample.size = 8000)

# Save rarefied phyloseq object
saveRDS(ps.rare,file.path("ps.rare.all.rds"))
