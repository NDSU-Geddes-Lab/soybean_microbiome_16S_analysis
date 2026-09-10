# Preparation of the greenhouse 16S rRNA phyloseq object
# Merges DADA2 output with metadata, assigns ASV IDs, and creates the rarefied dataset

# Load required packages
library(phyloseq)
library(Biostrings)

# Load greenhouse phyloseq object from DADA2 processing
physeq <- readRDS("IDCY1.GHrhizo.ps.rds")

# Load metadata
metadata <- import_qiime_sample_data("Metadata_GH.txt")

# Merge phyloseq object with metadata
ps <- merge_phyloseq(physeq, sample_data(metadata))


# Add DNA sequences and assign ASV IDs
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
final.ps <- merge_phyloseq(ps, dna)
taxa_names(final.ps) <- paste0("ASV", seq(ntaxa(final.ps)))

# Save greenhouse phyloseq object
saveRDS(final.ps, "Final.ps.GH.Rhizo.2021.rds")


# Rarefy to 6,000 reads per sample
set.seed(777)
greenhouse.ps.rare <- rarefy_even_depth(final.ps, sample.size = 6000)

# Save rarefied phyloseq object

saveRDS(greenhouse.ps.rare, "Final.ps.GH.Rhizo.2021.rarefied.rds")
