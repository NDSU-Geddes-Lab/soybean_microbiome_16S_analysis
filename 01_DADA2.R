# DADA2 processing of field and greenhouse 16S rRNA amplicon sequencing data
# The same sequence-processing workflow and parameters were applied to both dataset

# Load required packages
suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(argparser)
  library(logger)
})

# Parse command-line arguments
parser <- arg_parser(
  "16S rRNA amplicon sequence processing",
  hide.opts = TRUE
)

parser <- add_argument(
  parser, "fastq_dir",
  help = "Directory containing input FASTQ files (.fastq.gz)"
)

parser <- add_argument(
  parser, "--db",
  help = "Path to the SILVA taxonomy reference database",
  default = NULL
)

parser <- add_argument(
  parser, "--outdir",
  help = "Directory for output files",
  default = "./output"
)

argv <- parse_args(parser)

indir <- argv$fastq_dir
outdir <- argv$outdir
dbpath <- argv$db


# Create output directories
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "filtered"), showWarnings = FALSE, recursive = TRUE)

# Create PDF file for plots
pdf(file.path(outdir, "Rplots.pdf"))

# Determine number of available CPU cores
ncpus <- parallel::detectCores()

if (is.na(ncpus)) {
  ncpus <- 1
}

# Sort files to ensure forward and reverse reads are in the same order
fnFs <- sort(list.files(indir, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(indir, pattern = "_R2_001.fastq", full.names = TRUE))


# Stop if no forward reads are found
if (length(fnFs) == 0) {
  log_error(paste("No reads found in", indir))
  quit(save = "no", status = 1)
}


# If number of forward and reverse reads don't match, throw error and exit.
if (length(fnFs) != length(fnRs)) {
  log_error("Number of forward and reverse reads don't match.")
  quit(save="no", status=1)
}

# Extract sample names, adjusted to handle names with underscores from tutorial code
fnFs.names <- sapply(strsplit(basename(fnFs), "_R1_001.fastq.gz"), `[`, 1) 
fnRs.names <- sapply(strsplit(basename(fnRs), "_R2_001.fastq.gz"), `[`, 1) 


# Make sure sample names match for forward and reverse reads.
if (! all(fnFs.names == fnRs.names)) {
  log_error("Forward and reverse sample names don't match.")
  quit(save="no", status=1)
}

sample.names <- fnFs.names



# Plot quality profiles for the first few sequences, which should be enough to
# indicate the overall quality of the sequencing run.
plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])



# Place filtered files in filtered/ subdirectory
filtFs <- file.path(outdir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names



# Filter and trim reads based on quality profiles and primer lengths
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230),trimLeft = c(19,21),
                     maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=ncpus)
head(out)

write.csv(out, file = file.path(outdir, "filtered_summary.csv"), row.names = FALSE)
saveRDS(out, file = file.path(outdir, "filtered_summary.rds"))



#Plot again to see if the trimming worked
plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])



#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=ncpus)
errR <- learnErrors(filtRs, multithread=ncpus)

saveRDS(errF, file = file.path(outdir, "errF.rds"))
saveRDS(errR, file = file.path(outdir, "errR.rds"))


# Plot estimated error rates
plotErrors(errF, nominalQ=TRUE)


# Infer sequence variants
dadaFs <- dada(filtFs, err=errF, multithread=ncpus)
dadaRs <- dada(filtRs, err=errR, multithread=ncpus)

saveRDS(dadaFs, file = file.path(outdir, "dadaFs.rds"))
saveRDS(dadaRs, file = file.path(outdir, "dadaRs.rds"))


#Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]


# Merge paired-end reads 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

saveRDS(mergers, file = file.path(outdir, "mergers.rds"))



# Construct ASV sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

saveRDS(seqtab, file = file.path(outdir, "seqtab.rds"))
seqtab_df <- as.data.frame(seqtab)
write.csv(seqtab_df, file = file.path(outdir, "seqtab.csv"), row.names = TRUE)


# Inspect distribution of sequence lengths
size <- table(nchar(getSequences(seqtab)))
write.csv(size, file.path(outdir, "amplicon_size.csv"))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:253]


# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=ncpus, verbose=TRUE)

saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab_nochim.rds"))
seqtab_nochim_df <- as.data.frame(seqtab.nochim)
write.csv(seqtab_nochim_df, file = file.path(outdir, "seqtab_nochim.csv"), row.names = TRUE)



# Calculate the proportion of reads retained after chimera removal
sum(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Remove singleton ASVs
is1 <- colSums(seqtab.nochim) <= 1
seqtab.nochim1 <- seqtab.nochim[,!is1] 
sum(seqtab.nochim1)/sum(seqtab.nochim)


# Assign taxonomy
taxa.nochim1 <- assignTaxonomy(seqtab.nochim1, dbpath, multithread=ncpus)

saveRDS(taxa.nochim1, file = file.path(outdir, "taxa_nochim1.rds"))
taxa_nochim1_df <- as.data.frame(taxa.nochim1)
write.csv(taxa_nochim1_df, file = file.path(outdir, "taxa_nochim1.csv"), row.names = TRUE)



# Remove chloroplast reads
is.chloro <- taxa.nochim1[,"Order"] %in% "Chloroplast" 
seqtab.nochloro <- seqtab.nochim1[,!is.chloro]
taxa.nochloro <- assignTaxonomy(seqtab.nochloro, dbpath, multithread=ncpus)
sum(seqtab.nochloro)
sum(seqtab.nochloro)/sum(seqtab.nochim1)


# Save the  no-chloroplast sequence table and taxonomy
saveRDS(seqtab.nochloro, file = file.path(outdir, "seqtab_nochloro.rds"))
saveRDS(taxa.nochloro, file = file.path(outdir, "taxa_nochloro.rds"))


# Remove mitochondrial reads
is.mito <- taxa.nochloro[,"Family"] %in% "Mitochondria"
seqtab.nomito <- seqtab.nochloro[,!is.mito]
taxa.nomito <- assignTaxonomy(seqtab.nomito, dbpath, multithread=ncpus)
sum(seqtab.nomito)
sum(seqtab.nomito)/sum(seqtab.nochloro)
sum(seqtab.nomito)/sum(seqtab.nochim1)


# Save the  no-mitochondria sequence table and taxonomy
saveRDS(seqtab.nomito, file = file.path(outdir, "seqtab_nomito.rds"))
saveRDS(taxa.nomito, file = file.path(outdir, "taxa_nomito.rds"))

# Remove eukaryotic reads
is.euk <- taxa.nomito[,"Kingdom"] %in% "Eukaryota"
seqtab.noeuk <- seqtab.nomito[,!is.euk]
taxa.noeuk <- assignTaxonomy(seqtab.noeuk, dbpath, multithread=ncpus)
taxa.print <- taxa.noeuk
rownames(taxa.print) <- NULL

sum(seqtab.noeuk)
sum(seqtab.noeuk)/sum(seqtab.nomito)
sum(seqtab.noeuk)/sum(seqtab.nochim1)




#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track.noeuk <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim1), rowSums(seqtab.nochloro), rowSums(seqtab.nomito), rowSums(seqtab.noeuk))


colnames(track.noeuk) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim", "nosingleton", "nochloro", "nomito", "noeuk")
rownames(track.noeuk) <- sample.names
head(track.noeuk)
write.csv(track.noeuk, file.path(outdir, "track.sanity.lastparameters.noeuk.csv"))

# Save final objects for use in phyloseq
saveRDS(taxa.noeuk, file = file.path(outdir, "taxa.noeuk.rds"))
saveRDS(seqtab.noeuk, file = file.path(outdir, "seqtab.noeuk.rds"))


# Write sequence and taxonomy table
dtb <- cbind(t(seqtab.noeuk), taxa.noeuk)
write.csv(dtb, file.path(outdir, "dtb.seqtab.noeuk.csv"))

# Create phyloseq object
ps.noeuk <- phyloseq(otu_table(seqtab.noeuk, taxa_are_rows=F), 
                     tax_table(taxa.noeuk))

saveRDS(ps.noeuk, file = file.path(outdir, "ps.noeuk.rds"))

dev.off()
