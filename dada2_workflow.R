# R Script for dada2 pipeline to convert 16S rDNA amplicons to ASVs

# Load packages
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and 
# SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = 6))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = 6))

# Extract sample names. **THIS WILL DEPEND ON .FASTQ NAMES**
# E.g., if my samples had the names cruise_depth_plate_date_sample_R1.fastq, I would separate
# the fields into each individual piece and name samples with a subset
sample.fields <- lapply(strsplit(basename(fnFs), "_"), `[`, c(1,2,3,4,5))
sample.names <- sapply(sample.fields, function(fields) paste(fields[1], fields[2], fields[4], fields[5], sep="_"))

# Filter & Trimming Steps: 
# Assign the filenames for the filtered fastq.gz files.
# Place filtered files in ./filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Set standard filtering parameters; maxEE = max number of expected errors allowed 
# truncLen will determine the bp length of the F and R reads, respectively
# trimLeft will remove X bps from the beginning and the F and R reads (useful for trimming primers)
# Filtering and trimming processes will occur with these parameters: 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=6,
                     compress=TRUE, multithread=6, trimLeft = c(15,15)) 

# Inspect number of reads per samples retained after filter/trim
head(out)

# Learn the Error Rates: the DADA2 algorithm uses a model (err) to adjust 
# for the different set of error rates for different datasets, using function errF()
errF <- learnErrors(filtFs, multithread=6)
errR <- learnErrors(filtRs, multithread=6)

# Plot errors in PDF file and inspect
pdf("plotErrors.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

# Applying CORE SAMPLE INFERENCE algorithm to filtered/trimmed sequences to identify
# true biological variants
dadaFs <- dada(filtFs, err=errF, multithread=6)
dadaRs <- dada(filtRs, err=errR, multithread=6)

# Merging Steps:
# Merging paired reads (F and R together) to obtain full denoised sequences 
# Aligning denoised forward reads with reverse-complement of  the corresponding denoised 
# reverse reads to construct ontigs 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Constructing a sequence table of ASVs
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Chimera Removal Step using the bimera method
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=6, verbose=TRUE)

# Inspect dimensions of the chimera-removed df
dim(seqtab.nochim)

# Calculate proportion of non-chimeric merged sequence variants/total 
# merged sequence variants (proportion). This should be close to 1
sum(seqtab.nochim)/sum(seqtab)

#Track Reads through Pipeline
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Read counts are tracked and stored in the object "track" at each step
head(track)

# Taxonomoy Assignment steps:
# This is where the SILVA database formatted for dada2 lives on the shamwow2 server
taxa <- assignTaxonomy(seqtab.nochim, "/mnt/nfs/sharknado/LimsData/Hallam_Databases/formatted/Dada2/silva_database/silva_nr99_v138.1_train_set.fa.gz", multithread=6)
# Assigning species: 
# This optional step provides higher resolution database annotation if possible
taxa <- addSpecies(taxa, "/mnt/nfs/sharknado/LimsData/Hallam_Databases/formatted/Dada2/silva_species_assignment_v138.1fa.gz")

# Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Handoff to phyloseq to create an object for later use:
# Note that "samdf" will be needed. This is a dataframe with sample.names (matching the
# seqtab.nochim column names) in the first column and metadata variables in the next columns

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Get nucleotide DNA sequences of all ASVs into a new object
asv_dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(asv_dna) <- taxa_names(ps)

# Add sequence data to phyloseq object and rename rows as ASV# instead of sequence
ps <- merge_phyloseq(ps, asv_dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Remove zero abundance ASVs (in case samples were removed)
ps = filter_taxa(ps, function(x) var(x) > 0, TRUE)

# Keep absolute count phyloseq object separate from relative abundance object
ps_absolute = ps

# Create relative abundance phyloseq object for later use
ps = transform_sample_counts(ps, function(x) x / sum(x) )

# Save the phyloseq object and all R variables for reproducibility
save.image("dada2_phyloseq_object.rdata")
