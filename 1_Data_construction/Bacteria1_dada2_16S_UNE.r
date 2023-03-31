# DADA2 Pipeline Tutorial (1.8)
# See https://benjjneb.github.io/dada2/tutorial_1_8.html

#Getting ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

DATABASE = "~/R/Database/silva_nr_v138_train_set.fa.gz"
setwd("~/R/Analysis/2_UNE/16S")  ## CHANGE ME to the directory containing the fastq files.
list.files()

# Define which is forward fastq and reverse fastq
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(getwd(), pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(getwd(), pattern="_R2_001.fastq", full.names = TRUE))

#Identify primers
FWD <- "GTGCCAGCMGCCGCGGTAA"  ## 515f
REV <- "GGACTACHVGGGTWTCTAAT"  ## 806r

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(getwd(), "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(getwd(), "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
    
# See https://github.com/benjjneb/dada2/issues/977

#Remove Primers
cutadapt <- "/Users/chikae/opt/miniconda3/envs/cutadaptenv/bin/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(getwd(), "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
    
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_L")[[1]][1]
#get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Inspect read quality profiles
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])
# See https://github.com/benjjneb/dada2/issues/1316

#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#Inspect read quality profiles
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #plot

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
 mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #>90% is good?

# Track reads through the pipeline
# look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt")

# Assign taxonomy
# Install "silva_nr_v138_train_set.fa.gz" from: https://zenodo.org/record/1172783#.XUmvQ_ZFw2w
taxa <- assignTaxonomy(seqtab.nochim,DATABASE, multithread=TRUE,tryRC=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.table(taxa,file="taxonomy_withMitoChlo.txt")
write.table(seqtab.nochim,file="seqtabnochim_withMitoChlo.txt")

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa))
#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#To output ASV table that still have non-bacterial sequences
otu_table.t<-t(ps@otu_table)
ps.t<-cbind(otu_table.t,ps@tax_table)
write.table(ps.t,  file="ASV_table_withMitoChlo.txt")

# Remove 
ps_removed_1 = subset_taxa(ps,(Kingdom=="Bacteria"))
ps_removed_2 = subset_taxa(ps_removed_1,(Family  != "Mitochondria"|is.na(Family) &
                             Order   != "Chloroplast"|is.na(Order))
                             
# To output ASV table
otu_table.t<-t(ps_removed_2@otu_table)
ps.t<-cbind(otu_table.t,ps_removed_2@tax_table)
write.table(ps.t,  file="ASV_table.txt")
write.table(ps_removed_2@tax_table, file="taxonomy.txt")
write.table(ps_removed_2@refseq, file="refseq.txt")

# Rarefication
ps.rarefied = rarefy_even_depth(ps_removed_2, rngseed=1,
sample.size=min(sample_sums(ps_removed_2)), replace=F)
otu_table.t<-t(ps.rarefied@otu_table)
ps.t<-cbind(otu_table.t,ps.rarefied@tax_table)
write.table(ps.t,  file="rarefied_ASV_table.txt")
sum(as.numeric(ps.t[,1])) 

# Deseq2 (https://joey711.github.io/phyloseq-extensions/DESeq2.html)
library("DESeq2"); packageVersion("DESeq2")
design <- read.csv("~/R/Analysis/2_UNE/experimental_design.csv",header=T)
rownames(design) <- rownames(seqtab.nochim)
ps.deseq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa), sample_data(design))
                         
dna <- Biostrings::DNAStringSet(taxa_names(ps.deseq))
names(dna) <- taxa_names(ps.deseq)
ps.deseq <- merge_phyloseq(ps.deseq, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps.deseq)))
ps.deseq        

# DESeq2 conversion and call (see https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html)
diagdds = phyloseq_to_deseq2(ps.deseq, ~DFB*DFE)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# Investigate test results table
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

# Get ASV table (https://github.com/joey711/phyloseq/issues/283)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)
diagvst = diagvst-min(diagvst)

otu_table(ps.deseq) <- otu_table(diagvst, taxa_are_rows = TRUE)

# Save
otu_table<-diagvst
ps.t<-cbind(otu_table,ps.deseq@tax_table)
write.table(ps.t,  file="deseq_ASV_table.txt")

# ADONIS for DESeq
library (vegan)
Seq_Depth <- rowSums(seqtab.nochim)
comm<- t(diagvst)
sink(file = "adonis.seq_depthforDESeq.doc")
adonis(comm ~ Seq_Depth, permutations=5000)
sink()

# ADONIS for Rarefication
comm<-ps.rarefied@otu_table@.Data
sink(file = "adonis.seq_depth_forRarefication.doc")
adonis(comm ~ Seq_Depth, permutations=5000)
sink()
