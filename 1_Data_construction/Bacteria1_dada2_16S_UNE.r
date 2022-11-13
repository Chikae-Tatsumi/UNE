#DADA2 Pipeline Tutorial (1.12)
#See http://benjjneb.github.io/dada2/tutorial.html

#Getting ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

DATABASE = "~/R/Database/tax/silva_nr_v132_train_set.fa.gz"
setwd("~/R/Analysis/2_UNE/16S")  ## CHANGE ME to the directory containing the fastq files.
list.files()
# a set of Illumina-sequenced paired-end fastq files that have been split by sample and from which the barcodes/adapters have already been removed

# Define which is forward fastq and reverse fastq
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(getwd(), pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(getwd(), pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
# visualizing the quality profiles of the forward reads and reverse reads (only sample 1 and 2):
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
# You should judge how length is good quality, and truncate the bad quality part

# Filter and trim
filtFs <- file.path(getwd(), "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(getwd(), "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(215,215),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                     compress=TRUE, multithread=FALSE, trimRight=c(32,31)) # On Windows set multithread=FALSE
# R1 contains RC of 806r primer (20) + Linker (2) + Pad (10)
# R2 contains RC of 806r primer (19) + Linker (2) + Pad (10)
head(out)

# If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails
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
dadaRs <- dada(derepRs, err=errR, multithread=TRUE

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

#Assign taxonomy
#Install "silva_nr_v132_train_set.fa.gz" from: https://zenodo.org/record/1172783#.XUmvQ_ZFw2w
#Other taxonomic reference database: http://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim,DATABASE, multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. 

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

#To output OTU table for Tax4Fun2
otu_table.t<-t(ps@otu_table)
ps.t<-cbind(otu_table.t,ps@tax_table)
write.table(ps.t,  file="ASV_table_withMitoChlo.txt")

# Remove 
ps_removed = subset_taxa(ps,(
                             Family  != "Mitochondria"|is.na(Family) &
                             Class   != "Chloroplast"|is.na(Class)  &
                             Kingdom  != "Eukaryota" ))
                             
#To output OTU table
otu_table.t<-t(ps_removed@otu_table)
ps.t<-cbind(otu_table.t,ps_removed@tax_table)
write.table(ps.t,  file="ASV_table.txt")
write.table(ps_removed@tax_table, file="seqtab.nochim.txt")

# Rarefication
ps.rarefied = rarefy_even_depth(ps_removed, rngseed=1,
sample.size=min(sample_sums(ps_removed)), replace=F)
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

#ADONIS for DESeq
library (vegan)
Seq_Depth <- rowSums(seqtab.nochim)
comm<- t(diagvst)
sink(file = "adonis.seq_depthforDESeq.doc")
adonis(comm ~ Seq_Depth, permutations=5000)
sink()

#ADONIS for Rarefication
comm<-ps.rarefied@otu_table@.Data
sink(file = "adonis.seq_depth_forRarefication.doc")
adonis(comm ~ Seq_Depth, permutations=5000)
sink()
