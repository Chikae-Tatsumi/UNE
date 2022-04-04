#Getting ready

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

DATABASE = "~/R/Database/sh_general_release_dynamic_all_10.05.2021.fasta"
setwd("~/R/Analysis/2_UNE/ITS")  ## CHANGE ME to the directory containing the fastq files.

#Identify primers
FWD <- "GTGARTCATCGAATCTTTG"  ## fITS7
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS4

fnFs <- sort(list.files(getwd(), pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(getwd(), pattern = "_R2_001.fastq.gz", full.names = TRUE))
fnFs
fnRs

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
    
#Remove Primers
cutadapt <- "/miniconda2/bin/cutadapt" 
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
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Inspect read quality profiles
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) 

#Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, 
method="consensus", multithread=TRUE, verbose=TRUE) 
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
    getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
    "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt")

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, DATABASE, multithread = TRUE, tryRC = TRUE) 
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file="taxonomy.txt")
rownames(seqtab.nochim) <- sample.names
write.table(seqtab.nochim, file="seqtabnochim.txt")

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa))
#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Remove 
ps_removed = subset_taxa(ps,(Kingdom  == "k__Fungi" ))
                             
#To output OTU table
otu_table.t<-t(ps_removed@otu_table)
ps.t<-cbind(otu_table.t,ps_removed@tax_table)
write.table(ps.t,  file="ASV_table.txt")

# Rarefication
ps.rarefied = rarefy_even_depth(ps_removed, rngseed=1, sample.size=min(sample_sums(ps_removed)), replace=F)
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

#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps.deseq))
names(dna) <- taxa_names(ps.deseq)
ps.deseq <- merge_phyloseq(ps.deseq, dna)
taxa_names(ps.deseq) <- paste0("ASV", seq(ntaxa(ps.deseq)))
ps.deseq

# Remove 
ps.deseq = subset_taxa(ps.deseq,(Kingdom  == "k__Fungi" ))

# Save
ps<-cbind(ps.deseq@otu_table,ps.deseq@tax_table)
write.table(ps,  file="deseq_ASV_table.txt")

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
