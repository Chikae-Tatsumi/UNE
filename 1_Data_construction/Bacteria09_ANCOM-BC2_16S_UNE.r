library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ANCOMBC); packageVersion("ANCOMBC")

# Import files
design <- read.csv("~/R/Analysis/2_UNE/experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
seqtab.nochim <- read.table("seqtabnochim.txt",header=T)
seqtab.nochim <- as.matrix(seqtab.nochim)
taxa <- read.table("taxonomy.txt",header=T)
taxa <- as.matrix(taxa)

# Format data
rownames(design) <- rownames(seqtab.nochim)
ps.ancom <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),tax_table(taxa), sample_data(design))
dna <- Biostrings::DNAStringSet(taxa_names(ps.ancom))
names(dna) <- taxa_names(ps.ancom)
ps.ancom <- merge_phyloseq(ps.ancom, dna)
taxa_names(ps.ancom) <- paste0("ASV", seq(ntaxa(ps.ancom)))
ps.ancom

# Remove non-bacterial ASVs
ps.ancom = subset_taxa(ps.ancom,(Kingdom=="Bacteria"))
ps.ancom = subset_taxa(ps.ancom,(Family  != "Mitochondria"|is.na(Family) &
                             Order   != "Chloroplast"|is.na(Order)))

# Run ANCOM-BC2
# See http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
out = ancombc2(ps.ancom, fix_formula = "DFB*DFE")
res <- out$res
samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)

# Make ASV table
taxonomy <- ps.ancom@tax_table
tax.df <- data.frame(taxonomy)
tax.df$ASV_ID <- rownames(tax.df)

log_corr_abn.df <- data.frame(log_corr_abn)
log_corr_abn.df$ASV_ID <- rownames(log_corr_abn.df)

merge <- merge(log_corr_abn.df,tax.df,by="ASV_ID")
rownames(merge) <- merge$ASV_ID
merge <- merge[,-1]

# Save
write.table(merge, file="ancom_ASV_table.txt")