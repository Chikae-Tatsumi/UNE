library(Tax4Fun2)
library(seqinr)
library(ape)

setwd("~/R/Database/Tax4Fun2")
ASV.table <- read.table(file="ASV_table_withMitoChlo.txt",header=T,row.names=1)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)] # I've changed 7--> 6
ASV <- cbind (rownames(ASV),ASV)
ASV <- rbind (colnames(ASV),ASV)
ASV[1,1] <- "ID"
rownames(ASV) <- NULL 
colnames(ASV) <- NULL
write.table(ASV, "ASV.txt",sep="\t",col.names = F, row.names = F,quote=F)
taxonomy <- read.table(file="taxonomy.txt",header=T, row.names=1)
tax <- subset(taxonomy,  (Family  != "Mitochondria"|is.na(Family) &
                             Order   != "Chloroplast"|is.na(Order)  &
                             Kingdom  != "Eukaryota" ))
write.fasta (sequences = as.list(rownames(tax)),names = rownames(ASV.table),file.out="seqs.fasta")
dir.create("Tax4Fun2")

#Step 2: Generate your own reference datasets
#1. Extracting SSU seqeunces (16S rRNA and 18S rRNA)
extractSSU(genome_file = "OneProkaryoticGenome.fasta", file_extension = "fasta", path_to_reference_data ="Tax4Fun2_ReferenceData_v2")
#2. Assigning functions to prokayotic genomes
assignFunction(genome_file = "OneProkaryoticGenome.fasta", file_extension = "fasta", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", num_of_threads = 1, fast = TRUE)
#3. Generate the reference data
generateUserData(path_to_reference_data ="Tax4Fun2_ReferenceData_v2", path_to_user_data = ".", name_of_user_data = "User_Ref0", SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")

#Step 3: Making functional predictions
#1. Making functional predictions using the default reference data only
#1. Run the reference blast
runRefBlast(path_to_otus = "seqs.fasta" , path_to_reference_data ="Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", use_force = T, num_threads = 6)
# 2) Predicting functional profiles
# Remove the first row's # and the second row's # of "17_otu_97_table_taxonomy.txt" --> "17_otu_97_table_taxonomy_tax4fun.txt"
makeFunctionalPrediction(path_to_otu_table = "ASV.txt", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)
# note. normalize_pathways = FALSE will affiliate the rel. abundance of each KO to each pathway it belongs to. By setting it to true, the rel. abundance is equally distributed to all pathways it was assigned to.)

#Step 4: Calculating (multi-)functional redundancy indices (experimental)
calculateFunctionalRedundancy(path_to_otu_table = "ASV.txt", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", min_identity_to_reference = 0.97)

# Don't forget to move the output file (Tax4Fun2) for your project file