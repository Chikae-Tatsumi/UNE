conda activate picrust2
cd R/Analysis/2_UNE/16S/PICRUSt2
picrust2_pipeline.py -s rarefied_seqs.fasta -i rarefied_ASV.txt -o picrust2_out_pipeline_rarefied -p 1 --per_sequence_contrib --stratified
gunzip picrust2_out_pipeline_rarefied/EC_metagenome_out/*gz
gunzip picrust2_out_pipeline_rarefied/KO_metagenome_out/*gz
gunzip picrust2_out_pipeline_rarefied/pathways_out/*gz