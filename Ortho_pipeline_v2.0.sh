#!/bin/sh

#===============================================================================================
# Script to detect putative orthologous TEs, deletions and insertions between two close genomes
# 
# USAGE: 
#	 bash Ortho_pipeline_v2.0.sh flank_size 
#
# DEPENDENCIES:
#	- BBmap
#	- Bedtools
#	- Samtools
#================================================================================================



# ------------------------------------ Set this BEFORE running the pipeline ---------------------------------

# Define full paths to necessary files:

scripts_folder="/path/to/scripts"     # Without the last '/'

# Srain1  (here reference):

reference_genome_file="genome_1.fai" # bedtools genome file, NOT fasta genome file
reference_fasta_file="genome_1.fa"


# Species 2 (here query):

query_genome_file="genome_2.fai"  # bedtools genome file, NOT fasta genome file
query_fasta_file="genome_2.fa"
LTRharvest_annot_query="TE_annotation.gff"   # Important! Third field should contain: Species_TE_1, Species_TE_2, ...


#---------------------------------------- PROGRAM -------------------------------------------------------------

################### STEP1: Identification of OrthoLTR sites

# Extract flanks: $1 gff file

bedtools flank -i $LTRharvest_annot_query -g $query_genome_file -b $1 > LTR_flanks.gff;
printf 'Flanks extracted..\n'

# Rename LTR flank IDs (Flanks _1 and _2)

python $scripts_folder/rename.py LTR_flanks.gff > LTR_flanks_renamed_tmp.gff;
awk '{if ($5-$4 > 100) print $0}' LTR_flanks_renamed_tmp.gff > LTR_flanks_renamed.gff 
rm LTR_flanks_renamed_tmp.gff

printf 'Flanks renamed..........\n'

# Extractfasta

bedtools getfasta -fi $query_fasta_file -bed LTR_flanks_renamed.gff -name > reference_LTR_flanks.fasta;
printf 'Extract flanks in fasta format..\n'


################### STEP2:  Map as paired end from interleaved reads (rcompmate=f --> reverse complement read 2) 

bbmap.sh local=f interleaved=TRUE rcompmate=f ambiguous=best in=reference_LTR_flanks.fasta ref=$reference_fasta_file nodisk out=mapped.sam slow k=12;

# make histogram of insert sizes (remember substracting flank_size later)

reformat.sh in=mapped.sam ihist=histogram.txt 

# transform to UNSORTED BAM

samtools view -bS mapped.sam > mapped.bam;

# obtain only records with BOTH mapped reads

samtools view -F 12 mapped.bam -b > mapped_paired.bam; 

# transform to BED

bedtools bamtobed -i mapped_paired.bam  > mapped_sorted.bed;


grep '/1' mapped_sorted.bed > forward.bed; 
grep '/2' mapped_sorted.bed > reverse.bed; 

paste forward.bed reverse.bed > combined.bed;

cat forward.bed reverse.bed | sort -k1,1 -k4,4n > combined_flanks_IGV.bed;

awk '{{FS="\t"} if ($1 == $7) print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$9"\t"$10}' combined.bed > combined_pairs.bed;

# Print length files to identify partial deletions

awk '{print$3"\t"($5-$4)}' $LTRharvest_annot_query > lengths.txt 

# Updated to avoid problems with 2/1, 1/1 in the name etc..

awk '{{FS="\t"} if ($5-$3 < 16000 && $5-$3 > -50 ) print $4"\t"$5-$3}' combined_pairs.bed | tr '_' '\t' | awk '{print $1"_"$2"_"$3"\t"$5}' > insertsize.txt
#awk '{{FS="\t"} if ($5-$3 < 16000 && $5-$3 > -50 ) print $4"\t"$5-$3}' combined_pairs.bed | sed -e 's:_1/1::g' | sed 's:_2/2::g' > insertsize.txt;

python $scripts_folder/classifyLTR.py

mkdir intermediate_files
mv LTR_flanks.gff combined.bed forward.bed lengths.txt mapped_paired.bam reference_LTR_flanks.fasta LTR_flanks_renamed.gff combined_flanks_IGV.bed histogram.txt mapped.bam mapped_sorted.bed reverse.bed combined_pairs.bed insertsize.txt mapped.sam intermediate_files 


#awk '{{FS="\t"} if ($5-$3 < 15 && $5-$3 > -15 ) print $1"\t"$5"\t"$3"\t"$4"\t"$5-$3}' combined_pairs.bed | sort -k5,5n  > Putative_insertions;
#awk '{{FS="\t"} if ($5-$3 > 100 && $5-$3 < 800 ) print $1"\t"$5"\t"$3"\t"$4"\t"$5-$3}' combined_pairs.bed | sort -k5,5n  > All_conserved.txt
