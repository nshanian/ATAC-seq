#!/bin/sh
#
# Set the name of the job
#$ -N AlignAll_ATAC
#
# Set the maximum memory allowed
#$ -l h_vmem=9G
#
# Set the maximum run time
#$ -l h_rt=48:00:00
#
#Check that there's no errors in the options I've set
#$ -w e
#
# The number of threads we will require
#$ -pe shm 8
#                                                                                                                                           
# Request large mem queue                                                                                                                   
#$ -P large_mem  
#
# Set output and error log files
#$ -o /<path_to_working_directory>/ATAC/step1/Log/logOutAlignAll_ATAC.txt
#$ -e /<path_to_working_directory>/ATAC/step1/Log/logErrorAlignAll_ATAC.txt
#
# Set the working directory for this jobs at the current directory the script has been submitted from
#$ -cwd
#
#
########## BEGIN ACTUAL COMMANDS


#Required modules
module purge
module load java/latest
module load samtools/1.2
module load python/2.7
module load cutadapt/1.8.1
module load picard-tools/1.92
module load bowtie/2.3.1

export BOWTIE2_INDEXES=/<path_to_bowtie2_index_directory>

#Begin commands for control group
cutadapt -a ACATCTCCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG -A CACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt1_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt1_trimmed_PE2.fastq.gz /<path_to_working_directory>/ATAC/fastq/cnt1_1_pf.fastq.gz /<path_to_working_directory>/ATAC/fastq/cnt1_2_pf.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt1_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt1_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1  /<path_to_working_directory>/refgenomes/bowtie2/hg38 -1 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt1_trimmed_PE1.fastq  -2 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt1_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/ATAC/step1/RawAlign/cnt1_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/ATAC/step1/RawAlign/cnt1_raw_aln.bam /<path_to_working_directory>/ATAC/step1/SortedAlign/cnt1_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_working_directory>/picard-tools/1.92/MarkDuplicates.jar INPUT=/<path_to_working_directory>/ATAC/step1/SortedAlign/cnt1_sorted_aln.bam OUTPUT=/<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/cnt1_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/ATAC/step1/Temp/cnt1_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/ATAC/step1/Temp/cnt1_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/cnt1_sorted_dedup_aln.bam > /<path_to_working_directory>/ATAC/step1/FinalAlign/cnt1_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a ACATCTCCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG -A CACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt2_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt2_trimmed_PE2.fastq.gz /<path_to_working_directory>/ATAC/fastq/cnt2_1_pf.fastq.gz /<path_to_working_directory>/ATAC/fastq/cnt2_2_pf.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt2_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt2_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1  /<path_to_working_directory>/refgenomes/bowtie2/hg38 -1 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt2_trimmed_PE1.fastq  -2 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt2_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/ATAC/step1/RawAlign/cnt2_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/ATAC/step1/RawAlign/cnt2_raw_aln.bam /<path_to_working_directory>/ATAC/step1/SortedAlign/cnt2_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_working_directory>/picard-tools/1.92/MarkDuplicates.jar INPUT=/<path_to_working_directory>/ATAC/step1/SortedAlign/cnt2_sorted_aln.bam OUTPUT=/<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/cnt2_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/ATAC/step1/Temp/cnt2_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/ATAC/step1/Temp/cnt2_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/cnt2_sorted_dedup_aln.bam > /<path_to_working_directory>/ATAC/step1/FinalAlign/cnt2_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a ACATCTCCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG -A CACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt3_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt3_trimmed_PE2.fastq.gz /<path_to_working_directory>/ATAC/fastq/cnt3_1_pf.fastq.gz /<path_to_working_directory>/ATAC/fastq/cnt3_2_pf.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt3_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt3_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1  /<path_to_working_directory>/refgenomes/bowtie2/hg38 -1 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt3_trimmed_PE1.fastq  -2 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/cnt3_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/ATAC/step1/RawAlign/cnt3_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/ATAC/step1/RawAlign/cnt3_raw_aln.bam /<path_to_working_directory>/ATAC/step1/SortedAlign/cnt3_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_working_directory>/picard-tools/1.92/MarkDuplicates.jar INPUT=/<path_to_working_directory>/ATAC/step1/SortedAlign/cnt3_sorted_aln.bam OUTPUT=/<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/cnt3_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/ATAC/step1/Temp/cnt3_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/ATAC/step1/Temp/cnt3_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/cnt3_sorted_dedup_aln.bam > /<path_to_working_directory>/ATAC/step1/FinalAlign/cnt3_sorted_dedup_filterUnmap_filterChrM_aln.bam


#Begin commands for treated group
cutadapt -a ACATCTCCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG -A CACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop1_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop1_trimmed_PE2.fastq.gz /<path_to_working_directory>/ATAC/fastq/prop1_1_pf.fastq.gz /<path_to_working_directory>/ATAC/fastq/prop1_2_pf.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop1_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop1_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1  /<path_to_working_directory>/refgenomes/bowtie2/hg38 -1 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop1_trimmed_PE1.fastq  -2 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop1_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/ATAC/step1/RawAlign/prop1_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/ATAC/step1/RawAlign/prop1_raw_aln.bam /<path_to_working_directory>/ATAC/step1/SortedAlign/prop1_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_working_directory>/picard-tools/1.92/MarkDuplicates.jar INPUT=/<path_to_working_directory>/ATAC/step1/SortedAlign/prop1_sorted_aln.bam OUTPUT=/<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/prop1_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/ATAC/step1/Temp/prop1_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/ATAC/step1/Temp/prop1_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/prop1_sorted_dedup_aln.bam > /<path_to_working_directory>/ATAC/step1/FinalAlign/prop1_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a ACATCTCCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG -A CACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop2_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop2_trimmed_PE2.fastq.gz /<path_to_working_directory>/ATAC/fastq/prop2_1_pf.fastq.gz /<path_to_working_directory>/ATAC/fastq/prop2_2_pf.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop2_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop2_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1  /<path_to_working_directory>/refgenomes/bowtie2/hg38 -1 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop2_trimmed_PE1.fastq  -2 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop2_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/ATAC/step1/RawAlign/prop2_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/ATAC/step1/RawAlign/prop2_raw_aln.bam /<path_to_working_directory>/ATAC/step1/SortedAlign/prop2_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_working_directory>/picard-tools/1.92/MarkDuplicates.jar INPUT=/<path_to_working_directory>/ATAC/step1/SortedAlign/prop2_sorted_aln.bam OUTPUT=/<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/prop2_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/ATAC/step1/Temp/prop2_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/ATAC/step1/Temp/prop2_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/prop2_sorted_dedup_aln.bam > /<path_to_working_directory>/ATAC/step1/FinalAlign/prop2_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a ACATCTCCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG -A CACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop3_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop3_trimmed_PE2.fastq.gz /<path_to_working_directory>/ATAC/fastq/1prop3_1_pf.fastq.gz /<path_to_working_directory>/ATAC/fastq/1prop3_2_pf.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop3_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop3_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1  /<path_to_working_directory>/refgenomes/bowtie2/hg38 -1 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop3_trimmed_PE1.fastq  -2 /<path_to_working_directory>/ATAC/step1/TrimmedFASTQ/prop3_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/ATAC/step1/RawAlign/prop3_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/ATAC/step1/RawAlign/prop3_raw_aln.bam /<path_to_working_directory>/ATAC/step1/SortedAlign/prop3_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_working_directory>/picard-tools/1.92/MarkDuplicates.jar INPUT=/<path_to_working_directory>/ATAC/step1/SortedAlign/prop3_sorted_aln.bam OUTPUT=/<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/prop3_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/ATAC/step1/Temp/prop3_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/ATAC/step1/Temp/prop3_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/ATAC/step1/SortedAndDedupAlign/prop3_sorted_dedup_aln.bam > /<path_to_working_directory>/ATAC/step1/FinalAlign/prop3_sorted_dedup_filterUnmap_filterChrM_aln.bam


