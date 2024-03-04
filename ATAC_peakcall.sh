#!/bin/bash -l
#
# Set the name of the job
#SBATCH --job-name=ATACpeakcall
#
# Set the maximum memory allowed
#SBATCH --mem=30G
#
# Set the maximum run time
#SBATCH -t 48:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email
#
# The number of threads we will require
#SBATCH -n 1
#
# Set output and error log files
#SBATCH -o /<path_to_working_directory>/step2/Log/logOutPeakCall_ATAC.txt
#SBATCH -e /<path_to_working_directory>/step2/Log/logErrorPeakCall_ATAC.txt
#
# set the account for hpc cluster user
#SBATCH --account=userid
#
#SBATCH --export=ALL
#
########## BEGIN ACTUAL COMMANDS 


# Required modules
module load python/2.7
module load MACS2/2.1.0
module load samtools/1.2
module load java/latest
module load ucsc_tools/3.0.9


#Set necessary environmental variables for align2rawsignal
#Setting appropriate pathway variables
MCRROOT=/home/user/MCR/v714
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR
export MCR_CACHE_ROOT=/home/user/Temp/


# Begin commands for peak calling and signal track generation for control
macs2 callpeak -t /<path_to_working_directory>/ATAC/step1/FinalAlign/cnt1_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n cnt1 -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/cnt1_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt1.bedGraph" -of="bg" -n=5 -l=1 -mm=30
bedGraphToBigWig /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt1.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt1.bw
rm /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt1.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/cnt1_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackMAT/cnt1.mat" -of="mat" -n=5 -l=1 -mm=30

macs2 callpeak -t /<path_to_working_directory>/ATAC/step1/FinalAlign/cnt2_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n cnt2 -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/cnt2_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt2.bedGraph" -of="bg" -n=5 -l=1 -mm=30
bedGraphToBigWig /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt2.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt2.bw
rm /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt2.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/cnt2_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackMAT/cnt2.mat" -of="mat" -n=5 -l=1 -mm=30

macs2 callpeak -t /<path_to_working_directory>/ATAC/step1/FinalAlign/cnt3_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n cnt3 -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/cnt3_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt3.bedGraph" -of="bg" -n=5 -l=1 -mm=30
bedGraphToBigWig /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt3.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt3.bw
rm /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/cnt3.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/cnt3_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackMAT/cnt3.mat" -of="mat" -n=5 -l=1 -mm=30


# Begin commands for peak calling and signal track generation for treated
macs2 callpeak -t /<path_to_working_directory>/ATAC/step1/FinalAlign/prop1_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n prop1 -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/prop1_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop1.bedGraph" -of="bg" -n=5 -l=1 -mm=30
bedGraphToBigWig /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop1.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop1.bw
rm /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop1.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/prop1_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackMAT/prop1.mat" -of="mat" -n=5 -l=1 -mm=30

macs2 callpeak -t /<path_to_working_directory>/ATAC/step1/FinalAlign/prop2_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n prop2 -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/prop2_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop2.bedGraph" -of="bg" -n=5 -l=1 -mm=30
bedGraphToBigWig /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop2.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop2.bw
rm /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop2.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/prop2_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackMAT/prop2.mat" -of="mat" -n=5 -l=1 -mm=30

macs2 callpeak -t /<path_to_working_directory>/ATAC/step1/FinalAlign/prop3_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n prop3 -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/prop3_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop3.bedGraph" -of="bg" -n=5 -l=1 -mm=30
bedGraphToBigWig /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop3.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset/<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop3.bw
rm /<path_to_working_directory>/ATAC/step2/SignalTrackBigWig/prop3.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/ATAC/step1/FinalAlign/prop3_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/ATAC/step2/SignalTrackMAT/prop3.mat" -of="mat" -n=5 -l=1 -mm=30


