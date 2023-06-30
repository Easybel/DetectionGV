#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH --account=XX
#SBATCH --mail-type=END 
#SBATCH --mail-user name@uni-koeln.de 
#SBATCH --time=12:00:00
#SBATCH --array=1
i=$SLURM_ARRAY_TASK_ID

# Here, we map the reads hard to the donor reference to find insertions!

#Define input and output paths.
myDataTrim="path to trimmed data"
myDictPath="path to donor reference dictionary"
myDataPath="path to save output"

#Define folders where software is installed  
bwaFold="path to bwa-0.7.17" 
samFold="path to samtools-1.16.1"
bcfFold="path to bcftools-1.16" 
bedFold="path to bedtools2/bin" 

#Define variables 
dict="name of donor dictionary" 

# Define the name of data file or retrieve the name for every run i. Automatically define a name for the output files.
#ID="Wns1020"
ID=$(ls -1 $myDataTrim | grep "20_1P.fastq.gz" | sed -n ''$i'p'| cut -d"_" -f1)
IDout=$ID"_2Donor" 

# From here on, no further changes are required.

# Aligning the reads on the reference donor genome dictionary.
# --- here, we map with "HARD" criteria that allow for little mapping inaccuracies
cd $bwaFold
./bwa mem $myDictPath/$dict".fasta" $myDataTrim/$ID"_1P.fastq.gz" $myDataTrim/$ID"_2P.fastq.gz" > $myDataPath/$IDout".sam"

# Sorting of the aligned reads.
cd $samFold
./samtools sort $myDataPath/$IDout".sam" --threads 8 --reference $myDictPath/$dict.fasta -o $myDataPath/$IDout"_sort.bam"
./samtools index -b $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_sort.bam.bai"

# Obtain the coverage at every position on the reference genome.
module unload gnu/4.4.7
module load gnu/5.1.0
cd $bedFold
./genomeCoverageBed -d -ibam $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_coverage.txt"

exit 0

