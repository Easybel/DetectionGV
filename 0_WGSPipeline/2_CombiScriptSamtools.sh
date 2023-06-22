#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH --account=AG-Maier
#SBATCH --mail-type=END 
#SBATCH --mail-user name@uni-koeln.de 
#SBATCH --time=12:00:00
#SBATCH --array=1
i=$SLURM_ARRAY_TASK_ID

#Define input and output paths.
myDataTrim="path to trimmed data"
myDictPath="path to refernce dictionary"
myDataPath="path to save output"

#Define folders where software is installed  
bwaFold="path to bwa-0.7.17" 
samFold="path to samtools-1.16.1"
bcfFold="path to bcftools-1.16" 
bedFold="path to bedtools2/bin" 

#Define variables 
dict="name of dictionary" 

# Define the name of data file or retrieve the name for every run i. Automatically define a name for the output files.
#ID="Wns1020"
ID=$(ls -1 $myDataRaw | grep "Wns" | grep "20_1.fastq.gz" | sed -n ''$i'p'| cut -d"_" -f1)
IDout=$ID"_2Ref" 

# From here on, no further changes are required.

# Aligning the reads on the reference genome dictionary.
cd $bwaFold
./bwa mem -t 8 -B 1 -O 1 $myDictPath/$dict".fasta" $myDataTrim/$ID"_1P.fastq.gz" $myDataTrim/$ID"_2P.fastq.gz" > $myDataPath/$IDout".sam"

# Sorting of the aligned reads.
cd $samFold
./samtools sort $myDataPath/$IDout".sam" --threads 8 --reference $myDictPath/$dict.fasta -o $myDataPath/$IDout"_sort.bam"
./samtools index -b $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_sort.bam.bai"

# Perform mpileup: Alignments are piled up for every position on the reference genome.
# if you are using an old samtools version, then use the following command: 
# # cd $samFold
# # ./samtools mpileup -e 10 -t AD -F 0.00001 -h 80 -L 10000 -o 20 -f $myDictPath/$dict.fasta -uv $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout".vcf"
# otherwise use the following:
cd $bcfFold
./bcftools mpileup -e 10 -F 0.00001 -h 80 -L 10000 -o 20 -a FORMAT/AD -d 8000 -f $myDictPath/$dict.fasta $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_bcf.vcf"

# Variant calling - Here, only variants are called arr given as an output.
cd $bcfFold
./bcftools call -vc $myDataPath/$IDout"_bcf.vcf" > $myDataPath/$IDout"_bcfcall.vcf"

# Obtain the coverage at every position on the reference genome.
module unload gnu/4.4.7
module load gnu/5.1.0
cd $bedFold
./genomeCoverageBed -d -ibam $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_coverage.txt"

exit 0

