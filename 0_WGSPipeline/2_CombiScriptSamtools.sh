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

## use config.env to store your paths and filenames, they get read from there
## MUST BE RUN WITH bash, sh DOES NOT WORK
source $(dirname "$0")/config.env


if [ -z "$ID_align" ]; then # if override for $ID is not set, then crawl folder
	ID=$(ls -1 $myDataRaw | grep $search_term_align | grep "20_1.fastq.gz" | sed -n ''$i'p'| cut -d"_" -f1); else ID=$ID_align
fi

# Aligning the reads on the reference genome dictionary.
cd $bwaFold
./bwa mem -t $threads_to_use -B 1 -O 1 $myDictPath/$dict".fasta" $myDataTrim/$ID"_1P.fastq.gz" $myDataTrim/$ID"_2P.fastq.gz" > $myDataPath/$IDout_align".sam"

# Sorting of the aligned reads.
cd $samFold
./samtools sort $myDataPath/$IDout_align".sam" --threads $threads_to_use -m $RAM_per_thread"M" --reference $myDictPath/$dict".fasta" -o $myDataPath/$IDout_align"_sort.bam"
./samtools index -b $myDataPath/$IDout_align"_sort.bam" > $myDataPath/$IDout_align"_sort.bam.bai"

# Perform mpileup: Alignments are piled up for every position on the reference genome.
# if you are using an old samtools version, then use the following command: 
# # cd $samFold
# # ./samtools mpileup -e 10 -t AD -F 0.00001 -h 80 -L 10000 -o 20 -f $myDictPath/$dict.fasta -uv $myDataPath/$IDout_align"_sort.bam" > $myDataPath/$IDout_align".vcf"
# normally use the following:
cd $bcfFold
./bcftools mpileup -e 10 -F 0.00001 -h 80 -L 10000 -o 20 -a FORMAT/AD -d 8000 -f $myDictPath/$dict.fasta $myDataPath/$IDout_align"_sort.bam" > $myDataPath/$IDout_align"_bcf.vcf"

# Variant calling - Here, only variants are called arr given as an output.
cd $bcfFold
./bcftools call -vc $myDataPath/$IDout_align"_bcf.vcf" > $myDataPath/$IDout_align"_bcfcall.vcf"

# Obtain the coverage at every position on the reference genome.
module unload gnu/4.4.7
module load gnu/5.1.0
cd $bedFold
./genomeCoverageBed -d -ibam $myDataPath/$IDout_align"_sort.bam" > $myDataPath/$IDout_align"_coverage.txt"

if $play_sound ; then echo $'\a'; fi # play finish sound
exit 0
