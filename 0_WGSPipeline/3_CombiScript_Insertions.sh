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


if [ -z "$ID_insert" ]; then # if override for $ID is not set, then crawl folder
	ID=$(ls -1 $myDataTrim | grep "_1P.fastq.gz" | sed -n ''$i'p'| cut -d"_1P.fastq.gz" -f1); else ID=$ID_insert
fi


### Here, we map the reads hard to the donor reference to find insertions!


# Aligning the reads on the reference donor genome dictionary.
# --- here, we map with "HARD" criteria that allow for little mapping inaccuracies
cd $bwaFold
./bwa mem $myDictPath/$dict".fasta" $myDataTrim/$ID"_1P.fastq.gz" $myDataTrim/$ID"_2P.fastq.gz" > $myDataPath/$IDout_insert".sam"

# Sorting of the aligned reads.
cd $samFold
./samtools sort $myDataPath/$IDout_insert".sam" --threads $threads_to_use -m $RAM_per_thread"M" --reference $myDictPath/$dict".fasta" -o $myDataPath/$IDout_insert"_sort.bam"
./samtools index -b $myDataPath/$IDout_insert"_sort.bam" > $myDataPath/$IDout_insert"_sort.bam.bai"

# Obtain the coverage at every position on the reference genome.
module unload gnu/4.4.7
module load gnu/5.1.0
cd $bedFold
./genomeCoverageBed -d -ibam $myDataPath/$IDout_insert"_sort.bam" > $myDataPath/$IDout_insert"_coverage.txt"

if $play_sound ; then echo $'\a'; fi # play finish sound
exit 0

