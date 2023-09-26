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

## use config.env to store your paths and filenames, they get read from there
## MUST BE RUN WITH bash, sh DOES NOT WORK
source $(dirname "$0")/config.env


if [ -z "$ID_QC" ]; then # if override for $ID is not set, then crawl folder
	ID=$(ls -1 $myDataRaw | grep $search_term_QC | grep "20_1.fastq.gz" | sed -n ''$i'p'| cut -d"_" -f1); else ID=$ID_QC
fi

# Perform FastQC quality control on the raw reads.
cd $FastQCFold
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_1.fastq.gz" -t $threads_to_use
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_2.fastq.gz" -t $threads_to_use

# Trim the raw reads with the list of adapters provided in the folder $TrimmFold
cd $TrimmFold
java "-Xmx"$RAM_to_use_max"G" "-Xms"$RAM_to_use_min"G" -jar trimmomatic-0.36.jar PE -threads $threads_to_use -trimlog $ID.TrimLog $myDataRaw/$ID"_1.fastq.gz" $myDataRaw/$ID"_2.fastq.gz" $myDataTrim/$ID"_1P".fastq.gz $myDataTrim $ID"_1U".fastq.gz $myDataTrim/$ID"_2P".fastq.gz $myDataTrim/$ID"_2U".fastq.gz ILLUMINACLIP:$TrimmFold/adapters/TruSeq3-PE-2.fa:2:30:10:4 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Perform FastQC quality control again on trimmed data.
cd $FastQCFold
./fastqc -o $myQCTrim/ $myDataTrim/$ID"_1P.fastq.gz" -t $threads_to_use
./fastqc -o $myQCTrim/ $myDataTrim/$ID"_2P.fastq.gz" -t $threads_to_use

if $play_sound ; then echo $'\a'; fi # play finish sound
exit 0