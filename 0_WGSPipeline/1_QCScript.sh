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
myDataRaw ="path to raw data"
myQCRaw   ="path to save quality control of raw data"
myDataTrim="path to save trimmed data"
myQCTrim  ="path to save quality control of trimmed data"

#Define folders where software is installed.
FastQCFold = "path to FastQC"
TrimmFold  = "path to Trimmomatic-0.36"

# Define the name of data file or retrieve the name for every run i.
#ID="Wns1020"
ID=$(ls -1 $myDataRaw | grep "Wns" | grep "20_1.fastq.gz" | sed -n ''$i'p'| cut -d"_" -f1)

# Perform FastQC quality control on the raw reads.
cd $FastQCFold
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_1.fastq.gz" -t 8
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_2.fastq.gz" -t 8

# Trim the raw reads with the list of adapters provided in the folder $TrimmFold
cd $TrimmFold
java -Xmx30G -Xms24G -jar trimmomatic-0.36.jar PE -threads 8 -trimlog $ID.TrimLog $myDataRaw/$ID"_1.fastq.gz" $myDataRaw/$ID"_2.fastq.gz" $myDataTrim/$ID"_1P".fastq.gz $myDataTrim $ID"_1U".fastq.gz $myDataTrim/$ID"_2P".fastq.gz $myDataTrim/$ID"_2U".fastq.gz ILLUMINACLIP:$TrimmFold/adapters/TruSeq3-PE-2.fa:2:30:10:4 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Perform FastQC quality control again on trimmed data.
cd $FastQCFold
./fastqc -o $myQCTrim/ $myDataTrim/$ID"_1P.fastq.gz" -t 8
./fastqc -o $myQCTrim/ $myDataTrim/$ID"_2P.fastq.gz" -t 8

exit 0
  
