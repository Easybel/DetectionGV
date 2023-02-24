#specify the software that is needed, the path to the dict and the name if the dict
#Define input and output paths
#In the myQCPath there should be a folder called FastQC_AfterTrimm, where the trimmed QC of the results will be put

myDictPath="/home/irathman/scripts/scripts_git/test"


#Define folders where software is installed

bwaFold="/home/irathman/sw/bwa-0.7.17"
samFold="/home/irathman/sw/samtools-1.8"
picardFold="/home/irathman/sw"
#GATKFold="/home/irathman/sw/GATK3.5"
#bcfFold="/home/irathman/sw/bcftools-1.8"

# here you map against: .fasta
dict="BsubNC_000964wt"

cd $bwaFold
./bwa index $myDictPath/$dict.fasta

module unload openjdk/1.8.0_202
module load openjdk/1.8.0_60
cd $picardFold
java -Xms1g -Xmx3g -jar picard.jar CreateSequenceDictionary R=$myDictPath/$dict.fasta O=$myDictPath/$dict.dict

cd $samFold
./samtools faidx $myDictPath/$dict.fasta
