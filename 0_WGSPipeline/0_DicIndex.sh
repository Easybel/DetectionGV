#!/bin/bash -l

## use config.env to store your paths and filenames, they get read from there
## MUST BE RUN WITH bash, sh DOES NOT WORK
source $(dirname "$0")/config.env

cd $bwaFold
./bwa index $myDictPath/$dict".fasta"

module unload openjdk/1.8.0_202
module load openjdk/1.8.0_60
cd $picardFold
java "-Xmx"$RAM_to_use_max"G" "-Xms"$RAM_to_use_min"G" -jar picard.jar CreateSequenceDictionary -R $myDictPath/$dict".fasta" -O $myDictPath/$dict".dict"

cd $samFold
./samtools faidx $myDictPath/$dict".fasta"

if $play_sound ; then echo $'\a'; fi # play finish sound
exit 0