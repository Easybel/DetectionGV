####
### B0_BlastGeneFasta_2_DB.sh
####
# crated by Isabel

# What the script does:
# # -- it takes a fasta file as query and blasts it against a database, the subject
# # --> this database is created from a reference fasta

##########    SET VARIABLES  ##########################
# the paths where to find the data
myQueryPath=" path to the query" # for example, a list of gene sequences
mySubjectPath=" path to the subject" # the reference fasta
myOutPath=" path to output"

# set the variables 
# The genes that you would like to blast ...
queryFile="name of fasta"

# ...against this dictionary 
subjectFile="name of reference fasta"

# the outputs will have this name:
outName="Blast_"$queryFile"_2_"$subjectFile

##########   SET VARIABLES ##########################
## OPTIONAL ##
### if you wish to use a fasta as reference database, then run this command before to create the reference
makeblastdb -in $mySubjectPath/$subjectFile.fasta -dbtype nucl -parse_seqids

##########    SCRIPT CODE   ##########################
# blast the single genes
blastn -db $mySubjectPath/$subjectFile.fasta -query $myQueryPath/$queryFile.fasta -outfmt 7 -out $myOutPath/$outName".tab"

# get only the hit lines
grep -v "^#" $myOutPath/$outName".tab" > $myOutPath/$outName"_editted.tab"
