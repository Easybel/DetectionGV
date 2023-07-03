# The scripts in this folder:

- **0_DicIndex.sh**
  - The dictionaries needed for the alignment steps are indexed and sorted.
  - Programs used in the script
      - bwa https://github.com/lh3/bwa
      - samfold from https://github.com/samtools/samtools
      - Picard https://github.com/broadinstitute/picard (requires Java)

- **1_QCScript.sh**
  - Raw reads are trimmed and adapters are removed. The read quality is assessed before and after the trimming.
  - Programs used in the script
      - FastQC https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
      - Trimmomatic https://github.com/usadellab/Trimmomatic

- **2_CombiScriptSamtools.sh**
  - The trimmed reads are aligned to the reference and variants are called. The coverage at each position is evaluated.
  - Programs used in the script
      - bwa and samtools, see above
      - BCFtools https://github.com/samtools/bcftools
      - Bedtools 2 https://github.com/arq5x/bedtools2
   
- **3_CombiScript_Insertions.sh**
  - The trimmed reads are aligned to the donor reference (hard) and the coverage at each position is evaluated.
  - Programs used in the script
      - bwa, samtools and Bedtools 2, see above
