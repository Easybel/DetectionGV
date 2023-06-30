# DetectionGV

In this repository, all scripts needed to analyse whole genome sequencing data of transformation hybrids are collected. First, raw sequencing reads are processed and then different genomic variations are detected.
The pipeline roughly contains the following steps:

<img src="https://github.com/Easybel/DetectionGV/blob/main/Pipeline_WGS.png" width="600">


- **[0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline)**  
  - This contains all bioinformatic scripts needed to analyse whole genome, raw sequencing reads of bacterial transformation hybrids.
    
- **[1_Detection](https://github.com/Easybel/DetectionGV/tree/main/1_Detection)**
  - With the outputs from [0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline), different genomic variations in the hybrid's genomes can be detected. This includes
    - Orthologous recombinations/ indels/ denovo SNPs ->  [A1_SNP2CNP.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A1_SNP2CNP.m)
    - deletions/ duplications -> [A3_Cov2DelDup.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A3_Cov2DelDup.m)
    - insertions -> [A4_Cov2Ins.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A4_Cov2Ins.m)
  - for all genomic variations, the genes that were affected can be detected with [A2_Lists2Genes.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A2_Lists2Genes.m)
    
- **[2_GeneralScripts](https://github.com/Easybel/DetectionGV/tree/main/2_GeneralScripts)**
  - Here, general scripts are collected that are needed to sort and convert annotations as well as find gene orthologues with blast.
 
- **[dictionaries_Bacillus](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus)**
  - All extra files that are needed across the different project folders for Bacillus hybrids are collected here. This includes:
    - masterlists for different donors (in **[ml](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus/ml)** -- these lists are created with [A0b_MasterListFiltering.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0b_MasterListFiltering.m)).
    - accessory genome lists for different donors (in **[acc](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus/acc)** -- these lists are created with [https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0b_MasterListFiltering.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0c_AccessoryGenome.m)).
    - other files:
      - reference files (**.fasta** ) -- usually downloaded, e.g. from NCBI
      - annotation files: these are initially downloaded as **.gff3** files and then converted to **bed.mat / .bed.txt** with the script [2_GeneralScripts/Convert_gff3_to_bed.m](https://github.com/Easybel/DetectionGV/blob/main/2_GeneralScripts/Convert_gff3_to_bed.m)
      - recipient specific list of multimapping regions ([Bs166NCe_mm.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_mm.txt (created with [A0d_Multimapper.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0d_Multimapper.m))
      -  SNP artifacts ([Bs166NCe_mm.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_ArteSNPs.vcf)), coverage artifacts for deletions/ duplications [Bs166NCe_ArteCov.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_ArteCov.txt) and insertions (???)  that have to be excluded in the following analysis. These lists are created by running the according scripts that detect these genomic variations with mapping data between the recipient and its own reference. 
  
The original version of all scripts, as cited in the PhD thesis of Isabel Rathmann, is preserved on branch [PhDThesis_Rathmann2023](https://github.com/Easybel/DetectionGV/tree/PhDThesis_Rathmann2023).
