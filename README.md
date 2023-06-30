# DetectionGV

In this repository, all scripts needed to analyse whole genome sequencing data of transformation hybrids are collected. First, raw sequencing reads are processed and then different genomic variations are detected.
The pipeline roughly contains the following steps:

<img src="https://github.com/Easybel/DetectionGV/blob/main/Pipeline_WGS.png" width="600">


- **[0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline)**  
  - This contains all scripts needed to analyze  whole genome, raw sequencing reads of bacterial transformation hybrids 
- **[1_Detection](https://github.com/Easybel/DetectionGV/tree/main/1_Detection)**
  - With the outputs from [0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline), different genetic variations in the hybrid's genomes can be detected.
- **[2_GeneralScripts](https://github.com/Easybel/DetectionGV/tree/main/2_GeneralScripts)**
  - Here, general scripts are collected that are needed to sort and convert annotations as well as find gene orthologues with blast.
 
- **[dictionaries_Bacillus](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus)**
  - All extra files that are needed across the different project folders for Bacillus hybrids are collected here. This includes:
    - ml: masterlists for different donors
    - acc: accessory genome listst for different donors
    - other files:
      - reference files (**.fasta** )
      - annotation files: these are initially downloaded as **.gff3** files and then converted to **bed.mat / .bed.txt** with the script [2_GeneralScripts/Convert_gff3_to_bed.m](https://github.com/Easybel/DetectionGV/blob/main/2_GeneralScripts/Convert_gff3_to_bed.m)
      - recipient specific list of multimapping regions ([Bs166NCe_mm.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_mm.txt)), SNP artifacts ([Bs166NCe_mm.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_ArteSNPs.vcf)), and coverage artifacts [Bs166NCe_ArteCov.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_ArteCov.txt) that have to be excluded in the following analysis
  
The original version of all scripts, as cited in the PhD thesis of Isabel Rathmann, is preserved on branch [PhDThesis_Rathmann2023](https://github.com/Easybel/DetectionGV/tree/PhDThesis_Rathmann2023).
