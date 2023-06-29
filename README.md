# DetectionGV
In this repository, bioinformatic scripts are collected in the following folders:

- **0_WGSPipeline**  
  - This contains all scripts needed to analyze  whole genome, raw sequencing reads of bacterial transformation hybrids 
- **1_Detection**
  - With the outputs from [0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline), different genetic variations in the hybrid's genomes can be detected.
- **2_GeneralScripts**
  - Here, general scripts are collected that are needed to sort and convert annotations as well as find gene orthologues with blast.
 
- **dictionaries_Bacillus**
  - All extra files that are needed across the different project folders for Bacillus hybrids are collected here. This includes:
    - masterlists for different donors
    - accessory genome listst for different donors
    - annotation files: these are initially downloaded as **.gff3** files and then converted to **bed.mat / .bed.txt** with the script [2_GeneralScripts/Convert_gff3_to_bed.m](https://github.com/Easybel/DetectionGV/blob/main/2_GeneralScripts/Convert_gff3_to_bed.m)
  
The original version of all scripts, as cited in the PhD thesis of Isabel Rathmann, is preserved on branch [PhDThesis_Rathmann2023](https://github.com/Easybel/DetectionGV/tree/PhDThesis_Rathmann2023).
