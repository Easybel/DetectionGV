# Introduction

In this repository, all scripts needed to analyse whole genome sequencing data of transformation hybrids are collected. First, raw sequencing reads are processed in the "Raw reads analysis pipeline" and then different genomic variations are detected in "Further analysis".  
The pipeline is visualized in the following scheme:  

  
<img src="https://github.com/Easybel/DetectionGV/blob/main/Pipeline_WGS.png" width="900">

A more detailed diagram can be found on LucidChart <https://lucid.app/documents/view/fa9f52da-edf4-4fb2-9b6c-d01996f4768e> or under [/documentation](/documentation)

Preview:

<img src="https://github.com/Easybel/DetectionGV/blob/main/documentation/DetectionGV_preview.png" width="900">

# Overview of folders and scripts
- **[0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline)**  
  - This contains all bioinformatic scripts needed to analyse whole genome, raw sequencing reads of bacterial transformation hybrids.
    
- **[1_Detection](https://github.com/Easybel/DetectionGV/tree/main/1_Detection)**
  - With the outputs from [0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline), different genomic variations in the hybrid's genomes can be detected. This includes:
    - orthologous replacements, denovo SNPs & indels :arrow_right:  [A1_SNP2CNP.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A1_SNP2CNP.m)
    - deletions/ duplications :arrow_right: [A3_Cov2DelDup.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A3_Cov2DelDup.m)
    - insertions :arrow_right: [A4_Cov2Ins.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A4_Cov2Ins.m)
  - For all genomic variations, the affected genes can be detected with [A2_Lists2Genes.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A2_Lists2Genes.m)
    
- **[2_GeneralScripts](https://github.com/Easybel/DetectionGV/tree/main/2_GeneralScripts)**
  - Here, general scripts are collected that are needed to sort and convert annotations as well as find gene orthologues with blast.
 
- **[dictionaries_Bacillus](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus)**
  - All extra files that are needed across the different project folders for _Bacillus_ hybrids are collected here. This includes:
    - masterlists for different donors (in **[ml](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus/ml)**. These lists are created with [A0b_MasterListFiltering.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0b_MasterListFiltering.m)).
    - accessory genome lists for different donors (in **[acc](https://github.com/Easybel/DetectionGV/tree/main/dictionaries_Bacillus/acc)**. These lists are created with [A0c_AccessoryGenome.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0c_AccessoryGenome.m)).
    - Other files:
      - reference files (**.fasta** ), usually downloaded, e.g. from NCBI
      - annotation files: initially downloaded as **.gff3** files and then converted to **bed.mat / .bed.txt** with the script [2_GeneralScripts/Convert_gff3_to_bed.m](https://github.com/Easybel/DetectionGV/blob/main/2_GeneralScripts/Convert_gff3_to_bed.m)
      - recipient specific list of multimapping regions ([Bs166NCe_mm.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_mm.txt) created with [A0d_Multimapper.m](https://github.com/Easybel/DetectionGV/blob/main/1_Detection/A0d_Multimapper.m))
      -  SNP artefacts ([Bs166NCe_mm.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_ArteSNPs.vcf)), coverage artefacts for deletions/ duplications [Bs166NCe_ArteCov.txt](https://github.com/Easybel/DetectionGV/blob/main/dictionaries/Bs166NCe_ArteCov.txt) and insertions have to be excluded in the following analysis. These lists are created by running the according scripts that detect these genomic variations with mapping data between the recipient and its own reference.
     
- **[ToyData](https://github.com/Easybel/DetectionGV/tree/main/ToyData)**
  - In this folder, we provide a toy data set of a hybrid sample (recipient: _Bacillus subtilis_, donor: _Bacillus spizizenii_). The files are the outputs from the analysis of raw reads from folder [0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline) and can be used to test the scripts in the detection scripts in [1_Detection](https://github.com/Easybel/DetectionGV/tree/main/1_Detection).
  - If you want raw reads to test the scripts in [0_WGSPipeline](https://github.com/Easybel/DetectionGV/tree/main/0_WGSPipeline), then check out our publication at [DOI: 10.1038/s41396-023-01440-x](https://www.nature.com/articles/s41396-023-01440-x). Here, the sequencing data is made available.
      
# How to cite our work:
- The bioinformatic scripts and basic pipeline for the detection of genomic variations were first developed for the work published here:

Jeffrey J. Power, Fernanda Pinheiro, Simone Pompei, Viera Kovacova, Melih Yüksel, Isabel Rathmann, Mona Förster, Michael Lässig, and Berenike Maier. **Adaptive evolution of hybrid bacteria by horizontal gene transfer.** _Proceedings of the National Academy of Sciences_, 2021, [DOI:10.1073/pnas.2007873118 ](https://www.pnas.org/doi/full/10.1073/pnas.2007873118)

- In the form presented here, the scripts were published in:

Isabel Rathmann, Mona Förster, Melih Yüksel, Lucas Horst, Gabriela Petrungaro, Tobias Bollenbach, and Berenike Maier. **Distribution of fitness effects of cross-species transformation reveals potential for fast adaptive evolution.** _The ISME Journal_, 2022, [DOI: 10.1038/s41396-023-01440-x](https://www.nature.com/articles/s41396-023-01440-x)

- and extensively used and improved for: 

Mona Förster, Isabel Rathmann, Melih Yüksel, Jeffrey J. Power, and Berenike Maier. **Genome-wide transformation reveals extensive exchange across closely related Bacillus species.** _BioXriv_, 2023, [DOI: 10.1101/2023.07.03.547483](https://biorxiv.org/cgi/content/short/2023.07.03.547483v1).

#

The original version of all scripts, as cited in the PhD thesis of Isabel Rathmann, is preserved on branch [PhDThesis_Rathmann2023](https://github.com/Easybel/DetectionGV/tree/PhDThesis_Rathmann2023).
