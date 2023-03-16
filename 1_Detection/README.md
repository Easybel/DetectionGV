
---------------------------------------------------
## Overview of the scripts 
In this folder, the scripts are provided with which different genetic variants can be detected. All scripts are implemented in MATLAB. The numbers of the scripts indicate the order in which they should be run. Here, we give a brief overview.

**A0_VariantFiltering**
- This script takes the list of called variants and filters by quality and read depth. This script is also used to filter the masterlist, which contains all SNPs detected between a recipient and donor species.

**A0c_AccessoryGenome**
- Here, we detect the accessory genome parts of one genome compared to a second genome.

**A0d_Multimapper**
- This script allows us to find the regions on a genome to which multimapping reads align.

**A1_SNP2CNP** 
- With the filtered variant list and the masterlist, we detect the replacements of donor segments into the recipient. Here, we exclude accessory genome and multimapping regions.

**A2_Lists2Genes** (prior A2_CNP2Genes and A4_Regions2Genes)
- input
- output
  - here now you can also upload deldups or a randomn txt file 
  - you can upload a txt file with BSU names of genes that you are interested in to search for (if you dont want to go through the whole bed file with all 4422 genes ..)
  - you can also exclude the accessory genome now
 
**A3_Cov2DelDup** 
- input
- output

---------------------------------------------------
