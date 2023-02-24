# HOW TO USE SNP2CNP.m :trollface:

SNP2CNP finds orthologously recombinated clusters and denovo mutations of an evolved strain by analysing the detected single nucleotide polymorphisms (SNPs). To be part of a cluster, a SNP needs to be on the "masterlist". More than $cms consecutively missing masterlist SNPs (mlSNPs) end a cluster. Accessory parts of the recipient with respect to the donor end a cluster as well. Multi mapper regions where no change can be detected are excluded. Clusters with length 1 (only containing one mlSNP) are called SPI (single polymorphism insert) and are saved seperately. Supplementary, the minimum, maximum and average lengths of the detected clusters are calculated. 

### Input:
- Masterlist (**$masterlist**): list of all single nucleotide differences in highly identicals parts (core genome) of recipient and donor with the recipient dictionary as reference; to find this, map the donor on the recipient dictionary.

- Genome size of the recipient strain (**$recipsize**)
- Accessory genome of the recipient with respect to the donor strain (**$accgenome**): parts of the recipient that are not present in the donor and thus cannot undergo orthologous recombination. **Core genome + accessory genome = whole genome**.
- Artefacts (**$artefacts**): list in vcf-format, containing artefacts found by mapping the ancestor reads on the recipient dictionary used.
- Multi mapper region list (**$mmlist**): list of multi mapping regions of the recipient strain; with start, end and length.
- OLD: Individual Mutations List (IndvMutList) of the evolved strains (**[$filepath, $filenames, $filesuffix]**): List of single nucleotide polymorphisms (SNPs) found in the evoved strains after mapping the Illumina short reads on the recipient dictionary.
- NEW: SNPSummary (**$SNPpath + $SNPname**): contains all Individual Mutations Lists of an experiment *plus* all information about the filtering 
- Number of allowed missing SNPs (**$cms**): Clusters end, when $cms SNPs can be found on the masterlist but are not present on the IndvMutList

### Output: 
#### OLD output, now all combined in the CNPSummary (see below)
- **denovo** contains all positions of mutations that have not been found on the masterlist.
- **C** contains all detected clusters described by (1) start position, (2) end position, (3) expected number of mlSNPs in the cluster and (4) number of missing mlSNPs in the cluster; (3) - (4) = number of detected SNPs in the cluster.
- **SPI** is an artefact (empty cell), since we allow cluster with just one SNP. #work on this !
- **Cdist** contains the **c**onservative length counting the length from the first SNP in the cluster to the last.
- **Mdist** contains the **m**aximum length, counting the length from the first mlSNP before the cluster starts to the first mlSNP after the cluster ends; expection: if an accessory part is closer to the cluster than the next mlSNP, the accessory part ends maximum length; in column 2 and 3, you can find the start and the end position of the mdist cluster
- **Adist** contains the **a**verage length, calculated as the mean of the difference between **Cdist** and **Mdist** (Adist = Cdist + (Mdist - Cdist)/2), in column 2 and 3, the start and the end position of the adist cluster can be found
#### New Output
- **CNPSummary**: new output-format. contains all information produced in the script: C, Cdist, Adist, Mdist, denovo, SPI, ORI_Crossing and the corresponding filename). 





