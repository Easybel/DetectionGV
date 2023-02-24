old :
There are different functions for the different analysis tasks and they are called from a script called 
**Script_CAna.m**
which answers the most basic and interesting tasks regarding the detection of clusters and the statisics and visualisation around this.
**Attention**: Script_CAna.m is not up to date.

---------------------------------------------------
## TheMaster Plan
.... What we have to do

- [x] **A1_SNP2CNP** 
- [x] **A2_Lists2Genes** (prior A2_CNP2Genes and A4_Regions2Genes)
  - [x] here now you can also upload deldups or a randomn txt file 
  - [x] you can upload a txt file with BSU names of genes that you are interested in to search for (if you dont want to go through the whole bed file with all 4422 genes ..)
  - [x] you can also exclude the accessory genome now
- **A2b_Gene2Fcn** - Beta version available - (Mona)
- **A3_CNP2Operons** (Isa + Melih)
- [x] **A5_Cov2DelDup** 
- [ ] **A6_CNP2AnchorRegions** - Analysis of sequence identity of potential recombination starting points (Mona) 

---------------------------------------------------

- **Plot** - Overview -> del&dup dazu (Isabel)
  - [x] Plotscript (Mona)
  - [x] P_JJPBoxplot

---------------------------------------------------

- **S1_ClusterStat** (isa)

- **S2_GeneStat** (isa)

- **S3_OperonStat** (melih)

- **S4_CatStat** (Mona) - included in A2b - Gene2fcn.m

- **S5_DelDupStat**

- **S6_DenovoStat** 
  - deletions connected to denovos??
  - duplications connected to denovos?
  - recombinations connected to denovos

----------------------------------------------------
- **Mock genome** 
  - [x] Recombination W23 (Mona) 
  - [x] Recombination Bval (Mona) 
  - [x] Deletions (Isa) 
  - [x] Duplications (Mona) 

-----------------------------------------------------
- **W23 vergleichen mit Bval**
  - masterlists vergleichen, acc vergleichen
  -> die beiden aufeinander mappen: ml, acc
  - acc. genome von W23 gegen (Bsub168 && Bval)
  - acc. genome von Bval gegen (Bsub168 && W23)  

------------------------------------------------------

- pseudogene in Bsub168 suchen - sind die bei euns auch pseudogene?? 
  - [x] (Isa) Result is in sciebo pseudogenes_IR
