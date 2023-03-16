# About the scripts in this folder:

- **0_DicIndex.sh**
  - The dictionaries needed for the alignment steps are indexed and sorted.  

- **1_QCScript.sh**
  - Raw reads are trimmed and adapters are removed. The read quality is assessed before and after the trimming.

- **2_CombiScriptSamtools.sh**
  - The trimmed reads are aligned to the reference and variants are called. The coverage at each position is evaluated. 

