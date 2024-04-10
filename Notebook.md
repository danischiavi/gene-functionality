# PROJECT NOTEBOOK 

## TO DO LIST

- [X] Decide method for negative control 

- [X] Generate negative control dataset

- [ ] Run all features for all 1000 sequences 

- [ ] Manuscript

- [ ] Review copy-number justification as feature 


Estefi 

- [X] Negative control analysis 
- [X] Incorporate epigenetics marks analysis to the workflow 
- [ ] intrinsic feature: di-nucleotide analysis


Dani 

- [X] Review correlated features
- [ ] Finish readme files and project setup
- [ ] Conda env file 
- [ ] Add filter for protein-coding for high-express genes used for RNA:RNA interactions
- [X] Intrinsic feature: sequence low complexity analysis --> Good preliminary results --> Incude it as a intrinsic feature
- [X] change CPC2 for tcode --> tcode uses a 200 nt window --> 

## IDEAS

Dani

- [ ] snakemake script for analysis 


## PROGRESS 

* 2024/04/05: Divided sequence-processing.sh script into 3 (one for each type of RNA) so they can be run in parallel. 

* 2024/04/05: Short-ncRNA dataset had many repeated sequences due to duplicated, but different ID, sequences in RNAcentral dataset --> Removed them from RNAcentral dataset and extracted new sequences 

* 2024/04/09: 1000 functional sequences with negative controls ready to analyse  --> new order of the dataset: by chr, start, end   

* 2024/04/10: Remove Fickett score (CPC2): we're already calculating the coding potential with RNAcode --> DJ work showed better results with RNAcode 
              Remove Accessibility (RNApfold): feature showed perfect spearman correlation with MFE (rho=-1) --> doesn't meet feature criteria: redundant
              Remove phasCons score: redundant with phyloP
              Add low complexity sequence density feature to intrinsic-features.sh script 
              



