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
- [X] intrinsic feature: di-nucleotide analysis


Dani 

- [X] Review correlated features
- [ ] Finish readme files and project setup
- [ ] Conda env file 
- [ ] Add filter for protein-coding for high-express genes used for RNA:RNA interactions
- [X] Intrinsic feature: sequence low complexity analysis --> Good preliminary results --> Incude it as a intrinsic feature
- [X] change CPC2 for tcode --> tcode uses a 200 nt window --> removed
- [ ] add "distance" label to short-ncrna-negative-control 

## IDEAS

Dani

- [ ] snakemake script for analysis 
- Extract gnomad data directly from a server ? (902G)

## PROGRESS 

* 2024/04/05 - Dani: 

                    Divided sequence-processing.sh script into 3 (one for each type of RNA) so they can be run in parallel. 

* 2024/04/05 - Dani: 

                    Short-ncRNA dataset had many repeated sequences due to duplicated, but different ID, sequences in RNAcentral dataset --> Removed them from RNAcentral dataset and extracted new sequences 

* 2024/04/09 - Dani: 

                    1000 functional sequences with negative controls ready to analyse  --> new order of the dataset: by chr, start, end   

* 2024/04/10 - Dani: 

                    Remove Fickett score (CPC2): we're already calculating the coding potential with RNAcode --> DJ work showed better results with RNAcode 
                    
                    Remove Accessibility (RNApfold): feature showed perfect spearman correlation with MFE (rho=-1) --> doesn't meet feature criteria: redundant
                    
                    Remove phasCons score: redundant with phyloP
                    
                    Add low complexity sequence density feature to intrinsic-features.sh script 
              
* 2024/04/12 - Dani: 
                    
                    update population script to remove the vcf files for each sequences after using it, since it's the process of extracting them from the local database is quite quick and it's not worth to store them in case of needed again (already storing the database)

* 2024/04/19 - Dani: 

                    update the way the gene_complement was being calculated --> now Gencode includes not only genes. 

                                                                                Gencode and RNAcentral are being merged before finding the complement (instead of finding complement of both separately and join afterwards)
                    change the steps used to sample the negative control: 
                        
                        The distance between "negative control exons" now is set at 100 (arbitrary) instead of the distance between the "functional exons" --> This was done because many exons distances were quite big respect the size of the "complement regions" --> many of this regions were not used reducing considerably the number of negative samples 

* 2024/05/27 - Dani:

					change how to extract the sequences from the genome: before --> grep on .csv; now bedtools getfasta --> much faster. Adapted scripts. Now the strand is considered as well: takes inverse complementary for negative strand. 
					
					The genome sequence now is downloaded from UCSC --> fasta file header chr1 which matches to the chr column we were using --> requirement for bedtools getfasta 

* 2024/05/28 - Dani: 

					add filter to exclude ncRNA which are included in the curated database use to calculate the interaction RNA:RNA (lncRNA and shor-ncRNA sequences script)

                    


