# Method and Code Outline

### Generation of the dataset

1. Download ncRNA from [RNAcentral](https://rnacentral.org/search?q=TAXONOMY:%229606%22%20AND%20expert_db:%22HGNC%22), filtering for _Homo sapiens_ and HGNC database, in the FASTA format. Also download the corresponding [chromosome coordinates file](ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/). Don't use lncRNA, rRNA or precurser RNA to generate the training dataset. 

ncRNA:
https://rnacentral.org/search?q=HGNC%20AND%20NOT%20rna_type:%22lncRNA%22%20%20AND%20NOT%20rna_type:%22rRNA%22%20%20AND%20NOT%20rna_type:%22precursor%20RNA%22


ncRNA chromosomal coordinates:
https://rnacentral.org/search?q=HGNC%20AND%20NOT%20rna_type:%22lncRNA%22%20%20AND%20NOT%20rna_type:%22rRNA%22%20%20AND%20NOT%20rna_type:%22precursor%20RNA%22%20AND%20has_genomic_coordinates:%22True%22


2. Download the GRCh38 genome from NCBI in FASTA format, and their corresponding BED files. Then reformat to a csv so that it can be easily read in python. The scaffolds that aren't associated with a main chromosome should also be removed.

```
#Convert from FASTA to tabular
fasta_formatter -i GRCh38_latest_genomic.fna -o GRCh38_part1.csv -t

#Remove scaffolds
grep "NC_" GRCh38_part1.csv > GRCh38_genome.csv
```

3. Using the files downloaded from RNAcentral as input, run in **RNAcentral.sh** to obtain a CSV file with the sequences for the training and control datasets, as well as the CM and e-value scores from CMscan.

* [**RNAcentral.sh:**](RNAcentral.sh) 
  * **Arguments:** $1 is the downloaded FASTA file from RNAcentral and $2 is the name you want to give the dataset (eg: 181203dataset).
  * **Additional Scripts:** [coordinates2.py](coordinates2.py), [RfamCM.py](RfamCM.py) and [random_sequences.py](random_sequences.py).
  * **Additional Files:** Homo_sapiens.GRCh38.bed (chromosome coordinates) and GRCh38_genome.csv (reformatted genome). NB: These files need these names
  * **Additional Functions:** [fasta_formatter](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html) and [cmscan](http://eddylab.org/infernal/).
  * **Output:** $2_final.csv is the converted dataset for RNAcentral data and $2_spare is the null dataset based off the RNAcentral data.

### Conservation of Sequence

1. The file snp151_trimmed_sorted.bed is provided but the code for generating it is as follows: 
Download snp151.txt.gz from ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/

```
gunzip snp151.txt.gz
grep  -v chr[A-Z]_* I think >snp151_no_alts
grep -v Chromosome rnacentraloutput_FINAL.csv |  cut -f 3,4,5 -d "," | tr ',' '\t' >rnacentraloutput_FINAL.bed
cut -f 2,3,4 snp151_no_alts >snp151_trimmed.bed
sort -k1,1 -k2,2n snp151_trimmed.bed > snp151_trimmed_sorted.bed 
sort -k1,1 -k2,2n rnacentraloutput_FINAL.bed > rnacentraloutput_FINAL_sorted.bed 

awk '$1="chr"$1' rnacentraloutput_FINAL_sorted.bed | tr ' ' '\t' >rnacentraloutput_FINAL_sorted1.bed 
bedtools intersect -c -a rnacentraloutput_FINAL_sorted1.bed -b snp151_trimmed_sorted.bed -wa -sorted >rnacentraloutput_snp_intersection
```

2. To obtain the phastCons and phyloP conservation scores, use the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=701273921_eML3CJEnXm4rsnQmqAuYnNHkZVkV) and grab the columns of interest (maximum and average). **Note:** doing this is much faster and easier than downloading them via mysql.

```
#Reformat chromosome coordinates for input
grep -v 'Chromosome' file.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}'

#Go onto UCSC table browser and input the following
group: Comparitive Genomes
track: Conservation
table: Cons 100 Verts (phyloP100 way) or Cons 100 Verts (phastCons100way)
region: define regions > submit chromosome coordinates (note only 1000 at a time)
button to click: Summary/Statistics
```
This will not return all results, so use ##script name to join the returned SNP statitistics to their corresponding IDs
### Population Statistics

1. Use **extractVCF.sh** to download VCF files (if available) from Phase 3 of the 1000 Genomes Project and puts them in a folder named FASTA (folder should be called this to make PopGenome run easier in R).

* [**extractVCF.sh:**](extractVCF.sh)
  * **Argument:** $1 is the name of the dataset file (including location to the file). 
  * **Additional Functions:** [tabix](http://www.htslib.org/doc/tabix.html).
  * **Output:** VCF files for each ncRNA.

2. Go into R and obtain the population statistics using PopGenome.

```
library(PopGenome)
genome.class <- readData("FASTA", format="VCF")
genome.class <- neutrality.stats(genome.class)
dataset <- get.neutrality(genome.class)[[1]]
write.csv(dataset, "name.csv")     
``` 

3. Attach the PopGenome output to the current dataset file, as well as calculate the snp average.

```
Pretty sure I calculated this in the excel file so I haven't got command line code for this.
snp average is number of snp divided by the ncRNA sequence length.
All NAs that were generated by PopGenome should also be changed to zeros.
```

4. Use **filter1000genomes.sh** to calculate the number of transitions, transversions and relevant MAF in the dataset using the downloaded VCF files.

* [**filter1000genomes.sh:**](filter1000genomes.sh) 
  * **Arguments:** $1 is the dataset file and location and $2 is the name of the updated file (without file extension).
  * **Additional Files:** VCF files in FASTA.
  * **Additional Functions:** [vcftools](http://vcftools.sourceforge.net/)
  * Transition:Transversion Ratio = (transition +1)/(transition + transversion + 2)
  * **Output:** $2.csv is the original dataset with the new data added in.

### Conservation of Secondary Structure

1. Use **calculate_CM.sh** to obtain the multiple sequence alignment files from multiz100way (UCSC). These files are then used to calculate CM and HMM using cmbuild.

* [**calculate_CM.sh:**](calculate_CM.sh) 
  * **Arguments:** $1 is the name of the dataset file (including location to file) and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Scripts:** [conversion3.py](conversion3.py)
  * **Additional Functions:** [mafFetch](http://hgdownload.soe.ucsc.edu/admin/exe/), [RNAalifold](https://www.tbi.univie.ac.at/RNA/), [R-scape](http://eddylab.org/R-scape/) and [cmbuild](http://eddylab.org/infernal/).
  * **Output:** $2.csv is the original dataset with the new data added in.
  
### RNA:RNA Interactions

1. Using NCBI Gene, download the text file for the "most relevant" promoters for protein coding genes, "most relevant" mRNA and "most relevant" ncRNA (eg: type ncRNA into NCBI Gene and download the results, alternatively the links for these are in the script filtergenes.sh). Using **filtergenes.py**, extract the sequences for the first 250 sequences for each of these groups.
To geenerate the random sequences from your dataset, this python code can be run:

```
from Bio import SeqIO
from random import sample
import sys
with open(sys.argv[1]) as f:
    seqs = SeqIO.parse(f, "fasta")
    samps = ((seq.name, seq.seq) for seq in  sample(list(seqs),250))
    for samp in samps:
        print(">{}\n{}".format(*samp))
```

* [**filtergenes.sh:**](filtergenes.py) 
  * **Note:** this code uses a lot of memory.
  * **Arguments:** $1 is the downloaded text file from NCBI, $2 is the name of your output (fasta).
  * **Additional Files:** GRCh38_genome.csv.
  * **Output:** fasta file containing sequence for the first 250 genes.

2. Generate the test interaction database using the 750 genes from NCBI plus 250 random ncRNA from the test dataset. Will need to cat the four fasta files together into one before creating the interaction database.

```
cat promoter.fa protein-coding.fa ncrna-ncbi.fa ncrna-dataset.fa > interactionstest.fa 
./risearch2.x -c interactionstest.fa -o interaction.suf
```

3. Use **risearch.sh** to calculate the number of interactions between the ncRNA dataset and the test interaction database.

* [**risearch.sh:**](risearch.sh)
  * **Arguments:** $1 is the name of the dataset file (including location to file) and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Scripts:** [average.py](average.py)
  * **Additional Files:** interaction.suf
  * **Additional Functions:** [RIsearch2](https://rth.dk/resources/risearch/)
  * **Output:** $2.csv is the original dataset with the new data added in.

### Genomic Copy Number

1. Create human genome database that can be used by command line blastn. First you need to download [GCF_000001405.38_GRCh38.p12_genomic.fna.gz](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12) and unzip the file. Then format this file into a BLAST database.

```
makeblastdb -in GCF_000001405.38_GRCh38.p12_genomic.fna -dbtype nucl -parse_seqids -out human_genome
```

2. Use **filterblastn.sh** to take the ncRNA and blast them against the GRCh38.p12 reference genome.

* [**Filterblastn.sh:**](filterblastn.sh)
  * **Arguments:** $1 is the name of the dataset file (including location to file) and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Functions:** [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
  * **Additional Files:** human_genome.*
  * **Output:** $2.csv is the original dataset with the new data added in.

### Transcription

1. Download cell line and differentiated cell total RNA-seq data from ENCODE. Note that Dataset and BAM file IDs for the same sample are different.

* Cell Lines: H7-hESC ([ENCFF089EWC](https://www.encodeproject.org/files/ENCFF089EWC/)), HepG2 ([ENCFF067CVP](https://www.encodeproject.org/files/ENCFF067CVP/)), GM12878 ([ENCFF893HSY](https://www.encodeproject.org/files/ENCFF893HSY/)), K562 ([ENCFF796BVP](https://www.encodeproject.org/files/ENCFF796BVP/)).
* Differentiated Cells: smooth muscle cell from H9 ([ENCFF369DYD](https://www.encodeproject.org/files/ENCFF369DYD/)), hepatocyte from H9 ([ENCFF907AIK](https://www.encodeproject.org/files/ENCFF907AIK/)), neural progenitor cell from H9 ([ENCFF766ESM](https://www.encodeproject.org/files/ENCFF766ESM/)), myocyte from LHCN-M2 ([ENCFF722EAR](https://www.encodeproject.org/files/ENCFF722EAR/)), bipolar neuron from GM23338 ([ENCFF713UNS](https://www.encodeproject.org/files/ENCFF713UNS/)), myotube from skeletal muscle myoblast ([ENCFF178TTA](https://www.encodeproject.org/files/ENCFF178TTA/)), hematopoietic multipotent progenitor cell ([ENCFF065MVD](https://www.encodeproject.org/files/ENCFF065MVD/)) and cardiac muscle from RUES2 ([ENCFF475WLJ](https://www.encodeproject.org/files/ENCFF475WLJ/)).

2. Index the downloaded bam files.

```
samtools index ENCFF089EWC.bam
```

3. Use **expencode.sh** to extract the number of reads for the chromosome coordinates of each ncRNA.

* [**expencode.sh:**](expencode.sh)
  * **Arguments:** $1 is the dataset (includes file location), $2 is the name of the output. 
  * **Additional Files:** the bam.bai files that have been generated for each BAM file.
  * **Additional Functions:** [samtools](http://www.htslib.org/)
  * **Output:** $2.csv is the original dataset with the new data added in.

### RandomForest Analysis

1. Use **randomforest.R** to run RandomForest for x number of times and exports the importance, error and prediction values for each run.

* [**randomforestx.R:**](randomforest.R)
  * **Arguments:** location of original data, name of importance file, name of errors file, name of prediction file, number of time to repeat random forest.
  * **Output:** Three CSV files corresponding to the importance of predictors across x models, error rates for predicting functionality across x models and the predictions made summed for all x models.

2. Calculate the mean of each of these values (I did this in the original csv file but this could probably be incorporated into the R script) and input back into R to graph the results. 

3. Graph the results in R using ggplot2. Can also create a [correlation plot](prettyCorrelation.R) of all the raw predictor data.

```
#Example: Plot of mean importance values; geom_hline should correspond to the "Start" predictor which represents neutral. Make sure the data is ordered from smallest to largest so that the plot is ordered logically.

ggplot(RF1, aes(x=Importance, y=Predictor))+geom_point()+ylab('')+xlab('Mean Variable Importance for 100 runs')+theme(text=element_text(size=18))+geom_line()+geom_hline(yintercept=17,color="Red")
```

```
#Example: Plot of ncRNA for specific predictors as a boxplot and scatter, which is seperated into functional (blue) and non-functional (red) ncRNA.

allData <- read.csv(“190131nullRF1.csv, header=TRUE)
fun <-data.frame(allData$Functionality, stringsAsFactors=FALSE)
CMs <- scale(allData$CM.scan)
test1<-cbind(fun,CMs)
CMb <- scale(allData$CM.build)
test14 <- cbind(fun,CMb)
test1<-cbind(fun,CMs)
test1$newcolumn <- "CM (CMscan)"
test14$newcolumn <- "CM (CMbuild)"
colnames(test1) <- c("Functionality", "Zscore","Predictor")
colnames(test14) <- c("Functionality", "Zscore","Predictor")
CMcompare <- rbind(test1,test14)
CMcompare$Predictor <- factor(CMcompare$Predictor,  levels=unique(CMcompare$Predictor))

qplot(data=CMcompare, x=Predictor, y=Zscore, geom=c("boxplot"), outlier.shape = NA) + geom_jitter(shape = 1, aes(colour=Functionality))+scale_color_manual(values=c("indianred2","cyan3"))+theme(text=element_text(size=22))+ylim(low=-2.7,high=5)+theme(legend.position="none")
```

**OR**
  
1. Use **rftest.R** to run RandomForest for x number of times using only a subset of the dataset being used as the testData (ie: predict specific types of ncRNA).

* [**rftest.R:**](RFtest.R)
  * **Arguments:** location of original data that wasn’t part of the testData, location of data that is part of the testData (70% will also be included in the trainData), name of prediction file, number of times to repeat random forest.
  * **Output:** One CSV file corresponding to the summed predictions made for all x models.
