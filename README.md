# pileup_stranded

A common task in bioinformatics is pileup of alignment files for detecting muation, inserstion or deletions sites.

The current tool provide a site specific pileup with the following advantages comparing to other tools:

1. Using mutliproccesing to decreas the time of analysis.

2. Allow strand specific pileup for RNA-seq data.

3. Allow selection of specific nucleotides to include in the analysis in order to decres the runing time and size of output files.

4. Allow filtering of the results based on coverate, muation rate or number of mutated reads.

### Requirements:

* Python3
* pip
* samtools


#### python moduls (in requirements.txt):

* pysam
* pandas
* biopython

### Installation 

```
git clone https://github.com/BioInfoData/pileup_stranded
cd pileup_stranded
pip install -r requirements.txt
```

### Usage:

#### Non-stranded data

For pileup ignoring strand use the **pileup_multi.py** tool 
with following parameters:

```
python pileup_multi.py --help

options:                                                                                                                             
  -h, --help            show this help message and exit                                                                              
  --bam BAM             path to bam file                                                                                             
  --genome GENOME       fasta file of the genome                                                                                     
  --sites_file SITES_FILE                                                                                                            
                        list of sites to check                                                                                       
  --out_file OUT_FILE   name of output file                                                                                          
  --nt NT               nucleotides to include in pileup. Default is all nucleotides ("aAgGcCtT")                                    
  --min_coverage MIN_COVERAGE                                                                                                        
                        min coverage.Default is 0.                                                                                   
  --min_mutation MIN_MUTATION                                                                                                        
                        min mutation rate.Default is 0.                                                                              
  --min_nonRef MIN_NONREF
                        min number of mutated reads.Default is 0.

```

For example, to pileup the test bam file in the folder **test_pileup** with 5 processors
use the following command (*hg19.fa* file should be downloaded form UCSC):

```commandline
python pileup_multi.py --sites_file test_pileup/sites_file.txt --bam test_pileup/sampled.bam --genome *hg19.fa* --num_threads 5

```

To pileup the test bam file in the folder **test_pileup** but only reporting base G positions with min coverage of 10 reads, min mutation rate of 50% and min 3 mutated reads at the position use the following command:

```
python pileup_multi.py --sites_file test_pileup/sites_file.txt --bam test_pileup/sampled.bam --genome *hg19.fa* --out_file test_pileup/out_pileup_sites.txt --num_threads 5 --min_coverage 5 --min_nonRef 3 --min_mutation 0.5 --nt Gg

```
