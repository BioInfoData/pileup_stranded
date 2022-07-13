# pileup_stranded

A common task in bioinformatics is pileup of alignment files for detecting muation, inserstion or deletions sites.

The current tool provide a site specific pileup with the following advantages comparing to other tools:

1. Using mutliproccesing to decreas the time of analysis.

2. Allow strand specific pileup for RNA-seq data.

3. Allow selection of specific nucleotides to include in the analysis in order to decres the runing time and size of output files.

4. Allow filtering of the results based on coverate, muation rate or number of mutated reads.

### Requirements:

Python3

pysam

pandas

biopython

### Usage:

```
python pileup.py

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

