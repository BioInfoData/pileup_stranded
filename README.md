# pileup_stranded

A common task in bioinformatics is pileup of alignment files for detecting muation, inserstion or deletions sites.

The current tool provide a site specific pileup with the following advantages comparing to other tools:

1. Using mutliproccesing to decreas the time of analysis.

2. Allow strand specific pileup for RNA-seq data.

3. Allow selection of specific nucleotides to include in the analysis in order to decres the runing time and size of output files.

4. Allow filtering of the results based on coverate, muation rate or number of mutated reads.

6. Allow pileup only for list of sites provided as input file to the tool. 

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
usage: pileup_multi.py --bam BAM --genome GENOME --out_file OUT_FILE 

options:
  -h, --help            show this help message and exit
  --bam BAM             path to bam file
  --genome GENOME       fasta file of the genome
  --sites_file SITES_FILE
                        List of sites to check. If not provided, the script
                        will run on all genome
  --out_file OUT_FILE   name of output file
  --nt NT               nucleotides to include in pileup. Default is all
                        nucleotides ("aAgGcCtT")
  --min_coverage MIN_COVERAGE
                        min coverage.Default is 1.
  --min_mutation MIN_MUTATION
                        min mutation rate.Default is 0.
  --min_nonRef MIN_NONREF
                        min number of mutated reads.Default is 0.
  --num_threads NUM_THREADS
                        number of threads to use

```

For example, to pileup the test bam file in the folder **test_pileup** only on sites listed at **sites_file.txt** using 5 processors
use the following command (*hg19.fa* file should be downloaded form UCSC):

```commandline
python pileup_multi.py --sites_file test_pileup/sites_file.txt --bam test_pileup/sampled.bam \
--genome hg19.fa --out_file test_pileup/out_pileup_sites.txt --num_threads 5 

```

To pileup the test bam file in the folder **test_pileup** but only reporting base G positions with min coverage of 10 reads, min mutation rate of 50% and min 3 mutated reads at the position use the following command:

```
python pileup_multi.py --sites_file test_pileup/sites_file.txt --bam test_pileup/sampled.bam \
--genome hg19.fa --out_file test_pileup/out_pileup_sites.txt \
--num_threads 5 --min_coverage 5 --min_nonRef 3 --min_mutation 0.5 --nt Gg

```

#### Strand-specific pileup

For strand-specific pileup use the **stranded_pileup.py** tool 
with following parameters:

```
usage: stranded_pileup.py --bam BAM --genome GENOME --out_file OUT_FILE
                          --out_folder OUT_FOLDER --mode MODE

options:
  -h, --help            show this help message and exit
  --bam BAM             path to bam file
  --genome GENOME       fasta file of the genome
  --out_file OUT_FILE   output file name (WITHOUT path)
  --out_folder OUT_FOLDER
                        output folder path to write output files
  --mode MODE           should be "genome" or "sites"
  --pos_strand POS_STRAND
                        Reads on positive strand must be "R1" or "R2"Default
                        is "R2"
  --chr_list CHR_LIST   File of chr length. if mode == "genome".Not required.
                        Otherwise will calculate size
  --sites_file_pos SITES_FILE_POS
                        List of sites for pileup on positive strand. Only if
                        mode == "sites"
  --sites_file_neg SITES_FILE_NEG
                        List of sites for pileup on positive strand. Only if
                        mode == "sites"
  --make_bam MAKE_BAM   If split bam not exists set this param to "False".
                        Default is "True"
  --min_coverage MIN_COVERAGE
                        min coverage
  --min_mutation MIN_MUTATION
                        min mutation rate
  --min_nonRef MIN_NONREF
                        min number of mutated reads
  --num_threads NUM_THREADS
                        number of threads to use
  --nt_pos NT_POS       nucleotides to include in pileup - positive strand.
                        Default is all nucleotides ("aAgGcCtT")
  --nt_neg NT_NEG       nucleotides to include in pileup - negative strand.
                        Default is all nucleotides ("aAgGcCtT")


```

For example, to pileup the test bam file in the folder **test_pileup** using 5 processors, reporting only G bases in the positive strand and only C bases at the negative strand and filtering for sites with min mutation rate of 20%, min coverage of 10 reads and min mutated reads of 4 reads, use this command:

```
python stranded_pileup.py --bam test_pileup/sampled.bam --genome hg19.fa \
  --out_file out_pileup_stranded.txt --out_folder test_pileup/ --mode genome --min_coverage 10 --min_nonRef 4 \
  --min_mutation 0.2 --num_threads 5 --nt_pos Gg --nt_neg Cc
```

For pileup with the same parameters this time only using list of sites from **test_pileup/sites_file.txt** use the following command:

```
python stranded_pileup.py --bam test_pileup/sampled.bam --genome hg19.fa \
 --out_file out_pileup_stranded.txt --out_folder test_pileup/ --mode sites --sites_file_pos test_pileup/sites_file.txt \
 --sites_file_neg test_pileup/sites_file.txt --min_coverage 10 --min_nonRef 4 -\
 -min_mutation 0.2 --num_threads 5 --nt_pos Gg --nt_neg Cc
```
