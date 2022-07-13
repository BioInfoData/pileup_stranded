import sys
import pandas as pd
import general
import pileup_multi
import os
import logging
import argparse


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHeandeler = logging.StreamHandler()
logger.addHandler(streamHeandeler)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHeandeler.setFormatter(formatter)



"""

mode1:  This will pileup all genome. 

mode2: set specific_sites = True. Give site files for pos and ned. This will go over specific sites 
List of sites must be a single nt positions.
Expect tsv file with header. first col is chr second is position. For example:
################
# chr   mut_site
# chrY  3921370
# chrY  3921400
# chrX  5734067
# chr1  175493476
#################
"""

def strand_specific_pileup(num_threads, bam, genome, out_file, out_folder, writeOutFile,
                           min_coverage, min_mutation,  min_nonRef,
                           nt_pos="aA", nt_neg = "tT",pos_strand = "R2",
                           make_bam = True,
                           mode = "genome", # should be "genome" or "sites"
                           chr_list = None,  # mode1. This should be file with chr size. Not requried. Othrewise will calculate site
                           sites_file_pos = None, sites_file_neg = None): # mode2. need to include sites files


    name_bam = os.path.basename(bam).split(".")[0]

    pos_bam = "{}/{}_pos.bam".format(out_folder,name_bam)
    neg_bam = "{}/{}_neg.bam".format(out_folder,name_bam)

    if (make_bam):
        general.split_bam(bam,pos_bam, neg_bam,pos_strand)
        general.index_file(pos_bam)
        general.index_file(neg_bam)
    else:
        general.check_file(pos_bam)
        general.check_file(neg_bam)

    logger.info("###### positive strand ####")
    if mode == "genome": # this pileup all genome
        logger.info("All genome mode")
        pos_pileup = pileup_multi.main_make_pileup_genome(pos_bam, genome, "tmp/pos_pileup.csv",  nt_pos, num_threads,
                                                    min_coverage, min_mutation,  min_nonRef,
                                                    chr_list, writeToFile= False)
        pos_pileup["strand"] = "+"

        logger.info("###### negative strand ####")
        neg_pileup = pileup_multi.main_make_pileup_genome(neg_bam, genome, "tmp/neg_pileup.csv", nt_neg, num_threads,
                                                    min_coverage, min_mutation, min_nonRef,
                                                    chr_list, writeToFile=False)
        neg_pileup["strand"] = "-"



    elif mode == "sites": # this pileup specific sites.
        logger.info("List sites mode")
        logger.info("###### positive strand ####")
        logger.info("sites file pos {}".format(sites_file_pos))
        pos_pileup = pileup_multi.main_make_pileup_sites(pos_bam, genome,  sites_file_pos, nt_pos,
                                                   "", num_threads,
                                                   min_coverage, min_mutation, min_nonRef,
                                                   writeToFile=False)
        pos_pileup["strand"] = "+"

        logger.info("###### negative strand ####")
        logger.info("sites file neg {}".format(sites_file_neg))
        neg_pileup = pileup_multi.main_make_pileup_sites(neg_bam, genome,  sites_file_neg, nt_neg,
                                                   "", num_threads,
                                                   min_coverage, min_mutation, min_nonRef,
                                                   writeToFile=False)
        neg_pileup["strand"] = "-"

    else:
        logger.critical("mode must be 'genome' or 'sites'")
        sys.exit(1)



    both_strand_pileup = pd.concat([pos_pileup,neg_pileup])
    if writeOutFile:
        both_strand_pileup.to_csv("{}/{}".format(out_folder,out_file), index=False)
    return(both_strand_pileup , pos_bam, neg_bam)




def parse_user_data():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', help='path to bam file', required=True)
    parser.add_argument('--genome', help='fasta file of the genome', required=True)
    parser.add_argument('--out_file', help='output file name (WITHOUT path)', required=True)
    parser.add_argument('--out_folder', help='output folder path to write output files', required=True)
    parser.add_argument('--mode', help='should be "genome" or "sites"', required=True)

    parser.add_argument('--pos_strand', help='Reads on positive strand must be "R1" or "R2"'
                                             'Default is "R2"',  required=False, default="R2")

    parser.add_argument('--chr_list', help='File of chr length. if mode == "genome".'
                                           'Not required. Otherwise will calculate size', required=False, default= None)

    parser.add_argument('--sites_file_pos', help='List of sites for pileup on positive strand. '
                                                 'Only if mode == "sites"', required=False)
    parser.add_argument('--sites_file_neg', help='List of sites for pileup on positive strand. '
                                                 'Only if mode == "sites"', required=False)

    parser.add_argument('--make_bam', help='If split bam not exists set this param to "False". '
                                           'Default is "True"', required=False, default=True)

    parser.add_argument('--min_coverage', help='min coverage', required=False, default=0)
    parser.add_argument('--min_mutation', help='min mutation rate', required=False, default=0)
    parser.add_argument('--min_nonRef', help='min number of mutated reads', required=False, default=0)
    parser.add_argument('--num_threads', help='number of threads to use', required=False, default=1)

    parser.add_argument('--nt_pos', help='nucleotides to include in pileup - positive strand. '
                                     'Default is all nucleotides ("aAgGcCtT")', required=False, default="aAgGcCtT")
    parser.add_argument('--nt_neg', help='nucleotides to include in pileup - negative strand. '
                                     'Default is all nucleotides ("aAgGcCtT")', required=False, default="aAgGcCtT")


    args = parser.parse_args()


    writeOutFile = True

    strand_specific_pileup(args.num_threads,args.bam, args.genome, args.out_file, args.out_folder, writeOutFile,
                           int(args.min_coverage), float(args.min_mutation), int(args.min_nonRef),
                           args.nt_pos, args.nt_neg, args.pos_strand,
                           args.make_bam,
                           args.mode,  # should be "genome" or "sites"
                           args.chr_list, # mode1. This should be file with chr size. Not requried. Othrewise will calculate site
                           args.sites_file_pos, args.sites_file_neg)  # mode2.

    logger.info("Finished successfully. Results are at: {}/{}".format(args.out_folder,args.out_file))

if __name__ == '__main__':
    parse_user_data()

'''
############## test code ##############

#
python stranded_pileup.py --bam test_pileup/sampled.bam --genome /nadata/users/reut/reference/GRCh37/hg19.fa \
  --out_file out_pileup_stranded.txt --out_folder test_pileup/ --mode genome --min_coverage 10 --min_nonRef 4 \
  --min_mutation 0.2 --num_threads 5 --nt_pos Gg --nt_neg Cc


python stranded_pileup.py --bam test_pileup/sampled.bam --genome /nadata/users/reut/reference/GRCh37/hg19.fa \
 --out_file out_pileup_stranded.txt --out_folder test_pileup/ --mode sites --sites_file_pos test_pileup/sites_file.txt \
 --sites_file_neg test_pileup/sites_file.txt --min_coverage 10 --min_nonRef 4 -\
 -min_mutation 0.2 --num_threads 5 --nt_pos Gg --nt_neg Cc


#####################################
'''

'''

bam = "/nadata/users/reut/raw_data/ASO_runNOVA4/bam_sorted/Liora_SPH_19_S19_Aligned_onlyPaired_collapsed_sorted.bam"
genome = '/nadata/users/reut/reference/hg38/STAR_human/hg38_ucsc.fa'
sites_file_pos =  "/nadata/users/reut/raw_data/ASO_runNOVA4/test_ASE/SNP_sites.txt"
sites_file_neg = "/nadata/users/reut/raw_data/ASO_runNOVA4/test_ASE/SNP_sites.txt"
out_folder = "/nadata/users/reut/raw_data/ASO_runNOVA4/test_ASE/"
out_file = "strand_speciifc"
strand_specific_pileup(bam, genome, out_file, out_folder, writeOutFile = True,
                           min_coverage = 0, min_mutation = 0,  min_nonRef =0,
                           nt_pos="aA", nt_neg = "tT",pos_strand = "R2",
                           make_bam = False,
                           mode = "sites",
                           chr_list = None,
                           sites_file_pos = sites_file_pos, sites_file_neg = sites_file_neg)

'''
