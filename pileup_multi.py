
import sys
import pysam
import logging
import pandas as pd
from Bio import SeqIO
import multiprocessing
import argparse
import general


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHeandeler = logging.StreamHandler()
logger.addHandler(streamHeandeler )
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHeandeler.setFormatter(formatter)

"""
Initiate dict for writing the results of one position pileup data
"""
def init_dict(one_chr_data,current_pos, chr_seq,chr):
    if one_chr_data[current_pos] == None:
        one_chr_data[current_pos] = {'chr': chr, 'mut_site' : current_pos+1, 'nt' : chr_seq[current_pos].upper(),
                                     'A': 0, 'T': 0, 'G': 0, 'C': 0,
                                     'delet': 0, 'insert': 0, 'reference_alignment': 0,
                                     'non_reference_alignment' : 0}
    return one_chr_data

"""
Parse each char in the pileup line:
reference_alignment
non_reference_alignment
deletions
"""
def parse_one_char(char,one_chr_data,current_pos):
    if char == '>' or char == '<':
        return one_chr_data[current_pos]
    elif char == '.' or char == ',':
        one_chr_data[current_pos]['reference_alignment'] += 1
    elif char in 'AGCTNagctn':
        one_chr_data[current_pos][char.upper()] += 1
        one_chr_data[current_pos]['non_reference_alignment'] += 1
    elif char == "*":
        one_chr_data[current_pos]['delet'] += 1
    return one_chr_data[current_pos]

"""
Parse one line of the pileup data = get one position data
"""
def parse_one_line_pileup(line_pileup,one_chr_data,current_pos, chr_seq):
            for char in line_pileup:
                if len(char) == 1:
                    one_chr_data[current_pos] = parse_one_char(char, one_chr_data, current_pos)
                elif len(char) > 1: # this means insertion or deletion of few nt starting from nighbering nt
                    multi_char = char
                    char = multi_char[0] # the fist char is the current nt
                    one_chr_data[current_pos] = parse_one_char(char, one_chr_data, current_pos)
                    sing = multi_char[1]
                    if sing == "+":  # means insertion
                        one_chr_data[current_pos]['insert'] += 1
                    # TODO This is probably not good
                    #elif sing == "-":  # means deletions
                    #    parse_multi_char_del(multi_char, one_chr_data, current_pos, chr_seq, chr)

def convert_to_df(final_one_chr_data):
            df = pd.DataFrame.from_dict(final_one_chr_data)
            df['coverage'] = df['non_reference_alignment'] + df['reference_alignment']
            df['mutation_rate'] = round(df['non_reference_alignment'] / df['coverage'],3)
            df['A'] = round(df['A']/ df['non_reference_alignment'],3)
            df['T'] = round(df['T'] / df['non_reference_alignment'],3)
            df['G'] = round(df['G'] / df['non_reference_alignment'],3)
            df['C'] = round(df['C'] / df['non_reference_alignment'],3)

            return df
"""
Filter the mutation sites according to user's parameters
"""
def filter_df(df, min_mutation, min_coverage, min_nonRef):
    df = df.loc[df['coverage'] >= min_coverage]
    df = df.loc[df['mutation_rate'] >= min_mutation]
    df = df.loc[df['non_reference_alignment'] >= min_nonRef]
    return df

"""
Generates pileup data for one chromosome. 
"""
def make_one_pileup(samfile,chr, start, end,chr_seq,fastfile,nt, one_chr_data):
    iter = samfile.pileup(chr,start, end, truncate =True, fastafile = fastfile, max_depth = 100000)
    for pileupcolumn in iter :
        current_pos = pileupcolumn.reference_pos
        if chr_seq[current_pos] in nt:
            one_chr_data = init_dict(one_chr_data,current_pos, chr_seq, chr)
            line_pileup = pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True)
            parse_one_line_pileup(line_pileup, one_chr_data, current_pos, chr_seq)

    return one_chr_data

"""
Core pileup function. 
Generates pipeup data for one chromosome, filter the data and fill the final results dict.
This can be all nt in the chromosome of specific list of sites in the same chromosome. 

Can run in 3 modes:

1. for all genome: get only chr name
2. for list of positions (single nt) get chr and list_sites
3. for specific area: get chr, start, end (only for debugging)
"""

def pileup_bam(bam, genome, nt, chr, start = None, end = None,
               min_coverage = 1, min_mutation = 0, min_nonRef = 0, return_dict =None,
               list_sites = None, sema = None):


    logger.info("analysing {}".format(chr))
    samfile = pysam.AlignmentFile(bam, "rb")
    fastfile = pysam.FastaFile(genome)

    chr_seq = fastfile.fetch(chr)
    one_chr_len = len(chr_seq)
    one_chr_data = [None]*one_chr_len

    if (not start) and (not end)  and (not list_sites): # all cromosome

        end = len(chr_seq)
        start = 0
        one_chr_data = make_one_pileup(samfile, chr, start, end, chr_seq, fastfile, nt, one_chr_data)
    elif list_sites: #list of sites in one chromosome

        for i in range(len(list_sites)):

            end = int(list_sites[i])
            start = end-1
            one_chr_data = make_one_pileup(samfile,chr, start, end,chr_seq,fastfile,nt, one_chr_data)

    elif start and end: # specific area
        end = end
        start = start
        one_chr_data = make_one_pileup(samfile, chr, start, end, chr_seq, fastfile, nt, one_chr_data)

    else:
        logger.critical("must give start and end, or only chr, or list of sites")
        sys.exit(1)

    samfile.close()

    final_one_chr_data = list(filter(None, one_chr_data))

    if len(final_one_chr_data) > 0:
        df = convert_to_df(final_one_chr_data )
        df_final = filter_df(df, min_mutation, min_coverage, min_nonRef)
        logger.info("finished {}".format(chr))
        return_dict[chr] = df_final

    else:
        logger.info("No site found in {}".format(chr))


    sema.release()

def get_chr_list(genome):
    chr_list = []
    with open(genome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            chr_list.append(record.id)
    return chr_list

"""
Mode 1 - all genome.
Run over each chromosome using multiprocessing.
"""
def main_make_pileup_genome(bam, genome, out_file = None, nt =  "aAgGcCtT", num_threads = 1,
                            min_coverage = 1, min_mutation = 0, min_nonRef = 0, chr_list = None,
                            writeToFile = True):

    if not chr_list:
        logger.info("Get chr length for {}".format(genome))
        chr_list = get_chr_list(genome)

    logger.info('starting pileup {}'.format(bam))
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    sema = multiprocessing.Semaphore(int(num_threads)) # limit number of processors)
    for chr_name in chr_list:
        sema.acquire()
        p = multiprocessing.Process(target=pileup_bam,
                                    args =(bam, genome, nt, chr_name,None,None,
                                           min_coverage, min_mutation, min_nonRef, return_dict,
                                           None, sema))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    list_df = list(return_dict.values())
    if len(list_df) > 0:
        all_df = pd.concat(list_df)
    else:
        return
    if writeToFile:
        all_df.to_csv(out_file, index=False )
    return all_df



"""
Mode 2 - list of sites.
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

def main_make_pileup_sites(bam, genome, sites_file, nt ="CcGgaAtT", out_file = None, num_threads = 1,
                           min_coverage=1, min_mutation=0, min_nonRef=0,
                           writeToFile = True):

    general.check_file(bam)
    general.check_file(genome)
    general.check_file(sites_file)

    site_dict = {}
    sites = pd.read_csv(sites_file, sep="\t")
    sites = sites.iloc[:,[0,1]]
    for i,row in sites.iterrows():
        chr = row[0]
        end = int(row[1])
        if chr in site_dict.keys():
            site_dict[chr].append(end)
        else:
            site_dict[chr] = []
            site_dict[chr].append(end)

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    sema = multiprocessing.Semaphore(int(num_threads))
    jobs = []
    for chrom in site_dict.keys():
        sema.acquire()
        p = multiprocessing.Process(target=pileup_bam, args=(bam, genome, nt, chrom, None, None,
                                                             min_coverage, min_mutation, min_nonRef,
                                                             return_dict,
                                                             site_dict[chrom],sema))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    list_df = list(return_dict.values())
    if len(list_df) > 0:
        all_df = pd.concat(list_df)
    else:
        return
    if (out_file):
        all_df.to_csv(out_file,sep = "\t", index=False)
    return all_df

"""
Mode 3 - area in the genome. Expect: chr, start, end
only for debugging
"""
#
def main_make_pileup_area(bam, genome, nt, chr, start, end,
                           min_coverage=1, min_mutation=0, min_nonRef=0):
    dict_res = {}
    all_df = pileup_bam(bam, genome, nt, chr, start, end,
               min_coverage, min_mutation, min_nonRef,dict_res)
    return(all_df)



def parse_user_data():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', help='path to  bam file', required=True)
    parser.add_argument('--genome', help='fasta file of the genome', required=True)
    parser.add_argument('--sites_file', help='List of sites to check. If not provided, the script will run on all genome',
                        required=False, default=False)
    parser.add_argument('--out_file', help='name of output file', required=True)
    parser.add_argument('--nt', help='nucleotides to include in pileup. '
                                     'Default is all nucleotides ("aAgGcCtT")', required=False, default="aAgGcCtT")
    parser.add_argument('--min_coverage', help='min coverage.'
                                               'Default is 1.', required=False, default=1)
    parser.add_argument('--min_mutation', help='min mutation rate.'
                                               'Default is 0.', required=False, default=0)
    parser.add_argument('--min_nonRef', help='min number of mutated reads.'
                                             'Default is 0.', required=False, default=0)
    parser.add_argument('--num_threads', help='number of threads to use', required=False, default=1)
    args = parser.parse_args()

    general.check_file(args.bam)
    general.check_file(args.genome)


    if args.sites_file:
        general.check_file(args.sites_file)
        main_make_pileup_sites(args.bam, args.genome, args.sites_file, args.nt, args.out_file, args.num_threads,
                               int(args.min_coverage), float(args.min_mutation), int(args.min_nonRef),
                               writeToFile=True)
    else:
        main_make_pileup_genome(args.bam, args.genome, args.out_file, args.nt, args.num_threads,
                                int(args.min_coverage), float(args.min_mutation), int(args.min_nonRef),
                                chr_list=None,
                                writeToFile=True)

    logger.info("Finished successfully. Results are at: {}".format(args.out_file))


if __name__ == '__main__':
    parse_user_data()

