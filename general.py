import sys  
import os
import logging
import pysam


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHeandeler = logging.StreamHandler()
logger.addHandler(streamHeandeler )
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHeandeler.setFormatter(formatter)





def check_file(file_name):
    if os.path.isfile(file_name):
        return
    else:
        logger.critical("file {} not found".format(file_name))
        sys.exit(1)

def remove_file(file_name):
    if os.path.isfile(file_name):
        logger.info("removing {}".format(file_name))
        os.remove(file_name)

def split_bam(bam, out_pos_file, out_neg_file, pos_strand = "R2"):
    check_file(bam)

    samfile = pysam.AlignmentFile(bam, "rb")
    out_pos = pysam.AlignmentFile(out_pos_file, "wb", template=samfile)
    out_neg = pysam.AlignmentFile(out_neg_file, "wb", template=samfile)

    if  (pos_strand  == "R2"):

        logger.info("Splitting bam. Positive stand in R2")

        for read in samfile:
            pos_read1 = read.is_read1 and read.is_reverse
            pos_read2 = read.is_read2 and (not read.is_reverse)
            neg_read1 = read.is_read1 and (not read.is_reverse)
            neg_read2 = read.is_read2 and read.is_reverse
            if pos_read1 or pos_read2:
                out_pos.write(read)
            elif neg_read1 or neg_read2 :
                out_neg.write(read)

    elif (pos_strand  == "R1"):

        logger.info("Splitting bam. Positive stand in R1")

        for read in samfile:
            pos_read1 = read.is_read1 and (not read.is_reverse)
            pos_read2 = read.is_read2 and read.is_reverse
            neg_read1 = read.is_read1 and read.is_reverse
            neg_read2 = read.is_read2 and (not read.is_reverse)
            if pos_read1 or pos_read2:
                out_pos.write(read)
            elif neg_read1 or neg_read2:
                out_neg.write(read)

    else:
        logger.critical("pos_strand must be R1 or R2")
        sys.exit(1)

    logger.info("Bam splitted by strand")

def index_file(bam):
    logger.info("indexing file {}".format(bam))
    os.system("samtools index {}".format(bam))