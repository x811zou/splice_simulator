#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett
# =========================================================================

from curses.ascii import FF
from distutils.debug import DEBUG

# from doctest import
from xmlrpc.client import boolean
from xxlimited import Xxo
import argparse
from helper import (
    loop_transcript_for_fragLen,
    tabix_regions,
    variant_processor,
    quality_string_processor,
    read_length_processor,
    simRead_patmat,
    makeAltTranscript,
    run2bitBatch,
    constructTwoBitInput,
    annotate_transcripts,
    chunk_iter,
    printRead,  #
    pick_random_Transcript,
)
import helper
import importlib

#spliced_reads = importlib.import_module("splice_reads.sim-spliced-rna")
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import random
import gzip
import os
import copy
import pickle
from misc_tools.GffTranscriptReader import GffTranscriptReader
from misc_tools.Pipe import Pipe
from misc_tools.ConfigFile import ConfigFile
from misc_tools.Rex import Rex
from Bio.Seq import Seq
from pathlib import Path
from datetime import datetime

# ============================================
# DEPENDENCIES:
#   module load htslib
#   python version 3.x
#   VCF file must be bzipped and indexed with tabix.
#   VCF must contain only one sample (individual).
#   VCF must include all sites, including homozygous and heterozygous sites.
#
# EXAMPLE CONFIG FILE:
#   util-dir = /Users/scarlett/BEASTIE/BEASTIE/reference/two_bit
#   genome = /Users/scarlett/BEASTIE/BEASTIE/reference/two_bit/hg19.2bit
#   aligned-rna = /Users/scarlett/allenlab/BEASTIE_other_example/HG00108_50M/HG00108.sam.gz
#   fragment-lengths = /Users/scarlett/allenlab/BEASTIE_other_example/HG00108_50M/random_splice_simulator/fragment-lengths.txt
#   vcf = /Users/scarlett/allenlab/BEASTIE_other_example/HG00108_50M/HG00108.with_chr.content.SNPs.hets.header.vcf.gz
#   gff = /Users/scarlett/BEASTIE/BEASTIE/reference/two_bit/gencode.v19.annotation.filtered.chr5.gtf
#   original-read-len = 75
#   out_path = /Users/scarlett/allenlab/BEASTIE_other_example/HG00108_50M/random_splice_simulator
#
# EXAMPLE running command:
# python sim-spliced-rna.py /Users/scarlett/allenlab/BEASTIE_other_example/HG00108_50M/random_splice_simulator/allchr.config 100 --chr chr21 -r -v
# ============================================

rex = Rex()
# =========================================================================
# main()
# =========================================================================
# if len(sys.argv) < 8:
#     exit(
#         os.path.basename(sys.argv[0])
#         + " <config-file> <per-base-read-depth> <if_random> <if_print> <prefix> <out-read1.gz> <out-read2.gz>\n"
#     )
# (configFile, DEPTH, if_random, if_print, prefix, outFile1, outFile2) = sys.argv[1:]
parser = argparse.ArgumentParser()
parser.add_argument("config_file", help="path to config file")
parser.add_argument("read_depth", help="per-base-read-depth")
parser.add_argument(
    "--out1",
    help="output name for forward strand fastq reads, default: read1.fastq.gz",
    default="read1.fastq.gz",
)
parser.add_argument(
    "--out2",
    help="output name for reverse strand fastq reads, default: read2.fastq.gz",
    default="read2.fastq.gz",
)
parser.add_argument("--chr", help="specific chromosome to simulate")
parser.add_argument(
    "-r",
    "--random",
    action="store_true",
    help="generate reads randomly from either haplotype",
)
parser.add_argument("-v", "--verbose", action="store_true", help="print lots of info")
args = parser.parse_args()
if args.chr:
    print(f"single chromosome mode turned on : {args.chr}")
if args.verbose:
    print("printing mode turned on")

print(args)
chromosome = args.chr
DEPTH = int(args.read_depth)
outFile1 = args.out1
outFile2 = args.out2
if_random = args.random
if_print = args.verbose
configFile = args.config_file

# Load config file
configFile = ConfigFile(configFile)
twoBitDir = configFile.lookupOrDie("util-dir")
genome2bit = configFile.lookupOrDie("genome")
vcfFile = configFile.lookupOrDie("vcf")
samFile = configFile.lookupOrDie("aligned-rna")
gffFile = configFile.lookupOrDie("gff")
readLen = int(configFile.lookupOrDie("original-read-len"))
out_path = configFile.lookupOrDie("out_path")

if_debug = False
if_print = str(if_print) == "True"
if if_print:
    print(
        f"{datetime.now()} Yes, print out simulation logs...",
        file=sys.stderr,
        flush=True,
    )
else:
    print(
        f"{datetime.now()} Not print out simulation logs...",
        file=sys.stderr,
        flush=True,
    )

if_random = str(if_random) == "True"
if if_random:
    print(
        f"{datetime.now()} generate RANDOM pat/mat reads...",
        file=sys.stderr,
        flush=True,
    )
else:
    print(
        f"{datetime.now()} generate EQUAL pat/mat reads...", file=sys.stderr, flush=True
    )


# Load GFF and fragment lengths
gffReader = GffTranscriptReader()
print(f"{datetime.now()} reading GFF...", file=sys.stderr, flush=True)
genes = gffReader.loadGenes(gffFile)
genes.sort(key=lambda gene: (gene.getSubstrate(), gene.getBegin()))
print(f"{datetime.now()} done reading GFF...", file=sys.stderr, flush=True)


if chromosome is not None:
    print(f"Looking at specific chromosome: {chromosome}")
    genes = list(filter(lambda x: x.getSubstrate() == chromosome, genes))

#    out_path = out_path + "/" + chromosome
#else:
#    out_path = out_path + "/allchr"
out_path_folder = out_path

#Path(out_path).mkdir(parents=True, exist_ok=True)
Path(out_path_folder).mkdir(parents=True, exist_ok=True)

## DEBUG
# DEBUG_GENES = ["ENSG00000180530.5"]
# if DEBUG_GENES is not None:
#     print(f"Debugging specific genes: {DEBUG_GENES}")
#     genes = list(filter(lambda x: x.getId() in DEBUG_GENES, genes))

# Create output files
OUT1 = gzip.open(out_path_folder + "/" + outFile1, "wt")
OUT2 = gzip.open(out_path_folder + "/" + outFile2, "wt")

# Simulate
print("simulating...", file=sys.stderr, flush=True)
nextReadID = 0
IN = open(samFile, "rt")
num_gene_gtf = 0
n_break = 0
counter = 0

print(
    f"{datetime.now()} processing {len(genes)} genes for transcripts",
    file=sys.stderr,
    flush=True,
)

twoBitInputFile = "/tmp/twobitinput"
constructTwoBitInput(genes, twoBitInputFile)

print(
    f"{datetime.now()} running 2bit",
    file=sys.stderr,
    flush=True,
)
twobitId_to_seq = run2bitBatch(twoBitDir, twoBitInputFile, genome2bit)

print(
    f"{datetime.now()} done running 2bit {len(twobitId_to_seq)} ids",
    file=sys.stderr,
    flush=True,
)

annotate_transcripts(genes, twobitId_to_seq)

print(
    f"{datetime.now()} done annotating transcripts",
    file=sys.stderr,
    flush=True,
)
####### extract regions for genes from reference gencode GTF
# Build list of gene regions to extract
#######

regions = set()
for idx, gene in enumerate(genes):
    counter += 1
    num_gene_gtf += 1
    region_str = f"{gene.getSubstrate()}:{gene.getBegin()}-{gene.getEnd()}"
    regions.add(region_str)

print(
    f"{datetime.now()} Got {num_gene_gtf} genes from gencode gtf",
    file=sys.stderr,
    flush=True,
)
####### extract regions for genes from VCF file
# dict2: gene_to_variants {gene region}:{records from VCF}
#######
region_str_to_variants = tabix_regions(
    regions, variant_processor, vcfFile, comment_char="#"
)

####### extract quality string score
# dict3: gene_to_qualityStr {gene region}:{quality string score from sam.gz}
#######
region_str_to_quality_strings = tabix_regions(
    regions, quality_string_processor, samFile, comment_char="@"
)
# print(region_str_to_quality_strings)
####### extract fragment length
# dict4: gene_to_fragLen {gene region}:{quality string score from sam.gz}
#######
region_str_to_fragLen = tabix_regions(
    regions, read_length_processor, samFile, comment_char="@"
)
# print(region_str_to_fragLen)

#######
# for each gene, generate reads using quality string from matching genes in SAM.GZ
#######
print(
    f"{datetime.now()} Start simulation",
    file=sys.stderr,
    flush=True,
)
processed_genes = 0
recorded_genes = 0
not_in_sam = 0
in_sam_not_in_vcf = 0
list_fragLen = []
for gene in genes:
    list_start1 = []
    list_start2 = []
    list_end1 = []
    list_end2 = []
    if gene.getSubstrate() == chromosome:
        transcript = gene.longestTranscript()
        transcript.exons = transcript.getRawExons()
        transcript.recomputeBoundaries()
        region_str = f"{gene.getSubstrate()}:{gene.getBegin()}-{gene.getEnd()}"
        chrN = gene.getSubstrate()
        geneid = gene.getId()
        length = gene.longestTranscript().getLength()

        numReads = int(float(DEPTH / readLen) * length)
        processed_genes += 1

        if not region_str in region_str_to_quality_strings:
            # if if_print:
            #     print(f"{chrN},gene: {geneid}, no mapped reads in SAM, skip")
            not_in_sam += 1
            continue
        if not region_str in region_str_to_variants:
            # if if_print:
            #     print(f"{chrN},gene: {geneid}, no variants/records in VCF, skip")
            in_sam_not_in_vcf += 1
            continue
        if not region_str in region_str_to_fragLen:
            # if if_print:
            #     print(f"{chrN},gene: {geneid}, no fragment length in SAM, skip")
            continue

        variants = region_str_to_variants[region_str]
        qual_strs = region_str_to_quality_strings[region_str]
        pos1_tlen = region_str_to_fragLen[region_str]
        # print(pos1_tlen)

        if len(qual_strs) == 0 or len(pos1_tlen) == 0:
            continue

        recorded_genes += 1
        # if if_print:
        #     print(
        #         "%s,%s,#reads: %s,total #variants in VCF: %d"
        #         % (region_str, geneid, numReads, len(variants))
        #     )

        maternal, transcriptIdToBiSNPcount, transcriptIdToBiSNPpos = makeAltTranscript(
            gene, 1, variants
        )
        paternal, _, _ = makeAltTranscript(gene, 0, variants)
        qual_idx = 0
        pos_idx = 0
        # counter = 0
        list_fragLen = []

        for i in range(numReads):
            (
                matTranscript,
                patTranscript,
                th_transcript,
                transcriptID,
                bivariant,
            ) = pick_random_Transcript(maternal, paternal, transcriptIdToBiSNPpos)
            transcript_length = matTranscript.getLength()
            # extract quality string from qual_strs in order
            ##########
            fwd_qual = qual_strs[qual_idx]
            qual_idx = (qual_idx + 1) % len(qual_strs)
            rev_qual = qual_strs[qual_idx]
            qual_idx = (qual_idx + 1) % len(qual_strs)
            ##########
            #print(
            #    f">>>>>>>>>>>>>>>> {i}th reads - transcript length {transcript_length}"
            #)
            fragLen, pos_idx = loop_transcript_for_fragLen(
                i,
                gene,
                pos_idx,
                pos1_tlen,
                readLen,
                transcript_length,
                if_debug=if_debug,
            )
            pos_idx = (pos_idx + 1) % len(pos1_tlen)

            if fragLen is None:
                if if_debug:
                    print(f">> {geneid} {i} th read skipped!")
                n_break += 1
                break
            list_fragLen.append(fragLen)
            # simulate reads for paternal and maternal copy
            (
                patSeq,
                patSeq_rev,
                matSeq,
                matSeq_rev,
                fwd_LEN,
                rev_LEN,
                start1,
                end1,
                start2,
                end2,
            ) = simRead_patmat(
                patTranscript, matTranscript, fwd_qual, rev_qual, fragLen, readLen
            )
            list_fragLen.append(fragLen)
            list_start1.append(start1)
            list_end1.append(end1)
            list_start2.append(start2)
            list_end2.append(end2)
            if patSeq is None or matSeq is None:
                n_break += 1
                continue  # gene is shorter than fragment length
            ################## random haplotype simulator: randomly generate a mat or pat copy
            if if_random:
                random_prob = random.random()
                if_mat = False
                if random_prob >= 0.5:
                    if_mat = True
                identifier_random = "@SIM-" + str(nextReadID) + "-" + str(geneid)
                # if if_print:
                #     print(
                #         "%s,if maternal: %s,rec1 quality string length %s, forward strand length %s, rec2 quality strand length %s, reverse strand length %s!! "
                #         % (
                #             identifier_random,
                #             str(if_mat),
                #             len(fwd_qual),
                #             fwd_LEN,
                #             len(rev_qual),
                #             rev_LEN,
                #         )
                #     )
                # random decision on maternal / paternal
                if if_mat is True:
                    randomSeq = matSeq
                    randomSeq_rev = matSeq_rev
                    # if if_print:
                    #     print(
                    #         "%s,if mat:%s,MAT FWD: %s"
                    #         % (identifier_random, str(if_mat), matSeq)
                    #     )
                    #     print(
                    #         "%s,if mat:%s,MAT REV: %s"
                    #         % (identifier_random, str(if_mat), matSeq_rev)
                    #     )
                else:
                    randomSeq = patSeq
                    randomSeq_rev = patSeq_rev
                    # if if_print:
                    #     print(
                    #         "%s,if mat:%s,PAT FWD: %s"
                    #         % (identifier_random, str(if_mat), patSeq)
                    #     )
                    #     print(
                    #         "%s,if mat:%s,PAT REV: %s"
                    #         % (identifier_random, str(if_mat), patSeq_rev)
                    #     )
                # if if_print:
                #     print("%s,qual FWD: %s" % (identifier_random, fwd_qual))
                #     print("%s,qual REV: %s" % (identifier_random, rev_qual))
                #     print("")
                # write to file
                printRead(
                    str(identifier_random)
                    + ":if_maternal"
                    + str(if_mat)
                    + ":FWD"
                    + " /1",
                    randomSeq,
                    fwd_qual,
                    OUT1,
                )
                nextReadID += 1
                printRead(
                    str(identifier_random)
                    + ":if_maternal"
                    + str(if_mat)
                    + ":REV"
                    + " /1",
                    randomSeq_rev,
                    rev_qual,
                    OUT2,
                )
            ################## equal haplotype copy simulator: both mat and pat has same number of reads
            else:
                identifier_PAT = "@SIM-" + str(nextReadID) + "-" + str(geneid)
                identifier_MAT = "@SIM-" + str(nextReadID + 1) + "-" + str(geneid)
                # if if_print:
                #     print(
                #         "REF: %s,ALT: %s,rec1 quality string length %s, forward strand length %s, rec2 quality strand length %s, reverse strand length %s!! "
                #         % (
                #             identifier_PAT,
                #             identifier_MAT,
                #             len(fwd_qual),
                #             fwd_LEN,
                #             len(rev_qual),
                #             rev_LEN,
                #         )
                #     )
                #     print("%s,PAT FWD: %s" % (identifier_PAT, patSeq))
                #     print("%s,MAT FWD: %s" % (identifier_MAT, matSeq))
                #     print("%s,qual FWD: %s" % (identifier_PAT, fwd_qual))
                #     print("%s,PAT REV: %s" % (identifier_PAT, patSeq_rev))
                #     print("%s,MAT REV: %s" % (identifier_MAT, matSeq_rev))
                #     print("%s,qual REV: %s" % (identifier_MAT, rev_qual))
                #     print("")
                # write to file
                printRead(
                    identifier_PAT + ":PAT:FWD" + " /1",
                    patSeq,
                    fwd_qual,
                    OUT1,
                )
                printRead(
                    identifier_PAT + ":PAT:REV" + " /2",
                    patSeq_rev,
                    rev_qual,
                    OUT2,
                )
                nextReadID += 1
                printRead(
                    identifier_MAT + ":MAT:FWD" + " /1",
                    matSeq,
                    fwd_qual,
                    OUT1,
                )
                printRead(
                    identifier_MAT + ":MAT:REV" + " /2",
                    matSeq_rev,
                    rev_qual,
                    OUT2,
                )
            nextReadID += 1
        # if processed_genes <= 10:
        #     out = out_path + "/start_end_pos" + "/" + str(geneid)
        #     Path(out).mkdir(parents=True, exist_ok=True)
        #     with open(out + "/start1", "wb") as fp:
        #         pickle.dump(list_start1, fp)
        #     with open(out + "/end1", "wb") as fp:
        #         pickle.dump(list_end1, fp)
        #     with open(out + "/start2", "wb") as fp:
        #         pickle.dump(list_start2, fp)
        #     with open(ovi testut + "/end2", "wb") as fp:
        #         pickle.dump(list_end2, fp)

        if processed_genes % 500 == 0:
            print(
                f"{datetime.now()} ... processed {processed_genes} / {len(genes)} genes : {round(100*processed_genes/len(genes),2)} %",
                file=sys.stderr,
                flush=True,
            )
#Path(out_path + "/fragLen").mkdir(parents=True, exist_ok=True)
#with open(out_path + "/fragLen.txt", "wb") as fp:
#    pickle.dump(list_fragLen, fp)

print(
    f"{datetime.now()} DONE",
    file=sys.stderr,
    flush=True,
)
print("")
print(">> total geneID in gtf: %d" % (num_gene_gtf))
print(">> total ReadID processed: %d" % (nextReadID))
# print(f">> # matches: {matches}, # mismatches: {mismatches}")
print(">> total num breaks /gene is shorter than fragment length : %d" % (n_break))
print(
    f">> # recorded genes : {recorded_genes}, # processed_genes: {processed_genes}, percentage: {round(recorded_genes/processed_genes,2)*100}%"
)
print(
    f">> # gene not in SAM: {not_in_sam}, # genes in SAM but not in VCF: {in_sam_not_in_vcf}, total filered out genes: {not_in_sam+in_sam_not_in_vcf} ({round((not_in_sam+in_sam_not_in_vcf)/processed_genes,2)*100})%"
)
