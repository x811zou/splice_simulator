#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett
# =========================================================================

from concurrent.futures import process
from curses.ascii import FF
from distutils.debug import DEBUG
import string
import time

# from doctest import
from xmlrpc.client import boolean
from xxlimited import Xxo
import argparse
from helper import (
    posTlen_to_fragLen,
    sam_data_processor,
    tabix_regions,
    variant_processor,
    simRead_patmat,
    run2bitBatch,
    constructTwoBitInput,
    annotate_transcripts,
    chunk_iter,
    printRead,  #
)
import helper
import importlib

# spliced_reads = importlib.import_module("splice_reads.sim-spliced-rna")
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
from datetime import datetime, timedelta
from collections import defaultdict
from misc_tools.Translation import Translation

rex = Rex()
# GLOBAL VARIABLES:
matches = 0  # number of sites containing alt or ref allele in transcript
mismatches = 0  # number of sites having neither ref nor alt allele in transcript


def makeAltTranscript(gene, haplotype, variants, if_print=False):
    #########
    # this function is used to create a set of REF copy of a gene transcrips, or ALT copy
    # REF copy (haplotype == 0): replace the REF allele for variants in ref gencode filtered transcript to actual REF allele in VCF
    # ALT copy (haplotype == 1): replace the REF allele for variants in ref gencode filtered transcript to actual ALT allele in VCF
    # allele_in_vcf : allele in VCF
    # allele_in_ref : allele in reference transcript (gencode GTF hg19.v19)
    # match: allele in gencode ref transcript (reversed if in reverse strand) MATCH the allele in VCF
    # mismatch: ... NOT MATCH ...
    #########
    if haplotype == 0:
        global matches
        global mismatches
    altGene = copy.deepcopy(gene)
    transcriptIdToBiSNPcount = {}
    transcriptIdToBiSNPpos = defaultdict(set)
    transcript_num = 0
    # loop through each ref gencode transcript, create a REF/ALT copy based on VCF
    for transcript in altGene.transcripts:
        # print(f"transcript {transcript_num}len {len(transcript.sequence)}")
        array = list(transcript.sequence)
        transcript_num += 1
        num = 0
        # loop through each bi-allelic SNP from VCF
        for variant in variants:
            trans_pos = transcript.mapToTranscript(variant.genomicPos)
            if trans_pos < 0:
                continue
            if variant.genotype[0] != variant.genotype[1]:
                transcriptIdToBiSNPpos[transcript.getID()].add(trans_pos)
                num += 1
            # sanity check
            if haplotype == 0:
                ref_in_ref = array[trans_pos]
                allele_in_vcf = variant.ref
                if (
                    gene.getStrand() == "-"
                ):  # reverse the allele if it is in the reverse strand
                    allele_in_vcf = Translation.reverseComplement(allele_in_vcf)
                allele_check = allele_in_vcf
                # print(
                #     f">>> {gene.getStrand()} -- VCF: {variant.genotype[0]} | {variant.genotype[1]} - {variant.ref}|{variant.alt}; vcf {allele_check} vs reference: {ref_in_ref}"
                # )
                if allele_check == ref_in_ref:
                    # print("match")
                    matches += 1
                else:
                    # print("mismatch")
                    mismatches += 1

            # use REF/ALT allele in VCF to modify the reference transcript
            if variant.genotype[haplotype] > 0:
                allele_in_vcf = variant.alt
                if (
                    gene.getStrand() == "-"
                ):  # reverse the allele if it is in the reverse strand
                    allele_in_vcf = Translation.reverseComplement(allele_in_vcf)
                array[trans_pos] = allele_in_vcf

            ##########################
            # if if_print is not False:
            #     print(
            #         "        %dth transcript, %dth bi-allelic SNP"
            #         % (transcript_num, num)
            #     )
            #     print(
            #         "            > if reverse strand: %s, variant trans pos: %s, haplotype: %s, in ref genome: %s, allele check %s, ref: %s, alt: %s, write_in_sequence: %s"
            #         % (
            #             str(if_rev),
            #             trans_pos,
            #             variant.genotype,
            #             allele_in_ref,
            #             allele_check,
            #             variant.ref,
            #             variant.alt,
            #             array[trans_pos],
            #         )
            #     )
        transcript.sequence = "".join(array)
        transcriptIdToBiSNPcount[transcript.getID()] = num
        # print(transcriptIdToBiSNPpos)
        # if if_print:
        #     print(
        #         "        >>>>>> %dth transcript, in total %d bi-allelic SNP"
        #         % (transcript_num, num)
        #     )
        #     print(" ")
    return altGene, transcriptIdToBiSNPpos, transcript_num


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
parser.add_argument("twobit", help="full path to two bit")
parser.add_argument("genome", help="full path to hg19.2bit")
parser.add_argument("gff", help="full path to gencode gtf file")
parser.add_argument("samgz", help="full path to sam.gz")
parser.add_argument("vcf", help="full path to VCF file with chr")
parser.add_argument("readLen", help="original read length from reads")
parser.add_argument("out_path", help="output path")
parser.add_argument("read_depth", help="per-base-read-depth", type=int)
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
parser.add_argument("--gene", help="specific gene to simulate")
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
twoBitDir = args.twobit
genome2bit = args.genome
gffFile = args.gff
samFile = args.samgz
vcfFile = args.vcf
readLen = int(args.readLen)
out_path = args.out_path
chromosome = args.chr
target_gene = args.gene
DEPTH = args.read_depth
outFile1 = args.out1
outFile2 = args.out2
if_random = args.random
if_print = args.verbose
# configFile = args.config_file

# Load config file
# configFile = ConfigFile(configFile)
# twoBitDir = configFile.lookupOrDie("util-dir")
# genome2bit = configFile.lookupOrDie("genome")
# vcfFile = configFile.lookupOrDie("vcf")
# samFile = configFile.lookupOrDie("aligned-rna")
# gffFile = configFile.lookupOrDie("gff")
# readLen = int(configFile.lookupOrDie("original-read-len"))
# out_path = configFile.lookupOrDie("out_path")

if_debug = False
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
if target_gene is not None:
    print(f"looking at specific gene {target_gene}")
    genes = list(filter(lambda x: x.getId() == target_gene, genes))

if len(genes) == 0:
    print(f"No genes to process.  Exiting.")
    sys.exit()

#    out_path = out_path + "/" + chromosome
# else:
#    out_path = out_path + "/allchr"
out_path_folder = out_path + "/simulated_fastq"

# Path(out_path).mkdir(parents=True, exist_ok=True)
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

suffix = "".join([random.choice(string.ascii_letters) for i in range(6)])
twoBitInputFile = f"/tmp/twobitinput_{suffix}"
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

####### extract sam data
# dict3: gene_to_qualityStr {gene region}:{quality string score from sam.gz}
#######
region_str_to_sam_data = tabix_regions(
    regions, sam_data_processor, samFile, comment_char="@"
)

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
start_time_ns = time.perf_counter_ns()
# list_fragLen = []
for gene in genes:
    # list_start1 = []
    # list_start2 = []
    # list_end1 = []
    # list_end2 = []
    if gene.getSubstrate() == chromosome:
        transcript = gene.longestTranscript()
        transcript.exons = transcript.getRawExons()
        transcript.recomputeBoundaries()
        region_str = f"{gene.getSubstrate()}:{gene.getBegin()}-{gene.getEnd()}"
        chrN = gene.getSubstrate()
        geneid = gene.getId()
        length = gene.longestTranscript().getLength()

        numReads = int(float(DEPTH / readLen) * length)
        

        if processed_genes > 0 and processed_genes % 100 == 0:
            sec_per_gene = (time.perf_counter_ns() - start_time_ns) / processed_genes / 1e9
            estimated_seconds_remaining = round((len(genes) - processed_genes) * sec_per_gene)
            print(
                f"{datetime.now()} ...  processed {processed_genes} / {len(genes)} genes ({round(100*processed_genes/len(genes),2)}%) genes_per_sec: {round(1/sec_per_gene, 2)} ETA: {timedelta(seconds=estimated_seconds_remaining)}",
                file=sys.stderr,
                flush=True,
            )
        processed_genes += 1

        ###### filtering 1
        if not region_str in region_str_to_sam_data:
            # if if_print:
            #     print(f"{chrN},gene: {geneid}, no mapped reads in SAM, skip")
            not_in_sam += 1
            continue
        if not region_str in region_str_to_variants:
            # if if_print:
            #     print(f"{chrN},gene: {geneid}, no variants/records in VCF, skip")
            in_sam_not_in_vcf += 1
            continue

        variants = region_str_to_variants[region_str]
        sam_data = region_str_to_sam_data[region_str]
        pos1_tlen = [x[0] for x in sam_data]
        qual_strs = [x[1] for x in sam_data]
        if len(qual_strs) == 0 or len(pos1_tlen) == 0:
            continue

        ###### filtering 2
        # print(">>>>>>>>>>>>>>>> input pos1_tlen list")
        # print(pos1_tlen)
        pos1_tlen_to_count = {}
        for x in pos1_tlen:
            pos1_tlen_to_count[x] = pos1_tlen_to_count[x] + 1 if x in pos1_tlen_to_count else 1

        # if len(pos1_tlen) > 999:
        # print(f"{datetime.now()}\t{geneid}\ttranscripts: {gene.getNumTranscripts()}\treads: {len(pos1_tlen)}\tdeduped reads: {len(pos1_tlen_to_count)}\tunique_pos1: {len(set([x[0] for x in pos1_tlen]))}")
        transcript_to_fragLen = posTlen_to_fragLen(gene, pos1_tlen_to_count, readLen)
        # print(f"{' '.join([str(len(l)) for l in transcript_to_fragLen.values()])}")

        if len(transcript_to_fragLen) == 0:
            # print("empty fragLen list!")
            continue
        # print(">>>>>>>>>>>>>>>> output fragLen list")
        # print(fragLen_list)
        recorded_genes += 1
        # if if_print:
        #     print(
        #         "%s,%s,#reads: %s,total #variants in VCF: %d"
        #         % (region_str, geneid, numReads, len(variants))
        #     )

        maternal, transcriptIdToBiSNPpos, transcript_num = makeAltTranscript(
            gene, 1, variants
        )
        paternal, transcriptIdToBiSNPpos, _ = makeAltTranscript(gene, 0, variants)
        qual_idx = 0
        # pos_idx = 0
        # counter = 0
        # list_fragLen = []

        candidate_transcripts = list(transcript_to_fragLen.keys())
        candidate_transcript_pairs = [
            (x, next(filter(lambda y: y.getID() == x.getID(), maternal.transcripts))) for x in candidate_transcripts
        ]

        for i in range(numReads):
            patTranscript, matTranscript = random.choice(candidate_transcript_pairs)
            transcript_length = matTranscript.getLength()
            # extract quality string from qual_strs in order
            ##########
            fwd_qual = qual_strs[qual_idx]
            qual_idx = (qual_idx + 1) % len(qual_strs)
            rev_qual = qual_strs[qual_idx]
            qual_idx = (qual_idx + 1) % len(qual_strs)
            max_qual_len = max(len(fwd_qual), len(rev_qual))
            ##########
            # print(
            #    f">>>>>>>>>>>>>>>> {i}th reads - transcript length {transcript_length}"
            # )

            fragLen = random.choice(transcript_to_fragLen[patTranscript])

            # fragLen = pick_fragLen(fragLen_list, max_qual_len, transcript_length)
            # pos_idx = (pos_idx + 1) % len(pos1_tlen)

            if fragLen is None:
                if if_debug:
                    print(
                        f">> {geneid} {i} th read skipped! coudl not find appropriate fraglen"
                    )
                n_break += 1
                continue

            # list_fragLen.append(fragLen)
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
            # list_fragLen.append(fragLen)
            # list_start1.append(start1)
            # list_end1.append(end1)
            # list_start2.append(start2)
            # list_end2.append(end2)
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

# Path(out_path + "/fragLen").mkdir(parents=True, exist_ok=True)
# with open(out_path + "/fragLen.txt", "wb") as fp:
#     pickle.dump(list_fragLen, fp)

print(
    f"{datetime.now()} DONE",
    file=sys.stderr,
    flush=True,
)
print("")
print(">> total geneID in gtf: %d" % (num_gene_gtf))
print(">> total ReadID processed: %d" % (nextReadID))
print(
    f">> # matches: {matches}, # mismatches: {mismatches}, percentage of matches: {round(matches/(mismatches+matches),2)*100}%"
)
print(">> total num breaks /gene is shorter than fragment length : %d" % (n_break))
print(
    f">> # recorded genes : {recorded_genes}, # processed_genes: {processed_genes}, percentage: {round(recorded_genes/processed_genes,2)*100}%"
)
print(
    f">> # gene not in SAM: {not_in_sam}, # genes in SAM but not in VCF: {in_sam_not_in_vcf}, total filered out genes: {not_in_sam+in_sam_not_in_vcf} ({round((not_in_sam+in_sam_not_in_vcf)/processed_genes,2)*100})%"
)
