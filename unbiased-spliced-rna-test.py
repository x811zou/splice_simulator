#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

from concurrent.futures import process
from curses.ascii import FF
from distutils.debug import DEBUG
import re
import string
import tempfile
import time

# from doctest import
from xmlrpc.client import boolean
from xxlimited import Xxo
import argparse
from helper import (
    posTlen_to_fragLen,
    sam_data_processor,
    tabix_regions,
    variant_processor_SNPs,
    variant_processor_hets,
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

def modifyTranscript_test(gene, variants, if_print=False):
    #########
    # this function is used to create a set of paternal/maternal copy of a gene transcrips
    # Paternal haplotype sequence (haplotype == 0): replace the allele for variants in ref gencode filtered transcript to actual ALT allele in VCF if first genotype column [0] is 1
    # Maternal haplotype sequence (haplotype == 1): replace the allele for variants in ref gencode filtered transcript to actual ALT allele in VCF if second genotype column [1] is 1
    # allele_in_vcf : allele in VCF
    # refallele_in_ref : allele in reference transcript (gencode GTF hg19.v19)
    # match: allele in gencode ref transcript (reversed if in reverse strand) MATCH the allele in VCF
    # mismatch: ... NOT MATCH in ref allele of variant...
    #########
    # for pat transcript only
    global matches
    global mismatches
    transcriptIdToBiSNPcount = {}
    transcriptIdToBiSNPpos = defaultdict(set)
    transcript_num = 0
    debug_print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> modifying transcript: ")
    patGene = copy.deepcopy(gene)
    matGene = copy.deepcopy(gene)
    # loop through each ref gencode transcript, create a mat/pat haplotype based on VCF, genotype[0] for pat, genotype[1] for mat
    for patTranscript,matTranscript in zip(patGene.transcripts,matGene.transcripts):
        transcript_num += 1
        num = 0
        debug_print(f"transcript {transcript_num} - {patTranscript.getID()}: len {len(patTranscript.sequence)}")
        patSeq = list(patTranscript.sequence)
        matSeq = list(matTranscript.sequence)
        # loop through each bi-allelic SNP from VCF
        for variant in variants:
            trans_pos = patTranscript.mapToTranscript(variant.genomicPos)
            if trans_pos < 0:
                continue
            debug_print(f"transcript {transcript_num} - {patTranscript.getID()}: genomic pos {variant.genomicPos} - trans pos: {trans_pos} - {variant.genotype[0]} | {variant.genotype[1]} - ref {variant.ref} | alt {variant.alt}")
            if variant.genotype[0] != variant.genotype[1]:
                transcriptIdToBiSNPpos[patTranscript.getID()].add(variant.genomicPos)
                num += 1
                refallele_in_ref = patSeq[trans_pos]
                refallele_in_vcf = variant.ref
                altallele_in_vcf = variant.alt
                if (
                    gene.getStrand() == "-"
                ):  # reverse the allele if it is in the reverse strand
                    if_rev = True
                    refallele_in_vcf = Translation.reverseComplement(refallele_in_vcf)
                    altallele_in_vcf = Translation.reverseComplement(altallele_in_vcf)
                ### debugging start
                #debug_print(f">>> information {gene.getStrand()} -- VCF: {variant.genotype[0]} | {variant.genotype[1]} - ref:{refallele_in_vcf} | alt:{altallele_in_vcf}; vs reference: {refallele_in_ref}, if reverse strand: {if_rev}")
                if refallele_in_vcf == refallele_in_ref:
                    #debug_print(f"match")
                    matches += 1
                else:
                    #debug_print(f"mismatch")
                    mismatches += 1
                ### debugging end
                patSeq[trans_pos] = refallele_in_vcf
                matSeq[trans_pos] = altallele_in_vcf
                debug_print(f">>> simulation {gene.getStrand()} -- VCF: {variant.genotype[0]} | {variant.genotype[1]} - ref:{refallele_in_vcf} | alt:{altallele_in_vcf}; vs reference: {refallele_in_ref} - pat seq :{patSeq[trans_pos]} | mat seq :{matSeq[trans_pos]}")
                ##########################
        patTranscript.sequence = "".join(patSeq)
        matTranscript.sequence = "".join(matSeq)
        transcriptIdToBiSNPcount[patTranscript.getID()] = num
        debug_print(f">>>>>> %dth transcript, in total %d bi-allelic SNP"
                % (transcript_num, num)
            )
        debug_print(f"")
    debug_print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        )
    debug_print(transcriptIdToBiSNPpos)
    return patGene,matGene, transcriptIdToBiSNPpos, transcript_num

def modifyTranscript(gene, variants, if_print=False):
    #########
    # this function is used to create a set of paternal/maternal copy of a gene transcrips
    # Paternal haplotype sequence (haplotype == 0): replace the allele for variants in ref gencode filtered transcript to actual ALT allele in VCF if first genotype column [0] is 1
    # Maternal haplotype sequence (haplotype == 1): replace the allele for variants in ref gencode filtered transcript to actual ALT allele in VCF if second genotype column [1] is 1
    # allele_in_vcf : allele in VCF
    # refallele_in_ref : allele in reference transcript (gencode GTF hg19.v19)
    # match: allele in gencode ref transcript (reversed if in reverse strand) MATCH the allele in VCF
    # mismatch: ... NOT MATCH in ref allele of variant...
    #########
    # for pat transcript only
    global matches
    global mismatches
    transcriptIdToBiSNPcount = {}
    transcriptIdToBiSNPpos = defaultdict(set)
    transcript_num = 0
    debug_print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> modifying transcript: ")
    patGene = copy.deepcopy(gene)
    matGene = copy.deepcopy(gene)
    # loop through each ref gencode transcript, create a mat/pat haplotype based on VCF, genotype[0] for pat, genotype[1] for mat
    for patTranscript,matTranscript in zip(patGene.transcripts,matGene.transcripts):
        transcript_num += 1
        num = 0
        debug_print(f"transcript {transcript_num} - {patTranscript.getID()}: len {len(patTranscript.sequence)}")
        patSeq = list(patTranscript.sequence)
        matSeq = list(matTranscript.sequence)
        # loop through each bi-allelic SNP from VCF
        for variant in variants:
            trans_pos = patTranscript.mapToTranscript(variant.genomicPos)
            if trans_pos < 0:
                continue
            debug_print(f"transcript {transcript_num} - {patTranscript.getID()}: genomic pos {variant.genomicPos} - trans pos: {trans_pos} - {variant.genotype[0]} | {variant.genotype[1]} - ref {variant.ref} | alt {variant.alt}")
            if variant.genotype[0] != variant.genotype[1]:
                transcriptIdToBiSNPpos[patTranscript.getID()].add(variant.genomicPos)
                num += 1
            refallele_in_ref = patSeq[trans_pos]
            refallele_in_vcf = variant.ref
            altallele_in_vcf = variant.alt
            if (
                gene.getStrand() == "-"
            ):  # reverse the allele if it is in the reverse strand
                if_rev = True
                refallele_in_vcf = Translation.reverseComplement(refallele_in_vcf)
                altallele_in_vcf = Translation.reverseComplement(altallele_in_vcf)
            ### debugging start
            #debug_print(f">>> information {gene.getStrand()} -- VCF: {variant.genotype[0]} | {variant.genotype[1]} - ref:{refallele_in_vcf} | alt:{altallele_in_vcf}; vs reference: {refallele_in_ref}, if reverse strand: {if_rev}")
            if refallele_in_vcf == refallele_in_ref:
                #debug_print(f"match")
                matches += 1
            else:
                #debug_print(f"mismatch")
                mismatches += 1
            ### debugging end
            if variant.genotype[0] == 0:
                patSeq[trans_pos] = refallele_in_vcf
            else:
                patSeq[trans_pos] = altallele_in_vcf
            if variant.genotype[1] == 0:
                matSeq[trans_pos] = refallele_in_vcf
            else:
                matSeq[trans_pos] = altallele_in_vcf
            debug_print(f">>> simulation {gene.getStrand()} -- VCF: {variant.genotype[0]} | {variant.genotype[1]} - ref:{refallele_in_vcf} | alt:{altallele_in_vcf}; vs reference: {refallele_in_ref} - pat seq :{patSeq[trans_pos]} | mat seq :{matSeq[trans_pos]}")
            ##########################
        patTranscript.sequence = "".join(patSeq)
        matTranscript.sequence = "".join(matSeq)
        transcriptIdToBiSNPcount[patTranscript.getID()] = num
        debug_print(f">>>>>> %dth transcript, in total %d bi-allelic SNP"
                % (transcript_num, num)
            )
        debug_print(f"")
    debug_print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        )
    debug_print(transcriptIdToBiSNPpos)
    return patGene,matGene, transcriptIdToBiSNPpos, transcript_num

# ============================================
# DEPENDENCIES:
#   module load htslib
#   python version 3.x
#   VCF file must be bzipped and indexed with tabix.
#   VCF must contain only one sample (individual).
#   VCF must include all sites, including homozygous and heterozygous sites.
# module load htslib
# module load samtools/1.11-rhel8
######### 1000 Genome individuals
# N=2
# sample="NA12878"
# util_dir=/hpc/home/bmajoros/twobit
# genome=/datacommons/allenlab/hg19/hg19.2bit
# gff=/datacommons/allenlab/hg19/filter/gencode.v19.annotation.level12.gtf
# sam=/hpc/group/allenlab/scarlett/output/RNAseq/GIAB/$sample/tmp/simulation_SNPs_even_101/${sample}.sam.gz
# vcfgz=/hpc/group/allenlab/scarlett/output/RNAseq/GIAB/$sample/${sample}.no_chr.content.SNPs.filtered.vcf.gz
# out_path=/hpc/group/allenlab/scarlett/output/RNAseq/GIAB/$sample/tmp/simulation_test
# read_depth=100
######### GSD individuals
# sample="123375"
# sam=/hpc/group/allenlab/scarlett/output/RNAseq/GSD/$sample/tmp/simulation_SNPs_even_100/${sample}.sam.gz
# vcfgz=/hpc/group/allenlab/scarlett/output/RNAseq/GSD/$sample/${sample}.no_chr.content.SNPs.filtered.vcf.gz
# out_path=/hpc/group/allenlab/scarlett/output/RNAseq/GSD/$sample/tmp/simulation_test
# read_depth=100
#
# EXAMPLE running command:
# python /hpc/group/allenlab/scarlett/script/spliced_simulator/unbiased-spliced-rna-test.py $util_dir $genome $gff $sam $vcfgz $out_path $read_depth --out-prefix chr${N} --chr chr${N} --allSNPs --gene ENSG00000213626.7 -v
# ============================================

# =========================================================================
# main()
# =========================================================================

parser = argparse.ArgumentParser()
parser.add_argument("twobit", help="full path to two bit")
parser.add_argument("genome", help="full path to hg19.2bit")
parser.add_argument("gff", help="full path to gencode gtf file")
parser.add_argument("samgz", help="full path to sam.gz")
parser.add_argument("vcf", help="full path to VCF file without chr")
parser.add_argument("out_path", help="output path")
parser.add_argument("read_depth", help="per-base-read-depth", type=int)
parser.add_argument(
    "--out1",
    help="output name for forward strand fastq reads, default: read_1.fastq.gz",
    default="_1.fastq.gz",
)
parser.add_argument(
    "--out2",
    help="output name for reverse strand fastq reads, default: read_2.fastq.gz",
    default="_2.fastq.gz",
)
parser.add_argument(
    "--out-prefix", help="prefix applied to output file names", default=""
)
parser.add_argument("--chr", help="specific chromosome to simulate")
parser.add_argument("--gene", help="specific gene to simulate",default=None)
parser.add_argument(
    "-r",
    "--random",
    action="store_true",
    help="generate reads randomly from either haplotype",
)
parser.add_argument(
    "-t",
    "--test",
    action="store_true",
    help="generate reads into PAT - ref /MAT - alt haplotype",
)
parser.add_argument(
    "--seed",
    help="random seed for simulation",
    default=random.randrange(sys.maxsize),
    type=int,
)
parser.add_argument("-v", "--verbose", action="store_true", help="print lots of info")
parser.add_argument("-snp", "--allSNPs", action="store_true", help="include all SNPs")

args = parser.parse_args()
if_random = args.random
if if_random:
    print(
        f"{datetime.now()} RANDOMLY generate pat/mat reads for each read ...",
        file=sys.stderr,
        flush=True,
    )
else:
    print(
        f"{datetime.now()} EQUALLY generate pat/mat reads for each read...", file=sys.stderr, flush=True
    )
#
if args.chr:
    print(f"{datetime.now()} CHRomosome mode turned on for {args.chr}")
else:
    print(f"{datetime.now()} default: looking at all genes")
#
if_print = args.verbose
if if_print:
    print(
        f"{datetime.now()} Yes, print out simulation logs for debugging purpose...",
        file=sys.stderr,
        flush=True,
    )
else:
    print(
        f"{datetime.now()} Not printing out simulation logs...",
        file=sys.stderr,
        flush=True,
    )
if if_print:
	debug_print = print
else:
	debug_print = lambda *args: None
#
SNPs = args.allSNPs
if SNPs:
    variant_processor_chosen=variant_processor_SNPs
    print(
        f"{datetime.now()} simulation sites for all SNPs...",
        file=sys.stderr,
        flush=True,
    )
else:
    variant_processor_chosen=variant_processor_hets
    print(
        f"{datetime.now()} simulation sites for all hets...",
        file=sys.stderr,
        flush=True,
    )
#
if args.test:
    print(f"{datetime.now()} testing mode turned on, pat haplotype only has ref, mat haplotype only has alt")
else:
    print(f"{datetime.now()} default setting: each haplpotype has both ref/alt allele")
if_test = args.test
if if_test:
    makingTranscript = modifyTranscript_test
else:
    makingTranscript = modifyTranscript
#
print(args)
#
twoBitDir = args.twobit
genome2bit = args.genome
gffFile = args.gff
samFile = args.samgz
vcfFile = args.vcf
out_path = args.out_path
chromosome = args.chr
target_gene = args.gene
DEPTH = args.read_depth
outFile1 = args.out1
outFile2 = args.out2
outPrefix = args.out_prefix
random.seed(args.seed)

#
print(
    f"{datetime.now()} simulation seed : {args.seed}",
    file=sys.stderr,
    flush=True,
)

# Load GFF and fragment lengths
gffReader = GffTranscriptReader()
print(f"{datetime.now()} reading GFF...", file=sys.stderr, flush=True)

gffFilterRegex = ""
if chromosome is not None:
    gffFilterRegex = f"^{chromosome}\s"
if target_gene is not None:
    gffFilterRegex += f".*{target_gene}"

if len(gffFilterRegex) > 0:
    if if_print:
        print(f"Filtering gffFile with {gffFilterRegex}")
    with tempfile.NamedTemporaryFile(mode="w") as filteredGffFile, (
        gzip.open(gffFile, "rt") if gffFile.endswith(".gz") else open(gffFile, "r")
    ) as inputGff:
        regex = re.compile(gffFilterRegex)
        for line in inputGff:
            if regex.match(line):
                filteredGffFile.write(line)
        filteredGffFile.flush()

        genes = gffReader.loadGenes(filteredGffFile.name)
else:
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

if len(outPrefix) > 0:
    outFile1 = outPrefix + outFile1
    outFile2 = outPrefix + outFile2
else:
    if target_gene is not None:
        outFile1 = target_gene + outFile1
        outFile2 = target_gene + outFile2
    else:
        outFile1 = "read" + outFile1
        outFile2 = "read" + outFile2

if if_print:
    if target_gene is not None:
        print(f"output file names {outFile1} {outFile2}")
    else:
        print(f"output file names {outFile1} {outFile2}")

out_path_folder = out_path
Path(out_path_folder).mkdir(parents=True, exist_ok=True)

if target_gene is not None:
    gene_folder = out_path + "/" + target_gene
    Path(gene_folder).mkdir(parents=True, exist_ok=True)
    out_path_folder = gene_folder
    Path(out_path_folder).mkdir(parents=True, exist_ok=True)

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

with tempfile.NamedTemporaryFile(mode="w") as twoBitInputFile:
    constructTwoBitInput(genes, twoBitInputFile.name)
    print(
        f"{datetime.now()} running 2bit",
        file=sys.stderr,
        flush=True,
    )
    twobitId_to_seq = run2bitBatch(twoBitDir, twoBitInputFile.name, genome2bit)

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


def get_region_str(gene):
    return f"{gene.getSubstrate().strip('chr')}:{gene.getBegin()}-{gene.getEnd()}"


#######

regions = set()
for idx, gene in enumerate(genes):
    counter += 1
    num_gene_gtf += 1
    region_str = get_region_str(gene)
    regions.add(region_str)

print(
    f"{datetime.now()} Got {num_gene_gtf} genes from gencode gtf",
    file=sys.stderr,
    flush=True,
)
if target_gene is not None:
    debug_print(f"target gene region: {regions}")

####### extract regions for genes from VCF file
# dict2: gene_to_variants {gene region}:{records from VCF}
#######

region_str_to_variants = tabix_regions(
    regions, variant_processor_chosen, vcfFile, comment_char="#"
)
# debug_print(region_str_to_variants)

####### extract sam data
# dict3: gene_to_qualityStr {gene region}:{quality string score from sam.gz}
#######
region_str_to_sam_data = tabix_regions(
    regions, sam_data_processor, samFile, comment_char="@", region_prefix="chr"
)
# debug_print(region_str_to_sam_data)

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
if target_gene is not None:
    list_fragLen = []
list_ratio = []

for gene in genes:
    mat = 0
    pat = 0
    region_str = get_region_str(gene)
    chrN = gene.getSubstrate()
    geneid = gene.getId()
    length = gene.longestTranscript().getLength()
    # DEBUGGING start
    if target_gene is not None:
        list_start1 = []
        list_start2 = []
        list_end1 = []
        list_end2 = []
        transcript = gene.longestTranscript()
        transcript.exons = transcript.getRawExons()
        debug_print(
            f"DEBUG... {geneid}, gene region: {region_str}, longest transcript length: {length}"
        )
    # DEBUGGING end
    if processed_genes > 0 and processed_genes % 100 == 0:
        sec_per_gene = (time.perf_counter_ns() - start_time_ns) / processed_genes / 1e9
        estimated_seconds_remaining = round(
            (len(genes) - processed_genes) * sec_per_gene
        )
        print(
            f"{datetime.now()} ...  processed {processed_genes} / {len(genes)} genes ({round(100*processed_genes/len(genes),2)}%) genes_per_sec: {round(1/sec_per_gene, 2)} ETA: {timedelta(seconds=estimated_seconds_remaining)}",
            file=sys.stderr,
            flush=True,
        )
    processed_genes += 1

    if not region_str in region_str_to_sam_data:
        debug_print(f"{chrN},gene: {geneid}, no mapped reads in SAM, skip")
        not_in_sam += 1
        continue
    if not region_str in region_str_to_variants:
        debug_print(f"{chrN},gene: {geneid}, no variants/records in VCF, skip")
        in_sam_not_in_vcf += 1
        continue
    # gene_to_variants {gene region}:{records from VCF}
    variants = region_str_to_variants[region_str]
    if all([not x.isHet() for x in variants]):
        debug_print("all SNPs are homozygous")
        continue
    # gene_to_qualityStr {gene region}:{quality string score from sam.gz}
    sam_data = region_str_to_sam_data[region_str]
    # within each gene, we obtain the reads information from SAM file
    pos1_tlen = [x[0] for x in sam_data]  # start pos of reads, quality string length
    qual_strs = [x[1] for x in sam_data]  # quality string
    if len(qual_strs) == 0 or len(pos1_tlen) == 0:
        continue
    # summarize start pos of reads, quality string length for this gene: {(89635, 124): 1, (89665, 187): 1}
    pos1_tlen_to_count = {}
    for x in pos1_tlen:
        pos1_tlen_to_count[x] = (
            pos1_tlen_to_count[x] + 1 if x in pos1_tlen_to_count else 1
        )
    # calculate the minimum quality length of quality string for this gene, maybe 75
    minQualLen = min([len(x) for x in qual_strs])
    ######## filter1: loop through each transcript, to
    transcript_to_fragLen = posTlen_to_fragLen(
        gene, pos1_tlen_to_count, minQualLen, if_print
    )

    if len(transcript_to_fragLen) == 0:
        #print("empty fragLen list!")
        continue
    recorded_genes += 1
    ######## write pat/mat transcripts
    paternal,maternal,transcriptIdToBiSNPpos, transcript_num = makingTranscript(
        gene, variants, if_print
    )
    qual_idx = 0
    #print(transcript_to_fragLen)
    candidate_transcripts = list(transcript_to_fragLen.keys())
    #print(candidate_transcripts)
    if target_gene is not None:
        debug_print(
                f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {len(transcript_to_fragLen)} available transcripts:"
            )
        for trans in candidate_transcripts:
            print(
                f"transcript {trans.getID()} CDS region: {trans.getCDSbeginEnd()}; transcript region: {trans.getBegin()}-{trans.getEnd()}; valid fragLen {transcript_to_fragLen[trans]}"
            )
        for j in transcriptIdToBiSNPpos:
            print(f">>>>> transcript ID {j},hetSNPs: {transcriptIdToBiSNPpos[j]}")
    candidate_transcript_pairs = [
        (x, next(filter(lambda y: y.getID() == x.getID(), maternal.transcripts)))
        for x in candidate_transcripts
    ]
    # read coverage per SNP we want * length of the longest transcript / minimum quality string length (for the gene)
    if if_random:
        DEPTH=int(DEPTH/2)
    numReads = int(float(DEPTH / minQualLen) * length)
    
    debug_print(
            f">>>>> gene id: {geneid}, gene region: {region_str}, #reads: {numReads},total #variants in VCF: {len(variants)}"
        )
    ##########################################################################
    for i in range(numReads):
        patTranscript, matTranscript = random.choice(candidate_transcript_pairs)
        transcript_length = matTranscript.getLength()
        frag_lens = transcript_to_fragLen[patTranscript]
        min_frag_len = min(frag_lens)

        assert min_frag_len <= transcript_length
        fragLen = random.choice(frag_lens)

        candidate_quals = list(filter(lambda x: len(x) <= fragLen, qual_strs))
        assert len(candidate_quals) > 0

        if if_print:
            if target_gene is not None:
                print("")
                print(
                    ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> simulated candidate transcripts"
                )
                print(
                    f"{i}-th reads, randomly chosen transcript {patTranscript.getID()},transcript length: {transcript_length}, randomly chosen fragLen: {fragLen}, chosen quality string length: {len(candidate_quals)}"
                )

        fwd_qual = random.choice(candidate_quals)
        rev_qual = random.choice(candidate_quals)
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
            patTranscript,
            matTranscript,
            fwd_qual,
            rev_qual,
            fragLen,
            if_print,
        )
        if target_gene is not None:
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
                mat += 1
            else:
                pat += 1
            identifier_random = "@SIM-" + str(nextReadID) + "-" + str(geneid)
            if if_print:
                print(
                    "%s,if maternal: %s,rec1 quality string length %s, forward strand length %s, rec2 quality strand length %s, reverse strand length %s!! "
                    % (
                        identifier_random,
                        str(if_mat),
                        len(fwd_qual),
                        fwd_LEN,
                        len(rev_qual),
                        rev_LEN,
                    )
                )
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
                str(identifier_random) + ":if_maternal" + str(if_mat) + ":FWD" + " /1",
                randomSeq,
                fwd_qual,
                OUT1,
            )
            nextReadID += 1
            printRead(
                str(identifier_random) + ":if_maternal" + str(if_mat) + ":REV" + " /1",
                randomSeq_rev,
                rev_qual,
                OUT2,
            )
        ################## equal haplotype copy simulator: both mat and pat has same number of reads
        else:
            mat += 1
            pat += 1
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

        if target_gene is not None:
            out = gene_folder
            Path(out).mkdir(parents=True, exist_ok=True)
            with open(out + "/trans_start1", "wb") as fp:
                pickle.dump(list_start1, fp)
            with open(out + "/trans_end1", "wb") as fp:
                pickle.dump(list_end1, fp)
            with open(out + "/trans_start2", "wb") as fp:
                pickle.dump(list_start2, fp)
            with open(out + "/trans_end2", "wb") as fp:
                pickle.dump(list_end2, fp)
        else:
            out = out_path
    ratio = mat / (pat + mat)
    list_ratio.append((pat, mat, ratio))

if target_gene is not None:
    print(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ratio of #mat/#total"
    )
    print(list_ratio)

if chromosome is not None:
    out = out_path
    with open(out + "/mat_ratio_" + str(chromosome) + ".pkl", "wb") as fp:
        pickle.dump(list_ratio, fp)

if target_gene is not None:
    with open(gene_folder + "/fragLen", "wb") as fp:
        pickle.dump(list_fragLen, fp)

print(
    f"{datetime.now()} DONE",
    file=sys.stderr,
    flush=True,
)
print("")
print(">> total geneID in gtf: %d" % (num_gene_gtf))
print(">> total ReadID processed: %d" % (nextReadID))
if matches==0:
    print(
        f">> # matches: {matches}, # mismatches: {mismatches}"
    )
else:
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
