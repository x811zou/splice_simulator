#!/usr/bin/env python
# =========================================================================
# This version has latest modification by Scarlett at 02/22/2022
# =========================================================================

from distutils.debug import DEBUG
from xmlrpc.client import boolean
from xxlimited import Xxo
import argparse
from helper import (
    constructTwoBitInput,
    run2bitBatch,
    pickTranscript,
    makeAltCopy,
    annotate_transcripts,
    printRead,
    Variant,
    chunk_iter,
)

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


def tabix_regions(regions, line_processor, target_file_path=None, comment_char="#"):
    CHUNK_SIZE = 1000
    regions_processed = 0
    region_to_results = {}

    print(
        f"{datetime.now()} Start tabix extraction of {len(regions)} regions from file {target_file_path}"
    )

    for x in chunk_iter(iter(regions), CHUNK_SIZE):
        """
        x stands for gene level chr: start-end
        look up 1000 genes at a time
        """
        region_batch = " ".join(x)
        cmd = f"tabix --separate-regions {target_file_path} {region_batch}"
        output = Pipe.run(cmd)
        if len(output) == 0:
            continue

        lines = output.split("\n")
        records = []
        for line in lines:
            if len(line) == 0:
                continue
            # start accumulating new region
            if line.startswith(comment_char):
                region_str = line[1:]
                records = []
                region_to_results[region_str] = records
                continue

            result = line_processor(line)
            if result is not None:
                records.append(result)

        regions_processed += len(x)
        print(
            f"{datetime.now()} ... finished {regions_processed} / {len(regions)} regions"
        )

    print(f"{datetime.now()} Got {len(region_to_results)} regions with data")

    return region_to_results


def variant_processor(line):
    #########
    # this function is used to fetch variants from .vcf
    #########
    fields = line.rstrip().split()
    if len(fields) != 10:
        raise Exception("Expecting 10 fields in VCF file")

    variant = Variant(fields)
    if not variant.isOK():
        return None
    if not variant.isHet() and not variant.isHomozygousAlt():
        return None
    return variant


def quality_string_processor(line):
    #########
    # this function is used to fetch quality score from .sam.gz
    # reference: https://samtools.github.io/hts-specs/SAMv1.pdf
    # fields[8]: template length
    # fields[10]: quality string
    #########
    fields = line.rstrip().split()

    # skip reads with 0 length
    if int(fields[8]) == 0:
        return None
    # skip reads without quality string
    if len(fields) < 11:
        return None

    quality_string = fields[10]
    return quality_string


def read_length_processor(line):
    #########
    # this function is used to fetch read start pos and template length from .sam.gz
    # fields[3]: start POS
    # fields[8]: template length
    #########
    fields = line.rstrip().split()
    # skip reads whose template length <= 0
    # ['ERR188166.25392658', '147', 'chr6', '8558820', '255', '57M95982N18M', '=', '8436351', '-218526', 'GATGTCATTGAACTCTTAAGTGCAAGATGAAACAAGTCTTTCTGGGGTTCTAAGTAGAAGTGATCTACCTTTCTA', 'IJIIJJGHCDGHHGGF@BGHEHBCEFEHECGHHGGGGDEDHEEGEHHCHEHGHBAHEHGIHHH?FHHDDFFDC@@', 'MD:Z:75', 'PG:Z:MarkDuplicates', 'NH:i:1', 'HI:i:1', 'jI:B:i,8558877,8654858', 'NM:i:0', 'jM:B:c,0', 'nM:i:0', 'AS:i:140', 'XS:A:+']
    # print(fields)
    pos1 = int(fields[3])
    tlen = int(fields[8])
    # skip duplicated ones
    if tlen <= 0:
        return None
    result = tuple((pos1, abs(tlen)))
    # print(result)
    return result


def simRead_patmat(refTranscript, altTranscript, qual1, qual2, fragLen, readLen):
    #####################################################################
    # transcript length
    L = len(refTranscript.sequence)
    # transcript seq index start/end
    L_end = L - 1
    L_start = 0

    # fragLen: actual read length drawn from SAM
    if L_end < fragLen or L_end < readLen or L_end < len(qual1) or L_end < len(qual2):
        return (None, None, None, None, None, None)
    # transcript coord
    lastStart = min(L_end - fragLen, L_end - len(qual1), L_end - len(qual2))  # 20
    start1 = random.randrange(lastStart + 1)  # 10
    end1 = start1 + len(qual1)  # rec1.readLen  # 10+75 = 85
    LEN1 = abs(end1 - start1)
    end2 = start1 + fragLen  # 10+80 = 90
    start2 = end2 - len(qual2)  # rec2.readLen  # 90-75 = 15
    LEN2 = abs(end2 - start2)
    print(f"L{L}-Lstart:{L_start} - start2:{start2}")
    assert start1 >= L_start
    assert end1 <= L_end
    assert start2 >= L_start
    assert end2 <= L_end
    assert len(qual1) == LEN1
    assert len(qual2) == LEN2

    # print(
    #     f"qual1 {len(qual1)} qual2{len(qual2)} len1 {LEN1} len2 {LEN2} fwdSeq length{len(refSeq)} revSeq length {len(refSeq_rev)}"
    # )

    ######## forward strand, same sequence pos for mat/aptf fov
    refSeq = refTranscript.sequence[start1:end1]
    altSeq = altTranscript.sequence[start1:end1]
    ######## reverse strand, same sequence pos for mat/apt rev
    refSeq_rev = Seq(refTranscript.sequence[start2:end2]).reverse_complement()
    altSeq_rev = Seq(altTranscript.sequence[start2:end2]).reverse_complement()
    assert len(qual1) == len(refSeq)
    assert len(qual2) == len(refSeq_rev)
    return (refSeq, refSeq_rev, altSeq, altSeq_rev, LEN1, LEN2)


# GLOBAL VARIABLES:
# matches = 0  # number of sites containing alt or ref allele exactly as in gencode ref transcript
# mismatches = 0  # number of sites having either or neither ref nor alt allele as in gencode ref transcript
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


out_path_folder = out_path + "/" + "simulated_fastq"
Path(out_path).mkdir(parents=True, exist_ok=True)
Path(out_path_folder).mkdir(parents=True, exist_ok=True)

# Load GFF and fragment lengths
gffReader = GffTranscriptReader()
print(f"{datetime.now()} reading GFF...", file=sys.stderr, flush=True)
genes = gffReader.loadGenes(gffFile)
genes.sort(key=lambda gene: (gene.getSubstrate(), gene.getBegin()))
print(f"{datetime.now()} done reading GFF...", file=sys.stderr, flush=True)


if chromosome is not None:
    print(f"Looking at specific chromosome: {chromosome}")
    genes = list(filter(lambda x: x.getSubstrate() == chromosome, genes))


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

for gene in genes:
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
            if if_print:
                print(f"{chrN},gene: {geneid}, no mapped reads in SAM, skip")
            not_in_sam += 1
            continue
        if not region_str in region_str_to_variants:
            if if_print:
                print(f"{chrN},gene: {geneid}, no variants/records in VCF, skip")
            in_sam_not_in_vcf += 1
            continue
        if not region_str in region_str_to_fragLen:
            if if_print:
                print(f"{chrN},gene: {geneid}, no fragment length in SAM, skip")
            continue

        variants = region_str_to_variants[region_str]
        qual_strs = region_str_to_quality_strings[region_str]
        pos1_tlen = region_str_to_fragLen[region_str]
        # print(pos1_tlen)

        if len(qual_strs) == 0 or len(pos1_tlen) == 0:
            continue

        recorded_genes += 1
        if if_print:
            print(
                "%s,%s,#reads: %s,total #variants in VCF: %d"
                % (region_str, geneid, numReads, len(variants))
            )

        maternal, transcriptIdToBiSNPcount, transcriptIdToBiSNPpos = makeAltCopy(
            gene, 1, variants
        )
        paternal, _, _ = makeAltCopy(gene, 0, variants)
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
            ) = pickTranscript(maternal, paternal, transcriptIdToBiSNPpos)
            transcript_length = matTranscript.getLength()
            # extract quality string from qual_strs in order
            ##########
            # print(f"{qual_idx}-{len(qual_strs)}")
            fwd_qual = qual_strs[qual_idx]
            qual_idx = (qual_idx + 1) % len(qual_strs)
            rev_qual = qual_strs[qual_idx]
            qual_idx = (qual_idx + 1) % len(qual_strs)
            # print(f"{qual_idx}-{len(qual_strs)}")
            ##########
            # randomly draw actual-fragment-length from file
            # fragLen = fragLens[random.randrange(len(fragLens))]
            # print(pos_idx)
            # print(f"{geneid} - {len(pos1_tlen)} - {pos_idx}-{pos1_tlen[pos_idx]}")
            fragLen = None
            initial_pos_idx = pos_idx
            while True:
                pos1, tlen = pos1_tlen[pos_idx]
                pos_idx = (pos_idx + 1) % len(pos1_tlen)
                begin = patTranscript.mapToTranscript(pos1)
                end = patTranscript.mapToTranscript(pos1 + tlen)

                # begin = transcript.mapToTranscript(pos1)
                # end = transcript.mapToTranscript(pos1 + tlen)
                candidateFragLen = abs(end - begin)
                print(
                    f"longest transcript idx:{pos_idx} - fragLen:{candidateFragLen} - pos1:{pos1} - tlen:{tlen} - begin:{begin} - end:{end}"
                )
                if begin >= 0 and end >= 0 and candidateFragLen < transcript_length:
                    fragLen = candidateFragLen
                    break
                if pos_idx == initial_pos_idx:
                    break

            if fragLen is None:
                print("coudl not find appropriate fraglen in longest transcript")
                n = gene.getNumTranscripts()
                # i = random.randrange(n)
                for i in range(n):
                    select_transcript = gene.getIthTranscript(i)
                    initial_pos_idx = pos_idx
                    while True:
                        pos1, tlen = pos1_tlen[pos_idx]
                        pos_idx = (pos_idx + 1) % len(pos1_tlen)
                        begin = patTranscript.mapToTranscript(pos1)
                        end = patTranscript.mapToTranscript(pos1 + tlen)

                        # begin = transcript.mapToTranscript(pos1)
                        # end = transcript.mapToTranscript(pos1 + tlen)
                        candidateFragLen = abs(end - begin)
                        print(
                            f"transcript {i} idx:{pos_idx} - fragLen:{candidateFragLen} - pos1:{pos1} - tlen:{tlen} - begin:{begin} - end:{end}"
                        )
                        if (
                            begin >= 0
                            and end >= 0
                            and candidateFragLen < transcript_length
                        ):
                            fragLen = candidateFragLen
                            break
                        if pos_idx == initial_pos_idx:
                            break

            if fragLen is None:
                print("coudl not find appropriate fraglen in all transcript")
                continue

            # list_fragLen.append(fragLen)
            print(
                f"{chrN}:{geneid}-{i} reads: pos indx {pos_idx} transcript len {transcript_length} - frag len {fragLen}"
            )
            # print(f"{geneid} - {pos1_tlen[pos_idx]} - {pos_idx}")
            # print(f"begin {transcript.mapToTranscript(pos1)}")
            # print(f"end {transcript.mapToTranscript(pos1+tlen)}")
            # print(f"fragment length {fragLen} - {pos_idx} - {len(pos1_tlen)}")
            ##########
            # make sure transcript length >= randomly drawn actual-fragment-length
            # if transcript_length < fragLen:
            #     fragLen = transcript_length

            # simulate reads for paternal and maternal copy
            (patSeq, patSeq_rev, matSeq, matSeq_rev, fwd_LEN, rev_LEN) = simRead_patmat(
                patTranscript, matTranscript, fwd_qual, rev_qual, fragLen, readLen
            )

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

        if processed_genes % 500 == 0:
            print(
                f"{datetime.now()} ... processed {processed_genes} / {len(genes)} genes : {round(100*processed_genes/len(genes),2)} %",
                file=sys.stderr,
                flush=True,
            )

# with open(out_path + "/fragLen", "wb") as fp:
#     pickle.dump(list_fragLen, fp)
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
