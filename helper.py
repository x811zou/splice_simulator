#!/usr/bin/env python
#!/usr/bin/env python
# =========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
# =========================================================================
# This version has latest modification by Scarlett
# =========================================================================
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
    generators,
    nested_scopes,
    with_statement,
)
from builtins import (
    bytes,
    dict,
    int,
    list,
    object,
    range,
    str,
    ascii,
    chr,
    hex,
    input,
    next,
    oct,
    open,
    pow,
    round,
    super,
    filter,
    map,
    zip,
)

# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import random
import copy
import gzip
import os
from collections import defaultdict
from misc_tools.GffTranscriptReader import GffTranscriptReader
from misc_tools.Pipe import Pipe
from misc_tools.ConfigFile import ConfigFile
from misc_tools.Translation import Translation
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
#   util-dir = /hpc/home/bmajoros/twobit
#   genome = hg19.2bit
#   aligned-rna = HG00096.sam
#   fragment-lengths = fragmentlengths.txt
#   vcf = HG00096.vcf.gz
#   gff = gencode1000.gff
#   original-read-len = 75
# ============================================

rex = Rex()

# CLASSES:


class Variant:
    def __init__(self, fields):
        self.Chr = fields[0]
        self.genomicPos = int(fields[1]) - 1  # convert to 0-based
        self.ref = fields[3]
        self.alt = fields[4]
        genotype = fields[9]
        self.genotype = None
        if rex.find("(\d)\|(\d)", genotype):
            self.genotype = (int(rex[1]), int(rex[2]))
            if self.genotype[0] not in {0, 1} or self.genotype[1] not in {0, 1}:
                self.genotype = None

    def isOK(self):
        return self.genotype is not None

    def isHet(self):
        return self.genotype[0] != self.genotype[1]

    def isHomozygousAlt(self):
        return self.genotype[0] > 0 and self.genotype[1] > 0


def readFrom2bit(genome, path, Chr, start, end):
    cmd = (
        os.path.join(path, "twoBitToFa")
        + " -seq="
        + Chr
        + " -start="
        + str(start)
        + " -end="
        + str(end)
        + " "
        + genome
        + " stdout"
    )
    output = Pipe.run(cmd)
    fields = output.rstrip().split("\n")
    seq = ""
    for i in range(1, len(fields)):
        seq += fields[i]
    return seq


def twoBitID(chr, begin, end):
    return f"{chr}:{begin}-{end}"


def pickTranscript(refGene, altGene, transcriptIdToBiSNPcount):
    n = refGene.getNumTranscripts()
    i = random.randrange(n)
    num = transcriptIdToBiSNPcount[refGene.getIthTranscript(i).getID()]
    return (
        refGene.getIthTranscript(i),
        altGene.getIthTranscript(i),
        i,
        refGene.getIthTranscript(i).getID(),
        num,
    )


def printRead(header, seq, qual, FH):
    print(header + "\n" + seq + "\n+\n" + qual, file=FH)


def loadFragLens(filename):
    array = []
    with open(filename, "rt") as IN:
        for line in IN:
            L = int(line)
            array.append(L)
    return array


def loadTranscriptSeqs(gene, genome, path):
    Chr = gene.getSubstrate()
    strand = gene.getStrand()

    twobitReads = []
    for transcript in gene.transcripts:
        transcript.exons = transcript.getRawExons()
        for exon in transcript.rawExons:
            twobitReads.append(f"{twoBitID(Chr, exon.begin, exon.end)}\n")

    with open("/tmp/2bitinput", "w") as f:
        f.writelines(twobitReads)

    id_to_seq = run2bitBatch(path, "/tmp/2bitinput", genome)

    for transcript in gene.transcripts:
        transcript.sequence = ""
        transcript.exons = transcript.getRawExons()
        for exon in transcript.rawExons:
            exon.sequence = id_to_seq[twoBitID(Chr, exon.begin, exon.end)]
            exon.sequence = exon.sequence.upper()
            if strand == "-":
                exon.sequence = Translation.reverseComplement(exon.sequence)
            if "N" in exon.sequence:
                raise Exception(
                    "N FOUND: Chr="
                    + Chr
                    + "begin="
                    + str(exon.begin)
                    + "end="
                    + str(exon.end)
                )
            transcript.sequence += exon.sequence


################## modified function ##################
def makeAltCopy(gene, haplotype, variants, if_print=False):
    #########
    # this function is used to create a set of REF copy of a gene transcrips, or ALT copy
    # REF copy (haplotype == 0): replace the REF allele for variants in ref gencode filtered transcript to actual REF allele in VCF
    # ALT copy (haplotype == 1): replace the REF allele for variants in ref gencode filtered transcript to actual ALT allele in VCF
    # allele_in_vcf : allele in VCF
    # allele_in_ref : allele in reference transcript (gencode GTF hg19.v19)
    # match: allele in gencode ref transcript (reversed if in reverse strand) MATCH the allele in VCF
    # mismatch: ... NOT MATCH ...
    #########
    # global matches
    # global mismatches
    altGene = copy.deepcopy(gene)
    transcriptIdToBiSNPcount = {}
    transcriptIdToBiSNPpos = defaultdict(set)
    if if_print:
        if haplotype == 0:
            print("    ref/paternal")
        else:
            print("    alt/maternal")
    transcript_num = 0
    # loop through each ref gencode transcript, create a REF/ALT copy based on VCF
    for transcript in altGene.transcripts:
        # print(f"transcript {transcript_num}len {len(transcript.sequence)}")
        array = list(transcript.sequence)
        transcript_num += 1
        num = 0
        # loop through each bi-allelic SNP from VCF
        for variant in variants:
            if_rev = False
            trans_pos = transcript.mapToTranscript(variant.genomicPos)
            # genom_pos = transcript.mapToGenome(trans_pos)
            allele_in_ref = array[trans_pos]
            if len(variant.ref) == 1 and len(variant.alt) == 1 and trans_pos >= 0:
                if_bi = True
                if if_print:
                    print(
                        " genomic pos is %d, trans pos is %d"
                        % (variant.genomicPos, trans_pos)
                    )
                if variant.genotype[0] != variant.genotype[1]:
                    transcriptIdToBiSNPpos[transcript.getID()].add(trans_pos)
                    num += 1
                ########################## use REF/ALT allele in VCF to modify the reference transcript
                if variant.genotype[haplotype] != 0:
                    allele_in_vcf = variant.alt
                else:
                    allele_in_vcf = variant.ref
                ########################## reverse the allele if it is in the reverse strand
                if gene.getStrand() == "-":
                    if_rev = True
                    allele_in_ref = Translation.reverseComplement(allele_in_ref)
                    array[trans_pos] = Translation.reverseComplement(allele_in_vcf)
                    allele_check = array[trans_pos]
                else:
                    array[trans_pos] = allele_in_vcf
                    allele_check = array[trans_pos]
                ##########################  check match/mismatch
                # if allele_check == allele_in_ref:
                #    matches += 1
                # else:
                #    mismatches += 1
                ##########################
                if if_print is not False:
                    print(
                        "        %dth transcript, %dth bi-allelic SNP"
                        % (transcript_num, num)
                    )
                    print(
                        "            > if reverse strand: %s, variant trans pos: %s, haplotype: %s, in ref genome: %s, allele check %s, ref: %s, alt: %s, write_in_sequence: %s"
                        % (
                            str(if_rev),
                            trans_pos,
                            variant.genotype,
                            allele_in_ref,
                            allele_check,
                            variant.ref,
                            variant.alt,
                            array[trans_pos],
                        )
                    )
        transcript.sequence = "".join(array)
        transcriptIdToBiSNPcount[transcript.getID()] = num
        # print(transcriptIdToBiSNPpos)
        if if_print:
            print(
                "        >>>>>> %dth transcript, in total %d bi-allelic SNP"
                % (transcript_num, num)
            )
            print(" ")
    return (
        altGene,
        transcriptIdToBiSNPcount,
        transcriptIdToBiSNPpos,
    )


def run2bitBatch(path, inputFile, genome):
    cmd = f"{os.path.join(path, 'twoBitToFa')} -seqList={inputFile} {genome} stdout"

    print(
        f"{datetime.now()} twobit command {cmd}",
        file=sys.stderr,
        flush=True,
    )

    output = Pipe.run(cmd)
    lines = output.rstrip().split("\n")

    id_to_seq = {}
    id = None
    seq = None
    for line in lines:
        if line.startswith(">"):
            if id is not None:
                id_to_seq[id] = seq
            id = line.strip(">")
            seq = ""
            continue
        seq += line

    if id is not None:
        id_to_seq[id] = seq

    return id_to_seq


def constructTwoBitInput(genes, file):
    count = 0
    with open(file, "w") as f:
        for gene in genes:
            chr = gene.getSubstrate()
            for transcript in gene.transcripts:
                transcript.exons = transcript.getRawExons()
                for exon in transcript.rawExons:
                    f.write(f"{twoBitID(chr, exon.begin, exon.end)}\n")
                    count += 1

    print(
        f"{datetime.now()} generated two bit input of {count} lines",
        file=sys.stderr,
        flush=True,
    )


def annotate_transcripts(genes, id_to_seq):
    for gene in genes:
        Chr = gene.getSubstrate()
        strand = gene.getStrand()

        for transcript in gene.transcripts:
            transcript.sequence = ""
            transcript.exons = transcript.getRawExons()
            for exon in transcript.rawExons:
                exon.sequence = id_to_seq[twoBitID(Chr, exon.begin, exon.end)]
                exon.sequence = exon.sequence.upper()
                if strand == "-":
                    exon.sequence = Translation.reverseComplement(exon.sequence)
                if "N" in exon.sequence:
                    raise Exception(
                        "N FOUND: Chr="
                        + Chr
                        + "begin="
                        + str(exon.begin)
                        + "end="
                        + str(exon.end)
                    )
                transcript.sequence += exon.sequence


def chunk_iter(iter, n):
    """Yield successive n-sized chunks from iter."""
    res = []
    try:
        while True:
            while len(res) < n:
                v = next(iter)
                res.append(v)
            yield res
            res = []
    except StopIteration:
        if len(res) > 0:
            yield res
