#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
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

rex = Rex()
#######################################
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


def printRead(header, seq, qual, FH):
    print(header + "\n" + seq + "\n+\n" + qual, file=FH)


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
    return (
        refSeq,
        refSeq_rev,
        altSeq,
        altSeq_rev,
        LEN1,
        LEN2,
        start1,
        end1,
        start2,
        end2,
    )


if_print = False
def print_verbose(s):
    if if_print:
        print(s)


def loop_transcript_for_fragLen(
    th_read, gene, pos_idx, pos1_tlen_list, readLen, transcript_length, if_debug
):
    longest_transcript = gene.longestTranscript()
    fragLen, chosen_idx = pick_fragLen(
        th_read,
        pos1_tlen_list,
        pos_idx,
        longest_transcript,
        readLen,
        transcript_length,
        if_debug=if_debug,
    )
    if fragLen is not None:
        return fragLen, chosen_idx
    if if_debug:
        print(
            ">> {th_read} th read: could not find appropriate fraglen in longest transcript"
        )
    n = gene.getNumTranscripts()
    for i in range(n):
        new_idx = 0
        selected_transcript = gene.getIthTranscript(i)
        fragLen, _ = pick_fragLen(
            th_read,
            pos1_tlen_list,
            new_idx,
            selected_transcript,
            readLen,
            transcript_length,
            if_debug=if_debug,
        )
    if if_debug:
        if fragLen is None:
            print(
                ">> {th_read} th read: coudl not find appropriate fraglen in all transcript"
            )
    return fragLen, pos_idx


def pick_fragLen(
    th_read,
    pos1_tlen_list,
    pos_idx,
    transcript,
    readLen,
    transcript_length,
    if_debug=False,
):
    """
    output:
    fragLen
    pos_idx: position index of the tuple chosen for this fragLen
    """
    fragLen = None
    initial_pos_idx = pos_idx
    while True:
        pos1, tlen = pos1_tlen_list[pos_idx]
        pos_idx = (pos_idx + 1) % len(pos1_tlen_list)
        begin = transcript.mapToTranscript(pos1)
        end = transcript.mapToTranscript(pos1 + tlen)
        candidateFragLen = abs(end - begin)
        if if_debug:
            print(
                f"{th_read} th read: pos idx:{pos_idx} - fragLen:{candidateFragLen} - pos1:{pos1} - tlen:{tlen} - begin:{begin} - end:{end}"
            )
        if (
            begin >= 0
            and end >= 0
            and candidateFragLen >= readLen
            and candidateFragLen < transcript_length
        ):
            fragLen = candidateFragLen
            return fragLen, pos_idx
        if pos_idx == initial_pos_idx:
            return fragLen, pos_idx


def pick_random_Transcript(refGene, altGene, transcriptIdToBiSNPcount):
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


def twoBitID(chr, begin, end):
    return f"{chr}:{begin}-{end}"


def load_twoBit_TranscriptSeqs(gene, genome, path):
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
