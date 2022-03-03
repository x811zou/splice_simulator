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
from misc_tools.GffTranscriptReader import GffTranscriptReader
from misc_tools.Pipe import Pipe
from misc_tools.ConfigFile import ConfigFile
from misc_tools.Rex import Rex
from Bio.Seq import Seq
from pathlib import Path
from datetime import datetime
from misc_tools.Translation import Translation

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
    # not het and not homoalt
    if not (variant.isHet() and not variant.isHomozygousAlt()):
        return None
    if not (len(variant.ref) == 1 and len(variant.alt) == 1):
        return None
    return variant

def sam_data_processor(line):
    #########
    # this function is used to fetch data from .sam.gz
    # returns ( (start_pos, template_len), quality_string )
    # reference: https://samtools.github.io/hts-specs/SAMv1.pdf
    # fields[3]: start POS
    # fields[8]: template length
    # fields[10]: quality string
    #########
    fields = line.rstrip().split()
    if len(fields) < 11:
        return None

    pos1 = int(fields[3])
    tlen = int(fields[8])
    quality_string = fields[10]

    if tlen <= 0:
        return None

    result = ((pos1, tlen), quality_string)
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
    # print(f"L{L}-Lstart:{L_start} - start2:{start2}")
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


def posTlen_to_fragLen(gene, pos1_tlen_list, readLen):
    fragLen_list = []
    transcript_len = []
    for pos1, tlen in pos1_tlen_list:
        # print(f">>> {pos1}-{tlen}:")
        n = gene.getNumTranscripts()
        for i in range(n):
            selected_transcript = gene.getIthTranscript(i)
            begin = selected_transcript.mapToTranscript(pos1)
            end = selected_transcript.mapToTranscript(pos1 + tlen)
            candidateFragLen = abs(end - begin)
            transcript_len.append(selected_transcript.getLength())
            if begin >= 0 and end >= 0 and candidateFragLen >= readLen:
                fragLen_list.append(candidateFragLen)
                # print(f"   {i}th transcript - {candidateFragLen}")
            # else:
            #     print(
            #         f"    skipped --> {i}th transcript - {candidateFragLen} = {end} - {begin}"
            #     )
    if len(fragLen_list) > 0:
        if min(fragLen_list) > max(transcript_len):
            fragLen_list = []
    return fragLen_list


def pick_fragLen(fragLens, max_qual_len, transcript_length):
    # print(fragLens)
    random.shuffle(fragLens)
    for fragLen in fragLens:
        # print(fragLen)
        if fragLen < transcript_length and fragLen >= max_qual_len:
            return fragLen
    return None


def pick_fragLen_X(
    th_read,
    pos1_tlen_list,
    pos_idx,
    transcript,
    readLen,
    max_qual_len,
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
            and candidateFragLen >= max_qual_len
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
