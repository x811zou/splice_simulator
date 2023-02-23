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
    int,
    list,
    range,
    str,
    next,
    open,
    map,
)
import logging

# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import random
import os
import tempfile
from misc_tools.Pipe import Pipe
from misc_tools.Rex import Rex
from Bio.Seq import reverse_complement
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
        if rex.find("(\d)[\|\/](\d)", genotype):
            self.genotype = (int(rex[1]), int(rex[2]))
            if self.genotype[0] not in {0, 1} or self.genotype[1] not in {0, 1}:
                self.genotype = None
        # cat HG00096.no_chr.content.SNPs.filtered.vcf.gz | gunzip | awk '{print $10}' | sort | uniq -c
        # 74304603 0|0
        # 1071282 0|1
        # 1065885 1|0
        # 1376362 1|1
        # cat 123375.no_chr.content.SNPs.filtered.vcf.gz | gunzip | awk '{print $10}' | sort | uniq -c
        # 21194 ./.
        # 2917072 0/0
        # 2354362 0/1
        # cat NA12878.no_chr.content.SNPs.filtered.vcf.gz | gunzip | awk '{print $10}' | sort | uniq -c
        # 5947 0/1
        # 976007 0|1
        # 986082 1|0
        # 1289007 1|1

    def isOK(self):
        return self.genotype is not None

    def isHet(self):
        return self.genotype[0] != self.genotype[1]

    def isHomozygousAlt(self):
        return self.genotype[0] > 0 and self.genotype[1] > 0


def printRead(header, seq, qual, FH):
    print(header + "\n" + seq + "\n+\n" + qual, file=FH)


def tabix_regions(
    regions, line_processor, target_file_path, comment_char="#", region_prefix=""
):
    region_to_results = {}

    print(
        f"{datetime.now()} Start tabix extraction of {len(regions)} regions from file {target_file_path}"
    )

    regions = list(map(lambda x: region_prefix + x, regions))

    if len(regions) > 1000:
        with tempfile.NamedTemporaryFile(mode="w") as file:
            for region in regions:
                chr, rest = region.split(":")
                start, end = rest.split("-")
                file.write(f"{chr}\t{start}\t{end}\n")
            file.flush()
            command = (
                f"tabix --separate-regions {target_file_path} --regions {file.name}"
            )
            output = Pipe.run(command)
    else:
        region_batch = " ".join(regions)
        command = f"tabix --separate-regions {target_file_path} {region_batch}"
        output = Pipe.run(command)
    if len(output) == 0:
        return region_to_results

    lines = output.split("\n")
    records = []
    for line in lines:
        # print(line)
        if len(line) == 0:
            continue
        # start accumulating new region
        if line.startswith(comment_char):
            region_str = line[1:].strip(region_prefix)
            records = []
            region_to_results[region_str] = records
            continue

        result = line_processor(line)
        if result is not None:
            records.append(result)

    print(
        f"{datetime.now()} Got {len(region_to_results)} / {len(regions)} regions with data"
    )

    return region_to_results


def variant_processor_hets(line):
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


def variant_processor_SNPs(line):
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


def simRead_patmat(
    input_ref_seq,
    input_alt_seq,
    input_ref_seq_rev_comp,
    input_alt_seq_rev_comp,
    qual1,
    qual2,
    fragLen,
):
    #####################################################################
    # transcript length
    L = len(input_ref_seq)
    qual1_len = len(qual1)
    qual2_len = len(qual2)
    # transcript seq index start/end
    L_end = L - 1
    L_start = 0
    # fragLen: actual read length drawn from SAM
    if L_end < fragLen or L_end < qual1_len or L_end < qual2_len:
        return None
    # transcript coord
    lastStart = min(L_end - fragLen, L_end - qual1_len, L_end - qual2_len)  # 20
    start1 = random.randrange(lastStart + 1)  # 10
    end1 = start1 + qual1_len  # rec1.readLen  # 10+75 = 85
    LEN1 = abs(end1 - start1)
    end2 = start1 + fragLen  # 10+80 = 90
    start2 = end2 - qual2_len  # rec2.readLen  # 90-75 = 15
    LEN2 = abs(end2 - start2)

    assert start1 >= L_start
    assert end1 <= L_end
    assert start2 >= L_start
    assert end2 <= L_end
    assert len(qual1) == LEN1
    assert len(qual2) == LEN2

    ######## forward strand, same sequence pos for mat/aptf fov
    simulated_ref_seq = input_ref_seq[start1:end1]
    simulated_alt_seq = input_alt_seq[start1:end1]
    ######## reverse strand, same sequence pos for mat/apt rev
    simulated_ref_seq_rev = input_ref_seq_rev_comp[L - end2 : L - start2]
    simulated_alt_seq_rev = input_alt_seq_rev_comp[L - end2 : L - start2]

    assert len(qual1) == len(simulated_ref_seq)
    assert len(qual2) == len(simulated_ref_seq_rev)

    return (
        simulated_ref_seq,
        simulated_ref_seq_rev,
        simulated_alt_seq,
        simulated_alt_seq_rev,
    )


def posTlen_to_fragLen(gene, pos1_tlen_to_count, readLen):
    transcript_id_to_mapped_lengths = {}
    # logging.debug(
    #     ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  filtering : valid start/end pos of SAM records map to transcript"
    # )
    for i in range(gene.getNumTranscripts()):
        transcript = gene.getIthTranscript(i)
        mapped_lengths = []

        pos_cache = {}

        def mapPos(pos):
            nonlocal transcript
            nonlocal pos_cache
            if pos not in pos_cache:
                pos_cache[pos] = transcript.mapToTranscript(pos)
            return pos_cache[pos]

        for pos1_tlen, count in pos1_tlen_to_count.items():
            pos1, tlen = pos1_tlen
            # get transcript coordiante position
            begin = mapPos(pos1)
            # filter out sites whose both start/end could not map to transcript
            if begin < 0:
                continue
            end = mapPos(pos1 + tlen)
            if end < 0:
                continue
            mapped_length = abs(end - begin)
            # logging.debug(
            #     f"transcript {i+1} - {transcript.getID()}:{pos1_tlen} mapped start,end: ({begin},{end}) fragLen {mapped_length}"
            # )
            if mapped_length < readLen:
                continue

            mapped_lengths.extend([mapped_length for i in range(count)])

        if len(mapped_lengths) > 0:
            mapped_lengths.sort()
            transcript_id_to_mapped_lengths[transcript.getID()] = mapped_lengths

    return transcript_id_to_mapped_lengths


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
