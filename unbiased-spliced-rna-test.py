#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import logging
import re
import tempfile

import argparse
import tracemalloc

import numpy as np
import psutil
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
)

import sys
import random
import gzip
from misc_tools.GffTranscriptReader import GffTranscriptReader
from misc_tools.Translation import Translation


def loadGenes(gffFile, chromosome=None, gene=None):
    gffReader = GffTranscriptReader()

    gffFilterRegex = ""
    if chromosome is not None:
        gffFilterRegex = f"^{chromosome}\s"
    if gene is not None:
        gffFilterRegex += f".*{gene}"

    if len(gffFilterRegex) > 0:
        logging.debug(f"Filtering gffFile with /{gffFilterRegex}/")
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

    if chromosome is not None:
        genes = list(filter(lambda x: x.getSubstrate() == chromosome, genes))
    if gene is not None:
        genes = list(filter(lambda x: x.getId() == gene, genes))

    return genes


def loadRegionData(vcfFile, samFile, genes, variant_processor):
    regions = set([getRegionStr(gene) for gene in genes])

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
        regions, sam_data_processor, samFile, comment_char="@", region_prefix="chr"
    )

    return region_str_to_variants, region_str_to_sam_data


def getRegionStr(gene):
    return f"{gene.getSubstrate().strip('chr')}:{gene.getBegin()}-{gene.getEnd()}"


def getMemoryUsageMB():
    return psutil.Process().memory_info().rss / 1024**2


# returns list<(transcript_id, pat_seq, mat_seq)>
def makeModifiedSequences(gene, variants, transcript_ids):
    result = []
    for transcript_id in transcript_ids:
        transcript = next(x for x in gene.transcripts if x.getID() == transcript_id)

        pat_seq = list(transcript.sequence)
        mat_seq = list(transcript.sequence)

        for variant in variants:
            trans_pos = transcript.mapToTranscript(variant.genomicPos)
            if trans_pos < 0:
                continue

            refallele_in_vcf = variant.ref
            altallele_in_vcf = variant.alt
            if (
                gene.getStrand() == "-"
            ):  # reverse the allele if it is in the reverse strand
                refallele_in_vcf = Translation.reverseComplement(refallele_in_vcf)
                altallele_in_vcf = Translation.reverseComplement(altallele_in_vcf)

            pat_seq[trans_pos] = (
                refallele_in_vcf if variant.genotype[0] == 0 else altallele_in_vcf
            )
            mat_seq[trans_pos] = (
                refallele_in_vcf if variant.genotype[1] == 0 else altallele_in_vcf
            )

        result.append((transcript_id, "".join(pat_seq), "".join(mat_seq)))

    return result


def simulateRead(
    candidate_sequence_items,
    transcript_id_to_fragLen,
    qual_strs,
    qual_str_lens,
    use_random_reads,
):
    transcript_id, pat_seq, mat_seq = random.choice(candidate_sequence_items)
    transcript_length = len(pat_seq)
    frag_lens = transcript_id_to_fragLen[transcript_id]
    min_frag_len = frag_lens[0]

    assert min_frag_len <= transcript_length
    frag_len = random.choice(frag_lens)

    candidate_qual_str_count = np.searchsorted(qual_str_lens, frag_len, side="right")
    assert candidate_qual_str_count > 0
    candidate_quals = qual_strs[:candidate_qual_str_count]

    fwd_qual = random.choice(candidate_quals)
    rev_qual = random.choice(candidate_quals)
    # simulate reads for paternal and maternal copy
    simulated_sequences = simRead_patmat(pat_seq, mat_seq, fwd_qual, rev_qual, frag_len)
    if simulated_sequences is None:
        return None

    (
        simulated_ref_seq,
        simulated_ref_seq_rev,
        simulated_alt_seq,
        simulated_alt_seq_rev,
    ) = simulated_sequences

    ################## random haplotype simulator: randomly generate a mat or pat copy
    if use_random_reads:
        is_mat = random.random() >= 0.5
        if is_mat:
            fwd_seq = simulated_alt_seq
            rev_seq = simulated_alt_seq_rev
        else:
            fwd_seq = simulated_ref_seq
            rev_seq = simulated_ref_seq_rev

        return [(is_mat, fwd_seq, fwd_qual, rev_seq, rev_qual)]
    ################## equal haplotype copy simulator: both mat and pat has same number of reads
    return [
        (False, simulated_ref_seq, fwd_qual, simulated_ref_seq_rev, rev_qual),
        (True, simulated_alt_seq, fwd_qual, simulated_alt_seq_rev, rev_qual),
    ]


def runSimulation(
    genes,
    region_str_to_variants,
    region_str_to_sam_data,
    use_random_reads,
    read_depth,
    max_genes,
    write_fwd,
    write_rev,
):
    #######
    # for each gene, generate reads using quality string from matching genes in SAM.GZ
    #######
    n_break = 0
    recorded_genes = 0
    not_in_sam = 0
    in_sam_not_in_vcf = 0
    num_genes_with_all_transcripts_filtered = 0
    num_genes_with_no_qual_strs_or_pos_tlens = 0

    nextReadID = 0
    n_genes = len(genes)
    for i, gene in enumerate(genes):
        if i % 100 == 0:
            logging.info(f"processed {i} / {n_genes} genes")
            logging.info(f"memory usage {getMemoryUsageMB()} MB")
        if max_genes > 0 and i >= max_genes:
            break
        mat_reads = 0
        pat_reads = 0
        region_str = getRegionStr(gene)
        chrN = gene.getSubstrate()
        geneid = gene.getId()
        length = gene.longestTranscript().getLength()

        if not region_str in region_str_to_sam_data:
            logging.debug(f"{chrN},gene: {geneid}, no mapped reads in SAM, skip")
            not_in_sam += 1
            continue
        if not region_str in region_str_to_variants:
            logging.debug(f"{chrN},gene: {geneid}, no variants/records in VCF, skip")
            in_sam_not_in_vcf += 1
            continue

        # gene_to_variants {gene region}:{records from VCF}
        variants = region_str_to_variants[region_str]
        # gene_to_qualityStr {gene region}:{quality string score from sam.gz}
        sam_data = region_str_to_sam_data[region_str]

        if all([not x.isHet() for x in variants]):
            logging.debug("all SNPs are homozygous")
            continue

        # within each gene, we obtain the reads information from SAM file
        pos1_tlen = [
            x[0] for x in sam_data
        ]  # start pos of reads, quality string length
        qual_strs = [x[1] for x in sam_data]  # quality string
        if len(qual_strs) == 0 or len(pos1_tlen) == 0:
            num_genes_with_no_qual_strs_or_pos_tlens += 1
            continue

        # summarize start pos of reads, quality string length for this gene: {(89635, 124): 1, (89665, 187): 1}
        pos1_tlen_to_count = {}
        for x in pos1_tlen:
            pos1_tlen_to_count[x] = (
                pos1_tlen_to_count[x] + 1 if x in pos1_tlen_to_count else 1
            )
        # calculate the minimum quality length of quality string for this gene, maybe 75
        qual_strs.sort(key=lambda x: len(x))
        qual_strs = np.array(qual_strs)
        qual_str_lens = np.array([len(x) for x in qual_strs])
        min_qual_len = qual_str_lens[0]
        ######## filter1: loop through each transcript, to
        transcript_id_to_fragLen = posTlen_to_fragLen(
            gene, pos1_tlen_to_count, min_qual_len
        )
        if len(transcript_id_to_fragLen) == 0:
            num_genes_with_all_transcripts_filtered += 1
            continue

        recorded_genes += 1

        candidate_sequence_items = makeModifiedSequences(
            gene, variants, list(transcript_id_to_fragLen.keys())
        )

        # read coverage per SNP we want * length of the longest transcript / minimum quality string length (for the gene)
        if use_random_reads:
            read_depth = int(read_depth / 2)
        num_reads = int(float(read_depth / min_qual_len) * length)

        unique_qual_str_lengths = set([len(x) for x in qual_strs])

        logging.debug(
            f">>>>> gene id: {geneid}, gene region: {region_str}, #reads: {num_reads}, #qual_strs: {len(qual_strs)}, #qual_str_lens: {len(unique_qual_str_lengths)},total #variants in VCF: {len(variants)}"
        )
        for i in range(num_reads):
            read_results = simulateRead(
                candidate_sequence_items,
                transcript_id_to_fragLen,
                qual_strs,
                qual_str_lens,
                use_random_reads,
            )
            if read_results is None:
                logging.debug("break!")
                n_break += 1
                continue

            for read in read_results:
                is_mat, fwd_seq, fwd_qual, rev_seq, rev_qual = read

                if is_mat:
                    mat_reads += 1
                else:
                    pat_reads += 1

                id = f"@SIM-{nextReadID}-{geneid}:{'MAT' if is_mat else 'PAT'}"
                nextReadID += 1
                write_fwd(id, fwd_seq, fwd_qual)
                write_rev(id, rev_seq, rev_qual)
    # TODO return stats


def write_read(id, seq, qual, file):
    file.write(f"@{id}\n{seq}\n+\n{qual}\n")


def open_output_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "wt", compresslevel=3)
    else:
        return open(filename, "w")


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("twobit", help="full path to two bit")
    parser.add_argument("genome", help="full path to hg19.2bit")
    parser.add_argument("gff", help="full path to gencode gtf file")
    parser.add_argument("samgz", help="full path to sam.gz")
    parser.add_argument("vcf", help="full path to VCF file without chr")
    parser.add_argument("--read_depth", help="per-base-read-depth", type=int)
    parser.add_argument(
        "--out1",
        help="output path for forward strand fastq reads",
        default="./read_1.fastq.gz",
    )
    parser.add_argument(
        "--out2",
        help="output path for reverse strand fastq reads",
        default="./read_2.fastq.gz",
    )
    parser.add_argument("--chr", help="specific chromosome to simulate")
    parser.add_argument("--gene", help="specific gene to simulate", default=None)
    parser.add_argument(
        "-r",
        "--random",
        action="store_true",
        help="generate reads randomly from either haplotype",
    )
    parser.add_argument(
        "--seed",
        help="random seed for simulation",
        default=random.randrange(sys.maxsize),
        type=int,
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="print lots of info"
    )
    parser.add_argument("--all_snps", action="store_true", help="include all SNPs")
    parser.add_argument(
        "--max_genes", help="max number of genes to simulate", type=int, default=0
    )

    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    if args.random:
        logging.debug("RANDOMLY generate pat/mat reads for each read ...")
    else:
        logging.debug("EQUALLY generate pat/mat reads for each read...")

    if args.chr:
        logging.debug(f"{datetime.now()} CHRomosome mode turned on for {args.chr}")
    else:
        logging.debug(f"{datetime.now()} default: looking at all genes")

    if args.all_snps:
        logging.debug(f"{datetime.now()} all SNPs mode turned on")
    else:
        logging.debug(f"{datetime.now()} default: looking at all hets")

    if args.gene:
        logging.debug(f"{datetime.now()} target gene mode turned on for {args.gene}")

    logging.debug(args)

    logging.info(f"{datetime.now()} simulation seed : {args.seed}")
    random.seed(args.seed)

    target_chromosome = args.chr
    target_gene = args.gene

    outFile1 = args.out1
    outFile2 = args.out2
    logging.info(f"output files {outFile1} {outFile2}")

    # Load GFF and fragment lengths

    logging.info(f"{datetime.now()} reading GFF...")
    genes = loadGenes(args.gff, target_chromosome, target_gene)
    logging.info(f"{datetime.now()} done reading GFF...")

    logging.info(f"found {len(genes)} genes for transcripts")
    if len(genes) == 0:
        logging.info(f"{datetime.now()} no genes found in GFF")
        return

    # annotate transcripts
    twoBitDir = args.twobit
    genome2bit = args.genome
    logging.info(f"running 2bit {twoBitDir} {genome2bit}")
    with tempfile.NamedTemporaryFile(mode="w") as twoBitInputFile:
        constructTwoBitInput(genes, twoBitInputFile.name)
        twobitId_to_seq = run2bitBatch(twoBitDir, twoBitInputFile.name, genome2bit)
    logging.info(f"done running 2bit {len(twobitId_to_seq)} ids")

    global matches, mismatches
    matches = 0
    mismatches = 0
    annotate_transcripts(genes, twobitId_to_seq)
    logging.info(f"done annotating transcripts")

    variant_processor = (
        variant_processor_SNPs if args.all_snps else variant_processor_hets
    )
    region_str_to_variants, region_str_to_sam_data = loadRegionData(
        args.vcf, args.samgz, genes, variant_processor
    )

    # Simulate
    logging.info("Start simulation")

    if tracemalloc.is_tracing():
        snapshot1 = tracemalloc.take_snapshot()

    with open_output_file(outFile1) as OUT1, open_output_file(outFile2) as OUT2:
        runSimulation(
            genes,
            region_str_to_variants,
            region_str_to_sam_data,
            args.random,
            args.read_depth,
            args.max_genes,
            lambda id, seq, qual: write_read(id, seq, qual, OUT1),
            lambda id, seq, qual: write_read(id, seq, qual, OUT2),
        )

    logging.info(f"{datetime.now()} Finish simulation")

    if tracemalloc.is_tracing():
        snapshot2 = tracemalloc.take_snapshot()
        top_stats = snapshot2.compare_to(snapshot1, "lineno")
        print("[ Top 10 differences ]")
        for stat in top_stats[:10]:
            print(stat)

    # print stats
    # print("")
    # print(">> total geneID in gtf: %d" % (num_gene_gtf))
    # print(">> total ReadID processed: %d" % (nextReadID))
    # if matches == 0:
    #     print(f">> # matches: {matches}, # mismatches: {mismatches}")
    # else:
    #     print(
    #         f">> # matches: {matches}, # mismatches: {mismatches}, percentage of matches: {round(matches/(mismatches+matches),2)*100}%"
    #     )
    # print(">> total num breaks /gene is shorter than fragment length : %d" % (n_break))
    # print(
    #     f">> # recorded genes : {recorded_genes}, # processed_genes: {processed_genes}, percentage: {round(recorded_genes/processed_genes,2)*100}%"
    # )
    # print(
    #     f">> # gene not in SAM: {not_in_sam}, # genes in SAM but not in VCF: {in_sam_not_in_vcf}, total filered out genes: {not_in_sam+in_sam_not_in_vcf} ({round((not_in_sam+in_sam_not_in_vcf)/processed_genes,2)*100})%"
    # )


# if target_gene is not None:
#     print(
#         ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ratio of #mat/#total"
#     )
#     print(list_ratio)

# if chromosome is not None:
#     out = out_path
#     with open(out + "/mat_ratio_" + str(chromosome) + ".pkl", "wb") as fp:
#         pickle.dump(list_ratio, fp)

# if target_gene is not None:
#     with open(gene_folder + "/fragLen", "wb") as fp:
#         pickle.dump(list_fragLen, fp)

if __name__ == "__main__":
    main()
