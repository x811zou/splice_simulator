#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett
# =========================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("twobit", help="full path to two bit")
parser.add_argument("genome", help="full path to hg19.2bit")
parser.add_argument("gff", help="full path to gencode gtf file")
parser.add_argument("samgz", help="full path to sam.gz")
parser.add_argument("vcf", help="full path to VCF file with chr")
parser.add_argument("readLen", help="original read length from reads")
parser.add_argument("out_path", help="output path")

args = parser.parse_args()
twobit = args.twobit
genome = args.genome
gtf = args.gff
samgz = args.samgz
vcf = args.vcf
readLen = args.readLen
out = args.out_path

print(args)
FH = open(out+"/simulator.config", "w")
print(f"util-dir = {twobit}",file=FH)
print(f"genome = {genome}",file=FH)
print(f"gff  = {gtf}",file=FH)
print(f"aligned-rna  = {samgz}",file=FH)
print(f"vcf  = {vcf}",file=FH)
print(f"original-read-len  = {readLen}",file=FH)
print(f"out_path = {out}",file=FH)
