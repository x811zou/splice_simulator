#!/bin/zsh

reference_dir=~scarlett/data/reference
out_dir=$HOME/data/simulator_test
out_extension=".gz"

profile=""
profile="$profile -O"
# profile="$profile -m cProfile -o $HOME/data/simulator_test/profile"
# profile="$profile -X tracemalloc"

max_genes=200

out1=$out_dir/read_1.fastq${out_extension}
out2=$out_dir/read_2.fastq${out_extension}


python3 ${=profile} unbiased-spliced-rna-test.py \
  ~/data/two_bit \
  $reference_dir/two_bit/hg19.2bit \
  $reference_dir/gencode.v19.annotation.level12.gtf \
  ~/data/NA12878.sam.gz \
  ~/data/NA12878.no_chr.content.SNPs.filtered.vcf.gz \
  --out1 >(pigz -3 -c > $out1) \
  --out2 >(pigz -3 -c > $out2) \
  --max_genes $max_genes \
  --read_depth 100 \
  --seed 0 \
  --all_snps \
  --chr chr5 \

wait

  




  # --out1 $out1 --out2 $out2 \
  # --out1 >(pigz --fast -c > $out1) --out2 >(pigz --fast -c > $out2) \
  # --out1 >(bgzip -l3 -c -@4 > $out1) --out2 >(bgzip -l3 -c -@4 > $out2) \


