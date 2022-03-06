#!/bin/bash
#SBATCH --job-name=all
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=5G
#SBATCH --time=72:00:00
#SBATCH -c 5
#SBATCH --output=log/write-config-downsampled-%j.out
#SBATCH --error=log/write-config-downsampled-%j.out
#SBATCH -p scavenger
hostname 1>&2

# load modules
module load VCFtools/0.1.17
module load htslib
module load Parallel

util_dir=/hpc/home/bmajoros/twobit
genome=/datacommons/allenlab/hg19/hg19.2bit
gff=/datacommons/allenlab/hg19/filter/gencode.v19.annotation.level12.gtf
original_read_len=75

# specify paths
commands=""
for sample in "HG00099" "HG00108" "HG00251" "NA18498" "NA18502" "NA19213"
do
    now=$(date +"%T")
    echo ">> start -"${sample}":"${now}
    ### make 
    vcf_path=/hpc/group/allenlab/scarlett/output/RNAseq/$sample/VCF
    if [ -f "$vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz.tbi" ];then
            echo "$vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz.tbi exists."
    else
            bgzip -c $vcf_path/${sample}.with_chr.content.SNPs.hets.vcf > $vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz
            tabix -p vcf $vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz
    fi
    
    ###
    aligned_rna=/hpc/group/allenlab/scarlett/output/RNAseq/${sample}/star_2pass_WASP/mismatch10/${sample}.sam.gz
    vcfgz=/hpc/group/allenlab/scarlett/output/RNAseq/${sample}/VCF/${sample}.with_chr.content.SNPs.hets.vcf.gz
    out_path=/hpc/group/allenlab/scarlett/output/RNAseq/${sample}/spliced_reads
    mkdir -p $out_path
    command="python write_config.py ${util_dir} ${genome} ${gff} ${aligned_rna} ${vcfgz} ${original_read_len} ${out_path}"
    command+=" && echo \"DONE $N\""
    commands+="$command"
    commands+=$'\n'
done

echo "$commands"
echo '\n'

echo "$commands" | parallel -j 5
now=$(date +"%T")
echo ">> end -"${now}
