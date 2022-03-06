#!/bin/bash
#SBATCH --job-name=all
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH -c 25
#SBATCH --output=log/AFR-%j.out
#SBATCH --error=log/AFR-%j.out
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
depth=100

# specify paths
echo ">>>>>>>>>>>>>>>>>> step1"
now=$(date +"%T")
echo ">> start -"${now}
commands=""
for sample in "NA18498" "NA18502" "NA19213" #"HG00099" 
do
    now=$(date +"%T")
    echo "sample: "${sample}
    ### make directory
    vcf_path=/hpc/group/allenlab/scarlett/output/RNAseq/$sample/VCF
    if [ -f "$vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz.tbi" ];then
            echo "$vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz.tbi exists."
    else
            bgzip -c $vcf_path/${sample}.with_chr.content.SNPs.hets.vcf > $vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz
            tabix -p vcf $vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz
    fi

    aligned_rna=/hpc/group/allenlab/scarlett/output/RNAseq/${sample}/star_2pass_WASP/mismatch10/${sample}.sam.gz
    vcfgz=$vcf_path/${sample}.with_chr.content.SNPs.hets.vcf.gz
    out_path=/hpc/group/allenlab/scarlett/output/RNAseq/${sample}/spliced_reads
    mkdir -p $out_path
    for N in {1..22}
    do
        command="python unbiased-spliced-rna.py $util_dir $genome $gff $aligned_rna $vcfgz $out_path $depth --out-prefix chr${N}. --chr chr${N} -r"
        command+=" && echo \"DONE $N\""
        commands+="$command"
        commands+=$'\n'
    done
done
echo "$commands"
echo '\n'

echo "$commands" | parallel -j+0 --halt now,fail=1
now=$(date +"%T")
echo ">> end -"${now}



#echo ">>>>>>>>>>>>>>>>>> step2"
#for sample in "HG00099" #"HG00108" "HG00251" "NA18498" "NA18502" "NA19213"
#do
#    now=$(date +"%T")
#    echo ">> start -"${sample}":"${now}
    ### make
    #gzip -d file.gz

#    fq_path=/hpc/group/allenlab/scarlett/output/RNAseq/$sample/spliced_reads/simulated_fastq
#    if [ -f "$fq_path/${sample}.fwd.fastq" ];then
#            echo "$fq_path/${sample}.fwd.fastq exists."
#    else
#        zcat chr*fwd.fastq.gz | gzip > ${sample}.fwd.fastq.gz
#	zcat chr*rev.fastq.gz | gzip > ${sample}.rev.fastq.gz
        #cat chr1.fwd.fastq chr2.fwd.fastq chr3.fwd.fastq chr4.fwd.fastq chr5.fwd.fastq chr6.fwd.fastq chr7.fwd.fastq chr8.fwd.fastq chr9.fwd.fastq chr10.fwd.fastq chr11.fwd.fastq chr12.fwd.fastq chr13.fwd.fastq chr14.fwd.fastq chr15.fwd.fastq chr16.fwd.fastq chr17.fwd.fastq chr18.fwd.fastq chr19.fwd.fastq chr20.fwd.fastq chr21.fwd.fastq chr22.fwd.fastq > ${sample}.fwd.fastq
        #rm chr*.fwd.fastq
        #gzip ${sample}.fwd.fastq
        #cat chr1.rev.fastq chr2.rev.fastq chr3.rev.fastq chr4.rev.fastq chr5.rev.fastq chr6.rev.fastq chr7.rev.fastq chr8.rev.fastq chr9.rev.fastq chr10.rev.fastq chr11.rev.fastq chr12.rev.fastq chr13.rev.fastq chr14.rev.fastq chr15.rev.fastq chr16.rev.fastq chr17.rev.fastq chr18.rev.fastq chr19.rev.fastq chr20.rev.fastq chr21.rev.fastq chr22.rev.fastq > ${sample}.rev.fastq
        #rm chr*.rev.fastq
        #gzip ${sample}.rev.fastq
#    fi
#done
