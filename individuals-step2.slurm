#!/bin/bash
#SBATCH --job-name=all
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=10G
#SBATCH --time=72:00:00
#SBATCH -c 4
#SBATCH --output=log/6individuals-step2-%j.out
#SBATCH --error=log/6individuals-step2-%j.out
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

echo ">>>>>>>>>>>>>>>>>> step2"
for sample in "HG00099"  "HG00108" "HG00251" "NA18498" "NA18502" "NA19213"
do
    now=$(date +"%T")
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start -"${sample}":"${now}
    ### make
    #gzip -d file.gz

    fq_path=/hpc/group/allenlab/scarlett/output/RNAseq/$sample/spliced_reads/simulated_fastq
    if [ -f "$fq_path/${sample}.fwd.fastq" ];then
            echo "$fq_path/${sample}.fwd.fastq exists."
    else
        zcat $fq_path/chr*read1.fastq.gz | gzip > $fq_path/${sample}.fwd.fastq.gz
        zcat $fq_path/chr*read2.fastq.gz | gzip > $fq_path/${sample}.rev.fastq.gz
        rm $fq_path/chr*.fastq.gz
        #cat chr1.fwd.fastq chr2.fwd.fastq chr3.fwd.fastq chr4.fwd.fastq chr5.fwd.fastq chr6.fwd.fastq chr7.fwd.fastq chr8.fwd.fastq chr9.fwd.fastq chr10.fwd.fastq chr11.fwd.fastq chr12.fwd.fastq chr13.fwd.fastq chr14.fwd.fastq chr15.fwd.fastq chr16.fwd.fastq chr17.fwd.fastq chr18.fwd.fastq chr19.fwd.fastq chr20.fwd.fastq chr21.fwd.fastq chr22.fwd.fastq > ${sample}.fwd.fastq
        #rm chr*.fwd.fastq
        #gzip ${sample}.fwd.fastq
        #cat chr1.rev.fastq chr2.rev.fastq chr3.rev.fastq chr4.rev.fastq chr5.rev.fastq chr6.rev.fastq chr7.rev.fastq chr8.rev.fastq chr9.rev.fastq chr10.rev.fastq chr11.rev.fastq chr12.rev.fastq chr13.rev.fastq chr14.rev.fastq chr15.rev.fastq chr16.rev.fastq chr17.rev.fastq chr18.rev.fastq chr19.rev.fastq chr20.rev.fastq chr21.rev.fastq chr22.rev.fastq > ${sample}.rev.fastq
        #rm chr*.rev.fastq
        #gzip ${sample}.rev.fastq
    fi
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end -"${sample}":"${now}
done
