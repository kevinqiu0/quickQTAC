#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics             # Queue
#SBATCH -t 48:00:00             # Walltime/duration of the job                 # Number of Nodes
#SBATCH -N 1 
#SBATCH -n 12       
#SBATCH --export=NONE
#SBATCH -J qtac
#SBATCH --mail-user=kevin.qiu1@northwestern.edu
#SBATCH --mail-type=END


module load python/anaconda3.6
module load tophat/2.1.0
module load java
module load fastqc
module load STAR
module load samtools
module load pigz
PATH=$HOME/.local/bin:$PATH


sampleID=$1
echo "${sampleID}"
if [ -f "${sampleID}_2.fastq.gz" ]; then
	echo "Running paired-end trimming and alignment";
	/projects/b1036/kevinq/corr/TrimGalore-0.6.5/trim_galore --cores $nthread --paired --basename ${sampleID} --fastqc ${sampleID}_1.fastq.gz ${sampleID}_2.fastq.gz
	/projects/b1036/SOFTWARE/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR \
	--runThreadN $nthread \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outReadsUnmapped Fastx \
	--outFilterType BySJout \
	--genomeDir $genomedir \
	--sjdbGTFfile $gtf \
	--readFilesIn ${sampleID}_val_1.fq.gz ${sampleID}_val_2.fq.gz \
	--outFileNamePrefix ${sampleID}_paired_all_ \
	--outSAMstrandField intronMotif \
	--limitBAMsortRAM $mem
	
    
	htseq-count  \
	--format bam \
	--order pos \
	--strand no \
	-q ${sampleID}_paired_all_Aligned.sortedByCoord.out.bam $gtf > ${sampleID}_paired_all.ht.counts
else
	echo "Running single-end trimming and alignment";
	/projects/b1036/kevinq/corr/TrimGalore-0.6.5/trim_galore --cores $nthread --basename ${sampleID} --fastqc ${sampleID}_1.fastq.gz
	/projects/b1036/SOFTWARE/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR \
	--runThreadN $nthread \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outReadsUnmapped Fastx \
	--outFilterType BySJout \
	--genomeDir $genomedir \
	--sjdbGTFfile $gtf\
	--readFilesIn ${sampleID}_trimmed.fq.gz  \
	--outFileNamePrefix ${sampleID}_single_all_ \
	--outSAMstrandField intronMotif \
	--limitBAMsortRAM $mem
	module load python/anaconda3.6

	htseq-count  \
	--format bam \
	--order pos \
	--strand no \
	-q ${sampleID}_single_all_Aligned.sortedByCoord.out.bam $gtf > ${sampleID}_single_all.ht.counts
fi


rm ${sampleID}_val_1.fq.gz ${sampleID}_val_2.fq.gz

