#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics             # Queue
#SBATCH -t 1:00:00             # Walltime/duration of the job                 # Number of Nodes
#SBATCH -n 1       
#SBATCH --mem=8G               # Memory per node in GB needed for a job. Also see --mem-per-cp                        # Number of Cores (Processors)
#SBATCH --export=NONE
#SBATCH -J qtac




module load python/anaconda3.6


for i in `ls *_1.fastq.gz`;
do
	sampleID="${i%*_1.fastq.gz}"
	echo "${sampleID}"

	sbatch qtacSingle.sh ${sampleID}

done
echo $?
echo "SBATCH Files have been submitted!"
