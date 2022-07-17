#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics             # Queue
#SBATCH -t 1:00:00             # Walltime/duration of the job                 # Number of Nodes
#SBATCH -n 1       
#SBATCH --mem=8G               # Memory per node in GB needed for a job. Also see --mem-per-cp                        # Number of Cores (Processors)
#SBATCH --export=NONE
#SBATCH -J qtac




module load python/anaconda3.6 fastqc

if [ $# -eq 0 ]
  then
    echo "No arguments supplied, please specify either hg19 or hg38"
	exit 1
fi
genome=$1


if [ $genome = hg19 ];then
  then
	. ./hg19.config
elif [ $genome = hg38 ];
  then
	. ./hg38.config
else
   echo "Unknown argument given, please specify either hg19 or hg38"
   exit 1
fi
#######################################################################


for i in `ls *_1.fastq.gz`;
do
	sampleID="${i%*_1.fastq.gz}"
	sbatch qtacSingle.sh $sampleID $genome
	echo "Submitted job for: ${sampleID}"


done
echo $?
echo "SBATCH Files have been submitted!"
