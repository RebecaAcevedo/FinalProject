#!/bin/bash

### Don't change unless you know why
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=30gb

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
### Each task requires one cpu by default; include cpus-per-task if more are required

#SBATCH --time=6:00:00

### Change as you see fit
#SBATCH --account=416a516a
#SBATCH --mail-type=END


### STAR mapping job part
### set input file dir
courseData="/groups/guangyao/416a516a"
# F1="/home/u3/racevedobarboza/HW1.fastp/SRR8353399_1.fastq.gz"
# F2="/home/u3/racevedobarboza/HW1.fastp/SRR8353399_2.fastq.gz"

inputDir=/xdisk/guangyao/416a516a/Group15/Rebeca/Fastp
cd $inputDir
filelist=$(ls *_clean* | awk -F_ '{print $1}' | uniq)

### set output file dir
Fo="/xdisk/guangyao/416a516a/Group15/Rebeca/STAR" ## full path is needed (~/ won't work)
#mkdir ~/$Fo

module load star

for readfile in $filelist
do
date
STAR \
--runThreadN 2 \
--genomeDir /groups/guangyao/416a516a/humanIndex49/ \
--readFilesIn $readfile"_clean.fastq.gz" \
--readFilesCommand zcat \
--outFileNamePrefix $Fo/$readfile"_" \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif 
date
done


