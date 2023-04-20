#!/bin/bash

### Don't change unless you know why
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
### Each task requires one cpu by default; include cpus-per-task if more are required
#SBATCH --mem=30gb
#SBATCH --time=2:00:00

### Change as you see fit
#SBATCH --account=416a516a
#SBATCH --mail-type=END

### Featurecount job part
courseData="/groups/guangyao/416a516a"
MP="$courseData/modules/subread-2.0.3-Linux-x86_64/bin"
annotation="$courseData/hs.genome.annotation/gencode.v26.primary_assembly.annotation.gtf"
IP=$(ls /xdisk/guangyao/416a516a/Group15/Rebeca/STAR/*.bam)

### set output file dir
OP="/xdisk/guangyao/416a516a/Group15/Rebeca/FeatureCounts" ## replace OP with your own output dir
#mkdir $OP

date
/usr/bin/time $MP/featureCounts \
-T 2 \
-a $annotation \
-o $OP/counts.txt $IP
$IP
date


