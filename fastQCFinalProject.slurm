#!/bin/bash

### Don't change unless you know why
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=01:00:00

### Change as you see fit
#SBATCH --account=416a516a
#SBATCH --mail-type=END

### fastQC job part
# set a variable by giving the name (e.g., inputfiles)
# refer to (use) the variable by adding a $ (e.g., $inputfiles)
# note there is no space on either side of "=" when setting a variable
InputPath="/groups/guangyao/416a516a/example_fastq/REF.G0"
Read1="SRR8353411_1.fastq.gz"
Read2="SRR8353411_2.fastq.gz"

# set the dir to store output files
# mkdir ~/fastQC (remember to create the dir if you haven't done so)
cd /xdisk/guangyao/416a516a/Group15/Rebeca/FastQC

module load fastqc	### Load the HPC-installed module/program
# run your executable program with begin and end date and time output

for readfile in /xdisk/guangyao/416a516a/Group15/Rebeca/Rawreads/*_clean.fastq.gz
do
    date
    fastqc $readfile -o ./
    date
done


# -o ~/fastQC is important: fastqc otherwise generate the output files in the input file dir by default