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

### fastp job part
# set a variable by giving the name (e.g., inputfiles)
# refer to (use) the variable by adding a $ (e.g., $inputfiles)
# note there is no space on either side of "=" when setting a variable
# InputPath="/groups/guangyao/416a516a/example_fastq/REF.G0"
# Read1="SRR8353411_1.fastq.gz"
# Read2="SRR8353411_2.fastq.gz"

# set the dir to store output files
# mkdir ~/fastp (remember to create the dir if you haven't done so)
cd /xdisk/guangyao/416a516a/Group15/Rebeca/Fastp

# run your executable program with begin and end date and time output
#set the input and output files paths 
InputPath="/xdisk/guangyao/416a516a/Group15/Rebeca/Rawreads"
OutputPath="/xdisk/guangyao/416a516a/Group15/Rebeca/Fastp"

# Loop over input files
for Input_file in $InputPath/*.fastq.gz
do
    date
    # Set the output file path
    Output_file=${OutputPath}/$(basename ${Input_file} .fastq.gz)_clean.fastq.gz
    
    # Run fastp
    /groups/guangyao/416a516a/modules/fastp -i $Input_file -o $Output_file -5 -r
    date
done

# InputPath="/xdisk/guangyao/416a516a/Group15/Rebeca/Rawreads/"
# for readfile in /xdisk/guangyao/416a516a/Group15/Rebeca/Rawreads/*.fastq.gz
# do
#    date
#    /groups/guangyao/416a516a/modules/fastp -i $readfile -o $readfile -5 -r
#    date
# done
# more about fastp: https://github.com/OpenGene/fastp

# Output_file=${OutputPath}/$(basename ${Input_file} .fastq.gz)_clean.fastq.gz
# SRR349741
# SRR349742
# SRR349743
# SRR349744
# SRR349745
# SRR349746
# SRR349747
# SRR349748
# SRR349749
# SRR349750
# SRR349751
# SRR349752