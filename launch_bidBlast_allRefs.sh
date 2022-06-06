#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=bidBlast
#SBATCH --mem=10G
#SBATCH --time=00-15:00:00
#Define sdout path
#SBATCH --output=/home/guisana/scripts/bidirectBlastNwGenomes/output_bidBlast.txt
#Define sderr path
#SBATCH --error=/home/guisana/scripts/bidirectBlastNwGenomes/error_bidBlast.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=medium


module load python/3.8

# General variables for whole pipeline
#######################################################################################################

refGbk= # ------># Path to reference .GBK file

# Paths to gbk files to be blasted. Add as many lines for files that are desired to be analyzed. Add number to variable and add variable to targetArray
gbk1="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0107/PROKKA_04212022.gbk" # --------># Path to target .GBK file
gbk2="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0109/PROKKA_04212022.gbk" # --------># Path to target .GBK file
gbk3="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0110/PROKKA_04212022.gbk" # --------># Path to target .GBK file
gbk4="/storage/PGO/data/mtb/wrk_dataset/h37Rv_gbk/AL123456.3.gbk" # ----------------------------------------># Path to target .GBK file
gbk5="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0108/PROKKA_04212022.gbk" # --------># Path to target .GBK file

gbkArray=($gbk1 $gbk2 $gbk3 $gbk4 $gbk5)

whatBlast="both" # -----------------------------------------------------------------------------------------># What BLAST to run. For BLASTn write nt,
                                                                                                             # for BLASTp write pt and for both write both
covThrshld=0.95 # ------------------------------------------------------------------------------------------># Coverage threshold for considering a blast 
                                                                                                             # hit among the best hits.
pidThrshld=95 # --------------------------------------------------------------------------------------------># Threshold Percentage of Identity for being
                                                                                                             # considered an orthologue.
output="/storage/PGO/results/mtb/bidBlastNuSeqs_allRefs/" # ------------------------------------------------># Output directory where results will be stored

# Create output directory if it does not exist.
#######################################################################################################
echo 'Starting Genome comparison pipeline...'

if [ ! -d $output ]
then 
  mkdir $output
  echo "${output} created."
else
  echo "${output} already exists, using this location."
fi
 
# Do BLASTs for each GBK file against the each other as reference and parse results in orthologue 
# matrixes for each one of the references.
#######################################################################################################

for r in ${gbkArray[@]}
do
  delete=$r
  refArray=( ${gbkArray[@]/$delete} )
  for t in ${refArray[@]}
  do
    # Run bidirectional BLAST
    /usr/bin/time python3 doBiBlast.py --refGbk $r --targetGbk $t --whatBlast $whatBlast --output $output
    # Run BLASTn of each reference's gene against each target's full genome with the aim of detecting
    # unannotated genes or genes that have stop gain mutations by comparing these results to what is 
    # absent in the bidirectional BLAST.
    /usr/bin/time python3 fullGenomBlast.py --refGbk $r --targetGbk $t --covThrshld $covThrshld --pidThrshld $pidThrshld --output $output
  done
  # Parse the result of the BLASTs into a orthologue matrix for the given reference genome.
  /usr/bin/time python3 parseBBHs.py --refGbk $r --whatBlast $whatBlast --pidThrshld $pidThrshld --geneIDsMat --output $output
  # Obtain the genes of each genome that are conflictive in terms of having multiple matches for each target that satisfy 
  # coverage and PID requirements, and generate a table of the genes and the hits for each genome. 
  /usr/bin/time python3 getConfGenes.py --refGbk $r --whatBlast $whatBlast --pidThrshld $pidThrshld --covThrshld $covThrshld--output $output
done

