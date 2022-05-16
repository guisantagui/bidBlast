#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=bidBlast
#SBATCH --mem=10G
#SBATCH --time=00-00:30:00
#Define sdout path
#SBATCH --output=/home/guisana/scripts/bidirectBlastNwGenomes/output_bidBlast.txt
#Define sderr path
#SBATCH --error=/home/guisana/scripts/bidirectBlastNwGenomes/error_bidBlast.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=short


module load python/3.8

# General variables for whole pipeline
#######################################################################################################

refGbk="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0108/PROKKA_04212022.gbk" # ------># Path to reference .GBK file

# Paths to target files. Add as many lines for files that are desired to be analyzed. Add number to variable and add variable to targetArray
targetGbk1="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0107/PROKKA_04212022.gbk" # --># Path to target .GBK file
targetGbk2="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0109/PROKKA_04212022.gbk" # --># Path to target .GBK file
targetGbk3="/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0110/PROKKA_04212022.gbk" # --># Path to target .GBK file

targetArray=($targetGbk1 $targetGbk2 $targetGbk3)

whatBlast="both" # -----------------------------------------------------------------------------------------># What BLAST to run. For BLASTn write nt,
                                                                                                             # for BLASTp write pt and for both write both
covThrshld=0.80 # ------------------------------------------------------------------------------------------># Coverage threshold for considering a blast 
                                                                                                             # hit among the best hits.
pidThrshld=90 # --------------------------------------------------------------------------------------------># Threshold Percentage of Identity for being
                                                                                                             # considered an orthologue.
output="/storage/PGO/results/mtb/bidBlastNuSeqs/" # --------------------------------------------------------># Output directory where results will be stored


# Do BLASTs for each GBK file against the reference.
#######################################################################################################
for t in ${targetArray[@]}; do
  /usr/bin/time python3 doBiBlast.py --refGbk $refGbk --targetGbk $t --whatBlast $whatBlast
done

# Parse BLAST results into orthologue matrixes.
#######################################################################################################
/usr/bin/time python3 parseBBHs.py --refGbk $refGbk --whatBlast $whatBlast --pidThrshld $pidThrshld --geneIDsMat --output $output