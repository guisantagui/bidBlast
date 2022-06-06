import pandas as pd
import os
from glob import glob
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import copy as cp
import numpy as np

#
# Parser stuff
####################################################################################################

parser = argparse.ArgumentParser(description='Parses BLAST best hits files into orthologue matrixes.')

parser.add_argument('--refGbk', 
	                  type=str,
                    help='Path to .GBK file of annotated genome that is gonna be used as reference', 
                    default='/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0108/PROKKA_04212022.gbk')
                    
parser.add_argument('--whatBlast', 
	                  type=str,
                    help='The type of BLAST to be performed. Can be nucleotide (nt), protein (pt) or both.', 
                    default='both')
                    
parser.add_argument('--pidThrshld', 
	                  type=float,
                    help='Threshold percentage of identity for considering a best hit an orthologue.', 
                    default=90)
                    
parser.add_argument('--geneIDsMat', 
	                  action='store_true',
                    help='Exports matrix of locus tag equivalence between strains.')
                    
parser.add_argument('--output', 
	                  type=str,
                    help='Path to output directory, where BLAST results will be saved.', 
                    default='/storage/PGO/results/mtb/bidBlastNuSeqs/')

args = parser.parse_args()

refGbk = args.refGbk
whatBlast = args.whatBlast
pidThrshld = args.pidThrshld
geneIDsMat = args.geneIDsMat
outDir = args.output