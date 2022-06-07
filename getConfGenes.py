import os
import pandas as pd
from Bio import Entrez, SeqIO
import copy as cp
import numpy as np

#
# Parser stuff
####################################################################################################

parser = argparse.ArgumentParser(description='Determine the genes that are conflictive in terms of matching multiple genes of the target genomes.')

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
                    help='Percentage of Identity threshold for considering a blast hit among the best hits.', 
                    default=0.80)
                    
parser.add_argument('--covThrshld', 
	                  type=float,
                    help='Coverage threshold for considering a blast hit among the best hits.', 
                    default=0.80)
                    
parser.add_argument('--output', 
	                  type=str,
                    help='Path to output directory, where BLAST results will be saved.', 
                    default='/storage/PGO/results/mtb/bidBlastNuSeqs/')

args = parser.parse_args()
outDir = args.output
whatBlast = args.whatBlast
refGbk = args.refGbk
pidThrshld = args.pidThrshld
covThrshld = args.covThrshld

#
# Directory stuff
####################################################################################################

if whatBlast == 'pt':
    bbhDir = outDir + 'bbh_prot/'
elif whatBlast == 'nt':
    bbhDir = outDir + 'bbh_nucl/'

outDirPt = outDir + 'prot/'
outDirNt = outDir + 'nucl/'

#
# Functions
####################################################################################################

# Obtains ID from within a genbank file
def getGbkID(gbkFile):
    outID = ''
    for record in SeqIO.parse(gbkFile, 'genbank'):
        #print(record.features)
        outID = cp.copy(record.id)
    return outID
    
# Gets sequence length 
def get_gene_lens(queryID, type = 'prot'):
    if type == 'prot':
        in_folder = outDirPt
        file = os.path.normpath('%s%s.fa'%(in_folder, queryID))
    elif type == 'nucl':
        in_folder = outDirNt
        file = os.path.normpath('%s%s.fna'%(in_folder, queryID))
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out

#
# Determine genes that are conflictive in terms that have multiple blast hits that satisfy 
# requirements
####################################################################################################

# Obtain reference genome ID.
refID = getGbkID(refGbk)

# Load orthologue matrix.
orthoMat = pd.read_csv('%s%s_ortho_mat.csv'%(bbhDir, refID), index_col = 0)

targIDs = [x for x in orthoMat.columns if x not in ['locus_tag', 'name', 'product']]

# Determine what are the conflictive genes in the reference (i.e have more than one match to 
# the target that satisfy PID and COV requisites) and obtains a matrix of all the hits against 
# each target.
cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
if whatBlast == 'pt':
    refLens = get_gene_lens(refID, type='prot')
elif whatBlast == 'nt':
    refLens = get_gene_lens(refID, type='nucl')
confGenes = []
for ID in targIDs:
    bbh = pd.read_csv('%s%s_vs_%s_parsed.csv'%(bbhDir, refID, ID), index_col = 0)
    bbhConfGenes = bbh[bbh['nOrths_%s_vs_%s'%(refID, ID)] > 1]['gene']
    # Append conflictive genes to general list outside of the loop 
    [confGenes.append(g) for g in bbhConfGenes]
    
# Build the conflictive genes dataframe
confGenes = np.unique(confGenes)
confGenesDF = cp.copy(orthoMat.iloc[:, 0:3])
confGenesDF = confGenesDF[[x in confGenes for x in confGenesDF['locus_tag']]].reset_index(drop = True)

# Fill the dataframe with the blasts results
for ID in targIDs:
    idHits = []
    hitsDict = {}
    # Load blast result for the target ID
    blastRes = pd.read_csv('%s%s_vs_%s.txt'%(bbhDir, refID, ID), sep = '\t', names = cols)
    if whatBlast == 'pt':
        tarLens = get_gene_lens(ID, type='prot')
    elif whatBlast == 'nt':
        tarLens = get_gene_lens(ID, type='nucl')
    blastRes = pd.merge(blastRes, refLens)
    blastRes['COV'] = blastRes['alnLength']/blastRes['gene_length']
    for g in confGenesDF['locus_tag'].tolist():
        geneHits = blastRes.loc[blastRes['gene'] == g]
        geneHits = geneHits[[x >= pidThrshld and y >= covThrshld for x, y in zip(geneHits.PID.tolist(), geneHits.COV.tolist())]]
        if geneHits.shape[0] > 0:
            hits = ', '.join(geneHits['subject'].tolist())
        else:
            hits = np.nan
        idHits.append(hits)
    confGenesDF[ID] = idHits

confGenesDF.to_csv('%s%s_confGenes.csv'%(bbhDir, refID))
print('%s_confGenes.csv saved at %s.'%(refID, bbhDir))