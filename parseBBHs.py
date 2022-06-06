import numpy as np
import os
from Bio import Entrez, SeqIO
from glob import glob
import re
import pandas as pd
import copy as cp
import argparse

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

#
# Directory stuff
####################################################################################################

blastnDir = outDir + 'bbh_nucl/'
blastpDir = outDir + 'bbh_prot/'

#
# Functions
####################################################################################################

# Given a .GBK file generates a pandas dataframe with locus tags, gene names and product info.
def getGenesInfo(gbkFile, type = 'prot'):
    for record in SeqIO.parse(gbkFile, "genbank"):
        fTypes = []
        locusTags = []
        geneNames = []
        products = []
        for f in record.features:
            fType = f.type
            if 'locus_tag' in f.qualifiers.keys():
                locusTag = f.qualifiers['locus_tag'][0]
            else:
                locusTag = np.nan
            if 'gene' in f.qualifiers.keys():
                geneName = f.qualifiers['gene'][0]
            else:
                geneName = np.nan
            if 'product' in f.qualifiers.keys():
                product = f.qualifiers['product'][0]
            else:
                product = np.nan
            fTypes.append(fType)
            geneNames.append(geneName)
            locusTags.append(locusTag)
            products.append(product)
    outDF = pd.DataFrame(list(zip(fTypes, locusTags, geneNames, products)), columns = ['type', 'locus_tag', 'name', 'product'])
    # Remove the 2 NaNs that are in the DF: one is because of a feature describing the genome, and another because
    # a CRISPR region
    outDF = outDF[[isinstance(x, str) for x in outDF['locus_tag'].to_list()]]
    if type == 'prot':
        outDF = outDF[outDF['type'] == 'CDS']
    outDF = outDF.drop('type', axis = 1)
    outDF = outDF.reset_index(drop=True)
    return outDF

# Obtains ID from within a genbank file
def getGbkID(gbkFile):
    outID = ''
    for record in SeqIO.parse(gbkFile, 'genbank'):
        #print(record.features)
        outID = cp.copy(record.id)
    return outID

# Parses BLAST best hits files into orthologue matrixes of percentage of identity (PID) and/or 
# locus tag equivalence
def parseBBHFiles(type = 'prot', refGbk = refGbk, getGeneIDsMat = False):
    refID = getGbkID(refGbk)
    if type == 'prot':
        blastDir = blastpDir
        blastFiles = [x for x in os.listdir(blastDir) if '_parsed.csv' in x]
    elif type == 'nucl':
        blastDir = blastnDir
        blastFiles = [x for x in os.listdir(blastDir) if '_parsed.csv' in x]
    blastFiles = [x for x in blastFiles if '%s_vs_'%refID in x]
    targetIDs = [x.split('_vs_')[1] for x in blastFiles]
    targetIDs = [re.sub('_parsed.csv', '', x) for x in targetIDs]
    orthoMat = getGenesInfo(refGbk, type = type)
    if getGeneIDsMat == True:
        geneIDs_matrix = cp.copy(orthoMat)
    for id in targetIDs:
        orthoMat[id] = list(np.repeat(np.nan, orthoMat.shape[0]))
        if getGeneIDsMat == True:
            geneIDs_matrix[id] = list(np.repeat(np.nan, geneIDs_matrix.shape[0]))
    for blast in blastFiles:
        blastPath = blastDir + blast
        bbh=pd.read_csv(blastPath)
        listIDs=[]
        listPID=[]
        for r in orthoMat['locus_tag']:
            try:
                currentOrtholog=bbh[bbh['gene']==r].reset_index()
                listIDs.append(currentOrtholog.iloc[0]['subject'])
                listPID.append(currentOrtholog.iloc[0]['PID'])
            except:
                listIDs.append('None')
                listPID.append(0)
        for col in orthoMat.columns:
            if col in blast:
                orthoMat[col]=listPID
                if(getGeneIDsMat == True):
                    geneIDs_matrix[col]=listIDs
    if getGeneIDsMat == True:
        out = [orthoMat, geneIDs_matrix]
        return out
    else:
        return orthoMat


# Binarizes PID matrix
def binPIDs(PIDs, thrshld = 90):
    binPIDs = cp.copy(PIDs)
    IDs = [x for x in list(PIDs.columns) if x not in ['locus_tag', 'name', 'product']]
    for ID in IDs:
        binCol = [int(x) for x in PIDs[ID] >= thrshld]
        binPIDs[ID] = binCol
    return binPIDs

#
# Parsing
####################################################################################################

# Obtain reference genome ID
refID = getGbkID(refGbk)

# Obtain matrixes of orthologous genes.
if whatBlast == 'pt':
    if geneIDsMat == True:
        parsed = parseBBHFiles(type = 'prot', refGbk = refGbk, geneIDsMat = True)
        parsedPIDs = parsed[0]
        parsedIDs = parsed[1]
        parsedIDs.to_csv('%s%s_genIDs_mat'%(blastpDir, refID))
        print('%s_geneIDs_mat.csv saved at %s.'%(refID, blastpDir))
    else:
        parsedPIDs = parseBBHFiles(type = 'prot', refGbk = refGbk, getGeneIDsMat = False)
    binnedPIDs = binPIDs(parsedPIDs, thrshld = pidThrshld)
    parsedPIDs.to_csv('%s%s_ortho_mat.csv'%(blastpDir, refID))
    print('%s_ortho_mat.csv saved at %s.'%(refID, blastpDir))
    binnedPIDs.to_csv('%s%s_ortho_matBin.csv'%(blastpDir, refID))
    print('%s_ortho_matBin.csv saved at %s.'%(refID, blastnDir))
elif whatBlast == 'nt':
    if geneIDsMat == True:
        parsed = parseBBHFiles(type = 'nucl', refGbk = refGbk, getGeneIDsMat = True)
        parsedPIDs = parsed[0]
        parsedIDs = parsed[1]
        parsedIDs.to_csv('%s%s_genIDs_mat.csv'%(blastnDir, refID))
        print('%s_geneIDs_mat.csv saved at %s.'%(refID, blastnDir))
    else:
        parsedPIDs = parseBBHFiles(type = 'nucl', refGbk = refGbk)
    binnedPIDs = binPIDs(parsedPIDs, thrshld = pidThrshld)
    parsedPIDs.to_csv('%s%s_ortho_mat.csv'%(blastnDir, refID))
    print('%s_ortho_mat.csv saved at %s.'%(refID, blastnDir))
    binnedPIDs.to_csv('%s%s_ortho_matBin.csv'%(blastnDir, refID))
    print('%s_ortho_matBin.csv saved at %s.'%(refID, blastnDir))
elif whatBlast == 'both':
    if geneIDsMat == True:
        # Amino acid sequences...
        parsedPt = parseBBHFiles(type = 'prot', refGbk = refGbk, getGeneIDsMat = True)
        parsedPIDsPt = parsedPt[0]
        parsedIDsPt = parsedPt[1]
        parsedIDsPt.to_csv('%s%s_genIDs_mat.csv'%(blastpDir, refID))
        print('%s_geneIDs_mat.csv saved at %s.'%(refID, blastpDir))
        # Nucleotide sequences...
        parsedNt = parseBBHFiles(type = 'nucl', refGbk = refGbk, getGeneIDsMat = True)
        parsedPIDsNt = parsedNt[0]
        parsedIDsNt = parsedNt[1]
        parsedIDsNt.to_csv('%s%s_genIDs_mat.csv'%(blastnDir, refID))
        print('%s_geneIDs_mat.csv saved at %s.'%(refID, blastnDir))
    else:
        # Amino acid sequences...
        parsedPIDsPt = parseBBHFiles(type = 'prot', refGbk = refGbk)
        # Nucleotide sequences...
        parsedPIDsNt = parseBBHFiles(type = 'nucl', refGbk = refGbk)
    # Amino acid sequences...
    binnedPIDsPt = binPIDs(parsedPIDsPt, thrshld = pidThrshld)
    parsedPIDsPt.to_csv('%s%s_ortho_mat.csv'%(blastpDir, refID))
    print('%s_ortho_mat.csv saved at %s.'%(refID, blastpDir))
    binnedPIDsPt.to_csv('%s%s_ortho_matBin.csv'%(blastpDir, refID))
    print('%s_ortho_matBin.csv saved at %s.'%(refID, blastpDir))
    # Nucleotide sequences...
    binnedPIDsNt = binPIDs(parsedPIDsNt, thrshld = pidThrshld)
    parsedPIDsNt.to_csv('%s%s_ortho_mat.csv'%(blastnDir, refID))
    print('%s_ortho_mat.csv saved at %s.'%(refID, blastnDir))
    binnedPIDsNt.to_csv('%s%s_ortho_matBin.csv'%(blastnDir, refID))
    print('%s_ortho_matBin.csv saved at %s.'%(refID, blastnDir))
    
# Create versions of the matrixes accounting only for the genes that are missing in any of the target genomes, both for 
# Pt blasts and Nt blasts
# Create versions of the matrixes accounting only for the genes that are missing in any of the target genomes, both for 
# Pt blasts and Nt blasts
if whatBlast == 'both':
    # Protein BLASTs...
    trueLstPt = [sum(binnedPIDsPt.iloc[:, 3:binnedPIDsPt.shape[1]].loc[i]) != 4 for i in range(0, binnedPIDsPt.shape[0])]
    # Nucleotide BLASTs...
    trueLstNt = [sum(binnedPIDsNt.iloc[:, 3:binnedPIDsNt.shape[1]].loc[i]) != 4 for i in range(0, binnedPIDsNt.shape[0])]
    if geneIDsMat == True:
        # Pt
        genIDMatMissPt = parsedIDsPt[trueLstPt].reset_index(drop = True)
        genIDMatMissPt.to_csv('%s%s_geneIDs_mat_miss.csv'%(blastpDir, refID))
        print('%s_geneIDs_mat_miss.csv saved at %s.'%(refID, blastpDir))
        # Nt
        genIDMatMissNt = parsedIDsNt[trueLstNt].reset_index(drop = True)
        genIDMatMissNt.to_csv('%s%s_geneIDs_mat_miss.csv'%(blastnDir, refID))
        print('%s_geneIDs_mat_miss.csv saved at %s.'%(refID, blastnDir))
    # Pt
    orthoMatMissPt = binnedPIDsPt[trueLstPt].reset_index(drop = True)
    orthoMatMissPt.to_csv('%s%s_ortho_matBin_miss.csv'%(blastpDir, refID))
    print('%s_ortho_matBin_miss.csv saved at %s.'%(refID, blastpDir))
    # Nt
    orthoMatMissNt = binnedPIDsNt[[sum(binnedPIDsNt.iloc[:, 3:binnedPIDsNt.shape[1]].loc[i]) != 4 for i in range(0, binnedPIDsNt.shape[0])]].reset_index(drop = True)
    orthoMatMissNt.to_csv('%s%s_ortho_matBin_miss.csv'%(blastnDir, refID))
    print('%s_ortho_matBin_miss.csv saved at %s.'%(refID, blastnDir))
elif whatBlast == 'pt':
    # Protein BLASTs...
    trueLstPt = [sum(binnedPIDsPt.iloc[:, 3:binnedPIDsPt.shape[1]].loc[i]) != 4 for i in range(0, binnedPIDsPt.shape[0])]
    if geneIDsMat == True:
        # Pt
        genIDMatMissPt = parsedIDsPt[trueLstPt].reset_index(drop = True)
        genIDMatMissPt.to_csv('%s%s_geneIDs_mat_miss.csv'%(blastpDir, refID))
        print('%s_geneIDs_mat_miss.csv saved at %s.'%(refID, blastpDir))
    # Pt
    orthoMatMissPt = binnedPIDsPt[trueLstPt].reset_index(drop = True)
    orthoMatMissPt.to_csv('%s%s_ortho_matBin_miss.csv'%(blastpDir, refID))
    print('%s_ortho_matBin_miss.csv saved at %s.'%(refID, blastpDir))
elif whatBlast == 'nt':
    # Nucleotide BLASTs...
    trueLstNt = [sum(binnedPIDsNt.iloc[:, 3:binnedPIDsNt.shape[1]].loc[i]) != 4 for i in range(0, binnedPIDsNt.shape[0])]
    if geneIDsMat == True:
        # Nt
        genIDMatMissNt = parsedIDsNt[trueLstNt].reset_index(drop = True)
        genIDMatMissNt.to_csv('%s%s_geneIDs_mat_miss.csv'%(blastnDir, refID))
        print('%s_geneIDs_mat_miss.csv saved at %s.'%(refID, blastnDir))
    # Nt
    orthoMatMissNt = binnedPIDsNt[[sum(binnedPIDsNt.iloc[:, 3:binnedPIDsNt.shape[1]].loc[i]) != 4 for i in range(0, binnedPIDsNt.shape[0])]].reset_index(drop = True)
    orthoMatMissNt.to_csv('%s%s_ortho_matBin_miss.csv'%(blastnDir, refID))
    print('%s_ortho_matBin_miss.csv saved at %s.'%(refID, blastnDir))