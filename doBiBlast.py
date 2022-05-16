import argparse
import os
from Bio import Entrez, SeqIO
from glob import glob
import copy as cp
import pandas as pd

#
# Parser stuff
####################################################################################################

parser = argparse.ArgumentParser(description='Do bidirectional BLAST between two fasta genome files.')

parser.add_argument('--refGbk', 
	                  type=str,
                    help='Path to .GBK file of annotated genome that is gonna be used as reference', 
                    default='/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0108/PROKKA_04212022.gbk')
                    
parser.add_argument('--targetGbk', 
	                  type=str,
                    help='Path to .GBK file of annotated target genome', 
                    default='/storage/PGO/results/mtb/pacbioFISABIO/anotacion_flye/PROKKA_CCS-D0107/PROKKA_04212022.gbk')
                    
parser.add_argument('--whatBlast', 
	                  type=str,
                    help='The type of BLAST to be performed. Can be nucleotide (nt), protein (pt) or both.', 
                    default='both')
                    
parser.add_argument('--covThrshld', 
	                  type=float,
                    help='Coverage threshold for considering a blast hit among the best hits.', 
                    default=0.80)
                    
parser.add_argument('--output', 
	                  type=str,
                    help='Path to output directory, where BLAST results will be saved.', 
                    default='/storage/PGO/results/mtb/bidBlastNuSeqs/')

args = parser.parse_args()
refGbk = args.refGbk
targetGbk = args.targetGbk
whatBlast = args.whatBlast
covThrshld = args.covThrshld
outDir = args.output

#
# Directory stuff
####################################################################################################

# Create output protein and/or nucleotide output directories, depending on what 
# blasts are wanted to be performed. 
outDirPt = outDir + 'prot/'
outDirNt = outDir + 'nucl/'
bbhDirPt = outDir + 'bbh_prot/'
bbhDirNt = outDir + 'bbh_nucl/'

if whatBlast == 'nt' and not os.path.exists(outDirNt):
    os.mkdir(outDirNt)
    os.mkdir(bbhDirNt)
elif whatBlast == 'pt' and not os.path.exists(outDirPt):
    os.mkdir(outDirPt)
    os.mkdir(bbhDirPt)
elif whatBlast == 'both':
    if not os.path.exists(outDirNt):
        os.mkdir(outDirNt)
    if not os.path.exists(bbhDirNt):
        os.mkdir(bbhDirNt)
    if not os.path.exists(outDirPt):
        os.mkdir(outDirPt)
    if not os.path.exists(bbhDirPt):
        os.mkdir(bbhDirPt)

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

# Parses genbank genome file into protein or nucleotide fastas. If protein selected, only CDSs 
# will be included. If not all ORFs will be included (CDS, rRNA, tRNA, tmRNA, repeat_region...).
def parse_genome(file, type='prot', overwrite=1):
    file = glob(file)[0]
    outID = getGbkID(file)
    if type == 'prot':
        out_folder = outDirPt
        out_file= '%s%s.fa'%(out_folder, outID)
    elif type == 'nucl':
        out_folder = outDirNt
        out_file= '%s%s.fna'%(out_folder, outID)
    
    files =os.listdir(out_folder)
    
    if os.path.basename(out_file) in files and overwrite==0:
        print (out_file, 'already parsed')
        return
    else:
        print ('parsing %s'%outID)
    
    handle = open(file)
    
    fout = open(out_file,'w')
    x = 0
    
    records = SeqIO.parse(handle, "genbank")
    for record in records:
        for f in record.features:
            if type == 'prot':
                if f.type=='CDS':
                    seq=f.extract(record.seq)
                    seq=str(seq.translate())
                    
                    if 'locus_tag' in f.qualifiers.keys():
                        locus = f.qualifiers['locus_tag'][0]
                    elif 'gene' in f.qualifiers.keys():
                        locus = f.qualifiers['gene'][0]
                    else:
                        locus = 'gene_%i'%x
                        x+=1
                    fout.write('>%s\n%s\n'%(locus, seq))
            elif type == 'nucl':
                if f.type!='source':
                    seq=f.extract(record.seq)
                    seq=str(seq)
                    if 'locus_tag' in f.qualifiers.keys():
                        locus = f.qualifiers['locus_tag'][0]
                    elif 'gene' in f.qualifiers.keys():
                        locus = f.qualifiers['gene'][0]
                    else:
                        locus = 'gene_%i'%x
                        x+=1
                    fout.write('>%s\n%s\n'%(locus, seq))
    fout.close()
    
# Make database of a given genome.    
def make_blast_db(id, db_type='prot'):
    import subprocess
    if db_type =='prot':
        out_file =os.path.normpath('%s%s.fa.pin'%(outDirPt, id))
        folder = os.path.dirname(out_file)
        files = glob('%s/*.fa.pin'%folder)
    elif db_type == 'nucl':
        out_file =os.path.normpath('%s%s.fna.nin'%(outDirNt, id))
        folder = os.path.dirname(out_file)
        files = glob('%s/*.fna.nin'%folder)
    if out_file in files:
        print (id, 'already has a blast db')
        return
    if db_type=='nucl':
        ext='fna'
    else:
        ext='fa'
    
    cmd_line='makeblastdb -in %s/%s.%s -dbtype %s' %(folder, id, ext, db_type)
    
    print ('making blast db with following command line...')
    print (cmd_line)
    subprocess.call(cmd_line, shell = True)

# Runs BLAST
def run_blast(seqID, dbID, type = 'prot', outfmt=6, evalue=0.001, threads=1):
    import subprocess
    if type == 'prot':
        out_folder = bbhDirPt
        in_folder = outDirPt
        db = os.path.normpath('%s/%s.fa'%(in_folder, dbID))
        seq = os.path.normpath('%s/%s.fa'%(in_folder, seqID))
        blastType = 'blastp'
    elif type == 'nucl':
        out_folder = bbhDirNt
        in_folder = outDirNt
        db = os.path.normpath('%s/%s.fna'%(in_folder, dbID))
        seq = os.path.normpath('%s/%s.fna'%(in_folder, seqID))
        blastType = 'blastn'
    out=os.path.normpath('%s/%s_vs_%s.txt'%(out_folder, seqID, dbID))
    cmd_line='%s -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(blastType, db, seq, out, evalue, outfmt, threads)
    files =os.listdir(out_folder)
    if os.path.basename(out) in files:
        print ('%s already blasted against %s'%(seqID, dbID))
        return
    
    print ('blasting %s vs %s'%(seq, db))
    
    print ('running blast with following command line...')
    print (cmd_line)
    subprocess.call(cmd_line, shell = True)
    return out

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

# Gets Bi-Directional BLASTp Best Hits
def get_bbh(queryID, subjectID, type = 'prot', covThrshld = 0.80):  
    if type == 'prot':
        in_folder = bbhDirPt
    elif type == 'nucl':
        in_folder = bbhDirNt
    
    #Utilize the defined protein BLAST function
    run_blast(queryID, subjectID, type = type)
    run_blast(subjectID, queryID, type = type)
    
    query_lengths = get_gene_lens(queryID, type=type)
    subject_lengths = get_gene_lens(subjectID, type=type)
    
    #Define the output file of this BLAST
    out_file = os.path.normpath('%s/%s_vs_%s_parsed.csv'%(in_folder,queryID, subjectID))
    files=glob('%s/*_parsed.csv'%in_folder)
    
    #Combine the results of the protein BLAST into a dataframe
    print ('parsing BBHs for', queryID, subjectID)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,queryID, subjectID), sep='\t', names=cols)
    bbh = pd.merge(bbh, query_lengths) 
    bbh['COV'] = bbh['alnLength']/bbh['gene_length']
    
    bbh2=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,subjectID, queryID), sep='\t', names=cols)
    bbh2 = pd.merge(bbh2, subject_lengths) 
    bbh2['COV'] = bbh2['alnLength']/bbh2['gene_length']
    out = pd.DataFrame()
    
    # Filter the genes based on coverage
    bbh = bbh[bbh.COV >= covThrshld]
    bbh2 = bbh2[bbh2.COV >= covThrshld]
    
    #Delineate the best hits from the BLAST
    for g in bbh.gene.unique():
        res = bbh[bbh.gene==g]
        if len(res)==0:
            continue
        best_hit = res.loc[res.PID.idxmax()]
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene==best_gene]
        if len(res2)==0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        if g==best_gene2:
            best_hit['BBH'] = '<=>'
        else:
            best_hit['BBH'] = '->'
        out=pd.concat([out, pd.DataFrame(best_hit).transpose()])
    
    #Save the final file to a designated CSV file    
    out.to_csv(out_file)
    
#
# BLAST
####################################################################################################

# Get IDs from within genbank annotated genome files.
refGbkID = getGbkID(refGbk)
targetGbkID = getGbkID(targetGbk)

# Do BLASTs
if whatBlast == 'pt':
    # Parse genomes
    parse_genome(refGbk, type = 'prot', overwrite=0)
    parse_genome(targetGbk, type = 'prot', overwrite=0)
    # Do BLAST DBs
    make_blast_db(refGbkID, db_type = 'prot')
    make_blast_db(targetGbkID, db_type = 'prot')
    # Do blasts and get best hits files
    get_bbh(refGbkID, targetGbkID, type = 'prot', covThrshld = covThrshld)
elif whatBlast == 'nt':
    # Parse genomes
    parse_genome(refGbk, type = 'nucl', overwrite=0)
    parse_genome(targetGbk, type = 'nucl', overwrite=0)
    # Do BLAST DBs
    make_blast_db(refGbkID, db_type = 'nucl')
    make_blast_db(targetGbkID, db_type = 'nucl')
    # Do blasts and get best hits files
    get_bbh(refGbkID, targetGbkID, type = 'nucl', covThrshld = covThrshld)
elif whatBlast == 'both':
    # Parse genomes
    parse_genome(refGbk, type = 'prot', overwrite=0)
    parse_genome(targetGbk, type = 'prot', overwrite=0)
    parse_genome(refGbk, type = 'nucl', overwrite=0)
    parse_genome(targetGbk, type = 'nucl', overwrite=0)
    # Do BLAST DBs
    make_blast_db(refGbkID, db_type = 'prot')
    make_blast_db(targetGbkID, db_type = 'prot')
    make_blast_db(refGbkID, db_type = 'nucl')
    make_blast_db(targetGbkID, db_type = 'nucl')
    # Do blasts and get best hits files
    get_bbh(refGbkID, targetGbkID, type = 'prot', covThrshld = covThrshld)
    get_bbh(refGbkID, targetGbkID, type = 'nucl', covThrshld = covThrshld)