from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import copy as cp
import os
from glob import glob
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

parser.add_argument('--covThrshld', 
	                  type=float,
                    help='Coverage threshold for considering a blast hit among the best hits.', 
                    default=0.80)
                    
parser.add_argument('--pidThrshld', 
	                  type=float,
                    help='Percentage of Identity threshold for considering a blast hit among the best hits.', 
                    default=95)
                    
parser.add_argument('--output', 
	                  type=str,
                    help='Path to output directory, where BLAST results will be saved.', 
                    default='/storage/PGO/results/mtb/bidBlastNuSeqs/')

args = parser.parse_args()
refGbk = args.refGbk
targetGbk = args.targetGbk
covThrshld = args.covThrshld
pidThrshld = args.pidThrshld
outDir = args.output

#
# Directory stuff
####################################################################################################

# Set blast directories
nuclDir = outDir + 'nucl/'
protDir = outDir + 'prot/'
targetDir = os.path.dirname(targetGbk) + '/'
blastnDir = outDir + 'bbh_nucl/'
blastpDir = outDir + 'bbh_prot/'


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
    
# Generate FASTA from the GBK files
def gbk2fasta(gbk_filename):
    ID = getGbkID(gbk_filename)
    outPath = os.path.dirname(gbk_filename)
    faa_filename = '%s/%s.fasta'%(outPath, ID)
    print(faa_filename)
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print ("Converting GenBank record %s" % seq_record.id)
        output_handle.write(">%s %s\n%s\n" % (
               seq_record.id,
               seq_record.description,
               seq_record.seq))

    output_handle.close()
    input_handle.close()
    
# Make database of a given genome.    
def make_blast_db(id, db_type='fullNucl'):
    import subprocess
    if db_type =='prot':
        out_file =os.path.normpath('%s%s.fa.pin'%(outDirPt, id))
        folder = os.path.dirname(out_file)
        files = glob('%s/*.fa.pin'%folder)
    elif db_type == 'nucl':
        out_file =os.path.normpath('%s%s.fna.nin'%(outDirNt, id))
        folder = os.path.dirname(out_file)
        files = glob('%s/*.fna.nin'%folder)
    elif db_type == 'fullNucl':
        out_file = os.path.normpath('%s/%s.fasta.nin'%(os.path.dirname(targetGbk), id))
        folder = os.path.dirname(out_file)
        files = glob('%s/*.fasta.nin'%folder)
    if out_file in files:
        print (id, 'already has a blast db')
        return
    if db_type == 'nucl':
        ext='fna'
    elif db_type == 'prot':
        ext='fa'
    elif db_type == 'fullNucl':
        ext='fasta'
        db_type = 'nucl'
    cmd_line='makeblastdb -in %s/%s.%s -dbtype %s' %(folder, id, ext, db_type)
    
    print ('making blast db with following command line...')
    print (cmd_line)
    subprocess.call(cmd_line, shell = True)
    
# Runs the BLASTn
def run_blastn(seq, db,outfmt=6,evalue=0.001,threads=1):
    import subprocess
    # Do the blast
    out = os.path.normpath('%s%s_vs_%s.txt'%(nuclDir, seq, db))
    seq = os.path.normpath(nuclDir+seq+'.fa')
    db = os.path.normpath(targetDir+db+'.fasta')
    
    cmd_line='blastn -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastn with following command line...')
    print (cmd_line)
    subprocess.call(cmd_line, shell = True)
    return out
    
# Gets sequence length 
def get_gene_lens(queryID, type = 'prot'):
    if type == 'prot':
        in_folder = protDir
        file = os.path.normpath('%s%s.fa'%(in_folder, queryID))
    elif type == 'nucl':
        in_folder = nuclDir
        file = os.path.normpath('%s%s.fna'%(in_folder, queryID))
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out

# Parses nucleotide BLAST result into a CSV file of best hits (above PID and coverage thresholds)
def parse_nucl_blast(dbID, queryID, covThrshld = 0.95, pidThrshld = 95):
    infile = '%s%s_vs_%s.txt'%(nuclDir, queryID, dbID)
    queryLens = get_gene_lens(queryID, type = 'nucl') ## Edit
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(infile, sep='\t', names=cols)
    queryLens = queryLens[[x in data['gene'].tolist() for x in queryLens['gene']]]
    queryLens = queryLens.reset_index(drop = True)
    data = pd.merge(data, queryLens) ## Edit
    data['COV'] = data['alnLength']/data['gene_length']
    data = data[(data['PID']>pidThrshld) & (data['COV']>covThrshld)]
    data2=data.groupby('gene').first()
    return data2.reset_index()
    
    
#
# BLASTn
####################################################################################################

# Obtain IDs
refID = getGbkID(refGbk)
targetID = getGbkID(targetGbk)

# Generate fasta file of full genome sequence from target genbank file
gbk2fasta(targetGbk)

# Create BLAST database for target full genome fasta file
make_blast_db(targetID, db_type = 'fullNucl')

# Run nucleotide blast of each one of the reference genes against full genomic sequence of target. 
run_blastn(refID, db = targetID)

# Parse result into a .CSV file
parsedNuclBlast = parse_nucl_blast(targetID, refID, covThrshld = covThrshld, pidThrshld = pidThrshld)

# Save CSV file.
parsedNuclBlast.to_csv('%s%s_vs_%s_parsed.csv'%(nuclDir, refID, targetID))
print('Parsed nucleotide blast of %s vs %s saved at %s. Filename is %s_vs_%s_parsed.csv'%(refID, targetID, nuclDir, refID, targetID))