import sys,os
from Bio.PDB.PDBParser import PDBParser
from glob import glob
import seqtools as st
import numpy as np
from IPython import embed

three2one = dict([line.split(' ') for line in """Ala A
Arg R
Asn N
Asp D
Cys C
Glu E
Gln Q
Gly G
His H
Ile I
Leu L
Lys K
Met M
Phe F
Pro P
Ser S
Thr T
Trp W
Tyr Y
Val V""".upper().splitlines()])

aln_h, aln_s = st.read_fasta('PROTEIN/ALIGNED.fasta')
aln_h = [h[1:-1] for h in aln_h]  # Strip > and \n from '>SEQID_STR\n'
id2aln = dict(zip(aln_h,aln_s))

pdbs = glob('/mnt/hgfs/D/jon/Documents/2012/summer/pdbs/*')
pdbparser = PDBParser(PERMISSIVE=1)
AA_hdrs = []
AA_seqs = []  # Array of (ASA array of len(aligned sequence))
for fpath in pdbs:
	name = os.path.basename(fpath)
	protid = '_'.join(name.split('_')[1:-2]).upper()
	structure = pdbparser.get_structure('struct',fpath)
	models = structure.get_list()    # Hardcoded to get model 0, chain 0
	chain = models[0].get_list()[0]
	residues = chain.get_list()
	AAs = ''.join([three2one[r.get_resname()] for r in residues])
	AA_hdrs.append('>'+protid)
	AA_seqs.append(AAs)

st.write_fasta((AA_hdrs,AA_seqs),'pdb_seqs.faa')
