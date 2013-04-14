import sys,os
from Bio.PDB.PDBParser import PDBParser
from glob import glob
import seqtools as st
import numpy as np
from IPython import embed

aln_h, aln_s = st.read_fasta('pdb_seqs_aln.faa')
aln_h = [h[1:-1] for h in aln_h]  # Strip > and \n from '>SEQID_STR\n'
id2aln = dict(zip(aln_h,aln_s))

pdbs = glob('/mnt/hgfs/D/jon/Documents/2012/summer/pdbs/o_*')
pdbparser = PDBParser(PERMISSIVE=1)
ASA_ids = []
ASA_vals = []  # Array of (ASA array of len(aligned sequence))
for fpath in pdbs:
	name = os.path.basename(fpath)
	protid = '_'.join(name.split('_')[1:-2]).upper()
	structure = pdbparser.get_structure('struct',fpath)
	models = structure.get_list()    # Hardcoded to get model 0, chain 0
	chain = models[0].get_list()[0]
	residues = chain.get_list()
	ASAs = [r['CA'].get_bfactor() for r in residues]
	aa_alignment = id2aln[protid]
	aligned_ASAs = np.ones(len(aa_alignment))*np.nan
	asa_idx = 0
	for aln_idx in range(len(aa_alignment)):
		if aa_alignment[aln_idx]!='-':
			try:
				aligned_ASAs[aln_idx] = ASAs[asa_idx]
			except:
				embed()
			asa_idx += 1
	if asa_idx!=len(ASAs):
		print "Sequence length mismatch"
		embed()
	ASA_ids.append(protid)
	ASA_vals.append(aligned_ASAs)

of = open('ASAs.tsv','w')
of.write('\t'.join(ASA_ids))
of.write('\n')
m,n = len(ASA_vals[0]), len(ASA_vals)  # m=#rows, n=#cols
for i in range(m):  # One aligned residue per row
	for j in range(n):
		of.write(str(ASA_vals[j][i]))
		of.write("\t" if j!=n-1 else "\n")



