import sys,os, pickle
import seqtools as st
from shared import *
from glob import *
import networkx as nx
from IPython import embed

np_2_aa_seq = {}
for fname in glob('dat_nuccore/*.faa'):
	nc = os.path.splitext(fname)[0]
	for h,s in zip(*st.read_fasta(fname)):
		np = h[1:-1]
		np_2_aa_seq[np] = s

gene_graph = nx.Graph()
ncnc2rsds = {}  # maps (nc,nc) -> {(np,np)->rsd, ...}
for fname in glob('dat_rsd/*-*'):
	nc1,nc2 = os.path.basename(fname).split('-')
	nc1,nc2 = nc1,nc2 if nc1<nc2 else nc2,nc1
	rsds = {}
	ncnc2rsds[(nc1,nc2)] = rsds
	cum_rsd = 0
	cum_cnt = 1
	for l in open(fname):
		l = l.split('\t')
		if l[0] != 'OR': continue
		np1,np2,rsd = l[1:4]
		rsd = float(rsd)
		np1,np2 = np1,np2 if np1<np2 else np2,np1
		rsds[(np1,np2)] = rsd
		gene_graph.add_edge(np1,np2,weight=rsd)
		cum_rsd += rsd
		cum_cnt += 1

embed()
