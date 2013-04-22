import seqtools as st
import pandas as pd
import numpy as np
from math import log,pi
from collections import Counter
from IPython import embed

def count_by_site(letters, seqs):
    """ Returns an array the same length as any of seqs that
    counts the number of sequences in seqs containing one of
    <letters> at each site. """
    accum = np.zeros(len(seqs[0]))
    for i in xrange(len(seqs)):
        for j in xrange(len(seqs[0])):
            if seqs[i][j] in letters:
                accum[j] += 1
    return accum

def find_sites(letters, threshold, seqs):
    """ Returns the 0-based indexes of residues which are listed in the string
    <letters>. At lest the fraction <threshold> of sequences must have residues
    listed in <letters> for the site to count. """
    ret = []
    seqlen = len(seqs[0])
    for j in range(seqlen):
        matching_residues = 0
        for seq in seqs:
            if seq[j] in letters:
                matching_residues += 1
        if matching_residues*1.0/seqlen > threshold:
            ret.append(j)
    return ret


def entropy(seqs,traditional=False):
    """ Returns an array with each element corresponding to the sequence entropy
    at that position. S = log(N!/sum_i n_i!), approximated by Sterling """
    m = len(seqs)
    n = len(seqs[0])
    ret = np.ones(n)*(m*log(m)) - (0 if traditional else log(2*pi*m)/2)
    for j in range(n):
        ctr = Counter(seqs[i][j] for i in range(m))
        for chr,count in ctr.iteritems():
            ret[j] -= count*log(count) - (0 if traditional else log(2*pi*count)/2)
    return ret/m/log(2)

charged = 'RKHYCDE'
acidic = 'DE'
basic = 'KR'
polar = 'STNQ' #'STNQCY'
hydrophobic = 'AILFWF'

bnames="""NP_267306
NP_349605
NP_471321
NP_883118
YP_069170
YP_115129
YP_818239
YP_928696
YP_002505644
YP_002996559
YP_003068918
ZP_04443390
ZP_04617065
NP_217278
NP_297275
NP_641196
YP_342295
YP_001312015
YP_001454852
ZP_04674563
ZP_06040137""".splitlines()

rnames="""AAA22853
CAA58993
NP_246832
YP_315797
YP_785039
YP_857866
ZP_04783600
ZP_06186569
ABC47008
NP_747233
YP_203667
YP_341143
YP_804567
YP_002274348
YP_002932150
YP_003262789
YP_003518975
ZP_01544314
ZP_02464292
ZP_05320398
ZP_06350196""".splitlines()

alphabet = ['A','C','E','D','G','F','I','M','K','V','-','L','N','Q','P','S','R','T','W','H','Y']
alphabet3 = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'MET', 'LYS', 'VAL', '---', 'LEU', 'ASN', 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'HIS', 'TYR']
la = len(alphabet)
chr2idx = np.ones(256)*-1
for i in range(la):
	chr2idx[ord(alphabet[i])] = i
idx2chr = alphabet

aligned_hdrs, aligned_seqs = st.read_fasta(open('pdb_seqs_aln.faa'))
aligned_hdrs = [h[1:-1] for h in aligned_hdrs]
rah, ras = [],[]  # Red Aligned Headers, Red Aligned Sequences
bah, bas = [],[]
rahp, rasp = [],[] # Red Aligned Headers from PDB
bahp, basp = [],[]
for h,s in zip(aligned_hdrs,aligned_seqs):
    if h in rnames:
        rah.append(h.strip())
        ras.append(s.strip())
    elif h in bnames:
        bah.append(h.strip())
        bas.append(s.strip())
    elif h.startswith('RED'):
        rahp.append(h.strip())
        rasp.append(s.strip())
    elif h.startswith('BLK'):
        bahp.append(h.strip())
        basp.append(s.strip())
    else:
        print "WARNING: shared.py couldn't classify %s as red/black"%h.strip()
assert(len(rah)>0 and len(bah)>0)
asmap = dict(zip(rah+bah,ras+bas))
aspmap = dict(zip(rahp+bahp,rasp+basp))

asa_frame=pd.read_table('ASAs.tsv')
asa_frame_norb = asa_frame.copy()
asa_frame_norb.columns = [cname[4:] for cname in asa_frame.columns]
red_asa_frame=pd.DataFrame(asa_frame[[k for k in asa_frame if k.startswith('RED')]])
blk_asa_frame=pd.DataFrame(asa_frame[[k for k in asa_frame if k.startswith('BLK')]])
red_asa = red_asa_frame.mean(axis=1, skipna=True)
red_asa.ffill(inplace=True); red_asa.bfill(inplace=True)
red_asa = red_asa.values
red_asa_std = blk_asa_frame.std(axis=1, skipna=True).values
blk_asa = blk_asa_frame.mean(axis=1, skipna=True)
blk_asa.ffill(inplace=True); blk_asa.bfill(inplace=True)
blk_asa = blk_asa.values
blk_asa_std = blk_asa_frame.std(axis=1, skipna=True).values
n = len(blk_asa_frame[blk_asa_frame.keys()[0]])  # Convention: m=#rows, n=#cols
mr = len(rah)
mb = len(bah)

try:
    ASA_n8n = open('ASA_normalization.npy')
    asa_max_mat = np.load(ASA_n8n)
    asa_mean_mat = np.load(ASA_n8n)
    asa_std_mat = np.load(ASA_n8n)
    def asa_max(i,j=None,k=None):
        if j==None: i,j,k = i
        return asa_max_mat[chr2idx[ord(i)],chr2idx[ord(j)],chr2idx[ord(k)]]
    def asa_mean(i,j=None,k=None):
        if j==None: i,j,k = i
        return asa_mean_mat[chr2idx[ord(i)],chr2idx[ord(j)],chr2idx[ord(k)]]
    def asa_std(i,j=None,k=None):
        if j==None: i,j,k = i
        return asa_std_mat[chr2idx[ord(i)],chr2idx[ord(j)],chr2idx[ord(k)]]
    del ASA_n8n
    nasa_frame = asa_frame.copy()
    for seqname in nasa_frame:
        nasa = nasa_frame[seqname]
        seq = asmap[seqname[4:]]
        for j in range(1,len(seq)-1):
            nasa[j] /= asa_mean(seq[j-1:j+2])
    nasa_frame.ffill(inplace=True)
    nasa_frame.bfill(inplace=True)
    red_nasa_frame = nasa_frame[[k for k in nasa_frame if k.startswith('RED')]]
    blk_nasa_frame = nasa_frame[[k for k in nasa_frame if k.startswith('BLK')]]
    red_nasa = red_nasa_frame.mean(axis=1, skipna=True).values
    blk_nasa = blk_nasa_frame.mean(axis=1, skipna=True).values
except IOError:
    pass
