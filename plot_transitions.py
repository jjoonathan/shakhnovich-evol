import sys,os
import numpy as np
import pandas as pd
import seqtools as st
import matplotlib.pyplot as plt
from IPython import embed

aligned_hdrs, aligned_seqs = st.read_fasta(open('pdb_seqs_aln.faa'))
aligned_hdrs = [h[1:-1] for h in aligned_hdrs]
aa = zip(aligned_hdrs, aligned_seqs)
r_hdrs, r_seqs = zip(*[(h,s) for h,s in aa if h.startswith('RED')])
b_hdrs, b_seqs = zip(*[(h,s) for h,s in aa if h.startswith('BLK')])

#Compute transition alphabet (set of 'A'->'B' tuples that appear)
alphabet = list(set(''.join(r_seqs+b_seqs)))
la = len(alphabet)
chr2idx = np.ones(256)*-1
for i in range(la):
	chr2idx[ord(alphabet[i])] = i
idx2chr = alphabet

def xmat(sequences):
	ret = np.zeros((la,la))
	for seq in sequences:
		for i in range(len(seq)-1):
			frm = chr2idx[ord(seq[i])]
			to = chr2idx[ord(seq[i+1])]
			ret[frm,to] += 1
	ret /= np.sum(ret)
	dashidx = chr2idx[ord('-')]
	ret[dashidx,dashidx] = 0  # Otherwise it drows out the signal
	return ret

def normalize(a):
	return a/a.max()

xr = normalize(xmat(r_seqs))
xb = normalize(xmat(b_seqs))
xmax = max(xr.max(),xb.max())
xr /=xmax
xb /=xmax
xsum = xr+xb
xdif = xr-xb
xsum = (xsum-xsum.mean())/xsum.std()
xdif = (xdif-xdif.mean())/xdif.std()

def annotate(plt,ax,mat):
	plt.xticks(np.arange(la),list(alphabet))
	plt.yticks(np.arange(la),list(alphabet))
	plt.colorbar(ims)
	for x in range(la):
		for y in range(la):
			v = mat[x][y]
			if abs(v)<1: continue
			ax.annotate(('%.2g'%v).lstrip('0'), xy=(y, x),
						horizontalalignment='center',
						verticalalignment='center')

fig = plt.figure(figsize=(15,30))
ax = fig.add_subplot(211)
ims = plt.imshow(xsum, cmap=plt.cm.gray, interpolation='nearest')
annotate(plt,ax,xsum)
plt.title('Average Transition Map in Std. Devs. from Mean')

ax = fig.add_subplot(212)
ims=plt.imshow(xdif, cmap=plt.cm.gray, interpolation='nearest')
plt.title('Difference Transition Map in Std. Devs. from Mean')
annotate(plt,ax,xdif)

plt.savefig('transitions.pdf')
# plt.show()
