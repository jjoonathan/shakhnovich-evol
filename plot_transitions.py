import sys,os
import numpy as np
import pandas as pd
import seqtools as st
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from shared import *
from IPython import embed

r_hdrs, r_seqs = rah, ras
b_hdrs, b_seqs = bah, bas

def xmat(sequences):
	ret = np.zeros((la,la))
	for seq in sequences:
		for i in range(len(seq)-1):
			frm = chr2idx[ord(seq[i])]
			to = chr2idx[ord(seq[i+1])]
			ret[frm,to] += 1
	dashidx = chr2idx[ord('-')]
	ret[dashidx,dashidx] = 0  # Otherwise it drows out the signal
	return ret

def normalize(a):
	return a/a.max()

def mat2str(mat,hdr):
    hstr = ' \t'+'\t'.join(hdr)
    bodystrs = [hstr]
    for i in range(len(mat)):
        bodystrs.append(hdr[i]+'\t'+'\t'.join(str(int(x)) for x in mat[i]))
    return '\n'.join(bodystrs)


xr = xmat(r_seqs)
xb = xmat(b_seqs)
open('xr.tsv','w').write(mat2str(xr,alphabet))
open('xb.tsv','w').write(mat2str(xb,alphabet))
# xmax = max(xr.max(),xb.max())
# xr /=xmax
# xb /=xmax
xsum = xr+xb
xdif = xr-xb
# xsum = (xsum-xsum.mean())/xsum.std()
# xdif = (xdif-xdif.mean())/xdif.std()

def annotate(plt,ax,mat):
	plt.xticks(np.arange(la),list(alphabet))
	plt.yticks(np.arange(la),list(alphabet))
	#plt.colorbar(ims)
	for x in range(la):
		for y in range(la):
			v = mat[x][y]
			#if abs(v)<1: continue
			ax.annotate(('%.2g'%v), xy=(y, x),
						horizontalalignment='center',
						verticalalignment='center')

fig = plt.figure(figsize=(20,20))


ax = fig.add_subplot(221)
ims = ax.imshow(xr, cmap=plt.cm.Reds, interpolation='nearest')
annotate(plt,ax,xr)
plt.title('Transition Counts Summed Across All %i MM Sequences'%len(r_seqs))

ax = fig.add_subplot(222)
xr1 = xmat([asmap['NP_349605']])
ims = ax.imshow(xr1, cmap=plt.cm.Reds, interpolation='nearest')
annotate(plt,ax,xr1)
plt.title('Transition Counts for the NP_349605 MM Sequence (Q=1)')


ax = fig.add_subplot(223)
ims = ax.imshow(xb, cmap=plt.cm.Greys, interpolation='nearest')
annotate(plt,ax,xb)
plt.title('Transition Counts Summed Across All %i Outlier Sequences'%len(b_seqs))

ax = fig.add_subplot(224)
xb1 = xmat([asmap['NP_747233']])
ims = ax.imshow(xb1, cmap=plt.cm.Greys, interpolation='nearest')
annotate(plt,ax,xb1)
plt.title('Transition Counts for the NP_747233 Outlier Sequence')

fig.tight_layout(pad=3)
plt.savefig('transitions.pdf')
# plt.show()
