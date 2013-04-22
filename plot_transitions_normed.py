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
open('xr.tsv','w').write(mat2str(xr,alphabet))
xb = (xmat(b_seqs)*len(r_seqs)/len(b_seqs)).round()
open('xb.tsv','w').write(mat2str(xb,alphabet))
# xmax = max(xr.max(),xb.max())
# xr /=xmax
# xb /=xmax
xsum = xr+xb
xdif = 2*(xr-xb)/(xr+xb)
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
plt.title('Transition Map for Red (MM) Sequences')

ax = fig.add_subplot(224)
ims = ax.imshow(xb, cmap=plt.cm.Greys, interpolation='nearest')
annotate(plt,ax,xb)
plt.title('Renormalized Transition Map for Black (Outlier) Sequences')

ax = fig.add_subplot(222)
ims = ax.imshow(xdif, cmap=plt.cm.Purples, interpolation='nearest')
annotate(plt,ax,xdif)
plt.title('(Tred-Tblack)/mean(Tred,Tblack)')

xdif = xdif.flatten()
xdif = xdif[~np.isnan(xdif)]
xdif_kde = gaussian_kde(xdif)
x = np.linspace(-3,3,300)
ax = fig.add_subplot(223)
plt.plot(x, xdif_kde(x), 'purple', label='Estimated PDF of Differences Between Transition Counts')
plt.plot(xdif, np.random.rand(len(xdif))*.1, color='purple', marker='.', linestyle='none', label='Transition Counts (jittered in y direction for visibility only)')
plt.legend()



fig.tight_layout(pad=3)
plt.savefig('transitions_normed.pdf')
# plt.show()
