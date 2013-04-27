import sys,os
import numpy as np
import pandas as pd
import seqtools as st
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from shared import *
from random import Random
from IPython import embed
rand = Random()
rand.seed(0)

r_hdrs, r_seqs = rah, ras
b_hdrs, b_seqs = bah, bas
dashidx = chr2idx[ord('-')]

def xmat(sequences):
	ret = np.zeros((la,la))
	for seq in sequences:
		for i in range(len(seq)-1):
			frm = chr2idx[ord(seq[i])]
			to = chr2idx[ord(seq[i+1])]
			ret[frm,to] += 1
	ret[dashidx,dashidx] = 0  # Otherwise it drows out the signal
	return ret

def xmat_null(sequences, shufs=100):
    """ Computes the marginal transition matrix based on AA frequency"""
    xmats = []
    npseqs = np.array([np.array(s,'c') for s in sequences])
    npseqs_shape = npseqs.shape
    npseqs = npseqs.reshape(-1).view(np.uint8)
    for i,v in enumerate(npseqs):
        npseqs[i] = chr2idx[v]
    npseqs = npseqs.reshape(npseqs_shape)
    seq_idxs = range(len(npseqs))
    chr_idxs = range(len(npseqs[0])-1)
    for i in range(shufs):
        shufseqs = npseqs.flatten()
        np.random.shuffle(shufseqs)
        shufseqs = np.reshape(shufseqs,npseqs_shape)
        xmat = np.zeros((la,la))
        for seq in shufseqs:
            for k in chr_idxs:
                xmat[seq[k],seq[k+1]] += 1
        xmats.append(xmat)
    return np.mean(xmats,axis=0), np.std(xmats,axis=0)

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
xr_mu,xr_std = xmat_null(r_seqs)
xb_mu,xb_std = xmat_null(b_seqs)
open('xr.tsv','w').write(mat2str(xr,alphabet))
open('xb.tsv','w').write(mat2str(xb,alphabet))


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

# Z-Score Matrix of {Red,Black} Xitions | Null Model Built Off {Red,Blk} Xitions
m11 = (xr-xr_mu)/xr_std
m22 = (xb-xb_mu)/xb_std
m12 = (xr-xb_mu)/xb_std - m11
m21 = (xb-xr_mu)/xr_std - m22
for m in [m11,m12,m21,m22]:
    m[dashidx,:] = 0
    m[:,dashidx] = 0

ax = fig.add_subplot(221)
ims = ax.imshow(np.abs(m11), cmap=plt.cm.Reds, interpolation='nearest')
annotate(plt,ax,m11)
plt.title('Z-Scores of MM Transitions in MM Shuffled Null Model')

ax = fig.add_subplot(222)
ims = ax.imshow(np.abs(m11-m22), cmap=plt.cm.Purples, interpolation='nearest')
annotate(plt,ax,m11-m22)
plt.title('Red Z-Score Matrix - Black Z-Score Matrix')

ax = fig.add_subplot(223)
xr[dashidx,:] = 0
xr[:,dashidx] = 0
ims = ax.imshow(xr, cmap=plt.cm.Greens, interpolation='nearest')
annotate(plt,ax,xr)
plt.title('Transition Counts Summed Across All MM Sequences ($\propto$ stability)')

# ax = fig.add_subplot(223)
# ims = ax.imshow(np.abs(m12), cmap=plt.cm.Purples, interpolation='nearest')
# annotate(plt,ax,m12)
# plt.title('Z-Scores of MM Transitions in Outlier Null Model - Z-Score in MM Null Model')
#
# ax = fig.add_subplot(222)
# ims = ax.imshow(np.abs(m21), cmap=plt.cm.Purples, interpolation='nearest')
# annotate(plt,ax,m21)
# plt.title('Z-Scores of Outlier Transitions in MM Null Model - Z-Score in Outlier Null Model')

ax = fig.add_subplot(224)
ims = ax.imshow(np.abs(m22), cmap=plt.cm.Greys, interpolation='nearest')
annotate(plt,ax,m22)
plt.title('Z-Scores of Outlier Transitions in Outlier Shuffled Null Model')

fig.tight_layout(pad=3)
plt.savefig('transitions_normed.pdf')
# plt.show()
