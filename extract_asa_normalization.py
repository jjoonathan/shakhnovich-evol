import sys,os
import numpy as np
from shared import *
import itertools

la = len(alphabet)
asa_max = np.zeros((la,la,la))
asa_mean = np.zeros((la,la,la))
asa_std = np.zeros((la,la,la))
expected_aa = set(alphabet)-set('-')

for i in range(la):
    a3 = alphabet3[i]
    if a3=='---': continue
    curt1 = None
    for l in open('HOA_hist_data/HOA_'+a3+'.txt'):
        triplet, amax, amean, astd = l.split('\t')
        t1,_,t2,_,t3 = triplet
        idx = (chr2idx[ord(t1)], chr2idx[ord(t2)], chr2idx[ord(t3)])
        asa_max[idx] = amax
        asa_mean[idx] = amean
        asa_std[idx] = astd

dash = chr2idx[ord('-')]

def goodmean(nda):  # Mean over non-zero, non-nan, non-inf members of nda
    idx = np.logical_and(np.isfinite(nda),nda!=0)
    return nda[idx].mean()

# '-AA' marginals
for i,j in itertools.product(range(la),range(la)):
    if idx2chr[i]=='-' or idx2chr[j]=='-': continue
    asa_max[dash,i,j] = goodmean(asa_max[:,i,j])
    asa_mean[dash,i,j] = goodmean(asa_mean[:,i,j])
    asa_std[dash,i,j] = goodmean(asa_std[:,i,j])  # Yes, 1-mean is dubious
# 'AA-' marginals
for i,j in itertools.product(range(la),range(la)):
    if idx2chr[i]=='-' or idx2chr[j]=='-': continue
    asa_max[i,j,dash] = goodmean(asa_max[:,i,j])
    asa_mean[i,j,dash] = goodmean(asa_mean[:,i,j])
    asa_std[i,j,dash] = goodmean(asa_std[:,i,j])
# '-A-' marginals
for i in range(la):
    if idx2chr[i]=='-': continue
    asa_max[dash,i,dash] = goodmean(asa_max[:,i,:])
    asa_mean[dash,i,dash] = goodmean(asa_mean[:,i,:])
    asa_std[dash,i,dash] = goodmean(asa_std[:,i,:])
asa_max[:,dash,:] = np.nan
asa_mean[:,dash,:] = np.nan
asa_std[:,dash,:] = np.nan

def fix_zeros(nda):
    for i,j,k in np.transpose((nda==0).nonzero()):
        if nda[k,j,i]!=0:
            nda[i,j,k] = nda[k,j,i]
        else:
            nda[i,j,k] = goodmean(nda[:,j,:])
fix_zeros(asa_max)
fix_zeros(asa_mean)
fix_zeros(asa_std)

of = open('ASA_normalization.npy','w')
np.save(of,asa_max)
np.save(of,asa_mean)
np.save(of,asa_std)
