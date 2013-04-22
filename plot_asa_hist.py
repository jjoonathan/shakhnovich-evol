import os,sys
from shared import *
import numpy as np
from matplotlib import pyplot as plt
from IPython import embed
from scipy.stats import gaussian_kde

def ASAs_for_sites(needles, haystack_h, haystack_s):
    ret = []
    seqlen = len(haystack_s[0])
    for i in range(len(haystack_h)):
        seq = haystack_s[i]
        try:
            asas = asa_frame_norb[haystack_h[i]]
        except KeyError:
            print "Couldn't find ASAs for %s"%haystack_h[i]
            continue
        if len(asas)!=seqlen:
            print "Mismatch between ASA, aligned seq len."
            embed()
        for j in range(seqlen):
            if seq[j] in needles and not np.isnan(asas[j]):
                amax = asa_max(seq[j-1:j+2])
                amean = asa_mean(seq[j-1:j+2])
                if np.isnan(amax) or np.isnan(amean):
                    print "Triplet had NAN  or 0 asa_max or asa_mean: "+seq[j-1:j+2]
                    embed()
                if amax==0 or amean==0:
                    print "Triplet had NAN or 0 asa_max or asa_mean: "+seq[j-1:j+2]
                    embed()
                ret.append(asas[j]/amean)
    return ret

acid_r = gaussian_kde(ASAs_for_sites(acidic, rah, ras))
acid_b = gaussian_kde(ASAs_for_sites(acidic, bah, bas))
basic_r = gaussian_kde(ASAs_for_sites(basic, rah, ras))
basic_b = gaussian_kde(ASAs_for_sites(basic, bah, bas))
polar_r = gaussian_kde(ASAs_for_sites(polar, rah, ras))
polar_b = gaussian_kde(ASAs_for_sites(polar, bah, bas))

x = np.linspace(0,.4,500)
fig = plt.figure(figsize=(20,15))
ax = fig.add_subplot(221)
ax.plot(x,acid_r(x),'r', x,acid_b(x),'k')
ax.set_ylabel('Density of Nearby ASA Observations')
ax.set_xlabel('ASA / (Max ASA for triplet)')
ax.legend(('Michaelis-Menton','Outlier'))
ax.set_title('Continuous Histogram (KDE) for ASA of Acidic Residues')

ax = fig.add_subplot(222)
ax.plot(x,basic_r(x),'r', x,basic_b(x),'k')
ax.set_ylabel('Density of Nearby ASA Observations')
ax.set_xlabel('ASA / (Max ASA for triplet)')
ax.legend(('Michaelis-Menton','Outlier'))
ax.set_title('Continuous Histogram (KDE) for ASA of Basic Residues')

ax = fig.add_subplot(223)
ax.plot(x,polar_r(x),'r', x,polar_b(x),'k')
ax.set_ylabel('Density of Nearby ASA Observations')
ax.set_xlabel('ASA / (Max ASA for triplet)')
ax.legend(('Michaelis-Menton','Outlier'))
ax.set_title('Continuous Histogram (KDE) for ASA of Polar Residues')

fig.tight_layout()
plt.savefig('asa_hist.pdf')

