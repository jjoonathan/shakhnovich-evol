from collections import Counter
from math import log, pi
import seqtools as st
import numpy as np
import matplotlib.pyplot as plt
import plot_asa
import os
from IPython import embed
from shared import *

r_hdrs, r_seqs = rah, ras
b_hdrs, b_seqs = bah, bas

fig = plt.figure(figsize=(20,20))
rm, bm, n = len(r_seqs), len(b_seqs), len(r_seqs[0])
x = np.arange(n)

########### Gross copypasta, but this is just a one-off script...
sr = entropy(r_seqs) # (S_red, S_blk) are replaced with (S_red, ASA_red)
sb = plot_asa.red_asa/4
# sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
pr, pb = count_by_site(polar,r_seqs)/rm, count_by_site(polar,b_seqs)/bm
wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pax = fig.add_subplot(311)
plt_rp = pax.bar(x,-pr,1,color='r')
plt_bp = pax.bar(x+(1-.3)/2,-pb,.3,color='k')
l1 = pax.legend([plt_rp,plt_bp],['MM Polar/Total','Outlier Polar/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='g')
plt_sb, = pax.plot(x+.5,sb/3,color='b')
l3 = pax.legend([plt_sr,plt_sb],['ASA','S'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.ylabel('S(nats/sequence), ASA(A^2/4)')
plt.title('(Fraction Polar, Fraction WY, Entropy, ASA) Within (MM,Outlier) Populations')

# sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
# wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pr, pb = count_by_site('DE',r_seqs)/rm, count_by_site('DE',b_seqs)/bm
pax = fig.add_subplot(312)
plt_rp = pax.bar(x,-pr,1,color='r')
plt_bp = pax.bar(x+(1-.3)/2,-pb,.3,color='k')
l1 = pax.legend([plt_rp,plt_bp],['MM DE/Total','Outlier DE/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='g')
plt_sb, = pax.plot(x+.5,sb/3,color='b')
l3 = pax.legend([plt_sr,plt_sb],['ASA','S'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.title('(Fraction Acidic, Fraction WY, Entropy, ASA) Within (MM,Outlier) Populations')

# sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
# wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pr, pb = count_by_site('KR',r_seqs)/rm, count_by_site('KR',b_seqs)/bm
pax = fig.add_subplot(313)
plt_rp = pax.bar(x,-pr,1,color='r')
plt_bp = pax.bar(x+(1-.3)/2,-pb,.3,color='k')
l1 = pax.legend([plt_rp,plt_bp],['MM KR/Total','Outlier KR/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='g')
plt_sb, = pax.plot(x+.5,sb/3,color='b')
l3 = pax.legend([plt_sr,plt_sb],['ASA','S'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.xlabel('AA Position')
plt.title('(Fraction Basic, Fraction WY, Entropy, ASA) Within (MM,Outlier) Populations')

plt.savefig('AA+S+ASA.pdf')
