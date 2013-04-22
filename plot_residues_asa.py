from collections import Counter
from math import log, pi
import seqtools as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import plot_asa
import os
import autowrap
from IPython import embed
from shared import *

# Use these parameters if creating an 8.5x11 PDF
# font_size = 8
# mp.rcParams['font.size'] = font_size
# mp.rcParams['axes.labelsize'] = font_size
# mp.rcParams['axes.linewidth'] = .5
# mp.rcParams['lines.linewidth'] = .5
# mp.rcParams['patch.linewidth'] = .5
# mp.rcParams['axes.titlesize'] = font_size
# mp.rcParams['legend.fontsize'] = font_size
# mp.rcParams['xtick.labelsize'] = font_size
# mp.rcParams['ytick.labelsize'] = font_size

r_hdrs, r_seqs = rah, ras
b_hdrs, b_seqs = bah, bas

fig = plt.figure(figsize=(20,20))
rm, bm, n = len(r_seqs), len(b_seqs), len(r_seqs[0])
x = np.arange(n)


########### Gross copypasta, but this is just a one-off script...
sr = entropy(r_seqs)
sb = entropy(b_seqs)
ar = red_nasa
ab = blk_nasa
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
plt_ar, = pax.plot(x+.5,ar/.5,color='g')
plt_sr, = pax.plot(x+.5,sr/4,color='b')
l3 = pax.legend([plt_ar,plt_sr],['ASA','S'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
# plt.ylabel('S(nats/sequence), ASA(A^2/4)')
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
plt_ar, = pax.plot(x+.5,ar/.5,color='g')
plt_sr, = pax.plot(x+.5,sr/4,color='b')
l3 = pax.legend([plt_ar,plt_sr],['ASA','S'],loc=2) ###################
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
plt_ar, = pax.plot(x+.5,ar/.5,color='g')
plt_sr, = pax.plot(x+.5,sr/4,color='b')
l3 = pax.legend([plt_ar,plt_sr],['ASA','S'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.xlabel('AA Position')
plt.title('(Fraction Basic, Fraction WY, Entropy, ASA) Within (MM,Outlier) Populations')

fig.tight_layout()
autowrap.autowrap(fig)
plt.savefig('AA+S+ASA.pdf')

















fig = plt.figure(figsize=(20,20))
windowlen = 10
window = np.ones(windowlen)/windowlen


########### Gross copypasta, but this is just a one-off script...
# sr = entropy(r_seqs)
# sb = entropy(b_seqs)
# ar = plot_asa.red_asa/4
# ab = plot_asa.blk_asa/4
sr = np.convolve(sr,window,'same')
sb = np.convolve(sb,window,'same')
ar = np.convolve(ar,window,'same')
ab = np.convolve(ab,window,'same')
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
plt_sr, = pax.plot(x+.5,sr/4,'r--',linewidth=1)
plt_sb, = pax.plot(x+.5,sb/4,'k--',linewidth=1)
plt_ar, = pax.plot(x+.5,ar/.5,'r',linewidth=2)
plt_ab, = pax.plot(x+.5,ab/.5,'k',linewidth=2)
l3 = pax.legend([plt_sr,plt_sb,plt_ar,plt_ab],['S','S','ASA','ASA'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
# plt.ylabel('S(nats/sequence), ASA(A^2/4)')
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
plt_sr, = pax.plot(x+.5,sr/4,'r--',linewidth=1)
plt_sb, = pax.plot(x+.5,sb/4,'k--',linewidth=1)
plt_ar, = pax.plot(x+.5,ar/.5,'r',linewidth=2)
plt_ab, = pax.plot(x+.5,ab/.5,'k',linewidth=2)
l3 = pax.legend([plt_sr,plt_sb,plt_ar,plt_ab],['S','S','ASA','ASA'],loc=2) ###################
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
plt_sr, = pax.plot(x+.5,sr/4,'r--',linewidth=1)
plt_sb, = pax.plot(x+.5,sb/4,'k--',linewidth=1)
plt_ar, = pax.plot(x+.5,ar/.5,'r',linewidth=2)
plt_ab, = pax.plot(x+.5,ab/.5,'k',linewidth=2)
l3 = pax.legend([plt_sr,plt_sb,plt_ar,plt_ab],['S','S','ASA','ASA'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.xlabel('AA Position')
plt.title('(Fraction Basic, Fraction WY, Entropy, ASA) Within (MM,Outlier) Populations')

fig.tight_layout()
autowrap.autowrap(fig)
plt.savefig('AA+S+ASA_smooth.pdf')
