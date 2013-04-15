from collections import Counter
from math import log, pi
import seqtools as st
import numpy as np
import matplotlib.pyplot as plt
import plot_asa
import os
from IPython import embed
rnames="""NP_267306
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
bnames="""AAA22853
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

def count_by_site(letters,seqs):
    accum = np.zeros(len(seqs[0]))
    for i in xrange(len(seqs)):
        for j in xrange(len(seqs[0])):
            if seqs[i][j] in letters:
                accum[j] += 1
    return accum

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
    return ret/m

charged = 'RKHYCDE'
acidic = 'DE'
basic = 'KR'
polar = 'STNQ' #'NCQSTY'
hydrophobic = 'AILFWF'

aligned_hdrs, aligned_seqs = st.read_fasta(open('pdb_seqs_aln.faa'))
aligned_hdrs = [h[1:-1] for h in aligned_hdrs]
r_hdrs, r_seqs = [],[]
b_hdrs, b_seqs = [],[]
for h,s in zip(aligned_hdrs,aligned_seqs):
    if h in rnames:
        r_hdrs.append(h)
        r_seqs.append(s)
    elif h in bnames:
        b_hdrs.append(h)
        b_seqs.append(s)
assert(len(r_hdrs)>0 and len(b_hdrs)>0)
rb_hdrs = r_hdrs + b_hdrs
rb_seqs = r_seqs + b_seqs

fig = plt.figure(figsize=(20,20))
rm, bm, n = len(r_seqs), len(b_seqs), len(r_seqs[0])
x = np.arange(n)

########### Gross copypasta, but this is just a one-off script...
sr = entropy(r_seqs) # (S_red, S_blk) are replaced with (S_red, ASA_red)
sb = plot_asa.redm/4
# sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
pr, pb = count_by_site(polar,r_seqs)/rm, count_by_site(polar,b_seqs)/bm
wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pax = fig.add_subplot(311)
plt_rp = pax.bar(x,-pr,1,color='r')
plt_bp = pax.bar(x+(1-.3)/2,-pb,.3,color='k')
l1 = pax.legend([plt_rp,plt_bp],['Outlier Polar/Total','MM Polar/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['Outlier WY/Total','MM WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='g')
plt_sb, = pax.plot(x+.5,sb/3,color='b')
l3 = pax.legend([plt_sr,plt_sb],['S','ASA'],loc=2) ###################
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
l1 = pax.legend([plt_rp,plt_bp],['Outlier DE/Total','MM DE/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['Outlier WY/Total','MM WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='g')
plt_sb, = pax.plot(x+.5,sb/3,color='b')
l3 = pax.legend([plt_sr,plt_sb],['S','ASA'],loc=2) ###################
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
l1 = pax.legend([plt_rp,plt_bp],['Outlier KR/Total','MM KR/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['Outlier WY/Total','MM WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='g')
plt_sb, = pax.plot(x+.5,sb/3,color='b')
l3 = pax.legend([plt_sr,plt_sb],['S','ASA'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.xlabel('AA Position')
plt.title('(Fraction Basic, Fraction WY, Entropy, ASA) Within (MM,Outlier) Populations')

plt.savefig('AA+S+ASA.pdf')
