# Turn back, all ye who look for good code here.
# Hardcoded values, copy+paste+modify, and general sloppiness run amok!
from collections import Counter
from math import log, pi
import seqtools as st
import numpy as np
import matplotlib.pyplot as plt
import plot_asa
import os
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

aligned_hdrs, aligned_seqs = st.read_fasta(open('PROTEIN/ALIGNED.fasta'))
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
    else:
        print "UNCLASSIFIED HEADER: "+h
assert(len(r_hdrs)>0 and len(b_hdrs)>0)
rb_hdrs = r_hdrs + b_hdrs
rb_seqs = r_seqs + b_seqs

fig = plt.figure(figsize=(20,20))
rm, bm, n = len(r_seqs), len(b_seqs), len(r_seqs[0])
x = np.arange(n)

########### Gross copypasta, but this is just a one-off script...
sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
pr, pb = count_by_site(polar,r_seqs)/rm, count_by_site(polar,b_seqs)/bm
wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pax = fig.add_subplot(311)
plt_rp = pax.bar(x,-pr,1,color='r')
plt_bp = pax.bar(x+(1-.3)/2,-pb,.3,color='k')
l1 = pax.legend([plt_rp,plt_bp],['MM Polar/Total','Outlier Polar/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='r')
plt_sb, = pax.plot(x+.5,sb/3,color='k')
l3 = pax.legend([plt_sr,plt_sb],['S(MM)','S(Outlier)'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.title('(Fraction Polar, Fraction WY, Entropy) Within (MM,Outlier) Populations')

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
plt_sr, = pax.plot(x+.5,sr/3,color='r')
plt_sb, = pax.plot(x+.5,sb/3,color='k')
l3 = pax.legend([plt_sr,plt_sb],['S(MM)','S(Outlier)'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.title('(Fraction Acidic, Fraction WY, Entropy) Within (MM,Outlier) Populations')

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
plt_sr, = pax.plot(x+.5,sr/3,color='r')
plt_sb, = pax.plot(x+.5,sb/3,color='k')
l3 = pax.legend([plt_sr,plt_sb],['S(MM)','S(Outlier)'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.xlabel('AA Position')
plt.title('(Fraction Basic, Fraction WY, Entropy) Within (MM,Outlier) Populations')

fig.tight_layout()
plt.savefig('AA+S.pdf')














fig = plt.figure(figsize=(20,20))
windowlen = 10
window = np.ones(windowlen)/windowlen

sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
pr, pb = count_by_site(polar,r_seqs)/rm, count_by_site(polar,b_seqs)/bm
wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pr = np.convolve(pr,window,'same')
pb = np.convolve(pb,window,'same')
pax = fig.add_subplot(311)
plt_rp, = pax.plot(x,-pr,color='r',linewidth=3)
plt_bp, = pax.plot(x+(1-.3)/2,-pb,color='k',linewidth=3)
l1 = pax.legend([plt_rp,plt_bp],['MM Polar/Total','Outlier Polar/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='r')
plt_sb, = pax.plot(x+.5,sb/3,color='k')
l3 = pax.legend([plt_sr,plt_sb],['S(MM)','S(Outlier)'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.title('(Fraction Polar, Fraction WY, Entropy) Within (MM,Outlier) Populations')

# sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
# wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pr, pb = count_by_site('DE',r_seqs)/rm, count_by_site('DE',b_seqs)/bm
pr = np.convolve(pr,window,'same')
pb = np.convolve(pb,window,'same')
pax = fig.add_subplot(312)
plt_rp, = pax.plot(x,-pr,color='r',linewidth=3)
plt_bp, = pax.plot(x+(1-.3)/2,-pb,color='k',linewidth=3)
l1 = pax.legend([plt_rp,plt_bp],['MM DE/Total','Outlier DE/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='r')
plt_sb, = pax.plot(x+.5,sb/3,color='k')
l3 = pax.legend([plt_sr,plt_sb],['S(MM)','S(Outlier)'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.title('(Fraction Acidic, Fraction WY, Entropy) Within (MM,Outlier) Populations')

# sr, sb, srb = entropy(r_seqs), entropy(b_seqs), entropy(r_seqs+b_seqs)
# wyr, wyb = count_by_site('WY',r_seqs)/rm, count_by_site('WY',b_seqs)/bm
pr, pb = count_by_site('KR',r_seqs)/rm, count_by_site('KR',b_seqs)/bm
pr = np.convolve(pr,window,'same')
pb = np.convolve(pb,window,'same')
pax = fig.add_subplot(313)
plt_rp, = pax.plot(x,-pr,color='r',linewidth=3)
plt_bp, = pax.plot(x+(1-.3)/2,-pb,color='k',linewidth=3)
l1 = pax.legend([plt_rp,plt_bp],['MM KR/Total','Outlier KR/Total'],loc=4) #
plt_rwy = pax.bar(x,wyr,1,color='r')
plt_bwy = pax.bar(x+(1-.3)/2,wyb,.3,color='k')
l2 = pax.legend([plt_rwy,plt_bwy],['MM WY/Total','Outlier WY/Total'],loc=1) #####
plt_sr, = pax.plot(x+.5,sr/3,color='r')
plt_sb, = pax.plot(x+.5,sb/3,color='k')
l3 = pax.legend([plt_sr,plt_sb],['S(MM)','S(Outlier)'],loc=2) ###################
plt.gca().add_artist(l1)
plt.gca().add_artist(l2)
plt.xlim([0,n])
plt.ylim([-1,1])
plt.xlabel('AA Position')
plt.title('(Fraction Basic, Fraction WY, Entropy) Within (MM,Outlier) Populations')

plt.tight_layout()
plt.savefig('AA+S_smooth.pdf')


basic_r = count_by_site('KR',r_seqs)/rm
basic_b = count_by_site('KR',b_seqs)/bm
acid_r = count_by_site('DE',r_seqs)/rm
acid_b = count_by_site('DE',b_seqs)/bm
polar_r = count_by_site(polar,r_seqs)/rm
polar_b = count_by_site(polar,b_seqs)/bm
tsr, tsb = entropy(r_seqs,traditional=True), entropy(b_seqs,traditional=True)
# sr,sb
of=open('DHFR.tsv','w')
of.write('pos\tbasic_red\tbasic_black\tacidic_red\tacidic_black\tpolar_red\tpolar_black\tentropy_red\tentropy_black\ttraditional_entropy_red\ttraditional_entropy_black\n')
for i in xrange(n):
    row = [i, basic_r[i], basic_b[i], acid_r[i], acid_b[i], polar_r[i], polar_b[i], sr[i], sb[i]]
    row.extend([tsr[i],tsb[i]])
    row = [str(s) for s in row]
    of.write('\t'.join(row))
    of.write('\n')

