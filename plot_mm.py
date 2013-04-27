import sys,os
from shared import *
import matplotlib.pyplot as plt
import scipy.optimize
from IPython import embed

# Growth Data Table
growth_tbl = pd.read_csv('growth_rates.csv').set_index('NP')

# Scrubbing
grate = 'Growth rate 37C'
grate_dlon = 'Growth rate 37C.1'
grate_grOE = 'Growth rate 37C.2'  # Growth rate with GroEL overexpressed
kckm = 'kcat/KM (sec-1mM-1)'
charge = 'charge at pH7'
growth_tbl['dlon'] = growth_tbl[grate_dlon]-growth_tbl[grate]
r_growth_tbl = growth_tbl[growth_tbl['IsRed']==1].dropna(subset=[grate])
b_growth_tbl = growth_tbl[growth_tbl['IsRed']==0].dropna(subset=[grate])

# Fitting
kckms = r_growth_tbl['kcat/KM (sec-1mM-1)'].values
kcs = r_growth_tbl['kcat (sec-1)'].values
# Hmm, kckms and kcs seem to be highly correlated:
# >>> np.corrcoef(kcs,kckms)
# array([[ 1.        ,  0.92266796],
#        [ 0.92266796,  1.        ]])
grates = r_growth_tbl[grate].values
def mm(x, c1, c2):
    #kckm = x[0,:]
    #kc = x[1,:]
    kckm = x
    kc = 1
    ret = c1*kckm/(1+c2*kckm/kc)
    return ret
kckm_kc = np.array([kckms,kcs])
rparam,rcov = scipy.optimize.curve_fit(mm, kckms, grates, [.24,5])
# rparam = array([ 0.0730465 ,  3.15967153])
# c1 = [E][S] = .073
# c2 = [S]/Kcat = 3.16

# Plotting
plt.figure(figsize=(5,4))
# Michaelis-Menton curve
pltx = np.linspace(1,75,100)
p0, = plt.plot(pltx, mm(pltx,*rparam), 'r-')
p1, = plt.plot(r_growth_tbl[kckm], r_growth_tbl[grate], 'r.')
p2, = plt.plot(b_growth_tbl[kckm], b_growth_tbl[grate], 'k.')
plt.xlim(0,80)
plt.ylim(.006,.026)
plt.ylabel(r'Growth Rate @ 37$^\circ$C')
plt.xlabel('$K_{cat}/K_m (sec^{-1}mM^{-1})$')
plt.legend([p0,p1,p2],['MM Fit Line','MM Points','Outlier Points'],'lower right')
plt.tight_layout()
plt.savefig('mm.pdf')
