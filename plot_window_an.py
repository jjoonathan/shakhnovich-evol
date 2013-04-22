from shared import *
from matplotlib import pyplot as plt

fig = plt.figure(figsize=(20,20))

sr = entropy(ras)
sb = entropy(bas)
polar_sites = count_by_site(polar,ras)/mr > .25
acidic_sites = count_by_site(acidic,ras)/mr > .25
basic_sites = count_by_site(basic,ras)/mr > .25
wy_sites = count_by_site('WY',ras)/mr > .25
x = np.arange(len(ras[0]))
polar_x = x[polar_sites]
acidic_x = x[acidic_sites]
basic_x = x[basic_sites]
wy_x = x[wy_sites]

def make_subplot_for_windowlen(wlen):
	ax = plt.gca()
	w = np.ones(wlen)/wlen
	sr_smooth = np.convolve(sr,w,'same')
	sb_smooth = np.convolve(sb,w,'same')
	ar_smooth = np.convolve(red_nasa,w,'same')
	ab_smooth = np.convolve(blk_nasa,w,'same')
	plt_sr, = ax.plot(x, sr_smooth/4, 'r')
	plt_sb, = ax.plot(x, sb_smooth/4, 'k')
	plt_ar, = ax.plot(x, -ar_smooth/.5, 'r')
	plt_ab, = ax.plot(x, -ab_smooth/.5, 'k')
	c0, = ax.plot(polar_x, np.zeros(len(polar_x)), 'mo')
	c1, = ax.plot(acidic_x, np.zeros(len(acidic_x)), 'ro')
	c2, = ax.plot(basic_x, np.zeros(len(basic_x)), 'bo')
	c3, = ax.plot(wy_x, np.zeros(len(wy_x)), 'go')
	plt.xlim([0, len(ras[0])])
	plt.ylim([-1,1])
	plt.title('%i Residue Moving Average of S, ASA'%wlen)
	plt.xlabel('Residue #')
	top_legend = ax.legend([plt_sr, plt_sb], ['S (bits/4)', 'S (bits/4)'], loc=2)
	bot_legend = ax.legend([plt_ar, plt_ab], ['ASA (HOA/2)', 'ASA (HOA/2)'], loc=3)
	dot_legend = ax.legend([c0,c1,c2,c3], ['Polar','Acidic','Basic','WY'], loc=1)
	ax.add_artist(top_legend)
	ax.add_artist(bot_legend)
	ax.add_artist(dot_legend)

fig.add_subplot(311)
make_subplot_for_windowlen(5)
fig.add_subplot(312)
make_subplot_for_windowlen(10)
fig.add_subplot(313)
make_subplot_for_windowlen(20)

fig.tight_layout()
plt.savefig('window_an.pdf')
