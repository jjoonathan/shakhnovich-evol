import scipy.stats
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython import embed

asas=pd.read_table('ASAs.tsv')
red=pd.DataFrame(asas[[k for k in asas if k.startswith('RED')]])
blk=pd.DataFrame(asas[[k for k in asas if k.startswith('BLK')]])
redm = red.mean(axis=1, skipna=True).fillna(method='pad').values
redstd = red.std(axis=1, skipna=True).fillna(method='pad').values
blkm = blk.mean(axis=1, skipna=True).fillna(method='pad').values
blkstd = blk.std(axis=1, skipna=True).fillna(method='pad').values
n = len(red['RED_NP_217278'])  # Convention: m=#rows, n=#cols

def fill_axes(ax):
	ax.plot(redm, lw=2, label='Outliers', color='red')
	ax.plot(blkm, lw=2, label='Michaelis-Menton', color='black')
	ax.fill_between(np.arange(n), redm+redstd, redm-redstd, facecolor='red', alpha=.3)
	ax.fill_between(np.arange(n), blkm+blkstd, blkm-blkstd, facecolor='black', alpha=.3)
	ax.set_title('ASA vs Residue')
	ax.set_ylabel('ASA (Angstrom^2)')
	ax.set_xlabel('Residue #')

if __name__ == '__main__':
	fig,ax = plt.subplots(1)
	fill_axes(ax)
	plt.show()