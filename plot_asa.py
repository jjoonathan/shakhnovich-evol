import scipy.stats
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython import embed
from shared import *

def fill_axes(ax):
	ax.plot(red_asa, lw=2, label='Michaelis-Menton', color='red')
	ax.plot(blk_asa, lw=2, label='Outlier', color='black')
	ax.fill_between(np.arange(n), red_asa+red_asa_std, red_asa-red_asa_std, facecolor='red', alpha=.3)
	ax.fill_between(np.arange(n), blk_asa+blk_asa_std, blk_asa-blk_asa_std, facecolor='black', alpha=.3)
	ax.set_title('Average +/- StdDev of ASA Within MM, Outlier Protein Families')
	ax.set_ylabel('ASA (Angstrom^2)')
	ax.legend(loc=1)
	ax.set_xlim([0,201])
	ax.set_xlabel('Residue #')

if __name__ == '__main__':
	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(111)
	fill_axes(ax)
	plt.savefig('ASA.pdf')
