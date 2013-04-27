import sys,os
import pandas as pd
import matplotlib.pyplot as plt
from shared import *
from IPython import embed

cats = pd.read_table('groEL_categorization.tsv').dropna()
cats['Group'] = ['Red' if np in rnames else 'Black' for np in cats['gene']]
cat_counts = cats.pivot_table('gene','class','Group',len)
cat_counts.plot(kind='bar', color=['k','r'], figsize=(5,4))
plt.xlabel('Class')

plt.savefig('groEL_cat_hist.pdf')
