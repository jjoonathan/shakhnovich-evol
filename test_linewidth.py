import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt


x = np.arange(200)
bottom_red_bar = -np.random.random(200)
bottom_black_bar = np.random.random(200) * bottom_red_bar

# Tried:
# fig.set_size_inches(17,22)
# fig.savefig('boxes.pdf', dpi=100,150,300)
#   dpi parameter has no effect on linewidth, text size, or the PDF's dimensions
# 'markerscale=.5' on plt.legend and pax.legend
#   no effect
# mp.rcParams['font.size']=8
#   worked! made text smaller, now for linewidth...
# mp.rcParams['lines.linewidth']=5
#   no effect
# fig.set_linewidth(5)
#   no effect
# pax.axhline(linewidth=5)
#   only changes x axis not box surrounding subplot
# fig.set_size_inches(8.5,11) immediately before plt.savefig('boxes.pdf')
#   identical results to calling plt.figure(figsize=(8.5,11)) in the first place
# mp.rcParams['axes.linewidth']=5
#   Made the box surrounding the plot thicker
# fig = plt.figure(linewidth=5)
#   no effect

# I tried passing plt.figure(figsize=(17,22)) and swapping it to 8.5x11 using
# fig.set_size_inches right before saving, but the lines were thick and the text
# was large in the PDF, exactly as if I had set figsize=(8.5,11) to begin with
mp.rcParams['lines.linewidth']=10
mp.rcParams['lines.color']='g'
fig = plt.figure()
for subplotnum in [1,2,3]:
    pax = plt.subplot(310+subplotnum)
    pax.set_ylim([-1,1])
    bot_rb = pax.bar(x,         bottom_red_bar,1,color='r')
    bot_bb = pax.bar(x+(1-.3)/2,bottom_black_bar,.3,color='k')
    # Note: necessary to use pax.legend NOT plt.legend because in
    # the actual figure I want to show more than one legend
    pax.legend([bot_rb,bot_bb],['Filler Text 1','Filler Text 2'],loc=4)
fig.set_size_inches(8.5,11)
fig.savefig('boxes_bad.png',dpi=300)
fig.set_size_inches(17,22)
fig.savefig('boxes_good.png',dpi=150)

fig.set_size_inches(8.5,11)
plt.savefig('boxes.pdf')
