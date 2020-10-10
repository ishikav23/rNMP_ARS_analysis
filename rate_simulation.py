#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# draw plot
def draw(slwt, slpd, slpe, out):
    fig, ax = plt.subplots(figsize=(7,6))
    plt.subplots_adjust(top=1, right=0.95, bottom=0.1, left=0.15)
    ax.plot(np.arange(-100, length_max-100), slpe, color='#FF5555', linewidth=5, alpha=0.7)
    ax.plot(np.arange(-100, length_max-100), slpd, color='#5555FF', linewidth=5, alpha=0.7)
    ax.plot(np.arange(-100, length_max-100), slwt,'k--',linewidth=5, alpha=0.7)
    plt.xlim([-100, length_max - 100])
    plt.ylim([0, 0.0045])
    ax.plot([0,0],[0,1], 'k:', linewidth=3)
    sns.despine()
    plt.savefig(out)

# get value from a random ars
def get_rand_leading(pos, dev, pold, base):
    if pos - dev < 0:
        return base[21]
    elif pos - dev < length_pola:
        return base[1]
    elif pos - dev < pold:
        return base[21]
    else:
        return base[-1]

# get value from a random ars
def get_rand_lagging(pos, dev, base, pole):
    if pos - dev < 0:
        return pole
    elif (pos - dev)%(length_pola + length_pold) < length_pola:
        return base[1]
    else:
        return base[21]

# parameters
length_pola = 20
length_pold = 180
length_max = 1100
rate_pola = 1/625
rate_pold = 1/5000
rate_pole = 1/1250
mrate_pold = 1/5000*10
mrate_pole = 1/1250*5
sns.set(font_scale=2, style='ticks')
np.random.seed(1919)

# single leading
slwt = [rate_pola]*length_pola + [rate_pold]*length_pold\
        + [rate_pole]*(length_max-length_pold-length_pola)
slpd = [rate_pola]*length_pola + [mrate_pold]*length_pold\
        + [rate_pole]*(length_max-length_pold-length_pola)
slpe = [rate_pola]*length_pola + [rate_pold]*length_pold\
        + [mrate_pole]*(length_max-length_pold-length_pola)
draw(slwt, slpd, slpe, 'single_leading.png')

# single lagging
slwt = ([rate_pola]*length_pola + [rate_pold]*length_pold) * 10
slpd = ([rate_pola]*length_pola + [mrate_pold]*length_pold) * 10
draw(slwt[:length_max], slpd[:length_max], slwt[:length_max], 'single_lagging.png')

# combined random
nars = 400
stdev = 100
min_pold = 100
max_pold = 500
devs = np.random.randn(nars)*stdev
pold_lengths = np.random.rand(nars)*(max_pold-min_pold) + min_pold

# combined leading
base = [rate_pola]*length_pola + [rate_pold]*length_pold\
        + [rate_pole]*(length_max-length_pold-length_pola)
slwt = []
for i in range(length_max):
    slwt.append(np.mean([get_rand_leading(i, x[0], x[1], base) for x in zip(devs, pold_lengths)]))
base= [rate_pola]*length_pola + [mrate_pold]*length_pold\
        + [rate_pole]*(length_max-length_pold-length_pola)
slpd = []
for i in range(length_max):
    slpd.append(np.mean([get_rand_leading(i, x[0], x[1], base) for x in zip(devs, pold_lengths)]))
base = [rate_pola]*length_pola + [rate_pold]*length_pold\
        + [mrate_pole]*(length_max-length_pold-length_pola)
slpe = []
for i in range(length_max):
    slpe.append(np.mean([get_rand_leading(i, x[0], x[1], base) for x in zip(devs, pold_lengths)]))
draw(slwt, slpd, slpe, 'combined_leading.png')

# combined lagging
base = ([rate_pola]*length_pola + [rate_pold]*length_pold) * 10
slwt = []
for i in range(length_max):
    slwt.append(np.mean([get_rand_lagging(i, x, base, rate_pole) for x in devs]))
base = ([rate_pola]*length_pola + [rate_pold]*length_pold) * 10
slpe = []
for i in range(length_max):
    slpe.append(np.mean([get_rand_lagging(i, x, base, mrate_pole) for x in devs]))
base = ([rate_pola]*length_pola + [mrate_pold]*length_pold) * 10
slpd = []
for i in range(length_max):
    slpd.append(np.mean([get_rand_lagging(i, x, base, rate_pole) for x in devs]))
draw(slwt, slpd, slpe, 'combined_lagging.png')

