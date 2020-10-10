import numpy as np
from collections import defaultdict

# read ars


def read_ars(fr):
    ars_orders = defaultdict(list)
    arss = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        arss[ws[3]] = ARS(ws[3], ws[0], float(ws[4]), int(ws[1]), int(ws[2]))
        ars_orders[ws[0]].append(ws[3])
    return arss, ars_orders


# read chrom size
def read_faidx(fr):
    chrom_sizes = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        chrom_sizes[ws[0]] = int(ws[1])
    return chrom_sizes


# calc boundary
def calc_boundary(arss, ars_orders, chrom_sizes, speed, ribosomal, max_len=2**32):
    # only keep left part for ribosomal ARS
    if ribosomal:
        for v in arss.values():
            v.left_bound = v.pos - max_len
    else:
        for chrom, arsl in ars_orders.items():
            # first ars
            first = arss[arsl[0]]
            first.left_end_point = 0
            first.left_boundary = max(0, first.pos - max_len)
            # last ars
            last = arss[arsl[-1]]
            last.right_end_point = chrom_sizes[chrom]
            last.right_boundary = min(chrom_sizes[chrom], last.pos + max_len)
            # other
            for i in range(len(arsl) - 1):
                split_interval(arss[arsl[i]], arss[arsl[i+1]], speed, max_len)


# split an interval between two arss
def split_interval(ars1, ars2, speed, max_len):
    if ars1.chrom != ars2.chrom:
        raise ValueError(f'{ars1.name} and {ars2.name} are not in same chrom!')
    # calc split point
    dist = ars2.pos - ars1.pos
    ep = int(
        np.ceil((dist + (ars2.firing_time - ars1.firing_time) * speed) / 2)) + ars1.pos
    if ep < ars1.pos:
        ep = ars1.pos
    elif ep > ars2.pos:
        ep = ars2.pos
    # get boundary
    ars1.right_end_point = ep
    ars2.left_end_point = ep
    ars1.right_boundary = min(ep, ars1.pos + max_len)
    ars2.left_boundary = max(ep, ars2.pos - max_len)

# separate ars


def sep_ars(arss, ts):
    if ts == None:
        ts = [0]
    else:
        ts = [0] + ts
    arss_sep = {}
    for i in ts:
        arss_sep[i] = []
    cate = 0
    for ars in arss.values():
        for t in ts:
            if ars.firing_time > t:
                cate = t
            else:
                break
        arss_sep[cate].append(ars.name)
    return arss_sep


# output bins to file
def output_bins(bins, basename):
    for t, v in bins.items():
        for s, v1 in v.items():
            for pos, v2 in v1.items():
                with open(f'{basename}_{t}_{pos}_{s}.bed', 'w') as fw:
                    for l in v2:
                        fw.write('\t'.join([str(x) for x in l]) + '\n')


# output timed bins to file
def output_timed_bins(bins, basename):
    for t, v in bins.items():
        for s, v1 in v.items():
            with open(f'{basename}_{t}_{s}.bed', 'w') as fw:
                for l in v1:
                    fw.write('\t'.join([str(x) for x in l]) + '\n')


class ARS(object):
    def __init__(self, name, chrom, t, left, right):
        self.name = name
        self.chrom = chrom
        self.firing_time = t
        self.left = left
        self.right = right
        self.pos = int((left + right)/2) + 1
        # need to implement later
        self.left_end_point = None
        self.right_end_point = None
        self.left_boundary = None
        self.right_boundary = None

    # split into bins
    def add_bins(self, bins, binsize):
        # left
        nbins = int(np.ceil((self.pos - self.left_boundary)/binsize))
        for i in range(nbins):
            l = (i+1)*binsize
            e = self.pos - i * binsize
            s = max(self.left_boundary, e - binsize)
            bins['leading'][l].append([self.chrom, s, e, self.name, -l, '-'])
            bins['lagging'][l].append([self.chrom, s, e, self.name, -l, '+'])
        # right
        nbins = int(np.ceil((self.right_boundary - self.pos)/binsize))
        for i in range(nbins):
            l = (i+1)*binsize
            s = self.pos + i * binsize
            e = min(self.right_boundary, s + binsize)
            bins['leading'][l].append([self.chrom, s, e, self.name, l, '+'])
            bins['lagging'][l].append([self.chrom, s, e, self.name, l, '-'])

    # split into bins according to replication time
    def add_bins_time(self, bins, min_time, max_time, binsize, speed):
        left_boundary_time = (self.pos - self.left_boundary) / \
            speed + self.firing_time
        right_boundary_time = (self.right_boundary -
                               self.pos)/speed + self.firing_time
        # left
        binning_times = self.generate_binning_time(self.firing_time, left_boundary_time, min_time, max_time, binsize)
        for t in binning_times[::-1]:
            s = max(self.left_boundary, int(self.pos - (t - self.firing_time) * speed))
            e = min(self.pos, int(self.pos - (t - self.firing_time - binsize) * speed))
            bins['leading'][t].append([self.chrom, s, e, self.name, self.firing_time, '-'])
            bins['lagging'][t].append([self.chrom, s, e, self.name, self.firing_time, '+'])
        # right
        binning_times = self.generate_binning_time(self.firing_time, right_boundary_time, min_time, max_time, binsize)
        for t in binning_times:
            s = max(self.pos, int(self.pos + (t - self.firing_time - binsize) * speed))
            e = min(self.right_boundary, int(self.pos + (t - self.firing_time) * speed))
            bins['leading'][t].append([self.chrom, s, e, self.name, self.firing_time, '+'])
            bins['lagging'][t].append([self.chrom, s, e, self.name, self.firing_time, '-'])  


    # generate binning time
    def generate_binning_time(self, firing_time, end_time, min_time, max_time, binsize):
        binning_times = list(np.arange(min_time, max_time+binsize, binsize))
        return [x for x in binning_times if x > firing_time and x - binsize < end_time]
