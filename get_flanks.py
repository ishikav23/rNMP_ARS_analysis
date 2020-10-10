#!/usr/bin/env python3

import argparse
from collections import defaultdict
import numpy as np
from getFlankUtils import *


def main():
     # argparse
     parser = argparse.ArgumentParser(description='Generate ARS flanks')
     parser.add_argument('ars', type=argparse.FileType('r'), help='Bed file for ars region')
     parser.add_argument('index', type=argparse.FileType('r'), help='index file for background genome')
     parser.add_argument('-l', type=int, default=15000, help='Length of flank region, default=5000')
     parser.add_argument('-b', type=int, default=0, help='Bin size, default = flank length.')
     parser.add_argument('-t', type=float, default=None, nargs='+', help='separator of firing time')
     parser.add_argument('-v', type=int, default=1600, help='Fork speed, base per minute')
     parser.add_argument('-r', action='store_true',  help='Input is ribosomal DNA, only generate the left half.')
     parser.add_argument('-o', default='ars', help='Output file basename')
     args = parser.parse_args()

     if args.b == 0:
         args.b = args.l

     # get ars
     arss, ars_orders = read_ars(args.ars)

     # read chrom size
     chrom_sizes = read_faidx(args.index)

     # calculate ars boundaries
     calc_boundary(arss, ars_orders, chrom_sizes, args.v, args.r, args.l)

     # output with collected desired ARS
     arss_sep = sep_ars(arss, args.t)

     # generate bins
     bins = defaultdict(lambda : defaultdict(lambda : defaultdict(list)))
     for t, v in arss_sep.items():
          for name in v:
               arss[name].add_bins(bins[t], args.b)

     # output
     output_bins(bins, args.o)

print('Done!')


if __name__ == '__main__':
     main()




