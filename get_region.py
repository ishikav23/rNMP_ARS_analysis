#!/usr/bin/env python3
import argparse
import sys
from collections import OrderedDict

def main():
    parser = argparse.ArgumentParser(description='Get a paticular range from an ARS info file')
    parser.add_argument('info', type=argparse.FileType('r'), help='ARS info file')
    parser.add_argument('-s', type=int, default=0, help='Start postion, exclude. (0)')
    parser.add_argument('-e', type=int, default=2**32, help='End position, include. (2**32)')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('--col_num', type=int, default=5, help='Column number for the postion, start with 0. (5)')
    args = parser.parse_args()

    # header
    header = args.info.readline().rstrip('\n').split('\t')
    args.o.write('chrom\t' + '\t'.join(header[args.col_num + 3:]) + '\n')

    # get data
    data = OrderedDict()
    for l in args.info:
        ws = l.rstrip('\n').split('\t')
        # skip unselected entries
        pos = float(ws[args.col_num])
        if pos <= args.s or pos > args.e:
            continue
        # sum up all selected entries
        name = '-'.join(ws[:args.col_num] + [ws[args.col_num+1]])
        if name not in data:
            data[name] = [float(x) for x in ws[args.col_num + 3:]]
        else:
            for i in range(len(ws[args.col_num+3:])):
                data[name][i] += float(ws[args.col_num+3+i])

    # output
    for k,v in data.items():
        args.o.write('\t'.join([k] + [str(x) for x in v]) + '\n')

    print('Done!')





if __name__ == '__main__':
    main()
