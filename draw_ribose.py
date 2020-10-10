#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from collections import defaultdict

def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='y',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
    ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw))
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.10, 0.2),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)

def main():
    # argparse
    parser = argparse.ArgumentParser(description='PCA for dinucleotide data')
    parser_io = parser.add_argument_group('I/O')
    parser_io.add_argument('DATA',type=argparse.FileType('r'), default=sys.stdin, help='Normalized dinucleotide data')
    parser_io.add_argument('-o', default='', help='output basename')
    parser_io.add_argument('-b', type=argparse.FileType('r'), help='Select background file. If a file is selected, the background percentage is added to labels.')
    parser_io.add_argument('--nr', action='store_true', help='Input file is NR type dinucleotide frequency. Of which the second position is rNMP.')
    parser_io.add_argument('--mono', action='store_true', help='Input file is mononucleotide frequency.')
    parser_io.add_argument('--tri', type=int, default=0, help='Input file is trinucleotide frequency, with rNMP incorporated at TRI position')
    parser_io.add_argument('--background_chrom', default='chrM', help='Chromosome name of background file, default = chrM')
    parser_graphs = parser.add_argument_group('Graphs')
    parser_h = parser.add_argument_group('Heatmap')
    parser_h.add_argument('--legend_group', default=0, choices=[0,4,16], type=int, help='Number of lables of which the sum is 1. If 0 is selected, the sum of all labels will be 1. default = 0.')
    parser_h.add_argument('--no_annot', action='store_true', help='Hide percentage annotation in each cell')
    parser_h.add_argument('--cmax', type=float, default=0.5, help='Maximum value in color scale. Any preferency beyond that will show as the maximum color.')
    args = parser.parse_args()

    # argument relations
    if sum([args.nr, args.mono, args.tri!=0]) > 1:
        print('Arguments conflict, can only be mono, di or trinuc')
        sys.exit(1)

    if args.mono:
        args.legend_sum1 = True


    # build ndarray
    sample_raw = []
    da = []
    di = args.DATA.readline().rstrip('\n').split('\t')[1:]
    for l in args.DATA:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < len(di):
            break
        try:
            da.append([float(i) for i in ws[1:]])
        except ValueError:
            print(ws)
            sys.exit('Cannot convert to float')
        sample_raw.append(ws[0])
    mat = np.asarray(da)
    sample = [i.split('/')[-1] for i in sample_raw]

    # set color settings for graph
    sns.set(palette='muted', style='white')

    # heatmap
    if args.mono:
        fig, ax = plt.subplots(figsize=(mat.shape[0]*0.45+6,5))
        plt.subplots_adjust(left=0.02, right=1.24, bottom=0.3, top=1)
        sns.heatmap(mat.T, vmin=0, vmax=args.cmax,center=args.cmax*0.45, ax=ax, annot=(not args.no_annot), annot_kws={"size":12})
    elif args.tri:
        fig, ax = plt.subplots(figsize=(mat.shape[0]*0.45+8,16))
        plt.subplots_adjust(left=0.03, right=1, bottom=0.15, top=0.98)
        sns.heatmap(mat.T, vmin=0, vmax=args.cmax,center=args.cmax*0.45, ax=ax, yticklabels=1, annot=(not args.no_annot), annot_kws={"size":12})
    else:
        fig, ax = plt.subplots(figsize=(mat.shape[0]*0.45+6,10))
        plt.subplots_adjust(left=0.03, right=1.24, bottom=0.15, top=1)
        sns.heatmap(mat.T, vmin=0, vmax=args.cmax,center=args.cmax*0.45, ax=ax, annot=(not args.no_annot), annot_kws={"size":12})

    # set labels
    # get percentage from background
    labels=[]
    f = defaultdict(float)
    if args.b:
        for l in args.b:
            ws = l.split('\t')
            if ws[0] == args.background_chrom:
                ws = ws[1:]
                for i in range(len(ws)):
                    f[di[i]] += float(ws[i])
    else:
        for d in di:
            f[d] = 0
    for k,v in f.items():
        labels.append([k,v])

    # sort labels
    if args.mono:
        labels = sorted(labels, key=lambda x:x[0])
        labels[-1][0]='U'

    elif args.tri:
        base_order = [[0,1,2], [1,0,2],[2,1,0]]
        labels = sorted(labels, key=lambda x:[x[0][i] for i in base_order[args.tri-1]])
        # change T to U
        for i in range(len(labels)-16, len(labels)):
            labels[i][0] = labels[i][0][0:args.tri-1] + 'U' + labels[i][0][args.tri:]

    elif args.nr:
        labels = sorted(labels, key=lambda x:(x[0][1],x[0][0]))
        labels[-4][0]='AU'
        labels[-3][0]='CU'
        labels[-2][0]='GU'
        labels[-1][0]='TU'
    else:
        labels = sorted(labels, key=lambda x:(x[0][0],x[0][1]))
        labels[-4][0]='UA'
        labels[-3][0]='UC'
        labels[-2][0]='UG'
        labels[-1][0]='UT'

    if args.b:
        # sum to 1
        if args.legend_group == 0:
            group_len = len(labels)
        else:
            group_len = args.legend_group
        s = []
        for i in range(int(len(labels)/group_len)):
            s += [sum([x[1] for x in labels[i*group_len:i*group_len+group_len]])]*group_len
        for j in range(len(labels)):
            labels[j][1] /= (s[j]/100)

    # merge labels
    label_texts = []
    for i in labels:
        if args.b:
            label_texts.append('{:.2f}% {}'.format(i[1], i[0]))
        else:
            label_texts.append(i[0])


    ax.set_xticklabels(sample, rotation='vertical')
    ax.set_yticklabels(label_texts,rotation='horizontal')

    # change top of colorbar to '0.5-1'
    cax = plt.gcf().axes[-1]
    color_labels = cax.get_ymajorticklabels()
    color_labels_texts = [i.get_text() for i in color_labels]
    color_labels_texts[-1] += ' - 1'
    cax.set_yticklabels(color_labels_texts)

    # add percentage labels
    if args.b:
        ax.set_yticklabels(label_texts)

    # formating
    # correct x ticks
    xticks = [tick.label.get_text().split('-')[0] for tick in ax.xaxis.get_major_ticks()]
    ax.set_xticklabels(xticks, fontsize=26)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.yticks(fontsize=28)
    # show or save
    if args.o == '':
        plt.show()
    else:
        ax.get_figure().savefig(args.o + '_heatmap.png')
        print('Heatmap is saved to {}_heatmap.png'.format(args.o))

if __name__ == '__main__':
    main()

