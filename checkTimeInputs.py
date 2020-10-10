from collections import defaultdict
import pandas as pd

# read lib info
def read_libinfo(fr):
    libinfo = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws)!= 4:
            continue
        libinfo[ws[0]] = tuple(ws[1:])
    return libinfo


# read ars
def read_ars(fr):
    ars = defaultdict(list)
    for l in fr:
        ws = l.split('\t')
        if len(ws) < 4:
            continue
        ars[ws[0]].append((ws[0], int(ws[1]), int(ws[2]), float(ws[4])))
    return ars


# generate small windows
def generate_windows(ars, l, block_ribosomal):
    windows = []
    for chrom, v in ars.items():
        ep = 0
        for i in range(len(v)):
            # postions
            _, sc, ec, t = v[i]
            # first half
            if sc > ep:
                s = max(sc - l, ep)
                windows.append((chrom, s, sc, t, 'leading','-'))
            # second half
            try:
                if ec < v[i+1][1]:
                    e = min(ec + l, int((v[i+1][1] + ec)/2))
                    windows.append((chrom, ec, e, t, 'leading','+'))
            except IndexError:
                e = ec+l
                windows.append((chrom, ec, ec+l, t, 'leading','+'))
            # update ep
            ep = e
    windows.sort(key=lambda x:(x[0], x[1]) )
    if block_ribosomal:
        windows = remove_ribosomal(windows)
    return windows


# remove ribosomal region
def remove_ribosomal(windows):
    # ribosomal information
    chrom = 'chrXII'
    s = 451576
    e = 467570
    # block
    win = []
    for l in windows:
        if s >= l[2] or e <= l[1]:
            win.append(l)
        elif l[1] < s < l[2] <= e:
            win.append((l[0], l[1], s, l[3], l[4], l[5]))
        elif s <= l[1] <= e < l[2]:
            win.append((l[0], e, l[2], l[3], l[4], l[5]))
        elif l[1] < s < e < l[2]:
            win.append((l[0], l[1], s, l[3], l[4], l[5]))
            win.append((l[0], e, l[2], l[3], l[4], l[5]))
    return win

# read data from bed file
def read_data(ars, libs, folder):
    # initialization
    data = {}
    for a in ars:
        data[a] = {'leading':defaultdict(int), 'lagging':defaultdict(int)}
    # read bed pair
    for lib in libs:
        with open(folder + '/{}.bed'.format(lib)) as fr:
            for l in fr:
                ws = l.rstrip('\n').split('\t')
                if len(ws) != 6:
                    continue
                ws[1] = int(ws[1])
                ws[2] = int(ws[2])
                pos = find_pos(ws, ars, 0, len(ars)-1)
                if pos:
                    if ws[5] == pos[5]:
                        data[pos]['leading'][lib] += 1
                    else :
                        data[pos]['lagging'][lib] += 1
    return data


# binary search for get position
def find_pos(ws, poss, s, e):
    # found
    if s == e:
        pos = poss[s]
        if comp_pos(ws, pos) == 0:
            return pos
        else:
            return None
    # get mid
    mid = int((s+e)/2)
    status = comp_pos(ws, poss[mid])
    if status == 0:
        return poss[mid]
    elif status == None:
        return None
    elif status == 1:
        return find_pos(ws, poss, mid+1, e)
    elif status == -1:
        if mid == s:
            return None
        else:
            return find_pos(ws, poss, s, mid-1)


# comparison for binary search
# ws > pos :1, ws < pos : 2
# ws in pos: 0
# else None
def comp_pos(ws, pos):
    if ws[0] > pos[0]:
        return 1
    elif ws[0] < pos[0]:
        return -1
    else:
        if ws[1] >= pos[2]:
            return 1
        elif ws[2] <= pos[1]:
            return -1
        elif ws[1] >= pos[1] and ws[2] <= pos[2]:
            return 0
        else:
            return None


# convert to dataframe
def generate_df(data, libinfo):
    d = []
    for fs, info in libinfo.items():
        for window, count in data.items():
            d.append([fs] + list(info) + list(window[:4]) + [window[5], count['leading'][fs], count['lagging'][fs]])
    df = pd.DataFrame(d, columns=['Library', 'String','Genotype','RESet',\
            'Window_chr','Window_start','Window_end','Firing_time','Leading_pos',\
            'Leading','Lagging'])
    return df

