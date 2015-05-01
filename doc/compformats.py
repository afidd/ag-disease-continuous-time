'''
NAADSM and our code predict different probabilities for infection. Why?
We print here a comparison of some internal values, showing that the
distance calculation is responsible. The values are:

infection factor for farm a
infection factor for farm b
distance factor from a to b
probability of infection

Next line:
ax, ay, bx, by, distance
Note that our code is in latlong, and NAADSM is in a projection.
'''

import math as m
import numpy as np
from termcolor import colored, cprint


def read_file(filename):
    columns=list()
    line_cnt=0
    for line in open(filename):
        vals=line.strip().split()
        if len(vals)==12:
            while len(vals)>len(columns):
                columns.append([""]*line_cnt)
            for v, c in zip(vals, columns):
                c.append(v)
            line_cnt+=1

    typed_columns=list()
    for col_kind in columns:
        try:
            intcol=list()
            for entry in col_kind:
                intcol.append(int(entry))
            typed_columns.append(np.array(intcol, dtype=np.int))
        except ValueError:
            try:
                doublecol=list()
                for entry in col_kind:
                    doublecol.append(float(entry))
                typed_columns.append(np.array(doublecol, dtype=np.double))
            except ValueError:
                pass

    return typed_columns

def keyed(cols):
    entries=dict()
    for i in range(len(cols[0])):
        target=cols[0][i]
        source=cols[1][i]
        vals={ "afactor" : cols[2][i], "prevalence" : cols[3][i],
            "distfactor" : cols[4][i], "bfactor" : cols[5][i],
            "probability" : cols[6][i], "ax" : cols[7][i], "ay" : cols[8][i],
            "bx" : cols[9][i], "by" : cols[10][i], "distance" : cols[11][i]}
        if not target in entries:
            entries[target]=dict()
        entries[target][source]=vals
    return entries



def showdiff(target, source, a, b):
    print("{0}--{1}".format(target, source))
    lines=list([list(), list()])
    for v in ["afactor", "bfactor", "distfactor", "probability"]:
        aprint, bprint=errdisp(a[v], b[v], 'red')
        lines[0].append(aprint)
        lines[1].append(bprint)
    for v in ["ax", "ay", "bx", "by", "distance"]:
        aprint, bprint=errdisp(a[v], b[v], 'blue')
        lines[0].append(aprint)
        lines[1].append(bprint)
    for which in [0, 1]:
        print("  ", end="")
        for idx0 in range(4):
            print(lines[which][idx0], end="")
            print("\t", end="")
        print("")
    for which in [0, 1]:
        print("  ", end="")
        for idx0 in range(4, 9):
            print(lines[which][idx0], end="")
            print("\t", end="")
        print("")



def compare_keyed(a, b):
    for target, sources in a.items():
        if target in b:
            for source, avalues in a[target].items():
                if source in b[target]:
                    showdiff(target, source, avalues, b[target][source])


def errdisp(a, b, color):
    tolerance=1e-4
    if abs((a-b)/a)>tolerance:
        return [colored("{0:8.4f}".format(a), color),
            colored("{0:8.4f}".format(b), color)]
    else:
        return ["{0:8.4f}".format(a), "{0:8.4f}".format(b)]




def compare():
    col1=keyed(read_file("adct_out.txt"))
    col2=keyed(read_file("test_output2.txt"))
    compare_keyed(col1, col2)


compare()
