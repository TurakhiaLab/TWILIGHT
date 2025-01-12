from argparse import ArgumentParser
import numpy as np
import time

parser = ArgumentParser(description='Reduce alignment length to speedup tree inference process')
parser.add_argument('inaln', help='Input alignment')
parser.add_argument('outaln', help='Output alignment')
parser.add_argument('threshold', type=float, help='Minimum gap porpotion for a column be removed')
args = parser.parse_args()

st = time.time()

threshold = args.threshold
name = []
aln = []

with open(args.inaln, "r") as alnFile:
    inContent = alnFile.read().splitlines()
    for c in inContent:
        if c[0] == '>':
            name.append(c)
        else:
            aln.append(c)

allAln = np.array([list(a) for a in aln])
lb = len(allAln[0])
allAln = np.transpose(allAln)
stayedRows = []
rowID = 0
for row in allAln:
    num_gap = (row == '-').sum()
    if num_gap/len(name) <= threshold:
        stayedRows.append(rowID)
    rowID += 1
newAln = []
for r in stayedRows:
    newAln.append(allAln[r])
newAln = np.array(newAln)
newAln = np.transpose(newAln)
la = len(newAln[0])
outFile = []
with open(args.outaln, "w") as outFile:
    n = 0
    for a in newAln:
        outFile.write(name[n]+'\n')
        n += 1
        outFile.write("".join(a)+'\n')
en = time.time()

print("Masked gappy site. Length before/after: "+str(lb)+"/"+str(la)+". Total time: ", en-st, "seconds.")
