#!/usr/bin/env python3

import fileinput
from pandas_plink import read_plink
import pandas as pd
from pyplink import PyPlink
import sys


infile = sys.argv[1]
outfile = sys.argv[2]
chrom = sys.argv[3]

(bim, fam, bed) = read_plink(
        infile, 
        verbose = True
        )

print(bim.head(), flush=True)

bim['snp'] = range(1, 1+len(bim))
#bim['pos'] = 0
bim['chrom'] = chrom

print(bim.head(), flush=True)

bim.drop(
        "i", 
        axis=1, 
        inplace=True
        )

print(bim.head(), flush=True)

filedata = bim.to_string(
        header=False, 
        index=False
        )
print("SNP IDs IN PLACE")

filedata = filedata.replace("0.0","0")


with open(outfile,"w") as file:
        file.write(filedata)


exit()
