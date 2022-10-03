#!/usr/bin/env python


'''
> MAPQ won't be very useful for filtering poor matches. You should look at [alignment] score, identity and positive. -- https://github.com/lh3/miniprot/issues/7

ms ..

> Please use this tag to estimate mapping uniqueness.

cs:Z::290
ms:i:1498

'''


import argparse
import re
import sys

import pandas as pd


parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--aln', required=True, help='Protein alignment [paf]')
parser.add_argument(
    '--identity', required=True, type=float, default=0.8, metavar='[0-1]', help='Minimum sequence identity [0-1]')
args = parser.parse_args()


assert 0 <= args.identity <= 1, 'Choose an identity between 0 and 1'

try:
    df = pd.read_csv(args.aln, sep='\t', header=None)
    # https://github.com/lh3/miniasm/blob/master/PAF.md
    # df[[13, 17]]
except pd.errors.EmptyDataError:
    sys.stdout.write('0')
    sys.exit()


for i in df.itertuples():
    score = int(i._14.split(':')[-1])
    # TODO: Use this score?
    
    # https://www.drive5.com/usearch/manual/cigar.html
    # http://okko73313.blogspot.com/2012/04/using-regular-expressions-to-analyze.html
    # https://github.com/lh3/miniasm/blob/master/PAF.md
    cigar = i._18.split(':')[-1]
    matches = sum(map(int, re.findall(r'(\d+)M', cigar)))
    ident = round(matches / i._2, 4)
    # PAF column 2: query sequence len

    # print(score, ident)
    if ident >= args.identity:
        sys.stdout.write('1')
        # Any hit is ok
        sys.exit()
    else:
        continue
