#!/usr/bin/env python


'''
> MAPQ won't be very useful for filtering poor matches. You should look at [alignment] score, identity and positive. -- https://github.com/lh3/miniprot/issues/7

ms ..

> Please use this tag to estimate mapping uniqueness.

cs:Z::290
ms:i:1498

'''


import argparse
from collections import defaultdict
from glob import glob
from gzip import GzipFile
from io import StringIO
import os
import re
import sys
from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
import pandas as pd
from pyfaidx import Fasta
from screed import rc


def parse_paf(aln, min_ident=.8, min_cov=.8):
    '''
    Returns the (single!) best hit from a miniprot alignment if it passes all
    thresholds.

    Spec: 

    https://github.com/lh3/miniasm/blob/master/PAF.md

    Example:

    score, q, contig, start, end, strand = parse_paf(
        aln, min_ident=.8, min_cov=.8)
    '''
    try:
        df = pd.read_csv(aln, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        return None

    results = []
    for i in df.itertuples():
        # Parse
        contig = i._6
        qry = i._1
        score = int(i._14.split(':')[-1])
        # https://www.drive5.com/usearch/manual/cigar.html
        # http://okko73313.blogspot.com/2012/04/using-regular-expressions-to-analyze.html
        # https://github.com/lh3/miniasm/blob/master/PAF.md
        cigar = i._18.split(':')[-1]
        query_len = i._2
        target_len = i._7
        
        # matches = sum(map(int, re.findall(r'(\d+)M', cigar)))
        '''
        From the alignment of top hit sequences, I get the impression that a
        "match" does not mean a literal residue match, but maybe just from the
        same group, think Dayhoff. Not sure though I understand the format
        correctly.

        Yup, there is another field "residue matches" which is different from
        the Ms one can count in the cigar string. Also, need to divide by 3 to
        get protein len.
        '''
        matches = i._10 / 3
        strand = i._5
        
        target_start, target_end = i._8, i._9
        if target_start > target_end:
            target_start, target_end = target_end, target_start 
        delta = target_end - target_start
        try:
            assert delta % 3 == 0      # can be translated into protein seq
            delta_aa = int(delta / 3)  # protein len of match on target
        except AssertionError:
            continue
        
        # Filter

        # Identity relative to query
        ident = round(matches / query_len, 4)
        if ident < min_ident:
            continue

        # cov_query = delta / (query_len * 3)
        # cov_target = delta / (target_len * 3)
        # 3x bc/ the query is a protein

        # Reciprocal coverage
        cov_query = delta_aa / query_len
        cov_target = query_len / delta_aa

        if (cov_query >= min_cov) and (cov_target >= min_cov):
            # We'll later identify the best hit using the score.
            results.append([
                score, qry, contig, target_start, target_end, strand
            ])
            # print(qry, query_len, contig, target_start, target_end, strand, cigar, matches, ident, delta_aa, cov_query, cov_target, i._10/3, i._11/3)
        else:
            continue

    # Return top hit
    try:
        return sorted(results, key=lambda x: x[0], reverse=True)[0]
    except IndexError:
        # No valid results, all items filtered out
        return None


def get_sequence(fp, contig, start, end, strand):
    '''
    Extract hit sequence, samtools faidx-style; expects compressed genome.fna.gz
    '''
    with open(fp, 'rb') as file, NamedTemporaryFile() as tmp:
        f = GzipFile(fileobj=file, mode='rb')
        tmp.write(f.read())
        genome = Fasta(tmp.name)
        # print(genome.keys())
        seq = genome[contig][start:end]
        fna = seq.seq

        if strand == '-':
            fna = rc(fna)

        # assert len(seq) % 3 == 0
        residues = Seq(fna).translate().__str__()

    return seq.fancy_name, residues


parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--genome', required=True, help='Target genome [.fna.gz]')
parser.add_argument(
    '--aln', required=True, help='Protein alignment [.paf]')
parser.add_argument(
    '--identity', type=float, default=0.8, metavar='[0-1]', help='Minimum sequence identity [0-1]')
parser.add_argument(
    '--coverage', type=float, default=0.8, metavar='[0-1]', help='Minimum coverage [0-1]')
parser.add_argument(
    '--out', required=True, type=str, help='Write top hit here')
args = parser.parse_args()


try:
    *rest, contig, start, end, strand = parse_paf(
            args.aln, min_ident=args.identity, min_cov=args.coverage)
except TypeError:
    # No valid hits
    sys.exit()

header, res = get_sequence(args.genome, contig, start, end, strand)
with open(args.out, 'w+') as out:
    out.write(f'>{header}\n{res}\n')
sys.exit()
