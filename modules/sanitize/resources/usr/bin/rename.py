#!/usr/bin/env python


'''
Rename contigs to alleviate the pain in downstream workflow processes. We
simultaneously unzip genomes.
'''


import argparse
import hashlib
import json
from pathlib import Path
import warnings

import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--genome', required=True, help='Genome')
parser.add_argument(
    '--name', required=True, help='Name of file as specified by user')
parser.add_argument(
    '--strict', action='store_true', help='Fail if duplicate contigs detected (omitted by default)')
args = parser.parse_args()


name = args.name
p = Path(args.genome)

log = {}
log['name'] = name
log['path'] = p.resolve().__str__()
log['filename'] = p.name


d, e = {}, {}
# screed can open zipped files
with screed.open(args.genome) as file:
    for record in file:

        seq = str(record.sequence)
        hsh = hashlib.md5(seq.encode('utf-8')).hexdigest()
        
        if hsh in d and not args.strict:
            warnings.warn('Duplicate contig detected, will omit.')
            continue

        # Err out if strict mode
        assert hsh not in d, 'Duplicate contig detected, exit.'

        d[hsh] = seq
        e[hsh] = str(record.name)


with open(name + '.fna', 'w+') as out:
    for hsh, seq in d.items():
        out.write(f'>{name}__{hsh}\n{seq}\n')

# Example of a silly contig header
# NZ_UIRU01000009.1 Klebsiella pneumoniae strain EuSCAPE_GR046, whole genome shotgun sequence
# -- contains space, comma, underscore, so we will split on "::"

log['map'] = e

with open(name + '.json', 'w+') as out:
    json.dump(log, out, indent=4)


