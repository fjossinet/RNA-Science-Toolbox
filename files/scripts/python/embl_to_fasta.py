#!/usr/bin/env python

import os, sys
from pyrna import parsers

def convert(working_dir):

    genomic_sequences = []

    for f in os.listdir(working_dir):
        if f.endswith('embl'):
            with open(os.path.join(working_dir,f)) as h:
                genomic_sequence, features =  parsers.parse_embl(h.read())
                genomic_sequences.append(genomic_sequence)

    with open(os.path.join(working_dir,'scaffolds.fasta'), 'w') as h:
        h.write(parsers.to_fasta(genomic_sequences))

if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        print "Usage: embl_to_fasta.py dir_with_embl_files"
        sys.exit(-1)
    
    convert(sys.argv[1])