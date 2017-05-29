#!/usr/bin/env python

import sys

from pyrna.parsers import parse_fasta, to_fasta
from pyrna.db import NCBI

"""
This script parse a FASTA file and search for orgnism names
"""

def convert(fasta_file):
    with open(fasta_file) as f:
        molecules = parse_fasta(f.read())
        ids = []
        for m in molecules:
            ids.append(m.name.split('/')[0])
        ncbi = NCBI()
        response = ncbi.efetch(db='nucleotide', ids=ids, rettype='gb')
        #print response.split('\n')
        organisms = filter(lambda line: line.strip().startswith('ORGANISM'), response.split('\n'))
        print len(ids), len(organisms)
        for i in xrange(0, len(molecules)):
            molecules[i].name =  ' '.join(organisms[i].strip().split(' ')[1:]).strip()
        with open("%s_with_organism_names.fasta"%fasta_file.split('.fasta')[0], 'w') as f2:
            f2.write(to_fasta(molecules))

if __name__ == '__main__':
    file = None

    if "-f" in sys.argv:
        file = sys.argv[sys.argv.index("-f")+1]

    if not file:
        print "Usage: from_genbank_ids_2_organism_names -f fasta_file"
        sys.exit(-1)

    convert(file)
