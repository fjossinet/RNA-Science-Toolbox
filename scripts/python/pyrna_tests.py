#!/usr/bin/env python

"""
A script to test the installation of PyRNA
"""

from pyrna.db import PDB
from pyrna.parsers import parse_pdb, secondary_structure_to_base_pairs
from pyrna.computations import Rnafold, Rnaview

def test():
    print "Recovering entry 1EHZ from Protein Databank...\n"
    pdb = PDB()
    tertiary_structures = parse_pdb(pdb.get_entry('1EHZ'))

    print "## 3D annotation ##\n"

    print "List of base-pairs computed with RNAVIEW:\n"

    for ts in tertiary_structures:
        secondary_structure, tertiary_structure = Rnaview().annotate(ts)
        print secondary_structure_to_base_pairs(secondary_structure, keep_tertiaries = True)

    print "\n## 2D prediction ##\n"

    for ts in tertiary_structures:
        print "RNA sequence from 1EHZ:\n"
        print ts.rna.sequence
        print "\nList of base-pairs computed with RNAfold (RNA Vienna Package):\n"
        print Rnafold().fold(molecule=ts.rna)

if __name__ == '__main__':
    test()
