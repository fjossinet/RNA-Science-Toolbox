#!/usr/bin/env python

"""
A script to test PyRNA
"""

from pyrna.db import PDB
from pyrna.parsers import parse_pdb
from pyrna.computations import Rnafold

def test():
    pdb = PDB()
    tertiary_structures = parse_pdb(pdb.get_entry('1EHZ'))

    for ts in tertiary_structures:
        print "\nRNA sequence from 1EHZ:\n"
        print ts.rna.sequence
        print "\nList of base-pairs computed from RNAfold (RNA Vienna Package):\n"
        print Rnafold().fold(molecule=ts.rna)

if __name__ == '__main__':
    test()