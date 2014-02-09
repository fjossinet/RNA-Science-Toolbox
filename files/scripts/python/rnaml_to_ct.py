#!/usr/bin/env python

"""
This scripts transform RNAML files into CT files.

Usage: rnaml_to_ct.py [-keep_tertiaries] [-canonical_only] RNAML files

- keep_tertiaries: 
    You can decide to keep or not tertiary interactions. 
    If you keep them, the CT file could not be « standard », in the sense that you could have several base-pairs described for a single position in the sequence.
- canonical_only: 
    You can decide to allow only helices made with canonical base-pairs. 
    For example, the algorithm RNAVIEW, which produces the RNAML files for the NDB, defines helices as a set a stacked base-pairs (canonical or not). 
    If you use this option, non-canonical interactions stacked in helices become tertiary interactions. 
    So, if you combine this option with the previous one, you can produce a CT file describing a kind of « core » secondary structure, made only with helices containing at least two stacked canonical base-pairs.
"""

import sys
from pyrna.parsers import parse_rnaml, secondary_structure_to_base_pairs, to_ct

if __name__ == '__main__':
    keep_tertiaries = "-keep_tertiaries" in sys.argv
    canonical_only = "-canonical_only" in sys.argv

    rnaml_files = []
    for arg in sys.argv:
        if not arg.endswith("rnaml_to_ct.py") and not arg.startswith("-"):
            rnaml_files.append(arg)
    
    for rnaml_file in rnaml_files:
        h = open(rnaml_file, 'r')
        rnaml_content = h.read()
        h.close()

        for secondary_structure in parse_rnaml(rnaml_content, canonical_only = canonical_only):
            print to_ct(secondary_structure_to_base_pairs(secondary_structure, keep_tertiaries = keep_tertiaries), secondary_structure.rna) 




    