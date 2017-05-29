#!/usr/bin/env python

"""
This script searches for the inner loop of HIV-1 SL1 in 2Ds infered from 3D structures and 2Ds infered from RFAM alignments.
"""

import os, math, sys, re
from pymongo import MongoClient

from pyrna.features import RNA
from pyrna.db import Rfam
from pyrna.parsers import consensus2d_to_base_pairs, to_bn, base_pairs_to_secondary_structure, to_clustalw


def search(db_host = 'localhost', db_port = 27017):
    print "Search in annotated 3Ds"
    client = MongoClient(db_host, db_port)
    db = client['PDB']
    for junction in db['junctions'].find({'location':{ '$size': 2 }}):
        strands = junction['description'].split(' ')
        if strands[0] == "AGG" or strands[1] == "AGG":
            print junction['description']
        if len(strands[0]) == 3 and len(strands[1]) == 1 and not strands[1] == '-' or len(strands[0]) == 1 and not strands[0] == '-' and len(strands[1]) == 3:
            print junction['description']

    print "Search in RFAM alignments"

    rfam = Rfam()
    rfam.generate_seed_alignments()
    
    familiesDetails = rfam.get_families_details()

    i = 0
    hit = 0
    for id in range(100, 101):
        try:
            (rnas, organisms, consensus_2D) = rfam.get_entry(rfam_id = "RF%05u" % id, aln_type = 'full')
            print "Search in RF%05u"% id
            interesting_rnas = []
            for rna in rnas:
                if rna.name in ['M.mulatta.244', 'E.telfairi.194', 'S.araneus.581', 'Macaca_fascicularis_.302', 'Macaca_fascicularis_.413', 'M.mulatta.290'] or re.match('^.+sapiens.+$', rna.name):
                    interesting_rnas.append(rna)
                i += 1
                base_pairs = consensus2d_to_base_pairs(rna, consensus_2D)
                non_aligned_rna = RNA(name = rna.name, sequence = rna.sequence.replace('-',''))
                ss = base_pairs_to_secondary_structure(non_aligned_rna, base_pairs)
                ss.find_junctions()
                for junction in ss.junctions:
                    strands = junction['description'].split(' ')
                    if len(strands) == 2 and (strands[0] == "AGG" and strands[1] == "G" or strands[1] == "AGG" and strands[0] == "G" ):
                        print junction['description']
                        hit += 1
            print "%i sequences processed, %i hit found"%(i, hit)
            with open("RNA_7SK.aln", 'w') as h:
                h.write(to_clustalw(consensus_2D, interesting_rnas, curate=True))

        except Exception, e:
            print str(e)   

if __name__ == '__main__':
    db_host = 'localhost'
    db_port = 27017

    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])

    search(db_host = db_host, db_port = db_port)
