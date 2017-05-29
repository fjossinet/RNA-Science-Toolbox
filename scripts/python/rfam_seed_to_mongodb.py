#!/usr/bin/env python

"""
This script imports all the RFAM seed alignments into an MongoDB database
"""

from pyrna.db import Rfam
from pyrna.parsers import to_clustalw
from bson.objectid import ObjectId
import sys
from pymongo import MongoClient

def import_seed_alignments(db_host = 'localhost', db_port = 27017):

    client = MongoClient(db_host, db_port)
    db = client["RFAM_Seed_Alignments"]

    rfam = Rfam()

    familiesDetails = rfam.get_families_details()

    rfam.generate_seed_alignments()

    for id in range(1, len(familiesDetails)+1):
        print "Import RF%05u" % id

        try:

            (rnas, organisms, consensus2D) = rfam.get_entry(rfam_id = "RF%05u" % id, aln_type = 'seed')

            alignment = {
                '_id': str(ObjectId())
            }

            for rna in rnas:
                rna.organism = rna.name
                rna.name = rna._id #we need to have the _id of the rna in the alignment, not its species name
            
            alignment['alignment'] = to_clustalw(consensus2D, rnas)

            for rna in rnas:
                nse = organisms[rna.organism]
                genomicPositions = map(int, nse.split('/')[1].split('-'))

                ncRNA = {
                    '_id': rna._id,
                    'name': familiesDetails[familiesDetails['accession'] == "RF%05u" % id].iloc[0,4],#id
                    'sequence': rna.sequence,
                    'class': familiesDetails[familiesDetails['accession'] == "RF%05u" % id].iloc[0,2],#family
                    'source': rna.source,
                    'organism': rna.organism,
                    'genomeName': nse.split('/')[0],
                    'genomicPositions':[genomicPositions[0], genomicPositions[1]] if genomicPositions[0] < genomicPositions[1] else [genomicPositions[1], genomicPositions[0]],
                    'genomicStrand': '+' if genomicPositions[0] < genomicPositions[1] else '-',
                    'alignment': alignment['_id']+"@alignments" 
                }

                db["ncRNAs"].insert(ncRNA)

            db['alignments'].insert(alignment)

        except Exception, e:
            import traceback
            traceback.print_exc()

if __name__ == '__main__':
    db_host = 'localhost'
    db_port = 27017

    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])

    import_seed_alignments(db_host = db_host, db_port = db_port)  

    

