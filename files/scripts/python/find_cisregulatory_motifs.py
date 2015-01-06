#!/usr/bin/env python

"""
This script tries to identify cis-regulatory motifs using different approaches:
- CMfinder
- mlocarna
"""

from pymongo import MongoClient

def find(db_host = "localhost", db_port = 27017, databases_to_be_used):

    client = MongoClient(db_host, db_port)

    if not databases_to_be_used:
        databases_to_be_used = client.database_names()   

    protein_sequences = []

    for database_to_be_used in databases_to_be_used:
        db = client[database_to_be_used]
        for cds in db['annotations'].find({'class':'CDS'}):
            protein_sequences.append({
                    'translation': cds['translation'],
                    'species': species,
                    '_id': cds['_id']
                });

if __name__ == '__main__':
    db_host = "localhost"
    db_port = 27017
    databases_to_be_used = None
    find(db_host, db_port, species_to_be_used)

    