#!/usr/bin/env python
"""
This script aligns NGS reads from a FASTQ file to genomic sequences stored in a MongoDB instance. The aligned reads will then be clustered and the clusters will be stored in the collection "annotations" of a MongoDB database (whose name is one of the script arguments).
"""

from pymongo import MongoClient
from pyrna.parsers import to_fasta
from pyrna.computations import Bowtie
import sys
from pyrna.features import DNA
from bson.objectid import ObjectId
import datetime
import re

from pandas import DataFrame, isnull

def align(client, db):
    target_molecules = []
    genome_id_dic = {}

    organism = None
    inputs = []

    #the genomic sequences are recovered from the collection "genomes"
    for genome in db['genomes'].find():
        dna = DNA(name=genome['name'], sequence=genome['sequence'])
        target_molecules.append(dna)  
        genome_id_dic[genome['name']] = genome['_id']+"@genomes"
        organism = genome['organism']
        inputs.append(genome['_id']+"@genomes")
    
    clusters_of_reads = Bowtie().align(target_molecules, fastq_file)
    
    inputs.append(fastq_file)
    
    outputs = []

    for row in clusters_of_reads.iterrows():
        cluster_of_reads = row[1]
        annotation = { '_id': str(ObjectId()),
                   'source': "tool:bowtie:N.A.",
                   'organism': organism,
        }
        annotation['genomeName'] = cluster_of_reads['genomeName']
        annotation['genome'] = genome_id_dic[cluster_of_reads['genomeName']]
        annotation['genomicPositions'] = [cluster_of_reads['genomicStart'], cluster_of_reads['genomicEnd']]
        annotation['score'] = cluster_of_reads['annotations_count']

        db["annotations"].insert(annotation)

        outputs.append(annotation['_id']+"@annotations")

    computation = {
            '_id': str(ObjectId()),
            'date': str(datetime.datetime.now()),
            'inputs': inputs,
            'outputs': outputs
    }

    db["computations"].insert(computation)
    client.disconnect()

    print "Number of clusters : %i"%len(clusters_of_reads)    

if __name__ == '__main__':

    db_name = None
    fastq_file = None
    db_host = "localhost"
    db_port = 27017

    if "-db" in sys.argv:
        db_name = sys.argv[sys.argv.index("-db")+1]
    if "-fq" in sys.argv:
        fastq_file = sys.argv[sys.argv.index("-fq")+1]
    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])
    if "-sf" in sys.argv:
        sam_file = sys.argv[sys.argv.index("-sf")+1]

    if not db_name and not fastq_file or not db_name:
        print "Usage: align_ngs_reads.py -db database_name -fq fastq_file [-h database_host] [-p database_port]"
        sys.exit()

    client = MongoClient(db_host, db_port)
    db = client[db_name]

    align(client, db)    
    
