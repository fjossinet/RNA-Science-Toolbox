#!/usr/bin/env python

import sys, commands, subprocess
from pymongo import MongoClient
from pyrna.computations import Tophat, Samtools, Cufflinks, Gmorse
from pyrna.features import DNA
from bson.objectid import ObjectId

"""
This scripts produces new gene models from RNA-seq data and store them in a MongoDB instance. The reads will be aligned with Tophat and the gene models infered with Cufflinks and GMorse
"""

def find(fastq_file, db_name, host, port, gmorse_depth_treshold = 4):
    client = MongoClient(host, port)
    db = client[db_name]

    genomic_sequences = []

    for genome in db['genomes'].find():
        dna = DNA(name = genome['name'], sequence = genome['sequence'])
        dna._id = genome['_id']
        dna.organism = genome['organism']
        genomic_sequences.append(dna)

    #step1: reads alignment with tophat
    print "Tophat..."
    tophat = Tophat()
    sam_file = tophat.align(target_molecules = genomic_sequences, fastq_file = fastq_file, no_convert_bam = True)
    
    #step2: gene models with cufflinks
    print "Cufflinks..."
    samtools = Samtools(sam_file)
    indexed_and_sorted_bam_file = samtools.sort_and_index()
    cufflinks = Cufflinks()
    gene_models = cufflinks.predict_genes(db_name = db_name, bam_file = indexed_and_sorted_bam_file, db_host = host, db_port = port)
    for row in gene_models.iterrows():
        gene_model = row[1]
        for genome in genomic_sequences:
            if genome.name == gene_model['genome']:
                gene_model_descr = {
                    '_id': str(ObjectId()),
                    'source': 'tool:cufflinks:NA',
                    'organism': genome.organism,
                    'genome': genome._id+"@genomes",
                    'genomeName': genome.name,
                    'genomicStrand': '?',
                    'genomicPositions': gene_model['genomicPositions'],
                    'score': gene_model['score'],
                    'class': gene_model['class'],
                    'coverage': gene_model['coverage'],
                    'gene_names': gene_model['gene_names'].split(','),
                    'file': fastq_file
                }
                db['gene_models'].insert(gene_model_descr)
                break

    #steps3: gene models with GMorse
    print "GMorse..."
    gmorse = Gmorse()
    soap_file = "%s.soap"%sam_file.split('.sam')[0]
    subprocess.call("sam2soap.py %s > %s"%(sam_file, soap_file), shell = True)
    unmapped_reads_file = gmorse.extract_unmapped_reads(soap_file, fastq_file)
    depth_coverage_file = gmorse.calculate_depth_coverage(soap_file)
    covtigs_file = gmorse.build_covtigs(depth_coverage_file, gmorse_depth_treshold)
    gene_models = gmorse.make_model(unmapped_reads_file, covtigs_file, genomic_sequences)

    for row in gene_models.iterrows():
        gene_model = row[1]
        for genome in genomic_sequences:
            if genome.name == gene_model['genome']:
                gene_model_descr = {
                    '_id': str(ObjectId()),
                    'source': 'tool:gmorse:NA',
                    'organism': genome.organism,
                    'genome': genome._id+"@genomes",
                    'genomeName': genome.name,
                    'genomicStrand': gene_model['genomicStrand'],
                    'genomicPositions': gene_model['genomicPositions'],
                    'class': gene_model['class'],
                    'file': fastq_file
                }
                db['gene_models'].insert(gene_model_descr)
                break

    client.disconnect()

if __name__ == '__main__':
    if "-f" in sys.argv and "-db" in sys.argv:
        host = "localhost"
        port = 27017
        if "-h" in sys.argv:
            host = sys.argv[sys.argv.index("-h")+1]
        if "-p" in sys.argv:
            port = int(sys.argv[sys.argv.index("-p")+1])
        find(fastq_file = sys.argv[sys.argv.index("-f")+1], db_name = sys.argv[sys.argv.index("-db")+1], host = host, port = port)    
    else:
        print "Usage: find_gene_models.py -f fastq -db database_name [-h mongo_host (default: localhost)] [-p mongo_port (default: 27017)] "
        sys.exit(-1)