#!/usr/bin/env python

"""
Refine snoRNA clusters in a Mongo database to remove overlapping clusters
"""

import sys
from pymongo import MongoClient
from bson.objectid import ObjectId
from pyrna.utils import cluster_genomic_annotations


def refine_snoclusters(db_name, db_host = "localhost", db_port = 27017):
    client = MongoClient(db_host, db_port)
    db = client[db_name]

    #sort snoRNA clusters from the database per snoRNA class (2 different classes), per strand (2 different strands) and per genome (n different genomes). 
    genome_names = []
    for genome in db['genomes'].find():
        genome_names.append(genome['_id']+"@genomes")

    sorted_clusters = []
    for i in range(0, len(genome_names)*4): #if n genomes => n*2*2 sub-lists of clusters
        sorted_clusters.append([])
    for snocluster in db['ncRNAs'].find({'source': "tool:letitsno:NA"}):
        index = 0 if snocluster['class'] == 'Gene, snRNA, snoRNA, CD-box' else 1
        index += 0 if snocluster['genomicStrand'] == '+' else 2
        j = 0
        while j != len(genome_names):
            genome_name = genome_names[j]
            if snocluster['genome'] == genome_name:
                index += range(0, len(genome_names)*4, 4)[j] #if n genomes => range_list=[0, 4, 8, 12, ..., n-1*2*2]
            j += 1
        snocluster['genomicStart'] = snocluster['genomicPositions'][0]
        snocluster['genomicEnd'] = snocluster['genomicPositions'][1]
        sorted_clusters[index].append(snocluster)

    for cluster_annotations in sorted_clusters:
        if len(cluster_annotations):
            #clustering of snoRNA clusters whose genomic positions overlap 
            new_clusters = cluster_genomic_annotations(cluster_annotations, threshold = 2, fill_cluster_with_genomic_annotations = True)

            #storage of enlarged snoRNA clusters in database and removing of all "old" clusters constituting these "new" clusters
            if len(new_clusters):
                for new_cluster in new_clusters:
                    new_cluster['_id'] = str(ObjectId())
                    new_cluster['name'] = new_cluster['genomic_annotations'][0]['name']
                    new_cluster['source'] = "tool:letitsno:NA"
                    new_cluster['class'] = new_cluster['genomic_annotations'][0]['class']
                    new_cluster['organism'] = new_cluster['genomic_annotations'][0]['organism']
                    new_cluster['genome'] = new_cluster['genomic_annotations'][0]['genome']
                    new_cluster['genomicStrand'] = new_cluster['genomic_annotations'][0]['genomicStrand']
                    new_cluster['genomicPositions'] = [new_cluster['genomicStart'], new_cluster['genomicEnd']]
                    new_cluster['ids'] = new_cluster['genomic_annotations'][0]['ids']

                    for cluster in new_cluster['genomic_annotations'][1:]:
                        new_cluster['name'] += cluster['name'][7:]
                        new_cluster['ids'].extend(cluster['ids']) #no double of ncRNAs _id
                        db["ncRNAs"].remove({'_id': cluster['_id']}) #db["ncRNAs"].remove(cluster) don't remove cluster

                    db["ncRNAs"].remove({'_id': new_cluster['genomic_annotations'][0]['_id']})
                    del new_cluster['genomicStart']
                    del new_cluster['genomicEnd']
                    del new_cluster['genomic_annotations']
                    del new_cluster['annotations_count']
                    db["ncRNAs"].insert(new_cluster)

    client.disconnect()

if __name__ == '__main__':
    db_name = None
    db_host = "localhost"
    db_port = 27017

    if "-db" in sys.argv:
        db_name = sys.argv[sys.argv.index("-db")+1]
    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])

    if not db_name:
        print "Usage: refine_snoclusters.py -db db_name [-h database_host] [-p database_port]"
        sys.exit()

    refine_snoclusters(db_name = db_name, db_host = db_host, db_port = db_port)

