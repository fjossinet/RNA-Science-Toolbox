#!/usr/bin/env python

import sys, pickle
from collections import Counter, OrderedDict
from pyrna.features import DNA
from pyrna.parsers import parse_sam
from pymongo import MongoClient
from bson.objectid import ObjectId

def import_sam(host, port, db, sam_file, sample_name, reverse = False):
    client = MongoClient(db_host, db_port)
    db = client[db]

    genomic_sequences = []

    for genome in db['genomes'].find():
        dna = DNA(name = genome['name'], sequence = genome['sequence'])
        dna._id = genome['_id']
        dna.organism = genome['organism']
        genomic_sequences.append(dna)

    print len(genomic_sequences)

    aligned_reads, total_reads, tids  = parse_sam(sam_file)
    reads = []
    total_aligned_reads = 0
    for aligned_reads_per_genome in aligned_reads:
        total_aligned_reads += len(aligned_reads_per_genome)
        current_molecule = None
        for genomic_sequence in genomic_sequences:
            if genomic_sequence.name == tids[aligned_reads_per_genome[0]['tid']]:
                current_molecule = genomic_sequence
                break

        length = len(current_molecule.sequence)
        
        if reverse:
            counts_plus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '-')
            counts_minus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '+')
        else:
            counts_plus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '+')
            counts_minus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '-')
        
        counts_undefined = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '?')

        counts_plus = OrderedDict(sorted(counts_plus.items(), key=lambda count: count[0]))
        counts_minus = OrderedDict(sorted(counts_minus.items(), key=lambda count: count[0]))
        counts_undefined = OrderedDict(sorted(counts_undefined.items(), key=lambda count: count[0]))

        _counts_plus = []
        _counts_minus = []
        _counts_undefined = []

        for pos in xrange(0, length):
            _counts_plus.append(counts_plus.get(pos+1, 0))
            _counts_minus.append(counts_minus.get(pos+1, 0))
            _counts_undefined.append(counts_undefined.get(pos+1, 0))

        with open("%s_%s.pickle"%(current_molecule.name, '+'), 'wb') as f:
            pickle.dump(_counts_plus, f)
        with open("%s_%s.pickle"%(current_molecule.name, '-'), 'wb') as f:
            pickle.dump(_counts_minus, f)
        #with open("%s_%s.pickle"%(current_molecule.name, '?'), 'wb') as f:
        #    pickle.dump(_counts_undefined, f)

    print "%i reads processed"%total_aligned_reads

if __name__ == '__main__':
    db_name = None
    db_host = "localhost"
    db_port = 27017
    sam_file = None
    sample_name = None
    reverse = False

    if "-f" in sys.argv:
        sam_file = sys.argv[sys.argv.index("-f")+1]
    if "-s" in sys.argv:
        sample_name = sys.argv[sys.argv.index("-s")+1]
    if "-db" in sys.argv:
        db_name = sys.argv[sys.argv.index("-db")+1]
    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])
    if "-r" in sys.argv:
        reverse = True


    import_sam(db_host, db_port, db_name, sam_file, sample_name, reverse)
