#!/usr/bin/env python

import os, sys, commands, subprocess, pickle
from collections import Counter, OrderedDict
from pymongo import MongoClient
from pyrna.computations import Tophat, Samtools, Cufflinks, Gmorse
from pyrna.features import DNA
from pyrna.parsers import parse_sam, parse_fasta
from bson.objectid import ObjectId
from pyrna.db import NCBI

"""
This scripts produces new gene models from RNA-seq data and store them in a MongoDB instance. The reads will be aligned with Tophat and the gene models infered with Cufflinks and GMorse
"""

def find(fastq_file, host, port, db_name = None, genome_ids = None, reverse = False, sam_file = None, gmorse_depth_treshold = 2):
    client = MongoClient(host, port)
    db = None

    genomic_sequences = []

    if db_name:
        db = client[db_name]
        for genome in db['genomes'].find():
            dna = DNA(name = genome['name'], sequence = genome['sequence'])
            dna._id = genome['_id']
            dna.organism = genome['organism']
            genomic_sequences.append(dna)

    if genome_ids:
        ncbi = NCBI()
        for genome_id in genome_ids:
            fasta_content = ncbi.efetch(ids = [genome_id], db='nucleotide', rettype='fasta')
            dna = parse_fasta(fasta_content, type = 'DNA')[0]
            dna.name = genome_id.split('.')[0]
            print dna.name
            genomic_sequences.append(dna)

        for genomic_sequence in genomic_sequences:
            print '%s: %i nts'%(genomic_sequence.name, len(genomic_sequence.sequence))

    current_dir = None
    if fastq_file:
        current_dir = os.path.dirname(os.path.realpath(fastq_file))
    elif sam_file:
        current_dir = os.path.dirname(os.path.realpath(sam_file))
    else:
        print "Neither fastq nor SAM file provided!!"
        sys.exit(1)

    #step1: reads alignment with tophat
    if not sam_file:
        print "Processing %s..."%fastq_file
        print "Tophat..."
        tophat = Tophat(current_dir, user_defined_options=['--read-realign-edit-dist 0', '-p 8']) #the option --read-realign-edit-dist 0 is to handle the problem with the pseudogenes
        sam_file = tophat.align(target_molecules = genomic_sequences, fastq_file = fastq_file, no_convert_bam = True)
        samtools = Samtools(sam_file)
        indexed_and_sorted_bam_file = samtools.sort_and_index()

    #step2: coverages
    if False:
        print "Coverages..."
        aligned_reads, total_reads, tids  = parse_sam(sam_file)
        reads = []
        total_aligned_reads = 0
        total_genome_covered = 0
        total_genome_size = 0
        for aligned_reads_per_genome in aligned_reads:
            total_aligned_reads += len(aligned_reads_per_genome)
            current_molecule = None
            for genomic_sequence in genomic_sequences:
                if genomic_sequence.name == tids[aligned_reads_per_genome[0]['tid']] or genomic_sequence.name.startswith(tids[aligned_reads_per_genome[0]['tid']]):
                    current_molecule = genomic_sequence
                    break

            length = len(current_molecule.sequence)
            total_genome_size += length
            
            if reverse:
                counts_plus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '-')
                counts_minus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '+')
            else:
                counts_plus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '+')
                counts_minus = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '-')
            
            #counts_undefined = Counter(aligned_read['genomicStart'] for aligned_read in aligned_reads_per_genome if aligned_read['genomicStrand'] == '?')

            counts_plus = OrderedDict(sorted(counts_plus.items(), key=lambda count: count[0]))
            counts_minus = OrderedDict(sorted(counts_minus.items(), key=lambda count: count[0]))
            #counts_undefined = OrderedDict(sorted(counts_undefined.items(), key=lambda count: count[0]))

            _counts_plus = []
            for i in xrange(0, length):
                _counts_plus.append(0)   
            _counts_minus = []
            for i in xrange(0, length):
                _counts_minus.append(0)
            #_counts_undefined = []
            
            for pos in xrange(0, length-49):
                plus_c = counts_plus.get(pos+1, 0)
                for _pos in xrange(0, 50):
                    _counts_plus[pos+_pos] = _counts_plus[pos+_pos] + plus_c
                minus_c = counts_minus.get(pos+1, 0)
                for _pos in xrange(0, 50):
                    _counts_minus[pos+_pos] = _counts_minus[pos+_pos] + minus_c
                #_counts_undefined.append(counts_undefined.get(pos+1, 0))

            _counts_both_strands = [x+y for x,y in zip(_counts_plus, _counts_minus)]

            s = sum(x > 0 for x in _counts_both_strands)

            print "%i/%i positions covered in %s"%(s, length, current_molecule.name)
            
            total_genome_covered += s

            with open("%s_%s.pickle"%(current_molecule.name, '+'), 'wb') as f:
                pickle.dump(_counts_plus, f)
            with open("%s_%s.pickle"%(current_molecule.name, '-'), 'wb') as f:
                pickle.dump(_counts_minus, f)
            #with open("%s_%s.pickle"%(current_molecule.name, '?'), 'wb') as f:
            #    pickle.dump(_counts_undefined, f)

        print "%i reads aligned"%total_aligned_reads
        print "%i/%i positions covered"%(total_genome_covered, total_genome_size)
    
    #step3: gene models with cufflinks
    if False:
        print "Cufflinks..."
        cufflinks = Cufflinks(current_dir)
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
                        'sample': fastq_file.split('.fastq')[0]
                    }
                    db['gene_models'].insert(gene_model_descr)
                    break

    #steps4: gene models with GMorse
    if False:
        print "GMorse..."
        gmorse = Gmorse(current_dir)
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
                        'sample': fastq_file.split('.fastq')[0]
                    }
                    db['gene_models'].insert(gene_model_descr)
                    break

    client.disconnect()

if __name__ == '__main__':
    host = "localhost"
    port = 27017
    db_name = None
    reverse = "-r" in sys.argv
    sam_file = None
    fastq_file = None
    genome_ids = None
    if "-h" in sys.argv:
        host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        port = int(sys.argv[sys.argv.index("-p")+1])
    if "-sam" in sys.argv:
        sam_file =  sys.argv[sys.argv.index("-sam")+1]
    if "-f" in sys.argv:
        fastq_file =  sys.argv[sys.argv.index("-f")+1]
    if "-ids" in sys.argv:
        genome_ids =  sys.argv[sys.argv.index("-ids")+1].split(',')
    if "-db" in sys.argv:
        db_name = sys.argv[sys.argv.index("-db")+1]

    find(fastq_file, host, port, db_name, genome_ids, reverse, sam_file)