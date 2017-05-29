#!/usr/bin/env python
from pyrna.parsers import parse_genbank, to_fasta
from pyrna.db import NCBI
from bson.objectid import ObjectId

"""
This script use NCBI ids to get genomic entries and reformat them as FASTA and GFF3 files.
"""

genome_ids = ["NC_004354.4","NT_033779.5","NT_033778.4","NT_037436.4","NT_033777.3","NC_004353.4","NC_024512.1","NC_024511.2"] #drosophila genome as example

genomic_sequences = []
annotations = []

ncbi = NCBI()
for genome_id in genome_ids:
    print "Processing %s..."%genome_id
    gb_content = ncbi.efetch(ids = [genome_id], db='nucleotide', rettype='gbwithparts')
    dna, features = parse_genbank(gb_content)
    genomic_sequences.append(dna)
    for (index, row) in features.iterrows():
        annotation = {
            '_id': str(ObjectId()),
            'class': row['type'],
            'genomeName': dna.name,
            'genomicStrand': row['genomicStrand'],
            'genomicPositions': row['genomicPositions'],
            'product': row['product']
        }
        annotations.append(annotation)

with open('sequences.fasta', 'w') as f:
    f.write(to_fasta(genomic_sequences))

with open('annotations.gff3', 'w') as f:
    for annotation in annotations:
        f.write("%s\t.\t%s\t%i\t%i\t0.0\t%s\t.\tProduct=%s\n"%(annotation['genomeName'], annotation['class'], annotation['genomicPositions'][0], annotation['genomicPositions'][1], annotation['genomicStrand'], annotation['product']))


