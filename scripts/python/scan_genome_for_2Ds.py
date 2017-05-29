#!/usr/bin/env python

from pyrna.features import DNA
from pymongo import MongoClient
from pyrna.computations import Rnafold, Rnasubopt
from pyrna.parsers import base_pairs_to_secondary_structure

client = MongoClient()

db = client["Candida_hispaniensis_reannotated"]

rnafold = Rnafold()
rnasubopt = Rnasubopt()

for genome in db['genomes'].find():
    plus_apical_loops = []
    minus_apical_loops = []
    molecule = DNA(name=genome['name'], sequence=genome['sequence'])
    print len(molecule.sequence)
    print "Processing %s"%molecule.name
    window_size = 300
    sliding_window = 150
    i = 0
    while i <= len(molecule)-sliding_window:
        all_secondary_structures = []
        print "%i %i"%(i,i+window_size)
        dna = DNA(name = molecule.name, sequence = molecule[i:i+window_size])
        ss = rnafold.fold(dna)
        ss = base_pairs_to_secondary_structure(dna, ss)
        ss.find_junctions()
        all_secondary_structures.append(ss)
        for ss in rnasubopt.fold(dna, random_sample = 20):
            ss = base_pairs_to_secondary_structure(dna, ss)
            ss.find_junctions()
            all_secondary_structures.append(ss)
        #search for apical loops
        for secondary_structure in all_secondary_structures:
            for junction in secondary_structure.junctions:
                if len(junction['location']) == 1 and len(junction['description']) >= 15:
                    positions = [x+i for x in junction['location'][0]]
                    if positions not in plus_apical_loops:
                        plus_apical_loops.append(positions)
                        print junction['description']    
        print len(plus_apical_loops) 
        i += sliding_window
    #last window
    if i < len(molecule):
        print "%i %i"%(i,len(molecule))
        dna = DNA(name = molecule.name, sequence = molecule[i:i+window_size])
        ss = rnafold.fold(dna)
        ss = base_pairs_to_secondary_structure(dna, ss)
        ss.find_junctions()
        print ss.junctions