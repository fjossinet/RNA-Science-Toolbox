#!/usr/bin/env python
"""
This scripts uses the SVM approach to classify RNA stem-loops based on their 2D
"""

import sys
from pyrna.db import Rfam
from pyrna.parsers import consensus2d_to_base_pairs, base_pairs_to_secondary_structure, to_bn
from pyrna.features import RNA
from pandas import DataFrame

def do_svm(rfam_id):
    rfam = Rfam(use_website = True)
    rnas, organisms, consensus_2d = rfam.get_entry(rfam_id = 'RF%05u'%rfam_id)
    #a matrix for each stem-loop
    stem_loop_descriptions = []
    for i, rna in enumerate(rnas):
        #print to_bn(consensus_2d, len(rna))
        #print rna.sequence
        ss = base_pairs_to_secondary_structure(rna, consensus_2d)
        ss.find_junctions()
        ss.find_stem_loops()

        if i == 0:
            print ss.stem_loops    

        for index, stem_loop in enumerate(ss.stem_loops):
            stem_loop_description = None
            if index >= len(stem_loop_descriptions):
                stem_loop_description = {}
                stem_loop_descriptions.append(stem_loop_description)
            else:
                stem_loop_description =  stem_loop_descriptions[index]       
            #print stem_loop
            for helix in stem_loop['helices']:
                location = helix['location']
                #we extract the sequence for each strand and we remove the gaps
                strand_1 = rna.sequence[location[0][0]-1:location[0][1]].replace('-','')
                strand_2 = rna.sequence[location[1][0]-1:location[1][1]].replace('-','')
                #if len(strand_1) != len(strand_2):
                #    print "not fully conserved helix"
                #    print rna.sequence[location[0][0]-1:location[0][1]], rna.sequence[location[1][0]-1:location[1][1]]
                l = stem_loop_description.get(helix['name']+'_strand_1', [])
                l.append(len(strand_1))
                stem_loop_description[helix['name']+'_strand_1'] = l
                l = stem_loop_description.get(helix['name']+'_strand_2', [])
                l.append(len(strand_2))
                stem_loop_description[helix['name']+'_strand_2'] = l
            for inner_loop in stem_loop['inner_loops']:
                for single_strand in inner_loop['single_strands']:
                    #we extract the sequence for this single-strand and we remove the gaps
                    seq = rna.sequence[single_strand['location'][0]-1:single_strand['location'][1]].replace('-','')
                    #print single_strand['name']
                    l = stem_loop_description.get("inner_loop_%s"%single_strand['name'], [])
                    l.append(len(seq))
                    stem_loop_description["inner_loop_%s"%single_strand['name']] = l

            apical_loop = stem_loop['apical_loop']['single_strands'][0]
            seq = rna.sequence[apical_loop['location'][0]-1:apical_loop['location'][1]].replace('-','')
            l = stem_loop_description.get("apical_loop_%s"%apical_loop['name'], [])
            l.append(len(seq))
            stem_loop_description["apical_loop_%s"%apical_loop['name']] = l
    
    for stem_loop_description in stem_loop_descriptions:
        df = DataFrame(stem_loop_description)
        columns = df.columns
        print columns
        print df.as_matrix(columns)


if __name__ == '__main__':
    if "-id" in sys.argv:
        do_svm(int(sys.argv[sys.argv.index("-id")+1]))    
    else:
        print "Usage: svm.py -id the RFAM ID to process"
        sys.exit(-1)