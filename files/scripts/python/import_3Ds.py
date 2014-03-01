#!/usr/bin/env python

"""
This script delegates the work to the grid task ../../import_3Ds.task
It is a way to bypass the need to have access to a grid (but the task is longer)
"""

from pyrna.db import RNA3DHub, PDB, PDBQuery
import os, math, sys

def import_3Ds(db_host = 'localhost', db_port = 27017, rna3dhub = False, canonical_only = True, annotate = False):

    total_structures = 0
    if not rna3dhub:
        query ="""<orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <description>Chain Type: there is a Protein and a RNA chain but not any DNA or Hybrid</description>
    <containsProtein>Y</containsProtein>
    <containsDna>N</containsDna>
    <containsRna>Y</containsRna>
    <containsHybrid>N</containsHybrid>
  </orgPdbQuery>"""
        total_structures = len(PDB().query(query))
    else:
        total_structures = len(RNA3DHub().get_clusters())

    print "%i 3Ds to process"%total_structures

    total_jobs = int(math.floor(total_structures/100)+1) #100 structures per job

    for job_id in range (1, total_jobs+1):
        os.system("import_3Ds.task -h %s -p %i -rna3dhub %s -id %i -canonical_only %s %s"%(db_host, db_port, "Y" if rna3dhub else "N", job_id, "Y" if canonical_only else "N", "-annotate" if annotate else ""))

if __name__ == '__main__':
    db_host = 'localhost'
    db_port = 27017
    rna3dhub = False
    canonical_only = False

    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])
    rna3dhub =  "-rna3dhub" in sys.argv
    canonical_only =  "-canonical_only" in sys.argv
    annotate = "-annotate" in sys.argv

    import_3Ds(db_host = db_host, db_port = db_port, rna3dhub = rna3dhub, canonical_only = canonical_only, annotate = annotate)


