#!/usr/bin/env python
"""
This script recovers all the genbank accession numbers linked to a taxid. It uses a file to keep track of the accession numbers already downloaded.
"""

import os, sys
import xml.etree.ElementTree as ET
from pyrna.db import NCBI

def recover(file, taxid):
    ncbi = NCBI()
    retstart = 0
    iteration_step = 10000
    if os.path.exists('%s_tmp'%file):
        os.renames('%s_tmp'%file, file)
    with open('%s_tmp'%file, 'w') as h2:
        if os.path.exists(file):
            with open(file) as h:
                lines = h.readlines()
                it = 0
                for i in range(0, len(lines)):
                    if lines[i].startswith('retstart:'):
                        it += 1
                        retstart = int(lines[i].split(' ')[1])
                        count = int(lines[i].split(' ')[2])
                        if not count == iteration_step:
                            print "Don't have %i accession numbers for retstart %i. New attempt..."%(iteration_step,retstart) 
                            accession_numbers = []
                            result = ncbi.esearch(db = "nucleotide", term = "txid"+str(taxid)+"[Organism:exp]", retstart = retstart, retmax = iteration_step)
                            try:
                                result = ET.fromstring(result)
                                ids = []
                                if result.find('IdList') is not None:
                                    for id in result.find('IdList').findall('Id'):
                                        ids.append(id.text)
                                    result = ncbi.esummary(db = "nucleotide", ids = ids, retmax = iteration_step)
                                    result = ET.fromstring(result)
                                    for docsum in result.findall('DocSum'):
                                        for item in docsum.findall("Item[@Name='Caption']"):
                                            accession_numbers.append(item.text)
                            except Exception, e:
                                print e
                            h2.write('retstart: %i %i accession_numbers\n'%(retstart, len(accession_numbers)))
                            h2.write(','.join(accession_numbers))
                            h2.write('\n')
                            h2.flush()
                        else:
                            accession_numbers = lines[i+1].split(',')
                            h2.write('retstart: %i %i accession_numbers\n'%(retstart, len(accession_numbers)))
                            h2.write(','.join(accession_numbers))
                            h2.write('\n')
                            h2.flush()
            os.remove(file)            
        
            retstart = it*iteration_step

        while True:
            print 'retstart: %i'%retstart
            accession_numbers =[]
            result = ncbi.esearch(db = "nucleotide", term = "txid"+str(taxid)+"[Organism:exp]", retstart = retstart, retmax = iteration_step)
            try:
                result = ET.fromstring(result)
                ids = []
                if result.find('IdList') is not None:
                    for id in result.find('IdList').findall('Id'):
                        ids.append(id.text)

                    result = ncbi.esummary(db = "nucleotide", ids = ids, retmax = iteration_step)
                    result = ET.fromstring(result)
                    for docsum in result.findall('DocSum'):
                        for item in docsum.findall("Item[@Name='Caption']"):
                            accession_numbers.append(item.text)
                else:
                    break
                
            except Exception, e:
                print e
            print 'got %i accession numbers'%len(accession_numbers)
            h2.write('retstart: %i %i accession_numbers\n'%(retstart, len(accession_numbers)))
            h2.write(','.join(accession_numbers))
            h2.write('\n')
            h2.flush()
            retstart += iteration_step

        os.renames('%s_tmp'%file, file)

    print "End"

if __name__ == '__main__':
    file = None
    taxid = None

    if "-f" in sys.argv:
        file = sys.argv[sys.argv.index("-f")+1]
    if "-id" in sys.argv:
        taxid = int(sys.argv[sys.argv.index("-id")+1])

    if not file or not taxid:
        print "Usage: get_accession_numbers_from_taxid.py -f file -id taxid"
        sys.exit(-1)

    recover(file, taxid)