#!/usr/bin/env python

import os, sys, re
from pyrna.parsers import parse_fasta
from pyrna.computations import TrnaScanSE
from pyrna.db import NCBI

def scan():

    fasta_sequences = []

    ncbi = NCBI()

    plasmodium_falciparum = ['NC_004325.1', 'NC_000910.2', 'NC_000521.3', 'NC_004318.1', 'NC_004326.1', 'NC_004327.2', 'NC_004328.2', 'NC_004329.2', 'NC_004330.1', 'NC_004314.2', 'NC_004315.2', 'NC_004316.3', 'NC_004331.2', 'NC_004317.2', 'NC_002375.1']

    print re.split('^$',ncbi.efetch(db = 'nucleotide', ids = plasmodium_falciparum[0:2], rettype = 'fasta'))[1]

    #fasta_sequences +=  parse_fasta(ncbi.efetch(db = 'nucleotide', ids = ['NC_000913.3'], rettype = 'fasta'), 'DNA')

    tool = TrnaScanSE()
    print tool.scan(fasta_sequences)

if __name__ == '__main__':

    scan()
