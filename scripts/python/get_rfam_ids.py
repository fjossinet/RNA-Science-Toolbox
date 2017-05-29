#!/usr/bin/env python

"""
This scripts generates some Java code for Assemble2. 
"""

from pyrna.db import Rfam

families = Rfam().get_families_details()

fam = {}

ids = []

for index,row in families.iterrows():
    if fam.has_key(row['family']):
        fam[row['family']].append(row['id'])
    else:
        fam[row['family']] = []
        fam[row['family']].append(row['id'])

for key in sorted(fam.keys()):
    print '"%s",'%key

print 'List<String> ids = null;'

for key in fam:
    print 'ids = new ArrayList<String>();';
    print 'ncRNA_ids.put("%s", ids);'%key
    for id in sorted(fam[key]):
        print 'ids.add("%s");'%id;