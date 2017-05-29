#!/usr/bin/env python

from pyrna.db import PDB, Rfam
from pyrna.parsers import parse_pdb, to_bn
from pyrna.computations import Cmalign, Rnaview
from bson.objectid import ObjectId
import os

pdb = PDB()
cmalign = Cmalign()
rnaview = Rnaview()
rfam = Rfam(cache_dir = "/home/fjossinet/tmp/Rfam")

families_with_structures = rfam.get_families_with_structures()

for index, row in families_with_structures.iterrows():

    rfam_id = row['rfam_id']

    tertiary_structures = parse_pdb(pdb.get_entry(row['pdb_id']))

    reference_rna = None
    ts = None

    for tertiary_structure in tertiary_structures:
        if tertiary_structure.rna.name == row['chain_name']:
            ts = tertiary_structure
            reference_rna = ts.rna
            break

    if ts:
        secondary_structure, tertiary_structure = rnaview.annotate(tertiary_structure = ts)
        rnas, orgs, consensus_2d = cmalign.align([reference_rna], rfam_id = rfam_id, rfam = rfam)
        os.mkdir("/home/fjossinet/tmp/%s"%rfam_id)
        os.mkdir("/home/fjossinet/tmp/%s/Molecules"%rfam_id)
        os.mkdir("/home/fjossinet/tmp/%s/SecondaryStructures"%rfam_id)
        os.mkdir("/home/fjossinet/tmp/%s/TertiaryStructures"%rfam_id)
        os.mkdir("/home/fjossinet/tmp/%s/StructuralAlignments"%rfam_id)

        with open("/home/fjossinet/tmp/%s/Molecules/%s.rnaml"%(rfam_id,reference_rna._id), 'w') as h:
            h.write("""<?xml version="1.0" encoding="UTF-8"?>
            <rnaml>
              <molecule type="rna">
                <identity>
                  <name>%s</name>
                </identity>
                <sequence length="%i">
                  <seq-data>%s</seq-data>
                </sequence>
              </molecule>
            </rnaml>"""%(reference_rna.name,len(reference_rna),reference_rna.sequence))

        for rna in rnas:
            if rna.name != reference_rna.name:
                with open("/home/fjossinet/tmp/%s/Molecules/%s.rnaml"%(rfam_id,rna._id), 'w') as h:
                    h.write("""<?xml version="1.0" encoding="UTF-8"?>
                    <rnaml>
                      <molecule type="rna">
                        <identity>
                          <name>%s</name>
                        </identity>
                        <sequence length="%i">
                          <seq-data>%s</seq-data>
                        </sequence>
                      </molecule>
                    </rnaml>"""%(rna.name,len(rna.sequence.replace('-','')),rna.sequence.replace('-','')))

        with open("/home/fjossinet/tmp/%s/TertiaryStructures/%s.rnaml"%(rfam_id, ts._id), 'w') as h:
            h.write("""<?xml version="1.0" encoding="UTF-8"?>\n""")
            h.write("""<rnaml><tertiary-structure name="Tertiary Structure" molecule-ids="%s.rnaml">\n"""%reference_rna._id)
            last_position = None
            for index, row in ts.get_atoms().iterrows():
                if last_position != row['absolute position']:
                    if last_position:
                        h.write("""</base>\n""")
                    last_position = row['absolute position']
                    h.write("""<base position="%i" base-id="%s">\n"""%(row['absolute position'], row['position label']))
                h.write("""<atom type="%s" x="%f" y="%f" z="%f" />\n"""%(row['name'], row['x'], row['y'], row['z'])) 
            h.write("""</base>\n""")
            h.write("""</tertiary-structure></rnaml>""")

        with open("/home/fjossinet/tmp/%s/SecondaryStructures/%s.rnaml"%(rfam_id,secondary_structure._id), 'w') as h:
            h.write("""<?xml version="1.0" encoding="UTF-8"?>\n""")
            h.write("""<rnaml><structure-annotation name="Computed with RNAVIEW" molecule-ids="%s.rnaml" tertiary-structure-id="%s.rnaml">\n"""%(reference_rna._id, ts._id))
            for helix in secondary_structure.helices:
                h.write("""<helix name="%s" base5-id="%i" base3-id="%i" length="%i" />\n"""%(helix['name'], helix['location'][0][0], helix['location'][-1][-1], helix['length']))
            for single_strand in secondary_structure.single_strands:
                h.write("""<single-strand name="%s" base5-id="%i" base3-id="%i" />\n"""%(single_strand['name'], single_strand['location'][0], single_strand['location'][-1]))
            h.write("""</structure-annotation></rnaml>""")

        with open("/home/fjossinet/tmp/%s/StructuralAlignments/%s.rnaml"%(rfam_id,str(ObjectId())), 'w') as h:
            h.write("""<?xml version="1.0" encoding="UTF-8"?>\n""")
            h.write("""<rnaml><alignment name="Structural Alignment" structure-annotation-id="%s.rnaml">\n"""%secondary_structure._id)
            h.write("""<consensus2D>%s</consensus2D>\n"""%to_bn(consensus_2d, len(rnas[0])))
            i = 1
            for rna in rnas:
                if rna.name == reference_rna.name:
                    h.write("""<ali-sequence molecule-id="%s.rnaml" position="0">\n"""%reference_rna._id)    
                else:
                    h.write("""<ali-sequence molecule-id="%s.rnaml" position="%i">\n"""%(rna._id,i) )
                    i+=1
                gaps_positions =  [x+1 for x in rna.get_gaps_positions()]
                print gaps_positions
                location = ""
                length = 1
                last_pos = gaps_positions[0]
                for pos in gaps_positions[1:]:
                    if pos-1 == last_pos:
                        length += 1
                    else:
                        if length == 1:
                            location += "%i,"%last_pos
                        else:
                            location += "%i:%i,"%(last_pos-length+1, length)
                        length = 1
                    last_pos = pos
                if length == 1:
                    location += "%i"%last_pos
                else:
                    location += "%i:%i"%(last_pos-length+1, length)
                print location
                print rna.sequence        
                h.write("""<structural-identity start="1" end="%i">%s</structural-identity>\n"""%(len(rna.sequence.replace('-','')), location))
                h.write("""</ali-sequence>\n""")
            h.write("""</alignment></rnaml>""")



