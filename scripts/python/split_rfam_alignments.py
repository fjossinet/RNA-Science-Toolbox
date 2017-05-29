#!/usr/bin/env python
"""
This script splits the file Rfam.full or Rfam.seed into one stockholm alignment per RFAM ID. 
The files are created in the same directory.
"""

import sys, os

def split(file):
    output_dir = os.path.dirname(os.path.abspath(file))
    with open(file) as h:
        output = None
        
        for line in h:
            if line.startswith("# STOCKHOLM 1.0"):
                if output:
                    output.close()
                    output = None
                header = line
            elif line.startswith('#=GF AC   RF'):
                output = open("%s/%s_%s.sto"%(output_dir, line.split("#=GF AC")[-1].strip(), 'seed' if file.endswith('seed') else 'full'), 'w')
                output.write(header)
                output.write(line)
                header = ""
            elif output:
                output.write(line)
            elif header:
                header += line   

        output.close()

if __name__ == '__main__':
    file = None
    
    if "-f" in sys.argv:
        file = sys.argv[sys.argv.index("-f")+1]
    if not file:
        print "Usage: split_rfam_alignments.py -f rfam_file"
        sys.exit(-1)

    split(file)