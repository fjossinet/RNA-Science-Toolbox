#!/usr/bin/env python
"""
Script to manage gLite jobs
"""

import sys, os
from pyrna import glite
from pymongo import MongoClient

if __name__ == '__main__':

    if "-f" in sys.argv:
        vienna_file = sys.argv[sys.argv.index("-f")+1]
    if len(sys.argv) == 1 or not sys.argv[1] in ["init", "check", "resubmit",  "cancel", "recover", "rfam_families_to_resubmit"]:
        print "Usage: glite.py [init | check | resubmit | cancel | recover | rfam_families_to_resubmit] "
        sys.exit(-1)

    if sys.argv[1] == "init":
        glite.init_proxy(print_output=True)
    elif sys.argv[1] == "check":
        if not "-o" in sys.argv :
            print "Usage: glite.py check -o submission_result_file_name"
            sys.exit(-1)    
        glite.check_jobs_statuses(sys.argv[sys.argv.index("-o")+1], serialize=True)
    elif sys.argv[1] == "resubmit":
        if not "-f" in sys.argv or not "-o" in sys.argv:
            print "Usage: glite.py resubmit -f job_list_file -o submission_result_file_name"
            sys.exit(-1)   
        glite.resubmit_jobs(sys.argv[sys.argv.index("-f")+1], sys.argv[sys.argv.index("-o")+1])
    elif sys.argv[1] == "cancel":
        if not "-o" in sys.argv :
            print "Usage: glite.py cancel -o job_list_file"
            sys.exit(-1)
        glite.cancel_jobs(sys.argv[sys.argv.index("-o")+1])
    elif sys.argv[1] == "recover":
        if not "-d" in sys.argv or not "-o" in sys.argv:
            print "Usage: glite.py recover -o job_list_file -d output_dir"
            sys.exit(-1)
        glite.recover_job_outputs(sys.argv[sys.argv.index("-o")+1], sys.argv[sys.argv.index("-d")+1])
    elif sys.argv[1] == "rfam_families_to_resubmit":
        if len(sys.argv) <= 2:
            print "Usage: glite.py rfam_families_to_resubmit job_list_file_1 job_list_file_1 ..."
            sys.exit(-1)
        families_to_resubmit = []
        client = MongoClient("localhost", 27017)
        db = client['comparative_genomics']
        for job_list_file in sys.argv[2:]:
            print "Recovering from", job_list_file
            with open(job_list_file) as f:
                for line in f:
                    job_id = line.split(' ')[0].split('job_')[-1].split('.jdl')[0]
                    with open(os.path.dirname(os.path.realpath(job_list_file))+"/script_%s.sh"%job_id) as script:
                        for line in script:
                            if line.startswith("#Families processed:"):
                                line = line.split('#Families processed:')[1].strip()
                                if len(line): #we could have a script with a line like #Families processed:
                                    ids = line.split(',')
                                    for id in ids:
                                        if not id.startswith('RF'):
                                            in_house_alignment = db['alignments'].find_one({'_id':id})
                                            families_to_resubmit.append(int(in_house_alignment['source'].split('RF')[1]))         
                                        else:
                                            families_to_resubmit.append(int(id.split('RF')[1]))    
                                break
        client.close()

        print '%i families to resubmit.'%len(families_to_resubmit)
        print ','.join([str(x) for x in families_to_resubmit])             

