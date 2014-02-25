#!/usr/bin/env python
"""
Script to manage GLite jobs
"""

import sys
from pyrna import glite

if __name__ == '__main__':

    if "-f" in sys.argv:
        vienna_file = sys.argv[sys.argv.index("-f")+1]
    if len(sys.argv) == 1 or not sys.argv[1] in ["init", "check", "resubmit",  "cancel", "recover"]:
        print "Usage: glite.py [init | check | resubmit | cancel | recover] "
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
