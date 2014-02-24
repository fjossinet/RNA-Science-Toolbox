"""
This scripts check all the let_it_sno tasks: 
- if some jobs have the status aborted, done_but_failed, done_without_success or unknown_status, they're resubmitted automatically
- if all jobs have succeed, the last step is launched (the refinement for clusters)

Use cron to time-schedule this script, like: 0,30 * * * * YOUR_PATH/check_letitsno_tasks.sh > /dev/null 
"""

import os, re, sys, datetime
from pyrna import glite
from pymongo import MongoClient
from bson.objectid import ObjectId
from pyrna.utils import cluster_genomic_annotations
from pyrna.glite import resubmit_jobs

def refine_snoclusters(db_name, db_host = "localhost", db_port = 27017):
    client = MongoClient(db_host, db_port)
    db = client[db_name]

    #sort snoRNA clusters from the database per snoRNA class (2 different classes), per strand (2 different strands) and per genome (n different genomes). 
    genome_names = []
    for genome in db['genomes'].find():
        genome_names.append(genome['_id']+"@genomes")

    sorted_clusters = []
    for i in range(0, len(genome_names)*4): #if n genomes => n*2*2 sub-lists of clusters
        sorted_clusters.append([])
    for snocluster in db['ncRNAs'].find({'source': "tool:letitsno:NA"}):
        index = 0 if snocluster['class'] == 'Gene, snRNA, snoRNA, CD-box' else 1
        index += 0 if snocluster['genomicStrand'] == '+' else 2
        j = 0
        while j != len(genome_names):
            genome_name = genome_names[j]
            if snocluster['genome'] == genome_name:
                index += range(0, len(genome_names)*4, 4)[j] #if n genomes => range_list=[0, 4, 8, 12, ..., n-1*2*2]
            j += 1
        snocluster['genomicStart'] = snocluster['genomicPositions'][0]
        snocluster['genomicEnd'] = snocluster['genomicPositions'][1]
        sorted_clusters[index].append(snocluster)

    for cluster_annotations in sorted_clusters:
        if len(cluster_annotations):
            #clustering of snoRNA clusters whose genomic positions overlap 
            new_clusters = cluster_genomic_annotations(cluster_annotations, threshold = 2, fill_cluster_with_genomic_annotations = True)

            #storage of enlarged snoRNA clusters in database and removing of all "old" clusters constituting these "new" clusters
            if len(new_clusters):
                for new_cluster in new_clusters:
                    new_cluster['_id'] = str(ObjectId())
                    new_cluster['name'] = new_cluster['genomic_annotations'][0]['name']
                    new_cluster['source'] = "tool:letitsno:NA"
                    new_cluster['class'] = new_cluster['genomic_annotations'][0]['class']
                    new_cluster['organism'] = new_cluster['genomic_annotations'][0]['organism']
                    new_cluster['genome'] = new_cluster['genomic_annotations'][0]['genome']
                    new_cluster['genomicStrand'] = new_cluster['genomic_annotations'][0]['genomicStrand']
                    new_cluster['genomicPositions'] = [new_cluster['genomicStart'], new_cluster['genomicEnd']]
                    new_cluster['ids'] = new_cluster['genomic_annotations'][0]['ids']

                    for cluster in new_cluster['genomic_annotations'][1:]:
                        new_cluster['name'] += cluster['name'][7:]
                        new_cluster['ids'].extend(cluster['ids']) #no double of ncRNAs _id
                        db["ncRNAs"].remove({'_id': cluster['_id']}) #db["ncRNAs"].remove(cluster) don't remove cluster

                    db["ncRNAs"].remove({'_id': new_cluster['genomic_annotations'][0]['_id']})
                    del new_cluster['genomicStart']
                    del new_cluster['genomicEnd']
                    del new_cluster['genomic_annotations']
                    del new_cluster['annotations_count']
                    db["ncRNAs"].insert(new_cluster)

    #now we just need to record something to let know that this step is done
    computation = {
            '_id': str(ObjectId()),
            'date': str(datetime.datetime.now()),
            'source': "tool:check_letitsno_task:NA"
        }

    db["computations"].insert(computation)

def check(db_host = "localhost", db_port = 27017):
    dir = os.getenv("HOME")+"/tmp/"
    for f in os.listdir(dir):
        if f.startswith("jobs_let_it_sno.task_on"):
            db_name = f.split('jobs_let_it_sno.task_on_')[1]
            client = MongoClient(db_host, db_port)
            db = client[db_name]
            if not db['computations'].find({"source": {'$regex':'tool:check_letitsno_task:.+'} }).count():
                submission_files = [submission_file for submission_file in os.listdir("%s/%s"%(dir, f)) if submission_file.startswith("submission_result")]
                submission_files_count = len(submission_files)
                everything_done = True
                for submission_file in submission_files:
                    glite.check_jobs_statuses("%s/%s/%s"%(dir, f, submission_file), serialize=True, verbose = False)
                    
                    #if we have unknown_status, aborted, done_but_failed or done_without_success jobs, we remove the annoying jobs from the submission file and we resubmit them
                    
                    if os.path.getsize("%s/%s/%s"%(dir, f, "done_but_failed_jobs")):
                        h = open("%s/%s/%s"%(dir, f, "done_but_failed_jobs"), 'r')
                        lines_to_remove = h.readlines()
                        h.close()

                        lines_to_keep = []
                        h = open("%s/%s/%s"%(dir, f, submission_file), 'r')
                        for l in h.readlines():
                            if not l in lines_to_remove:
                                lines_to_keep.append(l)    
                        h.close()

                        h = open("%s/%s/%s"%(dir, f, submission_file), 'w')
                        h.write(''.join(lines_to_keep))
                        h.close()

                        submission_files_count += 1
                        resubmit_jobs("%s/%s/%s"%(dir, f, "done_but_failed_jobs"), "submission_result_%i"%submission_files_count)
                        everything_done = False

                    if os.path.getsize("%s/%s/%s"%(dir, f, "done_without_success_jobs")):
                        h = open("%s/%s/%s"%(dir, f, "done_without_success_jobs"), 'r')
                        lines_to_remove = h.readlines()
                        h.close()

                        lines_to_keep = []
                        h = open("%s/%s/%s"%(dir, f, submission_file), 'r')
                        for l in h.readlines():
                            if not l in lines_to_remove:
                                lines_to_keep.append(l)    
                        h.close()

                        h = open("%s/%s/%s"%(dir, f, submission_file), 'w')
                        h.write(''.join(lines_to_keep))
                        h.close()

                        submission_files_count += 1
                        resubmit_jobs("%s/%s/%s"%(dir, f, "done_without_success_jobs"), "submission_result_%i"%submission_files_count)
                        everything_done = False
                        
                    if os.path.getsize("%s/%s/%s"%(dir, f, "unknown_status_jobs")):
                        h = open("%s/%s/%s"%(dir, f, "unknown_status_jobs"), 'r')
                        lines_to_remove = h.readlines()
                        h.close()

                        lines_to_keep = []
                        h = open("%s/%s/%s"%(dir, f, submission_file), 'r')
                        for l in h.readlines():
                            if not l in lines_to_remove:
                                lines_to_keep.append(l)    
                        h.close()

                        h = open("%s/%s/%s"%(dir, f, submission_file), 'w')
                        h.write(''.join(lines_to_keep))
                        h.close()
                        
                        submission_files_count += 1
                        resubmit_jobs("%s/%s/%s"%(dir, f, "unknown_status_jobs"), "submission_result_%i"%submission_files_count)
                        everything_done = False

                    if os.path.getsize("%s/%s/%s"%(dir, f, "aborted_jobs")):
                        h = open("%s/%s/%s"%(dir, f, "aborted_jobs"), 'r')
                        lines_to_remove = h.readlines()
                        h.close()

                        lines_to_keep = []
                        h = open("%s/%s/%s"%(dir, f, submission_file), 'r')
                        for l in h.readlines():
                            if not l in lines_to_remove:
                                lines_to_keep.append(l)    
                        h.close()

                        h = open("%s/%s/%s"%(dir, f, submission_file), 'w')
                        h.write(''.join(lines_to_keep))
                        h.close()

                        submission_files_count += 1
                        resubmit_jobs("%s/%s/%s"%(dir, f, "aborted_jobs"), "submission_result_%i"%submission_files_count)
                        everything_done = False

                    #if not all the jobs have succeeded, it is not done
                    if os.path.getsize("%s/%s/%s"%(dir, f, "succeeded_jobs")) != os.path.getsize("%s/%s/%s"%(dir, f, submission_file)):
                        everything_done = False
                        print "%s: Not all the jobs have succeeded for %s"%(str(datetime.datetime.now()), f)

                if everything_done:
                    print "%s: all the jobs have succeeded for %s => cluster refinement"%(str(datetime.datetime.now()), f)
                    refine_snoclusters(db_name = db_name, db_host = db_host, db_port = db_port)
            else:
                print "%s: cluster refinement already done for %s"%(str(datetime.datetime.now()), f)

if __name__ == '__main__':

    db_host = "localhost"
    db_port = 27017

    if "-h" in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if "-p" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])

    check(db_host = db_host, db_port = db_port)