import glite
import shutil, os
from pymongo import MongoClient

class Task:

    def __init__(self, db_name, db_host = "localhost", db_port = 27017, endPoint = None, algorithms = None , python = None):
        self.jobsToSubmit = 0
        self.jobsSubmitted = 0
        self.db_name = db_name
        self.endPoint = endPoint
        self.client = MongoClient(db_host, db_port)
        self.algorithms = algorithms
        self.python = python
        if db_name:
            self.db = self.client[db_name]

    def getTotalJobsToSubmit(self, data):
        pass

    def getScriptContent(self, job_id):
        pass

    def storeData(self):
        pass

    def doTheJob(self, job_id = 1):
        pass 

    def submitJobs(self, script_name):
        totalJobs = self.getTotalJobsToSubmit(self.storeData())
        if not self.endPoint:
            self.doTheJob(None)
        else:
            #glite.init_proxy()
            outputDir = os.getenv("HOME")+"/tmp/jobs_%s_on_%s"%(script_name, self.db_name)
            if not os.path.exists(outputDir):
                shutil.os.mkdir(outputDir)
            info = open("%s/info"%outputDir, 'a+')
            info.write("total jobs: %i\n"%totalJobs)
            info.close()

            for i in range(1, totalJobs+1):
                glite.write_pyrna_job(outputDir, self.getScriptContent(i), i, [], self.endPoint, algorithms = self.algorithms, python = self.python)

            if False: #submit jobs as a collection
                for i in range(1, totalJobs+1):
                    glite.write_pyrna_job(outputDir, self.getScriptContent(i), i, [], self.endPoint, algorithms = self.algorithms, python = self.python, submit=False) #submit=False means that we just create the job files

                #now we submit the collection of jobs
                glite.submit_collection_of_jobs(directory=outputDir+"/jobs", endpoint=self.endPoint)

