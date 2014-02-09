import commands, re, time, os, shutil

def init_proxy(virtual_organisation="vo.grand-est.fr", proxy_server = "myproxy.cern.ch", print_output=False):
    """
    This function init your voms-proxy. You will have to type your password twice:
    - the first one for the voms-proxy-init command
    - the second one for the myproxy-init command
    If you don't want to use a proxy_server, call this function with proxy_server = None
    """
    output = commands.getoutput("voms-proxy-init --voms %s"%virtual_organisation)
    if print_output:
        print output 
    if proxy_server:
        output = commands.getoutput("myproxy-init -v -s %s -d -n -t 48 -c 720 "%proxy_server)
        if print_output:
            print output    

def write_pyrna_job(directory, bash_script_content, uid, input_files = [], endpoint="https://sbgwms1.in2p3.fr:7443/glite_wms_wmproxy_server", virtual_organisation="vo.grand-est.fr", proxy_server= "myproxy.cern.ch", algorithms = None , python = None, submit = True):
    """
    Write a glite job which launches a python script based on the pyrna API:
    - directory: the directory where to package and submit the job
    - bash_script_content: the content to append to the bash script that will launch the python script
    - uid: a unique id used to create the name of the bash script, python script and jdl file
    - input_files: the files to be used as inputs for the python script
    - proxy_server: the name of the proxy server to be used for jobs taking time
    - submit: if True (default), the job is submitted. If False, the job file is just created.
    This function produces a file named submission_result in the directory precised in the first argument.
    """
    if not python:
        print("You have to define the absolute path of the python interpreter to be used")
        sys.exit(1)
    if not algorithms:
        print("You have to define the absolute path of the directory containing the algorithms to be used (and containing a file named setmyenv)")
        sys.exit(1)
    directory = os.path.realpath(directory)
    if not os.path.exists(directory):
        os.mkdir(directory)

    jobs_dir = os.path.join(directory, "jobs")

    if not os.path.exists(jobs_dir):
        os.mkdir(jobs_dir)

    input_file_names = []

    #we copy the input files and extract their names for the input_sandbox
    for input_file in input_files:
        file_name = input_file.split("/")[-1]
        shutil.copy2(input_file, directory+"/"+file_name)
        input_file_names.append(file_name)
    
    bash_script_name = directory+"/script_"+str(uid)+".sh"
    jdl_file_name = directory+"/jobs/job_"+str(uid)+".jdl"

    f = open(bash_script_name ,'w')
    f.write("#!/bin/bash\n")
    f.write("source %s/setmyenv\n"%algorithms) #this will add to the PATH all the binaries installed on the grid with the script files/scripts/shell/install_algorithms.sh
    f.write("export PATH=%s/bin:$PATH:${VO_VO_GRAND_EST_FR_SW_DIR}/public/bin\n"%python)
    f.write("hg clone https://fjossinet@bitbucket.org/fjossinet/pyrna\n")
    f.write("export PYTHONPATH=$PWD/pyrna:$PYTHONPATH\n")
    f.write("cd pyrna\n")
    f.write(bash_script_content)
    f.close()

    bash_script_name = "script_"+str(uid)+".sh"
    input_sandbox = [bash_script_name]
    input_sandbox += input_file_names

    create_jdl_file(file_name=jdl_file_name, executable="/bin/bash", arguments=[bash_script_name], input_sandbox=input_sandbox, virtual_organisation=virtual_organisation, proxy_server=proxy_server)

    if submit:
        submit_glite_job(directory=directory, jdl_file_name=jdl_file_name, endpoint=endpoint)
    
def submit_glite_job(directory, jdl_file_name, endpoint="https://sbgwms1.in2p3.fr:7443/glite_wms_wmproxy_server", submission_result_file_name="submission_result"):
    """
    Submit a glite job already packaged:
    - directory: the path for the directory containing the jdl file
    - jdl_file_name: the name of the jdl file to submit
    - endpoint: the endpoint for submission
    - use_proxy_server: to use a proxy server for jobs taking time
    This function produces a file in the directory containing the jdl file and storing the list of the jobs ids submitted. The name of this file is precised with the argument "submission_result_file_name". 
    If the file already exists, the new jobs ids are appended at the end of this file.
    """
    if not os.path.exists(directory):
        print "Cannot find %s"%directory
        return

    directory = os.path.realpath(directory)
    submission_output = commands.getoutput("cd %s ; glite-wms-job-submit -e %s -a %s"%(directory, endpoint, jdl_file_name))
    submission_result = open(directory+"/"+submission_result_file_name,'a')

    if re.search("Success",submission_output):
        for line in submission_output.split('\n'):
            if line.startswith("https"):
                submission_result.write(jdl_file_name+" "+line+"\n")
    else:
        submission_result.write(jdl_file_name+" Failed to submit\n")    
    submission_result.close()

def submit_collection_of_jobs(directory, endpoint="https://sbgwms1.in2p3.fr:7443/glite_wms_wmproxy_server", submission_result_file_name="submission_result"):
    directory = os.path.realpath(directory)
    submission_output = commands.getoutput("cd %s/.. ; glite-wms-job-submit -e %s -a --collection %s"%(directory, endpoint, directory))
    submission_result = open(directory+"/"+submission_result_file_name,'a')
    if re.search("Success",submission_output):
        for line in submission_output.split('\n'):
            if line.startswith("https"):
                submission_result.write(directory+" "+line+"\n")
    else:
        submission_result.write(directory+" Failed to submit\n")    
    submission_result.close()
    
def check_jobs_statuses(submission_output, serialize=False, suffix="", verbose = True):
    """
    Check the job statuses:
    - submission_output: the file storing the output of the jobs submission
    - serialize: if True, all statuses are saved in a file according to the status category
    - suffix: a suffix that will be added to the default file names if the job statuses are serialized
    """
    if not os.path.exists(submission_output):
        print "Cannot find %s"%submission_output
        return
    submission_dir =os.path.dirname(os.path.realpath(submission_output))
    f = open(submission_output,'r')
    running_jobs = 0
    scheduled_jobs = 0
    failed_jobs = 0
    aborted_jobs = 0
    succeeded_jobs = 0
    done_without_success = 0
    waiting_jobs = 0
    cancelled_jobs = 0
    cleared_jobs = 0
    unknown_status = 0
    ready_jobs = 0
    if serialize:
        aborted_jobs_output = open(submission_dir+"/aborted_jobs"+suffix,'w')
        succeeded_jobs_output = open(submission_dir+"/succeeded_jobs"+suffix,'w')
        done_without_success_jobs_output = open(submission_dir+"/done_without_success_jobs"+suffix,'w')
        failed_jobs_output = open(submission_dir+"/done_but_failed_jobs"+suffix,'w')
        cancelled_jobs_output = open(submission_dir+"/cancelled_jobs"+suffix,'w')
        unknown_status_output = open(submission_dir+"/unknown_status_jobs"+suffix,'w')
        running_status_output = open(submission_dir+"/running_jobs"+suffix,'w')
        cleared_jobs_output = open(submission_dir+"/cleared_jobs"+suffix,'w')
        scheduled_jobs_output = open(submission_dir+"/scheduled_jobs"+suffix,'w')
    for line in f.readlines():
        tokens = line.split(" ")
        output = commands.getoutput("glite-wms-job-status "+tokens[1].strip())
        if re.search("Current Status:.+Success",output):
            succeeded_jobs +=1
            if verbose:
                print "Success: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                succeeded_jobs_output.write(line)
        elif re.search("Current Status:.+Aborted",output):
            aborted_jobs +=1
            if verbose:
                print "!!! Aborted !!!: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                aborted_jobs_output.write(line)
        elif re.search("Current Status:.+Running",output):
            running_jobs +=1
            if verbose:
                print "Running: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                running_status_output.write(line)   
        elif re.search("Current Status:.+Scheduled",output):
            scheduled_jobs +=1
            if verbose:
                print "Scheduled: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                scheduled_jobs_output.write(line)
        elif re.search("Current Status:.+!=0",output):
            done_without_success +=1
            if verbose:
                print "!!! Done without success !!!: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                done_without_success_jobs_output.write(line)
        elif re.search("Current Status:.+(Failed)",output):
            failed_jobs +=1
            if verbose:
                print "!!! Done failed (no recovering available)  !!!: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                failed_jobs_output.write(line)
        elif re.search("Current Status:.+Cleared",output):
            cleared_jobs +=1
            if verbose:
                print "Already Cleared: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                cleared_jobs_output.write(line)
        elif re.search("Current Status:.+Cancelled",output):
            cancelled_jobs +=1
            if verbose:
                print "Cancelled: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                cancelled_jobs_output.write(line)
        elif re.search("Current Status:.+Waiting",output):
            waiting_jobs +=1
            if verbose:
                print "Waiting: "+tokens[0]+" "+tokens[1].strip()
        elif re.search("Current Status:.+Ready",output):
            ready_jobs +=1
            if verbose:
                print "Ready: "+tokens[0]+" "+tokens[1].strip()
        else:
            unknown_status +=1
            if verbose:
                print "Unknown Status: "+tokens[0]+" "+tokens[1].strip()
            if serialize:
                unknown_status_output.write(line)
    if verbose:
        print "Running: %s\nScheduled: %s\nReady: %s\nDone Failed: %s\nDone Without Success: %s\nAborted: %s\nSucceeded: %s\nWaiting: %s\nCancelled: %s\nCleared: %s\nUnknown Status: %s\n"%(str(running_jobs),str(scheduled_jobs),str(ready_jobs),str(failed_jobs),str(done_without_success), str(aborted_jobs),str(succeeded_jobs),str(waiting_jobs),str(cancelled_jobs),str(cleared_jobs), str(unknown_status))          
    f.close()
    if serialize:
        failed_jobs_output.close()
        aborted_jobs_output.close()
        succeeded_jobs_output.close()
        done_without_success_jobs_output.close()
        cancelled_jobs_output.close()
        unknown_status_output.close()
        running_status_output.close()
        cleared_jobs_output.close()
        scheduled_jobs_output.close()
    
def cancel_jobs(jobs_list):
    """
    Cancel a running job:
    - jobs_list: the file containing the list of the jobs to cancel
    """
    if not os.path.exists(jobs_list):
        print "Cannot find %s"%jobs_list
        return
    i=0
    f = open(jobs_list,'r') 
    for line in f.readlines():
        tokens = line.split(" ")     
        commands.getoutput("glite-wms-job-cancel "+tokens[1].strip()+' <<< "y"')
        i+=1
    f.close()
    print "%s job(s) cancelled."%(str(i))               

def resubmit_jobs(jobs_list, submission_result_file_name):
    """
    Resubmit a job:
    - jobs_list: the file containing the list of the jobs to resubmit
    """
    if not os.path.exists(jobs_list):
        print "Cannot find %s"%jobs_list
        return
    i=0
    f = open(jobs_list,'r') 
    for line in f.readlines():
        tokens = line.split(" ")     
        endpoint = "https://sbgwms1.in2p3.fr:7443/glite_wms_wmproxy_server"
        if tokens[1].strip().startswith("https://sbgwms2"):
            endpoint="https://sbgwms2.in2p3.fr:7443/glite_wms_wmproxy_server"
        submit_glite_job(os.path.dirname(os.path.realpath(jobs_list)),tokens[0].strip(),endpoint,submission_result_file_name)
        i+=1
        if i%20 == 0:
            time.sleep(60)#we wait 60 seconds before to resubmit the next 20 jobs

    f.close()
    print "%s job(s) resubmitted."%(str(i))             
    
def recover_job_outputs(jobs_list, output_dir):
    """
    Recover the results of jobs:
    - jobs_list: the file containing the list of the jobs to recover
    - output_dir: the directory where to store the outputs
    """
    if not os.path.exists(jobs_list):
        print "Cannot find %s"%jobs_list
        return
    if not os.path.exists(output_dir):
        print "Cannot find %s"%output_dir
        return
    output_dir = os.path.realpath(output_dir)
    i=0
    f = open(jobs_list,'r') 
    for line in f.readlines():
        tokens = line.split(" ")
        id = tokens[1].strip().split("/")[-1]
        already_recovered = False
        for _f in os.listdir(output_dir):
            if _f.endswith(id):
                already_recovered = True
                print "Seems already recovered : %s"%line
                break
        if not already_recovered:    
            commands.getoutput("glite-wms-job-output --dir "+output_dir+" "+tokens[1].strip())
            for _f in os.listdir(output_dir):
                if _f.endswith(id):
                    i+=1
                    if i%50 == 0:
                        time.sleep(5)#we wait 5 seconds before to recover the next 50 jobs
                    break                               
    f.close()
    print "%s job output(s) recovered in %s."%(str(i),output_dir)               

def create_jdl_file(file_name, executable, arguments, input_sandbox, virtual_organisation, job_type="Normal", output_sandbox=["stdout.out","stderr.err"],  stdoutput="stdout.out", stderror="stderr.err", proxy_server=None):
    f = open(file_name,'w')
    f.write('JobType = "'+job_type+'";\n')
    f.write('Executable = "'+executable+'";\n')
    args = ""
    for argument in arguments:
        args += argument+" "
    f.write('Arguments = "'+args.strip()+'";\n')
    f.write('StdOutput = "'+stdoutput+'";\n')
    f.write('StdError = "'+stderror+'";\n')

    _input_sandbox = "{"
    for isb in input_sandbox:
        _input_sandbox += '"'+isb+'",'
    _input_sandbox = _input_sandbox[0:len(_input_sandbox)-1]+"}"    
    f.write('InputSandbox = '+_input_sandbox+';\n')

    _output_sandbox = "{"
    for osb in output_sandbox:
        _output_sandbox += '"'+osb+'",'
    _output_sandbox = _output_sandbox[0:len(_output_sandbox)-1]+"}" 
    f.write('OutputSandbox = '+_output_sandbox+';\n')
    f.write('VirtualOrganisation = "'+virtual_organisation+'";\n')
    if proxy_server:
        f.write('MyProxyServer = "%s";\n'%proxy_server)
    f.close()
    
