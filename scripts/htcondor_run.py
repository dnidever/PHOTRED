import htcondor
import classad
import sys

def err_msg_exit(msg):
  sys.exit("ERROR!! " + msg)

# Get the HTCondor command
try:
  oper     = sys.argv[1]
except:
  err_msg_exit("Missing first parameter with HTCondor command")

# Get the Schedd (we will not use try/except since we want to 
# receive error in stderr when problems getting the Schedd)
schedd = htcondor.Schedd()

#########################
# JOBS SUMISSION
#########################
if oper == 'condor_submit':

  # Get parameters: number of jobs and base name of jobs (extension .sh will be added)
  try:
    N             = int(sys.argv[2])
    executable    = sys.argv[3]
    htcondor_cmd  = sys.argv[4]
  except:
    err_msg_exit("Syntax is not valid. Use: %s %s <num_jobs> <executable> <initdir> <inputs> <outputs> <environment> <definevars>" 
                 % (sys.argv[0], sys.argv[1]))

  # Add the standard submission commands
  sub = htcondor.Submit()
  sub['FNAME'] = executable.replace("$(Process)","")
  sub['ID']   = "$(Cluster).$(Process)"

  sub['output']  = "$(FNAME).$(ID).out"
  sub['error']   = "$(FNAME).$(ID).err"
  sub['log']     = "$(FNAME).$(Cluster).log"

  sub['universe']                = 'vanilla'
  sub['should_transfer_files']   = 'YES' 
  sub['when_to_transfer_output'] = 'ON_EXIT'

  # Check if we have to include some HTCondor Submit commands and/or define some variables
  # We recieve them in format "name1|=|value1|;|name2|=|value2|;|...|;|nameN|=|valueN" 
  # Commands beginning with # or with no or several assignations (|=|) will be ignored
  if htcondor_cmd != "":
    cmdlist = htcondor_cmd.split("|;|")
    for cmd in cmdlist:
      cmd_data = cmd.split("|=|")
      cmd_name = cmd_data[0].strip()
      if cmd_name and len(cmd_data) == 2 and not cmd_name.startswith('#'):
        sub[cmd_name] = cmd_data[1].strip()

  # Set executable (mandatory!!)
  sub['executable'] = executable
  #print ("=== START ===\n%s\n=== END ===" % sub)
  

  # Submit!!
  try:
    with schedd.transaction() as txn:  # txn will now represent the transaction.
      clusterId = sub.queue(txn, N)
  except:
    clusterId = 0

  print (clusterId)
  # Check submission
  #if not clusterId:
  #  err_msg_exit("There was an error submitting jobs")


#########################
# CHECKING QUEUE
#########################
elif oper == 'condor_q':

  # Get parameters: clusterId
  try:
    clusterId = int(sys.argv[2])
  except:
    err_msg_exit("Syntax is not valid. Use: %s %s <clusterId>" % (sys.argv[0], sys.argv[1]))

  # Get all jobs of clusterID (we will only want JobStatus info)
  jobs = schedd.xquery(requirements = 'ClusterId == %d' % clusterId, projection = ['JobStatus'])

  # Get total number of jobs and number of jobs of each job status
  jobs_status = {}
  total_jobs = 0
  for j in jobs:
    total_jobs +=1
    if not j['JobStatus'] in jobs_status:
      jobs_status[j['JobStatus']]  = 1
    else:
      jobs_status[j['JobStatus']] += 1

  # Build output: total_jobs;statusCode:num_jobs;statusCode:num_jobs...
  out = str(total_jobs)
  for key in jobs_status:
    out += ";" + str(key) + ":" + str(jobs_status[key])
  print(out)


#########################
# REMOVE JOBS
#########################
elif oper == 'condor_rm':

  # Get parameters: clusterId
  try:
    clusterId = int(sys.argv[2])
  except:
    err_msg_exit("Syntax is not valid. Use: %s %s <clusterId>" % (sys.argv[0], sys.argv[1]))

  # Remove all jobs
  schedd.act(htcondor.JobAction.Remove, 'ClusterId == %d' % clusterId)
  

#########################
# NON-VALID OPERATION
#########################
else:
  err_msg_exit("Operation is NOT valid. Valid operations: condor_submit, condor_q, condor_rm")
    
    
