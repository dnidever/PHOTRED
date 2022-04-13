# This is the main PHOTRED python program

# Need something that checks what "mode" we are in.
# Need a parser


# We are actually running something

# Load the list of stages

# Load the setup file

# Make sure that all of the stages actually exist and work

# Stages Loop
for stagename in stagelist:

    # Instantiate the stage object

    # Instantiate the lists object for this stage
    lists = photredListsObject(stagename)
    
    # call the jobs_daemon with the list if there's something to do
    if len(lists.inlist):
        photred.jobdaemon(lists)


# The photred jobs daemon
def jobdaemon(lists):

    # Initialize the stage object
    
    # Initalize jobs structure

    # While loop

        # Check status of running jobs

            # If some finished
            #  -check if the outputs are okay
            #  -update the lists
    
        # Give current status

        # Submit new jobs

        # Check if we are done
    
# Create a script for the jobs daemon
def job_makescript():
    
# This checks the status of the jobs
def job_checkstat(jobs):

    # Loop through the jobs
    for i in jobs:
    
# This reports the current status
def job_reportstatus(jobs):

# This submits new jobs
def job_submitnew
