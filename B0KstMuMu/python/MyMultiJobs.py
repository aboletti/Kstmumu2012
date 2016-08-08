#!/usr/bin/env python

from collections import deque
import subprocess
import os
import time
import sys

loopy = 0
while loopy:
    haveRunningJobs = 0
    outpipe = subprocess.Popen(['Qstat | grep dinardo | wc -l'],
                               shell = True, stdout = subprocess.PIPE).stdout
    for line in outpipe.readlines():
        if int(line):
            haveRunningJobs = 1
    if not haveRunningJobs:
        break
    sys.stdout.write('Jobs still running, sleeping for few minutes.\n')
    sys.stdout.flush()
    time.sleep(600)

## Set the number of jobs to have in the batch queue at a time
## and the number of seconds to wait between submitting further
## jobs to the queue

nJobsToRun = 200
nSecToWait = 60
sample = 'BKMM_'

### GEN MC ###
from B0ToKstMuMu_GEN_NoFilter_MC_cff import readFiles
#from B0ToJPsiKst_GEN_NoFilter_MC_cff import readFiles
#from B0ToPsi2SKst_GEN_NoFilter_MC_cff import readFiles
files = readFiles
listStart = 0
listEnd = len(files)

## preamble is the Qsub command
preamble = "Qsub -l lnxfarm -e -o "
command = " cmsRun B0KstMuMu.py "

## List holding jobs. You can create multiple ones and then append
## their contents to the final queue, useful for doing multiple
## datasets that you may want to run together or separately

myjobs = []
job = listStart
for jobnum in files[job:listEnd]:
    jobstr = preamble + 'job' + str(job) + '.out -N ' + sample + str(job) + command + jobnum
    myjobs.append([jobstr, sample+str(job)])
    job += 1

alljobs = deque([])

for name in myjobs:
    alljobs.append(name)

jobq = deque([])

sys.stdout.write("Will submit jobs using the following commands:\n")
for job in alljobs:
    sys.stdout.write(job[0] + '\n')

sys.stdout.write('------\n')
sys.stdout.flush()
#sys.exit(1)

while len(alljobs):
    jobr = []
    for job in jobq:
        outpipe = subprocess.Popen(['Qstat | grep ' + job[1] + ' | wc -l'],
                                   shell = True, stdout = subprocess.PIPE).stdout
        for line in outpipe.readlines():
            if not int(line):
                jobr.append(job)

    for job in jobr:
        jobq.remove(job)

    if (len(jobq) < nJobsToRun) and len(alljobs):
        jobcmd = alljobs.popleft()
        sys.stdout.write("Submitting job: " + jobcmd[1] + "\n")
        sys.stdout.flush()
        outpipe = subprocess.Popen([jobcmd[0]], shell = True, stdout = subprocess.PIPE).stdout
        spl = ''
        for line in outpipe.readlines():
            spl = line.split(' ')
        jnum = spl[2]
        jobq.append([jobcmd[0], jnum, jobcmd[1]])

    time.sleep(nSecToWait)

sys.stdout.write("All jobs submitted.\n")
sys.stdout.flush()

print "Jobs currently running: "

loopy = 1
while loopy:
    haveRunningJobs = 0
    for job in jobq:
        outpipe = subprocess.Popen(['Qstat | grep ' + job[1] + ' | wc -l'],
                                   shell = True, stdout = subprocess.PIPE).stdout
        for line in outpipe.readlines():
            if int(line):
                haveRunningJobs = 1
    if not haveRunningJobs:
        loopy = 0
    time.sleep(120)

sys.stdout.write("All jobs finished.\n")
sys.stdout.flush()
