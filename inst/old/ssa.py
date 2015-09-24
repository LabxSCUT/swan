#!/usr/bin/env python
# submit jobs in $1 while total cores<300, total mem<300G
# current limit is 299
# current mem limit is 300G

import os, shutil, subprocess, argparse, fnmatch, time, re
import numpy as np
core_max=50
mem_max=200
sleep_time=10
short_time=5
uname="lixia"
qname="main"

def mem_size(mem):
  #need to match pattern here
  if mem[-2:]=='--':
    return 0
  if mem[-2:]=='mb': 
    return int(mem[:-2])/1000
  elif mem[-1:]=='m':
    return int(mem[:-1])/1000
  elif mem[-2:]=='gb':
    return int(mem[:-2])
  else:
    return 0

def core_num(core):
  if core[-2:]=='--':
    return 1
  else:
    return(int(core))

def peek_pbs(pbsfile):
  pbs_content=open(pbsfile).read() 
  pat=re.compile("nodes=(?P<nodes>\d+):ppn=(?P<ppn>\d+)")
  m=re.search(pat, pbs_content)
  core=int(m.group('ppn'))*int(m.group('nodes'))
  pat=re.compile("mem=(?P<mem>\d+[mg]b),vmem=(?P<vmem>\d+[mg]b),pmem=(?P<pmem>\d+[mg]b)")
  m=re.search(pat, pbs_content)
  mem=max(mem_size(m.group('mem')),mem_size(m.group('pmem')),mem_size(m.group('vmem')))
  return np.array([core,mem])

def peek_current(uname, core_max, mem_max):
  tmp=subprocess.Popen("qstat -u %s | awk '{print $5,$6,$8}'" % uname, shell=True, stdout=subprocess.PIPE).communicate()
  cores=tmp[0].split('\n')
  sessions=tmp[0].split('\n')
  hskip=5
  tskip=1
  raw_table=[ line for line in tmp[0].rstrip('\n').split('\n') ]
  run_table=raw_table[hskip+1:len(raw_table)-tskip]
  mems=[ mem_size(line.split(' ')[2]) for line in run_table ]
  cores=[ core_num(line.split(' ')[1]) for line in run_table ]
  jobs=[ line.split(' ')[0] for line in run_table ]
  return np.array([core_max-sum(cores), mem_max-sum(mems)])

def submit( pbsFile ):
  tmp=subprocess.Popen("qsub -q %s %s" % (qname, pbsFile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
  submitted=True
  if tmp[1] == '':
    print tmp[0]
    print pbsFile, "submitted"
  else:
    print tmp[1]
    print pbsFile, "error"
    submitted=False
  return submitted

def main():

  parser = argparse.ArgumentParser(description="MCB Queue Checking and Submission Tool")
  parser.add_argument("pbsFile", metavar="pbsFile", help="single/multiple pbs file to be submitted")
  arg_namespace = parser.parse_args()
  pbsFile = vars(arg_namespace)['pbsFile']
  pbsFiles = set(pbsFile.split(","))
  print len(pbsFiles)

  while( len(pbsFiles) != 0 ):
    open_up=peek_current(uname, core_max, mem_max)
    #print "found open_up", open_up
    #if not all(v>0 for v in open_up):  #full
    remain=open_up
    while( len(pbsFiles) != 0 ):
      pbsFile=pbsFiles.pop()
      use_up=peek_pbs(pbsFile)
      print "next pbs use", use_up
      print "current remain", remain
      if not all(v>=0 for v in remain-use_up):
        print "not enough open_up"
        pbsFiles.add(pbsFile) #impossible, put pbs back, continue peek and sleep
        time.sleep(sleep_time) #cool down
        break
      if not submit(pbsFile):
        pbsFiles.add(pbsFile) #possible, but submit failed, put pbs bak, continue peek and sleep
        time.sleep(sleep_time) #cool down
        break
      else:
        remain=remain-use_up  
        print "added job", pbsFile 
        print "remain open_up", open_up
        time.sleep(short_time)
    if( len(pbsFiles) != 0): #if some job remains, need to cool down, else quit the routine
      print "cool down", sleep_time
      time.sleep(sleep_time)

if __name__=="__main__":
  main()
