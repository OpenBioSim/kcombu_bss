#!/usr/bin/env python
import sys
import os
from datetime import datetime
import socket

LastModDate = "Apr 10, 2014"


def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


############
### MAIN ###
############

OPT = {}
OPT['A'] = 'F'
OPT['C'] = '1'

PID = os.getpid()
HOSTNAME = socket.gethostname()

if (len(sys.argv)<2):
  print "submit_cores.py [command_submit] <options>"
  print " for submit [C] nohup jobs with different -c [c] option."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -A       : Action (T or F) [%s]"%(OPT['A'])
  print " -C       : number of cores in one nodes [%s]"%(OPT['C'])
  sys.exit(1)

read_option(sys.argv,OPT)
command_submit = sys.argv[1] 

C = int(OPT['C'])

for c in range(C):
  command = 'nohup %s -c %d > %s_%d.log &'%(command_submit,c,HOSTNAME,c)
  if (OPT['A']=='T'):
    print "#action:%s"%(command)
    os.system(command)
  else:
    print "#no_action:%s"%(command)
