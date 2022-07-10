#!/usr/bin/env python
##
## <propana.py>
##


import sys
import os
import re
import math

LastModDate = "Mar 9, 2012"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['if'] = ''
OPT['c'] = '0'
OPT['qc'] = '0'
OPT['qs'] = ''
OPT['oql'] = ''

if (len(sys.argv)<2):
  print "propana.py <options>"
  print " for analyzing property data file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -if  : input property file (tab-splited)  [%s]"%(OPT['if']) 
  print " -A   : Action (T or F) [%s]"%(OPT['A']) 
  print " -c   : column number for focused property (1,2,...)[%s]"%(OPT['c'])
  print " -qc  : column number for query property (1,2,...)[%s]"%(OPT['qc'])
  print " -qs  : query string of property (1,2,...)[%s]"%(OPT['qs'])
  print " -oql : output file query column observing count [%s]"%(OPT['oql'])
  sys.exit(1)

read_option(sys.argv,OPT)
OPT['qc'] = int(OPT['qc']) - 1
OPT['c']  = int(OPT['c'])  - 1

S = SS = 0.0
N = 0
ifname = sys.argv[1]
Nqcount = {}
f = open(OPT['if'])
### (1) Read file line one by one ###
for line in f:
  line = line.rstrip('\n')
  if (len(line)>10) and (line.startswith('#')==0):
    field = line.split('\t')
    for i in range(len(field)):
      if (field[i].find(':')>=0):
        items = field[i].split(':')
        field[i] = items 
    #print field
    #print line
    ### (1) Check Query Pattern ###
    hit = 0
    if (OPT['qs']==''):
     hit = 1
    else:
      if (isinstance(field[OPT['qc']],str))  and (field[OPT['qc']]==OPT['qs']): 
        hit = 1
      if (isinstance(field[OPT['qc']],list)) and (field[OPT['qc']].count(OPT['qs'])>0): 
        hit = 1
        if (hit==1) and (OPT['oql']!=''):
          for x in (field[OPT['qc']]):
            Nqcount[x] = Nqcount.get(x,0) + 1
    ### (2) Check focused propert ##
    if (hit==1):
      #  print line 
      val = float(field[OPT['c']])
      S  += val
      SS += val*val
      N += 1
f.close()


### OUTPUT ###
print "S %f SS %f"%(S,SS)
M = SD = 0.0
if (N>0):
  M = S/N
  SD = math.sqrt(SS/N - M * M)
print "#N %d M %f SD %f"%(N,M,SD)


if (OPT['oql']!=''):
 print "#write_query_string_histogram -->'%s'"%(OPT['oql'])
 of = open(OPT['oql'],'w')
 items = sorted(Nqcount.keys(),lambda y,x:cmp(Nqcount[x],Nqcount[y]))
 for x in (items):
   of.write("%d %s\n"%(Nqcount[x],x))
 of.close()

