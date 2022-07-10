#!/usr/bin/env python
##
## <propsearch.py>
##


import sys
import os
import re
import math

LastModDate = "Aug 6, 2012"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def isfloat(str):
  for i in range(len(str)):
    if (str[i].isdigit()):
      return(1)
  return(0)


###############
#### MAIN #####
###############

OPT = {}
OPT['if'] = ''
OPT['qn0'] = ''
OPT['qi0'] = ''
OPT['qx0'] = ''
OPT['qs0'] = ''

if (len(sys.argv)<2):
  print "propsearch.py [property_file] <options>"
  print " for analyzing property data file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -if  : input property file (tab-splited)  [%s]"%(OPT['if']) 
  print " -qn0 : property id number (1,2,...) [%s]"%(OPT['qn0']) 
  print " -qi0 : property minimum value [%s]"%(OPT['qi0']) 
  print " -qx0 : property maximum value [%s]"%(OPT['qx0']) 
  print " -qs0 : property substring [%s]"%(OPT['qs0']) 
  sys.exit(1)

read_option(sys.argv,OPT)
ipropfile = sys.argv[1]
qn0 = int(OPT['qn0']) - 1
ifname = sys.argv[1]
f = open(ipropfile)
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
    hit = 1
    #print "qn0 '%s' len %d\n"%(OPT['qn0'],len(field))
    if (isfloat(field[qn0])==0):
      hit = 0
    if (OPT['qn0'] !='') and (len(field)>qn0) and (isfloat(field[qn0])==1):
      val = field[qn0]
      #print "val %s %s %s"%(val,OPT['pi0'],OPT['px0']) 
      if (OPT['qi0'] != '') and (float(val)<float(OPT['qi0'])):
        hit = 0 
      if (OPT['qx0'] != '') and (float(val)>float(OPT['qx0'])):
        hit = 0 
      if (OPT['qs0'] != ''):
#        print "'%s' '%s' %d\n"%(val,OPT['ps0'],val.find(OPT['ps0']))
        if (isinstance(val,str))  and (val.find(OPT['qs0'])==-1): 
          hit = 0
        if (isinstance(val,list)):
          for x in (val):
            if (x.find(OPT['qs0'])==-1):
              hit = 0 
    if (hit==1):
      print line            
f.close()

