#!/usr/bin/env python

##
## <extSameEquivLig.py>
##

import sys
import os
import random
from datetime import datetime
import re
import math
import subprocess

LastModDate = "Sep 27, 2012"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        opt_dic[argv[i][1:]] = argv[i+1]



def read_annotation_file(ifname,dat):
#>> FILE FORMAT EXAMPLE <<
#>02454613
#PAU
#>02454944
#TMP T DT DRT 0DT
#>02454534
#DTD D1D

  print "#read_annotation_file('%s')"%(ifname)

  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING:can't open sdf file '%s'"%(ifname)
    return(0)
  f = open(ifname)
  id = ''
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('>')):
      id = line[1:]
    elif ((line.startswith('#')==0) and (len(line)>1) and (id != '')):
      dat[id] = line
  f.close()
  pass


##############
#### MAIN #####
###############

OPT = {}
OPT['ianA'] = ''
OPT['ianB'] = ''
OPT['oan'] = 'idlist.out'
if (len(sys.argv)<3):
  print "extSameEquivLig.py <options>"
  print " for extracting the same equivalent ligand molecules."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ianA    : input ID_A-vs-equivalent ID_P annotation file [%s]"%(OPT['ianA'])
  print " -ianB    : input ID_B-vs-equivalent ID_P annotation file [%s]"%(OPT['ianB'])
  sys.exit(1)

read_option(sys.argv,OPT)

### READ ANNOTATION FILE ###
print "#COMMAND %s"%(OPT['COMMAND'])
print "#DATE    %s"%(OPT['START_DATE'])

annotA = {}
read_annotation_file(OPT['ianA'],annotA)

annotB = {}
read_annotation_file(OPT['ianB'],annotB)

inv_annotA = {}
for x in (annotA.keys()):
  inv_annotA[annotA[x]] = x

inv_annotB = {}
for x in (annotB.keys()):
  inv_annotB[annotB[x]] = x

print "#[idA] [idB] [SameEquivLigands]"
for x in (annotA.keys()):
  if (inv_annotB.has_key(annotA[x])):
     print "%-12s %-12s    %s"%(x,inv_annotB[annotA[x]],annotA[x])


