#!/usr/bin/env python

##
## <sim2anno.py>
##

import sys
import os
import random
from datetime import datetime
import re
import math
import subprocess

LastModDate = "July 25, 2012"

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


def read_similar_pairs_file(ifname,thre,simlistA,simlistB,RemoveHeadDir,RemoveTail):
  if (not os.access(ifname,os.R_OK)):
    print "#ERROR:Can't oepn similar_pair_file '%s'"%(ifname)
    sys.exit(1)  

#>> FORMAT EXAMPLE <<
#376 QTR 6612 D03703.mol 1.000000
#376 QTR 6746 D00001.mol 1.000000
#378 TAR 3180 D00103.mol 1.000000
#380 AG7 6341 D05776.mol 1.000000
#390 QPS 5236 D00216.mol 1.000000
#:
#:
#7 3NA.sdf 3251145 D00000000/D03803.sdf 1.000000
#9 QBT.sdf 11136 02450000/02454944.sdf 1.000000
#9 QBT.sdf 2213840 03305000/03305531.sdf 1.000000
#11 GVE.sdf 3169428 00090000/00094821.sdf 1.000000
#19 N5N.sdf 3250116 D00000000/D03728.sdf 1.000000
#20 QUE.sdf 486640 D00005000/D08112.sdf 1.000000
#20 QUE.sdf 493247 00385000/00386087.sdf 1.000000
#:

  print "#read_similar_pairs_file('%s')"%(ifname)
  f = open(ifname)
  for line in f:
    if (line.startswith('#')==0) and (len(line)>5):
      line = line.rstrip('\n')
      field = line.split() 
      nameA = field[1]
      nameB = field[3]

      if (RemoveTail== 'T') and (nameA.find('.')>=0):
        (head,tail) = nameA.split('.')
        nameA = head

      if (RemoveTail== 'T') and (nameB.find('.')>=0):
        (head,tail) = nameB.split('.')
        nameB = head

      if (RemoveHeadDir== 'T') and (nameA.find('/')>=0):
        (head,tail) = nameA.split('/')
        nameA = tail 

      if (RemoveHeadDir== 'T') and (nameB.find('/')>=0):
        (head,tail) = nameB.split('/')
        nameB = tail 

      sc    = float(field[4])
      if (sc >= thre):
        if (simlistA.has_key(nameA)==0):
          simlistA[nameA] = [] 
        simlistA[nameA].append(nameB)
        if (simlistB.has_key(nameB)==0):
          simlistB[nameB] = [] 
        simlistB[nameB].append(nameA)

  f.close()



##############
#### MAIN #####
###############

OPT = {}
OPT['ism'] = ''
OPT['th'] = '1.0'
OPT['oanA'] = 'annoA.out'
OPT['oanB'] = 'annoB.out'
OPT['rmtail'] = 'T'
OPT['rmhead'] = 'T'

if (len(sys.argv)<3):
  print "sim2anno.py <options>"
  print " for adding new annotations into SDF files."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ism    : input similar compound pair list file [%s]"%(OPT['ism'])
  print " -th     : threshold for similarity [%s]"%(OPT['th'])
  print " -oanA   : output annotation file for group A[%s]"%(OPT['oanA']) 
  print " -oanB   : output annotation file for group B[%s]"%(OPT['oanB']) 
  print " -rmhead : remove head of molecular name, such as 'D00005000/', '01985000/'(T or F)[%s]"%(OPT['rmhead']) 
  print " -rmtail : remove tail of molecular name, such as '.mol2','.sdf' (T or F)[%s]"%(OPT['rmtail']) 
  sys.exit(1)

read_option(sys.argv,OPT)

### READ ANNOTATION FILE ##
simlistA = {}
simlistB = {}
read_similar_pairs_file(OPT['ism'],float(OPT['th']),simlistA,simlistB,OPT['rmhead'],OPT['rmtail'])


### OUTPUT ANNOTATION FILE ###

if (OPT['oanA'] != ''):
  print "#write_annotation_for_A()-->'%s'"%(OPT['oanA'])
  of = open(OPT['oanA'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#COMMAND %s\n"%(OPT['START_DATE']))
  for A in (simlistA.keys()):
    of.write(">%s\n"%(A))
    init = 1 
    for B in (simlistA[A]):
      if (init==0):
        of.write(" ")
      of.write("%s"%(B))
      init = 0
    of.write("\n")
  of.close() 


if (OPT['oanB'] != ''):
  print "#write_annotation_for_B()-->'%s'"%(OPT['oanB'])
  of = open(OPT['oanB'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#COMMAND %s\n"%(OPT['START_DATE']))
  for B in (simlistB.keys()):
    of.write(">%s\n"%(B))
    init = 1 
    for A in (simlistB[B]):
      if (init==0):
        of.write(" ")
      of.write("%s"%(A))
      init = 0
    of.write("\n")
  of.close() 


