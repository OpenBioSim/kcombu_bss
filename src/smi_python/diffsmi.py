#!/usr/bin/python
import sys
import os

LastModDate = "Sep 28, 2012"

if (len(sys.argv)<2):
  print "diffsmi.py [smaller_smi] [larger_smi] > [smi_not_in_smaller]"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  sys.exit(1)

ifileS = sys.argv[1]
ifileL = sys.argv[2]
 

### (1) read smaller_smi ##
if (os.access(ifileS,os.R_OK)==0):
  print "#ERROR:Can't open smaller_smi '%s'"(ifileS)
  sys.exit(1)
f = open(ifileS)
SmallerID = {}
for line in f:
  if (line.startswith('#')==0) and (len(line)>1):
    line = line.rstrip('\n')
    #print line
    fields = line.split() 
    smiles = fields[0]
    id     = fields[1]
    SmallerID[id] = len(smiles)
f.close()


#for id in (SmallerID.keys()): 
#  print id,SmallerID[id]

### (2) read larger_smi and output only smi not in smaller_smi##
if (os.access(ifileL,os.R_OK)==0):
  print "#ERROR:Can't open larger_smi '%s'"(ifileL)
  sys.exit(1)
f = open(ifileL)
for line in f:
  if (line.startswith('#')==0) and (len(line)>1):
    line = line.rstrip('\n')
    fields = line.split() 
    smiles = fields[0]
    id     = fields[1]
    if (SmallerID.has_key(id)==0):
      print line 
f.close()
