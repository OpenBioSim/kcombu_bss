#!/usr/bin/env python
##
## <propsearch.py>
##


import sys
import os
import re
import math

LastModDate = "Apr 14, 2011"

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


def read_smiles_file(ifname,SMILESlist,IDlist):
#>> FILE EXAMPLE <<
#CCOC(=O)C2(c1ccc(c(c1)))CCN(C(C2))CCc3ccc(c(c3))N D02942.sdf
#CCCCN(CCCC)CCCOc(c(c1))ccc1C(=O)C(=C(O2)CCCC)c3cc(cc(c32))NS(=O)(=O)C D03914.sdf
#OCC(C(C2(F)F)O)OC2N(C1=O)C=CC(=N1)N D02368.sdf
#CC(O)(CC(C52))C5(C)CCC3C4(C)CC(=C1C(C4C(C(C32))))C=N(N1) D00444.sdf
  print "#read_smiles_file('%s')"%(ifname)
  f = open(ifname)
  for line in f:
    if (line.startswith('#')==0) and (len(line)>5):
      line = line.rstrip('\n')
      (smiles,id) = line.split(' ')
      if (id.endswith('.sdf')):
        id = id.replace('.sdf','')
      SMILESlist.append(smiles)
      IDlist.append(id)
  f.close()



###############
#### MAIN #####
###############

OPT = {}
OPT['if'] = ''
OPT['M']  = 'A'

if (len(sys.argv)<2):
  print "smisearch.py <options>"
  print " for searching SMILES strings"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -M   : 'A'll-vs-all 'S'earch     [%s]"%(OPT['M']) 
  print " -if  : input SMILES string file  [%s]"%(OPT['if']) 
  sys.exit(1)

read_option(sys.argv,OPT)

SMILESlist = []
IDlist = []

read_smiles_file(OPT['if'],SMILESlist,IDlist)
Nmol = len(SMILESlist)
print "#Nmolecule %d"%(Nmol)


if (OPT['M']=='a'): ## naive all-vs-all
  SAME_SMILES_IDs = {}
  for i in range(Nmol):
    for j in range(i+1,Nmol):
      if (SMILESlist[i]==SMILESlist[j]):
        #print "%s %s %s %s"%(SMILESlist[i],SMILESlist[j],IDlist[i],IDlist[j])
        IDi = IDlist[i]
        IDj = IDlist[j]
        if (SAME_SMILES_IDs.has_key(IDi)==0):
          SAME_SMILES_IDs[IDi] = []
        if (SAME_SMILES_IDs.has_key(IDj)==0):
          SAME_SMILES_IDs[IDj] = []
        SAME_SMILES_IDs[IDi].append(IDj)
        SAME_SMILES_IDs[IDj].append(IDi)
           
  for ID in (SAME_SMILES_IDs.keys()):
    print ">%s"%(ID)
    for i in range(len(SAME_SMILES_IDs[ID])):
      if (i>0):
        sys.stdout.write(" ")
      sys.stdout.write("%s"%(SAME_SMILES_IDs[ID][i]))
    sys.stdout.write("\n")

if (OPT['M']=='A'):
  LEN_SMILES_HASH = {}
  for i in range(Nmol):
    smi = SMILESlist[i]
    id  = IDlist[i]
    L = len(smi)
    if (LEN_SMILES_HASH.has_key(L)==0):
      LEN_SMILES_HASH[L] = []
    LEN_SMILES_HASH[L].append(i)
 
  SAME_SMILES_IDs = {}
  for i in range(Nmol):
    SMi = SMILESlist[i]
    IDi = IDlist[i]
    for j in (LEN_SMILES_HASH[len(SMi)]):
      if (i!=j):
        SMj = SMILESlist[j]
        IDj = IDlist[j]
        if (SMi==SMj):
          #print "%s %s %s %s"%(SMILESlist[i],SMILESlist[j],IDlist[i],IDlist[j])
          if (SAME_SMILES_IDs.has_key(IDi)==0):
            SAME_SMILES_IDs[IDi] = []
          SAME_SMILES_IDs[IDi].append(IDj)
           
  for ID in (SAME_SMILES_IDs.keys()):
    print ">%s"%(ID)
    for i in range(len(SAME_SMILES_IDs[ID])):
      if (i>0):
        sys.stdout.write(" ")
      sys.stdout.write("%s"%(SAME_SMILES_IDs[ID][i]))
    sys.stdout.write("\n")

