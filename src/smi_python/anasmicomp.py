#!/usr/bin/env python
##
## <anasmicomp.py>
##


import sys
import os
import re
import math

LastModDate = "July 10, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_smiles_file(ifname,IDlist):
#>> FILE EXAMPLE <<
#CCOC(=O)C2(c1ccc(c(c1)))CCN(C(C2))CCc3ccc(c(c3))N D02942.sdf
#CCCCN(CCCC)CCCOc(c(c1))ccc1C(=O)C(=C(O2)CCCC)c3cc(cc(c32))NS(=O)(=O)C D03914.sdf
#OCC(C(C2(F)F)O)OC2N(C1=O)C=CC(=N1)N D02368.sdf
#CC(O)(CC(C52))C5(C)CCC3C4(C)CC(=C1C(C4C(C(C32))))C=N(N1) D00444.sdf
  print "#read_smiles_file('%s')"%(ifname)
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open smiles file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    if (line.startswith('#')==0) and (len(line)>5):
      line = line.rstrip('\n')
      (smiles,id) = line.split(' ')
      if (id.endswith('.sdf')):
        id = id.replace('.sdf','')
     # SMILESlist.append(smiles)
      IDlist.append(id)
  f.close()


def read_corresponding_ID_file(ifname,simID_A,simID_B):
#>> FILE EXAMPLE <<
# >04272959
# ZINC12113654
# >04274705
# ZINC19235523 ZINC19235522
# >04274873
# ZINC19563680 ZINC19563677 ZINC19563680 ZINC19563677

  print "#read_corresponding_ID_file('%s')"%(ifname)
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open corresponding ID file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  idA = '' 
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0) and (len(line)>5):
      if (line.startswith('>')):
        idA = line[1:]
      elif (idA != ''):
        field = line.split(' ')
        if (simID_A.has_key(idA)==0):
          simID_A[idA] = [] 

        for idB in (field):
          simID_A[idA].append(idB)
          if (simID_B.has_key(idB)==0):
            simID_B[idB] = [] 
          simID_B[idB].append(idB)

  f.close()




###############
#### MAIN #####
###############

OPT = {}
OPT['isA'] = ''
OPT['isB'] = ''
OPT['ofA'] = 'A.out'
OPT['ofB'] = 'B.out'

if (len(sys.argv)<2):
  print "anasmicomp.py <options>"
  print " for analyzing SMILES comparison"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -isA  : input SMILES string file for molsA  [%s]"%(OPT['isA'])
  print " -isB  : input SMILES string file for molsB  [%s]"%(OPT['isB'])
  print " -iA2B : input corresponding ID for molsA to molsB  [%s]"%(OPT['isB'])
  print " -ofA  : summary file for molsA [%s]"%(OPT['ofA'])
  print " -ofB  : summary file for molsA [%s]"%(OPT['ofB'])
  sys.exit(1)

read_option(sys.argv,OPT)

IDlistA = []
IDlistB = []

read_smiles_file(OPT['isA'],IDlistA)
NidA = len(IDlistA)
print "#len(IDlistA):%d"%(NidA)
read_smiles_file(OPT['isB'],IDlistB)
NidB = len(IDlistB)
print "#len(IDlistB):%d"%(NidB)

simID_A = {}
simID_B = {}
read_corresponding_ID_file(OPT['iA2B'],simID_A,simID_B)

### CALCULATE Number of corresponding IDs ###
NidA_corr =  {} 
idA_corr_exp = {}
for idA in (IDlistA):
  if (simID_A.has_key(idA)):
    n = len(simID_A[idA])
  else:
    n = 0
    print idA
  NidA_corr[n] = NidA_corr.get(n,0) + 1 
  if (idA_corr_exp.has_key(n)==0):
    idA_corr_exp[n] = []
  if (NidA_corr[n]<=5):
    idA_corr_exp[n].append(idA)

NidB_corr =  {} 
idB_corr_exp = {}
for idB in (IDlistB):
  if (simID_B.has_key(idB)):
    n = len(simID_B[idB])
  else:
    n = 0
  NidB_corr[n] = NidB_corr.get(n,0) + 1 
  if (idB_corr_exp.has_key(n)==0):
    idB_corr_exp[n] = []
  if (NidB_corr[n]<=5):
    idB_corr_exp[n].append(idB)

### OUTPUT SUMMARY for molsA ###
nlist = sorted(NidA_corr.keys(),lambda x,y:cmp(x,y))
of = open(OPT['ofA'],'w')
print "#write_summary() -->'%s'"%(OPT['ofA'])
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#NidA %d NidB %d\n"%(NidA,NidB))
for n in (nlist):
  of.write("%-2d %-8d %7.3f %%"%(n,NidA_corr[n],100.0*float(NidA_corr[n])/NidA))
  for x in (idA_corr_exp[n]):
    of.write(" %s"%(x))
  of.write("\n")
of.close()

### OUTPUT SUMMARY for molsB ###

nlist = sorted(NidB_corr.keys(),lambda x,y:cmp(x,y))
of = open(OPT['ofB'],'w')
print "#write_summary() -->'%s'"%(OPT['ofB'])
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#NidA %d NidB %d\n"%(NidA,NidB))
for n in (nlist):
  of.write("%-2d %-8d %7.3f %%"%(n,NidB_corr[n],100.0*float(NidB_corr[n])/NidB))
  for x in (idB_corr_exp[n]):
    of.write(" %s"%(x))
  of.write("\n")
of.close()
of.close()
