#!/usr/bin/env python
##
## <chkCanSMILES.py>
##


import sys
import os
import re
import math

LastModDate = "Apr 5, 2012"

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
OPT['idir'] = 'tmpout'
OPT['tdir'] = 'tmpout'
OPT['nr']   = 10 
OPT['ofu'] = 'uniq.out'
OPT['ofs'] = 'smiles.out'
OPT['A']    =  'F' 

if (len(sys.argv)<2):
  print "chkCanSMILES.py <options>"
  print " for checking canonical SMILES by randomizing atom orders."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir : input property file (tab-splited)  [%s]"%(OPT['idir']) 
  print " -tdir : output directory for temporary files [%s]"%(OPT['tdir']) 
  print " -ofu  : output file for uniquness of SMILES[%s]"%(OPT['ofu']) 
  print " -ofs  : output file for SMILES string[%s]"%(OPT['ofs']) 
  print " -nr   : number of random repeat [%d]"%(OPT['nr']) 
  print " -A    : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)

Nrand = int(OPT['nr'])
molfile_list = os.listdir(OPT['idir'])

if (OPT['A']=='T'):
  ofu = open(OPT['ofu'],'w')
  ofs = open(OPT['ofs'],'w')

for file in (molfile_list):
  file_full = OPT['idir'] + '/' + file
  print ">%s"%(file_full)
  ### (1) output nr-different SMILES by random reordering ###
  for r in range(Nrand):
    comstr = "ckcombu -A %s -ratm T -srd %d -omA %s/%d.smi > tmp.out"%(file_full,r,OPT['tdir'],r) 
    if (r==0):
      print "#%s"%(comstr)
    if (OPT['A']=='T'):
      os.system(comstr)
    else:
      print "#%s"%(comstr)
  ### (2) check consistencies among nr-different SMILES ###
  if (OPT['A']=='T'):
     smiles = ['' for r in range(Nrand)]
     Nobs_smiles = {}
     Nmiss = 0
     for r in range(Nrand):
       smifile = "%s/%d.smi"%(OPT['tdir'],r)
       if (os.access(smifile,os.R_OK)==1):
         f = open(smifile)
         for line in f:
           line = line.rstrip('\n')
           smiles[r]  = line
           Nobs_smiles[line] = Nobs_smiles.get(line,0) + 1 
       else:
         Nmiss += 1
     Nuniq_smiles = len(Nobs_smiles.keys()) 
     ofu.write("%s Nuniq_smiles %2d Nmiss %2d\n"%(file,Nuniq_smiles,Nmiss))
     print "%s Nuniq_smiles %d Nmiss %d\n"%(file,Nuniq_smiles,Nmiss)
     for s in (Nobs_smiles.keys()):
       ofs.write("%s %2d %s\n"%(file,Nobs_smiles[s],s)) 
  
     for r in range(Nrand):
       smifile = "%s/%d.smi"%(OPT['tdir'],r)
       os.system("rm %s"%(smifile))

if (OPT['A']=='T'):
  print "#write_uniqueness_file -->'%s'"%(OPT['ofu'])
  print "#write_SMILESs_file -->'%s'"%(OPT['ofs'])
  ofu.close()
  ofs.close()
