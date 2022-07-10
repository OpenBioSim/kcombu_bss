#!/usr/bin/env python
import sys
import os
import math
from datetime import datetime

LastModDate = "Apr 22, 2014"


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



def read_list_file(ifname,list):
  print "#read_list_file('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     list.append(fields[0])
  f.close()



def perform_kcombu_and_get_key_values(command,valkey):
#>>Example of -KV output of pkcombu<<
# NHEAVYATOM_A 28
# NHEAVYATOM_B 26
# NATOM_MATCH  15
# TANIMOTO     0.384615
# RMSD         0.194683
  f = os.popen(command)
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#')==0) and (len(field)>=2):
      key   = field[0]
      value = field[1]
      valkey[key] = value
    pass
  f.close()


############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'
OPT['A'] = 'F'
OPT['idtars']  = 'SDF_3D'
OPT['idrefp'] = 'SupLIG'
OPT['of']     = '-'


if (len(sys.argv)<2):
  print "sum_molprop.py <options>"
  print " for making summary of accuracy of predicted simimilariy-based 3D model."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idtars : input dir for target molecules in SDF[%s]"%(OPT['idtars'])
  print " -idrefp : input dir for reference 3D model in PDB[%s]"%(OPT['idrefp'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print " -of     : output summary file [%s]"%(OPT['of']) 
  sys.exit(1)


### [1] read option ###
PID = os.getpid()
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_summary_file() --> '%s'"%(OPT['of'])

of.write("#COMMAND '%s'\n"%(OPT['COMMAND']))
of.write("#DATE    '%s'\n"%(OPT['START_DATE']))
of.write("#[ligname3:1] [Nheavyatom:2] [Nring:3] [Nringblock:4] [molform:5]\n")
liglist = []
read_list_file(OPT['L'],liglist)
print "#PID %d"%(PID)

for ligA in (liglist):
## G20_2qwf_1_G_1.pdb
  fieldA = ligA.split('_')
  ligname3A = fieldA[0]
  (ligheadA,ligtailA) = ligA.split('.')

  tarfile_sdf  = "%s/%s.sdf"%(OPT['idtars'],ligname3A)
  tarfile_pdb  = "%s/%s"%(OPT['idrefp'],ligA)
  
  command_ckcombu    = "ckcombu -A %s -KV T"%(tarfile_sdf)
  Cdat = {}
  perform_kcombu_and_get_key_values(command_ckcombu,Cdat)
  of.write("%s %3d %2d %2d %s\n"%(ligname3A,int(Cdat['NHEAVYATOM']),int(Cdat['NRING']),int(Cdat['NRINGBLOCK']),Cdat['MOLFORM']))

now = datetime.now()
OPT['END_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

if (OPT['of'] != '-'):
  of.write("#END_DATE %s\n"%(OPT['END_DATE']))
  of.close()


