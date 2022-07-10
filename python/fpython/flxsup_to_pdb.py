#!/usr/bin/env python
import sys
import os
from datetime import datetime
import time

LastModDate = "Apr 18, 2014"


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

def Action(comstr,action):
  if (action=='T'):
    os.system(comstr)  
  else:
    print "#%s"%(comstr)




############
### MAIN ###
############

OPT = {}
OPT['L'] = 'ligand.list'
OPT['A'] = 'F'
OPT['idsdf']  = ''
OPT['idpdb']  = ''
OPT['idsmi']  = ''
OPT['odmol2'] = ''
OPT['odsdf']  = ''
OPT['div']  = '0/1'

if (len(sys.argv)<2):
  print "flxsup_to_pdb.py <options>"
  print " for flexible superposition of mol2 file onto pdb file."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options for input/output>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idmol2 : input directory for 3D MOL2 molecules[%s]"%(OPT['odmol2'])
  print " -idpdb  : input directory for PDB molecules[%s]"%(OPT['idpdb'])
  print " -odmol2 : output directory for 3D MOL2 molecules[%s]"%(OPT['odmol2'])
  print " -div    : [bunshi]/[bunbo] [%s]"%(OPT['div']) 
  print " -A    : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)


### [1] read option ###
PID = os.getpid()
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

### [2] read liglist (-L) ###
liglist = []
read_list_file(OPT['L'],liglist)
print "#PID %d"%(PID)

Nlig = len(liglist)
Nsta = Nlig*bunshi/bunbo
Nend = len(liglist)*(bunshi+1)/bunbo
print "#bunshi %d bunbo %d Npair %d Nsta %d Nend %d"%(bunshi,bunbo,Nlig,Nsta,Nend)
if (Nend>Nlig):
  Nend = Nlig

for i in range(Nsta,Nend):
  lig = liglist[i]
  field = lig.split('_')
  ligname3 = field[0]
  (ligfilecore,tail) = lig.split('.')
  #(ligname3,pdbch) = lig.split('_')
  #print "%s '%s' '%s'"%(lig,ligname3,pdbch)
  imol2file = "%s/%s.mol2"%(OPT['idmol2'],ligname3)
  ipdbfile  = "%s/%s"%(OPT['idpdb'],lig)
  omol2file = "%s/%s.mol2"%(OPT['odmol2'],ligname3)
  
  command = "fkcombu -T %s -fT 2 -R %s -fR P -omol2T %s -hyop T"%(imol2file,ipdbfile,omol2file) 
  Action(command,OPT['A'])

