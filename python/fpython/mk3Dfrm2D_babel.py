#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Mar 22, 2012"

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
    print "#'%s'"%(comstr)


############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['idsdf']  = '/DB/LIGAND-EXPO/SDF'
OPT['odmol2'] = 'tmpout'
OPT['3D'] = 'g'

if (len(sys.argv)<2):
  print "mk3Dfrm2D_babel.py <options>"
  print " for making 3D mol2 file from 2D sdf file using OpenBabel"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "**note: 'babel --gen3d' sometimes fails for all-zero-XYZ coordinate file as its input."
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idsdf  : input directory for 2D SDF molecules[%s]"%(OPT['idsdf'])
  print " -odmol2 : output directory for 3D MOL2 molecules[%s]"%(OPT['odmol2'])
  print " -div  : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -3D   : 3D generation.  'g':babel --gen3d [%s]"%(OPT['3D'])
  print "       :                 '3': three commands(babel -h;obrotamer;obconformer)"
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
  if (lig.find('-')>=0):
    (ligname3,pdbch) = lig.split('_')
  else:
    ligname3 = lig
    pdbch = ''
  print "%s '%s' '%s'"%(lig,ligname3,pdbch)

  if (OPT['3D']=='g'):  
    comstr = "babel --gen3d -isdf %s/%s -omol2 %s/%s.mol2"%(OPT['idsdf'],ligname3,OPT['odmol2'],ligname3)
    Action(comstr,OPT['A'])

  if (OPT['3D']=='3'):  
    comstr = "babel -h -isdf %s/%s -omol2 %d-1.mol2"%(OPT['idsdf'],ligname3,PID)
    Action(comstr,OPT['A'])

    comstr = "obrotamer %d-1.mol2 > %d-2.mol2"%(PID,PID)
    Action(comstr,OPT['A'])

    comstr = "obconformer 250 100 %d-2.mol2 > %s/%s.mol2"%(PID,OPT['odmol2'],ligname3)
    Action(comstr,OPT['A'])
  
    comstr = "rm %d-1.mol2"%(PID)
    Action(comstr,OPT['A'])

    comstr = "rm %d-2.mol2"%(PID)
    Action(comstr,OPT['A'])
