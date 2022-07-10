#!/usr/bin/env python
import sys
import os

LastModDate = "Feb 2, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
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

############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['idp'] = 'SupLIG'
OPT['ids'] = '/DB/LIGAND-EXPO/SDF2d'
OPT['ods'] = ''
OPT['odp'] = ''
OPT['A'] = 'F'
OPT['cls'] = ''

if (len(sys.argv)<2):
  print "mk_paste3D.py [matchfileA] [matchfileB] <options>"
  print " for making pasting 3D xyz from PDB into SDF files"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L    : list of library molecules[%s]"%(OPT['L'])
  print " -idp  : input  dir for molecules in PDB with 3D xyz[%s]"%(OPT['idp'])
  print " -ids  : input  dir for molecules in SDF with 2D xyz[%s]"%(OPT['ids'])
  print " -ods  : output dir for molecules in SDF with 3D xyz[%s]"%(OPT['ods'])
  print " -odp  : output dir for molecules in PDB with 3D xyz[%s]"%(OPT['odp'])
  print " -cls  : class of molecule(such as 'NEU','EGFR') (optional)[%s]"%(OPT['cls'])
  print " -A    : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)

### [1] read option ###
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
liglist = []
read_list_file(OPT['L'],liglist)

command_common = "kcombu -fA S -fB P -alg B -con C -at K -nk 40 -pas3 T"
for i in range(len(liglist)):
  lig = liglist[i]
  (lig3,pdb) = lig.split('_')
  if (OPT['cls']!=''):
    osdffile = lig + '_' + OPT['cls']
  else:
    osdffile = lig
  if (OPT['ods'] != ''):
    command = "%s -A %s/%s -B %s/%s -osA %s/%s"%(command_common,OPT['ids'],lig3,OPT['idp'],lig,OPT['ods'],osdffile)

  if (OPT['odp'] != ''):
    command = "%s -A %s/%s -B %s/%s -opA %s/%s"%(command_common,OPT['ids'],lig3,OPT['idp'],lig,OPT['odp'],osdffile)

  print "#%s"%(command) 
  if (OPT['A']=='T'):
    os.system(command) 
