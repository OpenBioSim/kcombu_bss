#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Oct 16, 2013"


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

def read_ligand_protein_list_file(ifname,list):
  print "#read_ligand_protein_list_file('%s')"%(ifname)
#>> FILE EXAMPLE 
# #[ligand_file] [protein_file]
# G20_2qwf_1_G_1.pdb 2qwf_1_A_1.pdb
# GNA_2qwe_1_G_1.pdb 2qwe_1_A_1.pdb
# DAN_2qwc_1_G_1.pdb 2qwc_1_A_1.pdb
# 49A_1f8e_1_H_1.pdb 1f8e_1_A_1.pdb
# 9AM_1f8d_1_H_1.pdb 1f8d_1_A_1.pdb
# G39_2qwh_1_E_1.pdb 2qwh_1_A_1.pdb
# SIA_2c4a_1_B_1.pdb 2c4a_1_A_1.pdb

  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR(read_ligand_protein_list_file): Can't open '%s'."%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0) and (len(line)>1):
      fields = line.split()
      ligand  = fields[0]
      protein = fields[1]
      dat = {}
      dat['lig'] = ligand
      dat['pro'] = protein
      list.append(dat)
  f.close()


def Action(command,action):
  if (action=='T'):
    os.system(command)
    print "#done '%s'"%(command)
  else:
    print "#%s"%(command)



############
### MAIN ###
############

OPT = {}
OPT['L']   = 'ligand.list'

OPT['A'] = 'F'
OPT['odtar'] = 'tmpout'
OPT['shaep'] = 'F'

if (len(sys.argv)<2):
  print "mk_subdir_odtar.py <options>"
  print " for making similariy-based 3D model using 'fkcombu' program."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -odtar  : output directory for predicted (transformed) 3D model in PDB format [%s]"%(OPT['odtar'])
  print " -shaep  : ShaEP output directory. Copy *.sdf and *.txt file (T or F) [%s]"%(OPT['A'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)


### [1] read option ###
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)

npair = 0

### [4] do fkcombu ###
  
for T in (ligpro_list):   ## target ##
  ligT = T['lig']
  fieldT = ligT.split('_')
  ligname3T = fieldT[0]
  subdir     = "%s/%s"%(OPT['odtar'],ligname3T)

  if (os.path.isdir(subdir)==0):
    print "subdir '%s' does not exist."%(subdir)
    mkdir_command = "mkdir %s"%(subdir)
    Action(mkdir_command,OPT['A'])

  for R in (ligpro_list): ## reference ##
    ligR = R['lig']
    fieldR = ligR.split('_')
    ligname3R = fieldR[0]
    (ligheadR,ligtailR) = ligR.split('.') 
    origpdbfile = "%s/%s_%s"%(OPT['odtar'],ligname3T,ligR)
    outpdbfile = "%s_%s"%(ligname3T,ligR)
    if (os.path.isfile(origpdbfile)):
      mv_command = "mv %s %s/%s"%(origpdbfile,subdir,outpdbfile)
      Action(mv_command,OPT['A'])

    if (OPT['shaep']=='T'):
      origsdffile = "%s/%s_%s.sdf"%(OPT['odtar'],ligname3T,ligheadR)
      origtxtfile = "%s/%s_%s_hits.txt"%(OPT['odtar'],ligname3T,ligheadR)
      origscorefile = "%s/%s_%s"%(OPT['odtar'],ligname3T,ligheadR)
      if (os.path.isfile(origsdffile)):
        outsdffile = "%s/%s_%s.sdf"%(subdir,ligname3T,ligheadR)
        mv_command = "mv %s %s"%(origsdffile,outsdffile)
        Action(mv_command,OPT['A'])
      if (os.path.isfile(origtxtfile)):
        outtxtfile = "%s/%s_%s_hits.txt"%(subdir,ligname3T,ligheadR)
        mv_command = "mv %s %s"%(origtxtfile,outtxtfile)
        Action(mv_command,OPT['A'])
      if (os.path.isfile(origscorefile)):
        outscorefile = "%s/%s_%s"%(subdir,ligname3T,ligheadR)
        mv_command = "mv %s %s"%(origscorefile,outscorefile)
        Action(mv_command,OPT['A'])

