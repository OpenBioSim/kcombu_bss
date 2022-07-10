#!/usr/bin/env python

import sys
import os
import math
import re
from datetime import datetime
import time
import glob

LastModDate = "May 8, 2014"

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
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open '%s'."%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
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



##############
#### MAIN ####
##############

OPT = {}
OPT['L'] = ''
OPT['A'] = 'F'

OPT['isdf3Di']  = '/home/DB/mmCIF_LIGAND/SDF_3D_IDEAL'
OPT['osdf3Di']  = 'SDF_3D'
OPT['ofail'] = '-'

if (len(sys.argv)<2):
  print "cp_sdf3d.py <options>"
  print "  for copying 3D SDF file."
  print "  coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -L      : input file for ligand-protein list [%s]"%(OPT['L'])
  print " -isdf3Di: input dir for 3D sdf ideal coordinates [%s]"%(OPT['isdf3Di'])
  print " -osdf3Di: output dir for 3D sdf ideal coordinates [%s]"%(OPT['osdf3Di'])
  print " -ofail  : output file for failed ligands [%s]"%(OPT['ofail'])
  print " -A         : Action ('T' or 'F') [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)


ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)

fail_liglist = []
for ligpro in (ligpro_list):   ## target ##
# #[ligand_file] [protein_file]
# G20_2qwf_1_G_1.pdb 2qwf_1_A_1.pdb
# GNA_2qwe_1_G_1.pdb 2qwe_1_A_1.pdb

  ligfile = ligpro['lig']
  (lighead,tail) = ligfile.split('.')
  profile = ligpro['pro']
  (prohead,tail) = profile.split('.')
  print ligfile,profile
  (ligname3,pdbL,assemblyL,asymL,operL) = lighead.split('_')
  (pdbP,assemblyP,asymP,operP) = prohead.split('_')

  iligsdffile = OPT['isdf3Di'] + '/' + ligname3[0] + '/' + ligname3 + '.sdf'
  oligsdffile = OPT['osdf3Di'] + '/' + ligname3 + '.sdf'
  if (os.access(iligsdffile,os.R_OK)):
    command = "cp %s %s"%(iligsdffile,oligsdffile)
    print "#%s"%(command)
    if (OPT['A']=='T'):
      os.system(command)
  else:
   fail_liglist.append(ligfile)

### OUTPUT Failed Ligand list ##

if (OPT['ofail']!='-'):
  of = open(OPT['ofail'],'w')
  print "#write_fail_ligands() --> '%s'"%(OPT['ofail'])
else:
  of = sys.stdout
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))
for lig in (fail_liglist):
  of.write("FAIL_LIGAND %s\n"%(lig))
if (OPT['ofail']!='-'):
  of.close() 


