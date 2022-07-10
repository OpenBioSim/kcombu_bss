#!/usr/bin/env python

import sys
import os
import math
import re
from datetime import datetime
import time
import glob

LastModDate = "May 7, 2014"

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
OPT['refpro'] = ''
OPT['idpdb'] = '/home/DB/mmCIF_MOL/ASMBLMOL_PDB'
OPT['idpdb_old'] = '/home/takawaba/etc/Sutherland07/lost_pdb/ASMBLMOL_PDB/'
OPT['odpro']  = 'SupPRO'
OPT['odlig']  = 'SupLIG'



if (len(sys.argv)<2):
  print "mksup_matras.py <options>"
  print "  for superimpose proteins and ligands using the MATRAS program."
  print "  coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -L         : input file for ligand-protein list [%s]"%(OPT['L'])
  print " -idpdb     : input directory for PDBs [%s]"%(OPT['idpdb'])
  print " -idpdb_old : input directory for old PDBs [%s]"%(OPT['idpdb_old'])
  print " -refpro    : reference protein PDB file  [%s]"%(OPT['refpro'])
  print " -odpro     : output dir for rotated protein PDBs [%s]"%(OPT['odpro'])
  print " -odlig     : output dir for rotated ligand  PDBs [%s]"%(OPT['odlig'])
  print " -A         : Action ('T' or 'F') [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)


ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)

if (OPT['refpro']==''):
  OPT['refpro'] = ligpro_list[0]['pro']

(refhead,tail) = OPT['refpro'].split('.')
(pdbR,assemblyR,asymR,operR) = refhead.split('_')
ref_profile_full = OPT['idpdb'] + '/' + pdbR[1:3] + '/' + OPT['refpro']


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

  profile_full = OPT['idpdb'] + '/' + pdbP[1:3] + '/' + profile
  ligfile_full = OPT['idpdb'] + '/' + pdbL[1:3] + '/' + pdbL + '_' + assemblyL + '_' + asymL + '_' + operL + '.pdb' 

  if (os.access(profile_full,os.R_OK)==0):
    profile_full_old = OPT['idpdb_old'] + '/' + pdbP[1:3] + '/' + profile
    if (os.access(profile_full_old,os.R_OK)==0):
      print "#WARNING:Can't open profile_full '%s' and '%s'"%(profile_full,profile_full_old)
    else:
      profile_full = profile_full_old

  if (os.access(ligfile_full,os.R_OK)==0):
    ligfile_full_old = OPT['idpdb_old'] + '/' + pdbL[1:3] + '/' + pdbL + '_' + assemblyL + '_' + asymL + '_' + operL + '.pdb' 
    if (os.access(ligfile_full_old,os.R_OK)==0):
      print "#WARNING:Can't open ligfile_full '%s' and '%s'"%(ligfile_full,ligfile_full_old)
    else:
      ligfile_full = ligfile_full_old

  oprofile = OPT['odpro'] + '/' + profile
  oligfile = OPT['odlig'] + '/' + ligfile

  command = "Matras P -A %s -B %s -opdbA %s -ilgA %s -olgA %s "%(profile_full,ref_profile_full,oprofile,ligfile_full,oligfile)
  print command 





  if (OPT['A']=='T'):
    os.system(command)  
