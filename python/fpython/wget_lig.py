#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Mar 14, 2014"


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


def write_chimera_command_file_for_ligand(ofname,isdffile,omol2file):
#open SupPRO/1bji_1_A_1.pdb
#addh
#addcharge
#write format mol2 0 out.mol2
  print "#write_chimera_command_file_ligand() -->'%s'"%(ofname)
  of = open(ofname,'w')
  of.write("open %s\n"%(isdffile))
#  of.write("addh\n")
  of.write("addcharge nonstd :#0 0 method am1\n")
  of.write("write format mol2 0 %s\n"%(omol2file))
  of.write("stop\n")
  of.close()



############
### MAIN ###
############

OPT = {}
OPT['L']   = 'ligand.list'

OPT['A'] = 'F'
OPT['idsdf']   = 'SDF_3D'
OPT['odcif']   = 'CIF_3D'
OPT['ligexpo']   = 'http://ligand-expo.rcsb.org/reports'


PID = os.getpid()

if (len(sys.argv)<2):
  print "wget_lig.py <options>"
  print " for getting ligand files from LIGAND-EXPO server."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L       : list of library molecules[%s]"%(OPT['L'])
  print " -ligexpo : web address of LIGAND-EXPO [%s]"%(OPT['ligexpo'])
  print " -odcif   : output directory for mmCIF files[%s]"%(OPT['odcif'])
  print "           (One of '-idsdf' and '-isdmol2' should be assigned.)"
  print " -A       : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)

### [1] read option ###
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)

ocomfile = '%d.com'%(PID)

for ligpro in (ligpro_list):
# G39_2qwh_1_E_1.pdb 2qwh_1_A_1.pdb
# SIA_2c4a_1_B_1.pdb 2c4a_1_A_1.pdb
  #print ligpro
  ligT = ligpro['lig']
  (ligname3T,pdbT,assembly_idT,asym_idT,operT) = ligT.split('_')
  ociffile  = OPT['odcif'] + '/' + ligname3T + '.cif'
  subdir = ligname3T[0] 
  wgetcommand = "wget %s/%s/%s/%s.cif -O %s"%(OPT['ligexpo'],subdir,ligname3T,ligname3T,ociffile)
  if (OPT['A']=='T'):
    os.system(wgetcommand)
  else:
    print "#WGET '%s'"%(wgetcommand)

