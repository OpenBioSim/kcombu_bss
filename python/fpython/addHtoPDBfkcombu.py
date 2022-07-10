#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Sep 1, 2013"

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
OPT['L'] = ''
OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['idpdb']    = 'SupLIG'
OPT['idwithH']  = 'SDF_3D'
OPT['odwithH']  = 'tmpout'
OPT['odmol2'] = 'tmpout'
OPT['iftype'] = 'S'
OPT['oftype'] = 'S'
OPT['oligname3'] = 'F'
OPT['olog'] = 'out.log'
OPT['amste3D'] = 'T'

if (len(sys.argv)<2):
  print "addHtoPDBfkcombu.py <options>"
  print " for making 3D molecule with H from 3D PDB file and with H file using 'fkcombu'."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L        : list of library molecules[%s]"%(OPT['L'])
  print " -idpdb    : input directory for 3D PDB molecules[%s]"%(OPT['idpdb'])
  print " -idwithH  : input directory for molecular files with hydrogen[%s]"%(OPT['idwithH'])
  print " -odwithH  : output directory for 3D molecular files with hydrogen[%s]"%(OPT['odwithH'])
  print " -iftype   : output molecular file type. 'S':sdf '2':mol2, 'P':PDB [%s]"%(OPT['iftype'])
  print " -oftype   : output molecular file type. 'S':sdf '2':mol2, 'P':PDB [%s]"%(OPT['oftype'])
  print " -oligname3: output file name is just ligname3 or not 'T' or 'F' [%s]"%(OPT['oligname3'])
  print " -amste3D  : adjust atom match to agree 3D chirality by permutations (T or F)[%s]"%(OPT['amste3D'])
  print " -div      : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -olog     : output log file [%s]"%(OPT['olog'])
  print " -A        : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)


### [1] read option ###
PID = os.getpid()
print "#PID %d"%(PID)
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

### [2] read liglist (-L) ###
#liglist = []
#read_list_file(OPT['L'],liglist)

ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig = len(ligpro_list)
Nsta = Nlig*bunshi/bunbo
Nend = Nlig*(bunshi+1)/bunbo
print "#bunshi %d bunbo %d Npair %d Nsta %d Nend %d"%(bunshi,bunbo,Nlig,Nsta,Nend)
if (Nend>Nlig):
  Nend = Nlig

if (OPT['A']=='T'):
  oflog = open(OPT['olog'],'w')
  print "#write_log_file() -->'%s'"%(OPT['olog'])
  oflog.write("#COMMAND    %s\n"%(OPT['COMMAND']))
  oflog.write("#START_DATE %s\n"%(OPT['START_DATE']))
  oflog.write("#[infile:1] [outfile:2] [tanimoto:3] [rmsd_match:4]\n")

for i in range(Nsta,Nend):
  lig = ligpro_list[i]['lig']

  #(ligname3,pdbch) = lig.split('_')
  field = lig.split('_')
  ligname3 = field[0]
  (ligfilecore,tail) = lig.split('.')

  #print "%s '%s' '%s'"%(lig,ligname3,pdbch)
  ipdbfile     = "%s/%s"%(OPT['idpdb'],lig)

  iwithHfile   = "%s/%s.sdf"%(OPT['idwithH'],ligname3)

  if (OPT['oligname3']=='T'):
    owithHfile   = "%s/%s"%(OPT['odwithH'],ligname3)
  else:
    owithHfile   = "%s/%s"%(OPT['odwithH'],ligfilecore)

  if (OPT['oftype']=='2'):
    owithHfile   = "%s.mol2"%(owithHfile)
  elif (OPT['oftype']=='P'):
    owithHfile   = "%s"%(owithHfile)
  elif (OPT['oftype']=='S'):
    owithHfile   = "%s.sdf"%(owithHfile)

  #comstr = "fkcombu -con T -mtd 1 -hyop T -fT %s -fR P -E A -S X -stp T -gra T -chr F -rms T"%(OPT['iftype'])
  comstr = "fkcombu -con T -mtd 1 -hyop T -fT %s -fR P -E A -S F -amste3D %s"%(OPT['iftype'],OPT['amste3D'])

  comstr += " -T %s -R %s "%(iwithHfile,ipdbfile)

  if (OPT['oftype']=='2'):
    comstr += ' -omol2T %s'%(owithHfile)
  if (OPT['oftype']=='P'):
    comstr += ' -opdbT %s'%(owithHfile)
  if (OPT['oftype']=='S'):
    comstr += ' -osdfT %s'%(owithHfile)

  valkey = {}
  if (OPT['A']=='T'):
    perform_kcombu_and_get_key_values(comstr,valkey)
    oflog.write("%s %s %s %s\n"%(iwithHfile,ipdbfile,valkey['TANIMOTO'],valkey['RMSD_MATCH']))
  else:
    print "#%s"%(comstr)


if (OPT['A']=='T'):
  print "#write_log_file() -->'%s' done."%(OPT['olog'])
  oflog.close()
