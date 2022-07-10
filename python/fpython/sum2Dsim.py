#!/usr/bin/env python
import sys
import os
import math
from datetime import datetime

LastModDate = "Sep 16, 2013"


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

def perform_openbabel_for_fingerprint(command):
#>>Example of babel -ofpt <<
# babel SDF_3D/G39.sdf SDF_3D/ZMR.sdf -ofpt -xfFP2
# >G39
# >ZMR   Tanimoto from G39 = 0.422111
# 2 molecules converted
# 42 audit log messages
  f = os.popen(command_babel)
  tanimoto_fpt = -1.0
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (len(field)>=5) and (field[1]=='Tanimoto') and (field[2]=='from'):
      tanimoto_fpt = float(field[5])
  f.close()
  return(tanimoto_fpt)


############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'
OPT['A'] = 'F'
OPT['idtars']  = 'SDF_3D'
OPT['idrefp'] = 'SupLIG'
OPT['self']   = 'F'
OPT['of']     = '-'


if (len(sys.argv)<2):
  print "sum2Dsim.py <options>"
  print " for making summary of accuracy of predicted simimilariy-based 3D model."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idtars : input dir for target molecules in SDF[%s]"%(OPT['idtars'])
  print " -idrefp : input dir for reference 3D model in PDB[%s]"%(OPT['idrefp'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print " -self   : only self modeling (T or F)[%s]"%(OPT['self'])
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
of.write("#[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [NheavyB-NheavyA:5]\n")
of.write("#[tanimoto(C_MCS_K):6] [tanimoto(TD_MCS1_K):7] [tanimoto(C_MCS_E):8] [tanimoto(TD_MCS1_E):9]\n")
of.write("#[tanimoto_fingerprint:10] [tanimoto_volume:11]\n")
of.write("#[Nmatch(C_MCS_K):12] [Nmatch(TD_MCS1_K):13] [Nmatch(C_MCS_E):14] [Nmatch(TD_MCS1_E):15]\n")
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
  
  for ligB in (liglist):
    fieldB = ligB.split('_')
    ligname3B = fieldB[0]
    (ligheadB,ligtailB) = ligB.split('.')

    reffile_pdb  = "%s/%s"%(OPT['idrefp'],ligB)
    reffile_sdf  = "%s/%s.sdf"%(OPT['idtars'],ligname3B)

    if ((OPT['self']!='T') or (ligA == ligB)):
      command_CK      = "pkcombu -A %s -B %s -fA S -fB P -con C        -at K -nk 40  -KV T"%(tarfile_sdf,reffile_pdb)
      command_TD1K    = "pkcombu -A %s -B %s -fA S -fB P -con T -mtd 1 -at K -nk 40  -KV T"%(tarfile_sdf,reffile_pdb)
      command_CE      = "pkcombu -A %s -B %s -fA S -fB P -con C        -at E -nk 100 -KV T"%(tarfile_sdf,reffile_pdb)
      command_TD1E    = "pkcombu -A %s -B %s -fA S -fB P -con T -mtd 1 -at E -nk 100 -KV T"%(tarfile_sdf,reffile_pdb)
      command_fkcombu = "fkcombu -T %s -R %s -fT P -fR P -S N -con C -at K -nk 100 -KV T"%(tarfile_pdb,reffile_pdb)
      print "#'%s'"%(command_CK)
      print "#'%s'"%(command_TD1K)
      print "#'%s'"%(command_CE)
      print "#'%s'"%(command_TD1E)
      print "#'%s'"%(command_fkcombu)
      command_babel = "babel %s %s -ofpt -xfFP2"%(tarfile_sdf,reffile_sdf)
      print "#'%s'"%(command_babel)

      if (OPT['A']=='T'):
        CKdat = {}
        perform_kcombu_and_get_key_values(command_CK,CKdat)
        TD1Kdat = {}
        perform_kcombu_and_get_key_values(command_TD1K,TD1Kdat)
        CEdat = {}
        perform_kcombu_and_get_key_values(command_CE,CEdat)
        TD1Edat = {}
        perform_kcombu_and_get_key_values(command_TD1E,TD1Edat)
        FKCOMBUdat = {}
        perform_kcombu_and_get_key_values(command_fkcombu,FKCOMBUdat)

        tanimoto_fpt = perform_openbabel_for_fingerprint(command_babel)
        diffBA = int(CKdat['NHEAVYATOM_B'])-int(CKdat['NHEAVYATOM_A'])
        of.write("%s %s %s %s %2d"%(ligA,ligB,CKdat['NHEAVYATOM_A'],CKdat['NHEAVYATOM_B'],diffBA))
        of.write(" %s %s"%(CKdat['TANIMOTO'],TD1Kdat['TANIMOTO']))
        of.write(" %s %s"%(CEdat['TANIMOTO'],TD1Edat['TANIMOTO']))
        of.write(" %f %s"%(tanimoto_fpt, FKCOMBUdat['TANIMOTO_VOL']))
        of.write(" %s %s"%(CKdat['NATOM_MATCH'],TD1Kdat['NATOM_MATCH']))
        of.write(" %s %s"%(CEdat['NATOM_MATCH'],TD1Edat['NATOM_MATCH']))
        of.write("\n")

now = datetime.now()
OPT['END_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

if (OPT['of'] != '-'):
  of.write("#END_DATE %s\n"%(OPT['END_DATE']))
  of.close()


