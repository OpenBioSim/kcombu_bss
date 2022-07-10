#!/usr/bin/env python
import sys
import os
import math
from datetime import datetime

LastModDate = "Sep 4, 2013"


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



def perform_kcombu_and_get_key_values(commandp,valkey):
#>>Example of -KV output of pkcombu<<
# NHEAVYATOM_A 28
# NHEAVYATOM_B 26
# NATOM_MATCH  15
# TANIMOTO     0.384615
# RMSD         0.194683
  f = os.popen(commandp)
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#')==0) and (len(field)>=2):
      key   = field[0]
      value = field[1]
      valkey[key] = value
    pass
  f.close()


def read_element_composition_in_pdb_file(ipdbfile,Nelement):
  if (os.access(ipdbfile,os.R_OK)==0):
    print "#WARNING:Can't open '%s'."%(ipdbfile)
    pass 

  f = open(ipdbfile)
#HETATM 2128  S   SRX A   1      15.287   2.681  16.972  1.00  0.00  C   263  S
#HETATM 2129  CL  SRX A   1      20.075  11.815   9.439  0.50  0.00  C   263 CL
#HETATM 2130  RU  SRX A   1      19.495  11.763  11.832  0.50  0.00  C   263 RU
#HETATM 2131  C1  SRX A   1      15.807   4.002  16.167  1.00  0.00  C   263  C
#HETATM 2132  N1  SRX A   1      15.860   1.420  16.291  1.00  0.00  C   263  N
#HETATM 2133  O1  SRX A   1      15.793   2.807  18.314  1.00  0.00  C   263  O
#HETATM 2134  C2  SRX A   1      17.028   4.096  15.516  1.00  0.00  C   263  C
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#          1         2         3         4         5         6         7
  for line in f:
    if line.startswith('HETATM') or line.startswith('ATOM  '):
      line = line.rstrip('\n')
      element = line[76:78].replace(' ','')
      #print "'%s'->'%s'"%(line,element)
      Nelement[element] = Nelement.get(element,0) + 1  

  f.close()



############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'
OPT['A'] = 'F'
OPT['idsdf']  = 'SDF_3D'
OPT['idpdb'] = 'SupLIG'
OPT['self']     = 'F'
OPT['of']       = '-'
OPT['pkcombu']  = '/home/takawaba/bin/pkcombu'
OPT['fkcombu']  = '/home/takawaba/bin/fkcombu'
OPT['chste3D']  = 'F'
OPT['amste3D']  = 'T'
OPT['chkrng']   = 'T'

if (len(sys.argv)<2):
  print "cmpSDFwithPDB.py <options>"
  print " for comparison of ideal SDF file with actual PDB files."
  print " Also check metal ions."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idsdf  : input dir for molecules in SDF[%s]"%(OPT['idsdf'])
  print " -idpdb  : input dir for moleculee in PDB[%s]"%(OPT['idpdb'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print " -of     : output summary file [%s]"%(OPT['of']) 
  print " -pkcombu: commandp for pkcombu [%s]"%(OPT['pkcombu']) 
  print " -fkcombu: commandp for fkcombu [%s]"%(OPT['fkcombu']) 
  print " -chste3D: checking chirarity agreement between two 3D conformations (T or F)[%s]"%(OPT['chste3D'])
  print " -amste3D: adjust atom match to agree 3D chirality by permutations (T or F)[%s]"%(OPT['amste3D'])
  print " -chkrng : checking ring 3D conformations (T or F)[%s]"%(OPT['chkrng'])
  sys.exit(1)

METALS = ['LI','NA','K','RB','CS','FR','BE','MG','CA','SR','BA','SC','TI','V','CR','MN','FE','CO','NI','CU','ZN','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT','AU']



### [1] read option ###
PID = os.getpid()
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
if (OPT['of']=='-') or (OPT['A']!='T'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_summary_file() --> '%s'"%(OPT['of'])

of.write("#COMMAND '%s'\n"%(OPT['COMMAND']))
of.write("#DATE    '%s'\n"%(OPT['START_DATE']))
of.write("#[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [tanimoto:5] [Nmetal:6] [Ndiff_chiral_sp3:7] [Ndiff_ring_conf:8]\n")
liglist = []
read_list_file(OPT['L'],liglist)
print "#PID %d"%(PID)

for ligA in (liglist):
## G20_2qwf_1_G_1.pdb
  (ligheadA,ligtailA) = ligA.split('.')
  fieldA = ligA.split('_')
  ligname3A = fieldA[0]

  molfile_sdf  = "%s/%s.sdf"%(OPT['idsdf'],ligname3A)
  molfile_pdb  = "%s/%s.pdb"%(OPT['idpdb'],ligheadA)
  commandp   = "%s -A %s -B %s -fA S -fB P -con C -mtd -1 -KV T"%(OPT['pkcombu'],molfile_sdf,molfile_pdb)
  commandf   = "%s -T %s -R %s -fT S -fR P -con C -mtd -1 -KV T -chste3D %s -amste3D %s -chkrng %s -S N"%(OPT['fkcombu'],molfile_sdf,molfile_pdb,OPT['chste3D'],OPT['amste3D'],OPT['chkrng'])
  pdat = {}
  fdat = {}
  if (OPT['A']=='T'): 
    perform_kcombu_and_get_key_values(commandp,pdat)
    perform_kcombu_and_get_key_values(commandf,fdat)
    Nelement = {}
    read_element_composition_in_pdb_file(molfile_pdb,Nelement)
    NheavyA  = int(pdat['NHEAVYATOM_A'])
    NheavyB  = int(pdat['NHEAVYATOM_B'])
    tanimoto = float(pdat['TANIMOTO'])
    Ndiff_chiral_sp3 = int(fdat.get('NDIFF_CHIRAL_SP3','0'))
    Ndiff_ring_conf  = int(fdat.get('NDIFF_RING_CONF','0'))
    Nmetal = 0 
    metal_str = ''
    for m in (METALS):
      if (Nelement.get(m,0)>0):
        metal_str += " %s %d"%(m,Nelement[m])
        Nmetal += 1 
    of.write("%s %s %3d %3d %8.6f %3d %d %d"%(molfile_sdf,molfile_pdb,NheavyA,NheavyB,tanimoto,Nmetal,Ndiff_chiral_sp3,Ndiff_ring_conf))
    of.write(" %s"%(metal_str))
    if (NheavyA > NheavyB):
      of.write(" DELETE_FOR_LACKED_ATOMS")
    elif (NheavyA > 100):
      of.write(" DELETE_FOR_TOO_LARGE_NHEAVY")
    elif (Nmetal > 0):
      of.write(" DELETE_FOR_METAL")
    elif (tanimoto < 1.0):
      of.write(" DELETE_FOR_LESS_ONE_TANIMOTO")
    elif (Ndiff_chiral_sp3 > 0):
      of.write(" DELETE_DIFF_CHIRALITY")
    of.write("\n")
  else:
    print "#%s"%(commandp)
    print "#%s"%(commandf)

if (OPT['of'] != '-') and (OPT['A']=='T'): 
    of.close()


