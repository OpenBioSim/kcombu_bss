#!/usr/bin/env python
import sys
import os
import math
from datetime import datetime

LastModDate = "Apr 2, 2014"

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




def read_XYZ_from_pdb_file(ipdbfile,X,Y,Z):
#          1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#HETATM   20  O   MOL A   1      27.169  14.191  62.446  0.00  0.00  O1  O x - A 1 0 0 0 0 0   1   3   4   0   0
#HETATM   21  H   MOL A   1      22.638  17.272  63.699  0.00  0.00      H x - H 1 0 0 0 0 0   0   0   0   0   0
#HETATM 4008  C10 EQP A 500      27.445  14.640  61.896  0.95  0.00
#HETATM 4010  N5  EQP A 500      26.556  15.627  61.837  0.77  0.00
#HETATM 4013  O3P EQP A 500      22.655  20.364  65.637  1.00  0.00
#HETATM 4019  O10 EQP A 500      28.023  14.307  62.949  0.83  0.00
#HETATM 4020  H2  EQP A 500      24.274  19.039  63.113  1.00  0.00
#HETATM 4021  H31 EQP A 500      23.256  17.174  64.524  1.00  0.00
#HETATM 4022  H32 EQP A 500      24.901  16.973  65.199  1.00  0.00
#HETATM 4030 H111 EQP A 500      28.343  14.492  59.955  1.00  0.00
#HETATM 4031 H112 EQP A 500      26.751  13.647  60.120  1.00  0.00
#HETATM 4032 H113 EQP A 500      28.237  12.939  60.857  1.00  0.00

  if (os.access(ipdbfile,os.R_OK)==0):
    print "#ERROR:Can't pdb file open '%s'."%(ipdbfile)
    sys.exit(1)
  f = open(ipdbfile)
  Nheavyatom = 0;
  for line in f:
    if (line.startswith('HETATM') or line.startswith('ATOM')):
      atomstr = line[12:16].replace(' ','')
      if (atomstr[0]!='H'):
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        X.append(x)
        Y.append(y)
        Z.append(z)
      pass

  f.close()


def max_closest_distance(Ax,Ay,Az,Bx,By,Bz):
  NatomA = len(Ax)
  NatomB = len(Bx)
  maxDclosest = -10.0
  for a in range(NatomA):
    Dclosest = -10.0
    for b in range(NatomB):
      D = math.sqrt((Ax[a]-Bx[b])*(Ax[a]-Bx[b]) +(Ay[a]-By[b])*(Ay[a]-By[b]) +(Az[a]-Bz[b])*(Az[a]-Bz[b]))
      if (Dclosest < 0.0) or (D<Dclosest):
        Dclosest = D
    if (maxDclosest<0.0) or (Dclosest>maxDclosest):
      maxDclosest = Dclosest
  return(maxDclosest)


def cal_gcenter(Ax,Ay,Az):
  NatomA = len(Ax)
  GxA = 0.0
  GyA = 0.0
  GzA = 0.0
  for a in range(NatomA):
    GxA += Ax[a]
    GyA += Ay[a]
    GzA += Az[a]
  GxA /= NatomA
  GyA /= NatomA
  GzA /= NatomA
  return((GxA,GyA,GzA)) 





############
### MAIN ###
############

OPT = {}
OPT['L']   = 'ligand.list'
OPT['idref'] = 'SupLIG'

OPT['A'] = 'F'
OPT['dligpdb'] = 'SupLIG'
PID = os.getpid()

if (len(sys.argv)<2):
  print "dis_bwn_ligands.py <options>"
  print " for checking distance between ligands."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print "<options for MODE='L'>"
  print " -dligpdb  : input dir. for 3D ligand PDB file [%s]"%(OPT['dligpdb'])
  sys.exit(1)

### [1] read option ###
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)


MAX_maxDclose = -1.0 
MIN_maxDclose = -1.0
AVE_maxDclose = 0.0
VAR_maxDclose = 0.0
MIN_Dg = -1.0
MAX_Dg = -1.0
AVE_Dg = 0.0
VAR_Dg = 0.0
Npair = 0
maxpair = ''
for a in range(Nlig):
  ligproA = ligpro_list[a]
  ligA = ligproA['lig']
  ligAfile = OPT['dligpdb'] + '/' + ligA
  Ax = []
  Ay = []
  Az = []
  read_XYZ_from_pdb_file(ligAfile,Ax,Ay,Az)
  NatomA = len(Ax)
  #print "#%s Natom %d"%(ligAfile,NatomA)
  #for i in range(NatomA):
  #  print Ax[i],Ay[i],Az[i]
  (GxA,GyA,GzA) = cal_gcenter(Ax,Ay,Az)
  for b in range(Nlig):
    if (a != b):
      ligproB = ligpro_list[b]
      ligB = ligproB['lig']
      ligBfile = OPT['dligpdb'] + '/' + ligB
      Bx = []
      By = []
      Bz = []
      read_XYZ_from_pdb_file(ligBfile,Bx,By,Bz)
      NatomB = len(Bx)
      (GxB,GyB,GzB) = cal_gcenter(Bx,By,Bz)
      maxDclose = max_closest_distance(Ax,Ay,Az,Bx,By,Bz)
      Dg = math.sqrt((GxA-GxB)*(GxA-GxB) +(GyA-GyB)*(GyA-GyB) +(GzA-GzB)*(GzA-GzB))


      print "%s %s maxDclose %f Dg %f"%(ligAfile,ligBfile,maxDclose,Dg)
      Npair += 1
      if (MAX_maxDclose<0.0) or (maxDclose > MAX_maxDclose):
        MAX_maxDclose = maxDclose 
        maxpair = ligA + '_' + ligB
      if (MIN_maxDclose<0.0) or (maxDclose < MIN_maxDclose):
        MIN_maxDclose = maxDclose 
      AVE_maxDclose += maxDclose
      VAR_maxDclose += maxDclose*maxDclose
      if (MIN_Dg<0.0) or (Dg < MIN_Dg):
        MIN_Dg = Dg
      if (MAX_Dg<0.0) or (Dg > MAX_Dg):
        MAX_Dg = Dg
      AVE_Dg += Dg
      VAR_Dg += Dg*Dg

AVE_maxDclose = AVE_maxDclose/Npair
AVE_Dg        = AVE_Dg/Npair
SD_Dg        = math.sqrt(VAR_Dg/Npair - AVE_Dg * AVE_Dg)
SD_maxDclose = math.sqrt(VAR_maxDclose/Npair - AVE_maxDclose * AVE_maxDclose)
print "#%s maxDclose ave %.2f sd %.2f min %.2f max %.2f %s Dg ave %.2f sd %.2f min %.2f max %.2f"%(OPT['L'],AVE_maxDclose,SD_maxDclose,MIN_maxDclose,MAX_maxDclose,maxpair,AVE_Dg,SD_Dg,MIN_Dg,MAX_Dg) 
  
