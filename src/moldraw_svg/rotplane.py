#!/usr/bin/env python
##
## <rotplane.py>
##


import sys
import os
import random
from datetime import datetime
import molecule
import pca

LastModDate = "Oct 22, 2011"

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



###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/DB/LIGAND-EXPO/SDF'
OPT['odir'] = 'tmpout'
OPT['tail'] = ''
OPT['A'] = 'F'
OPT['avow'] = 'F'

if (len(sys.argv)<3):
  print "rotplane.py <options>"
  print " for rotating molecules under '-idir', and write them to '-odir'."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir : input molecule directory [%s]"%(OPT['idir']) 
  print " -odir : output molecule directory for rotated SDF file [%s]"%(OPT['odir']) 
  print " -tail : tail string required for file [%s]"%(OPT['tail'])
  print " -avow : avoid overwriting to existing files (T or F)[%s]"%(OPT['avow']) 
  print " -A    : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['idir'])

filelist = []

for x in (dirlist):
  xfull = OPT['idir'] + '/' + x
  if (os.path.isdir(xfull)==0):
    if (OPT['tail']=='') or (x.endswith(OPT['tail'])):
      filelist.append(x)


### (1) Getting filelist[] in the director 'idir' ###
for x in (filelist):
  ifile = OPT['idir'] + '/' + x
  print "#%s"%(ifile)
  mol = molecule.Molecule()

  if (ifile.endswith(".sdf")):
      mol.read_in_sdf(ifile)
  elif (ifile.endswith(".pdb")):
    mol.read_in_pdb(ifile)
    mol.guess_bonds_from_atom_xyz()
  elif (ifile.endswith(".mol2")):
    mol.read_in_mol2(ifile)
  else:
      mol.read_in_sdf(ifile)

 # (Rmat,gpos) = pca.cal_Rmat_and_gpos_molecule_XY_flat_plane_by_PCA(mol)

 # print "#Rmat0 %f %f %f"%(Rmat[0][0], Rmat[0][1], Rmat[0][2])
 # print "#Rmat1 %f %f %f"%(Rmat[1][0], Rmat[1][1], Rmat[1][2])
 # print "#Rmat2 %f %f %f"%(Rmat[2][0], Rmat[2][1], Rmat[2][2])
 # print "#gpos %f %f %f"%(gpos[0], gpos[1], gpos[2])

  pca.rotate_molecule_XY_flat_plane_by_PCA(mol)

  ofile = OPT['odir'] + '/' + x
  if (OPT['avow']=='T') and (os.access(ofile,os.F_OK)):
    print "#WARNING:file '%s' already exists."%(ofile)
  else:
    if (OPT['A']=='T'):
      mol.write_in_sdf(ofile)
    else:
      print "#output_file '%s'"%(ofile)

