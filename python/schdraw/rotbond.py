#!/usr/bin/env python

##
## <rotbond.py>
##

import sys
import os 
import math
from pyx import *
import molecule 
import molfunc

LastModDate = "Feb 2, 2011"


def make_rotational_axis(axis,mol,natom_sta,natom_ori):
  axis[0] = mol.atoms[natom_ori-1].x - mol.atoms[natom_sta-1].x
  axis[1] = mol.atoms[natom_ori-1].y - mol.atoms[natom_sta-1].y
  axis[2] = mol.atoms[natom_ori-1].z - mol.atoms[natom_sta-1].z
  len = math.sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2])
  axis[0] /= len
  axis[1] /= len
  axis[2] /= len


def make_rotational_matrix_from_axis_and_angle(R,axis,angle):
  ux = axis[0]
  uy = axis[1]
  uz = axis[2]
  sint   = math.sin(math.pi*angle/180.0)
  cost   = math.cos(math.pi*angle/180.0)
  one_cost   = 1.0-cost

  R[0][0] = cost + ux*ux*one_cost
  R[0][1] = ux*uy*one_cost - uz*sint
  R[0][2] = ux*uz*one_cost + uy*sint

  R[1][0] = uy*ux*one_cost + uz*sint
  R[1][1] = cost + uy*uy*one_cost
  R[1][2] = uy*uz*one_cost - ux*sint

  R[2][0] = uz*ux*one_cost - uy*sint
  R[2][1] = uz*uy*one_cost + ux*sint
  R[2][2] = cost + uz*uz*one_cost

def rotate_part_of_molecule(mol,ratom_list,orig_natom,R):
  ratomdic = {} 
  for r in (ratom_list):
    ratomdic[int(r)-1] = 1
  ori  = [0.0,0.0,0.0]
  pos  = [0.0,0.0,0.0]
  npos = [0.0,0.0,0.0]
  ori[0] = mol.atoms[orig_natom-1].x 
  ori[1] = mol.atoms[orig_natom-1].y 
  ori[2] = mol.atoms[orig_natom-1].z 
  
  for a in range(mol.Natom):
    if (ratomdic.has_key(a)):
      pos[0] = mol.atoms[a].x
      pos[1] = mol.atoms[a].y
      pos[2] = mol.atoms[a].z
      for i in range(3):
        npos[i] = 0.0
        for j in range(3):
          npos[i] += R[i][j] * (pos[j]-ori[j]) 
        npos[i] += ori[i] 
      print "(%f %f %f) -> (%f %f %f)"%(pos[0],pos[1],pos[2],npos[0],npos[1],npos[2])
      mol.atoms[a].x = npos[0]
      mol.atoms[a].y = npos[1]
      mol.atoms[a].z = npos[2]

  pass

#############
#### MAIN ###
#############

OPT = {}

 
## read options ##
OPT['isdf']  = '' 
OPT['ipdb']  = '' 
OPT['bond']  = ''
OPT['ratm']  = ''
OPT['xyz'] = '0:0:0'
OPT['ang'] = '0.0'
OPT['osdf'] = ''
OPT['orig'] = ''

if (len(sys.argv)<2):
  print "rotbond.py <options>"
  print " for rotating atoms around a given bond."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options for single molecule>"
  print " -isdf : input SDF file [%s]"%(OPT['isdf']) 
  print " -ipdb : input PDB file [%s]"%(OPT['ipdb']) 
  print " -osdf : output SDF file [%s]"%(OPT['osdf']) 
  print " -bond : two atom number for bond atom   [natm_start]:[natm_origin] [%s]"%(OPT['bond']) 
  print " -ratm : atom numbers for rotating atoms [natm1]:[natm2]:... [%s]"%(OPT['ratm']) 
  print " -ang  : rotational angle around bond axis (degree) [%s]"%(OPT['ang']) 
  print " -xyz  : (rasmol) XYZ-rotational angle around origin th_x:th_y:th_z:(degree)[%s]"%(OPT['xyz'])
  print " -orig : origin position ox:oy:oz:[%s]"%(OPT['orig'])
  sys.exit(1)

if (len(sys.argv)>2):
  molfunc.read_option(sys.argv,OPT)




molA = molecule.Molecule()
  
if (OPT['isdf'] != ''):
  molA.read_in_sdf(OPT['isdf'])
if (OPT['ipdb'] != ''):
  molA.read_in_pdb(OPT['ipdb'])
  molA.guess_bonds_from_atom_xyz()
  
if (OPT['xyz'] != '0:0:0'):
  thxyz = OPT['xyz'].split(':')
  R = molfunc.make_rotation_matrix_by_xyz_angle(float(thxyz[0]),float(thxyz[1]),float(thxyz[2]))
  if (OPT['orig'] == ''):
    gcen = molA.cal_gcenter()
  else:
    olist = OPT['orig'].split(':')
    gcen= [0.0,0.0,0.0]
    gcen[0] = float(olist[0])  
    gcen[1] = float(olist[1])  
    gcen[2] = float(olist[2])  
  print "#gcen %f:%f:%f"%(gcen[0],gcen[1],gcen[2])
  molA.rotate_gcenter(R,gcen)

if (OPT['bond'] != '') and (OPT['ratm'] != ''):
  (natom_sta,natom_ori) = OPT['bond'].split(':')
  ratmlist = []
  ratmlist = OPT['ratm'].split(':')
  natom_sta = int(natom_sta)
  natom_ori = int(natom_ori)
  axis = [0.0,0.0,0.0] 
  make_rotational_axis(axis,molA,natom_sta,natom_ori)
  print "natom_sta %d natom_ori %d"%(natom_sta,natom_ori)
  print "axis %f %f %f"%(axis[0],axis[1],axis[2])
  R = [[0.0 for j in range(3)] for i in range(3)]
  make_rotational_matrix_from_axis_and_angle(R,axis,float(OPT['ang'])) 
  print "R %f %f %f"%(R[0][0],R[0][1],R[0][2])
  print "  %f %f %f"%(R[1][0],R[1][1],R[1][2])
  print "  %f %f %f"%(R[2][0],R[2][1],R[2][2])
  rotate_part_of_molecule(molA,ratmlist,natom_ori,R)


if (OPT['osdf'] != ''):
  molA.write_in_sdf(OPT['osdf'])
