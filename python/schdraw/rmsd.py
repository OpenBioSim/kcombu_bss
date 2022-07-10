#!/usr/bin/python
#
# <rmsd.py>
#
#  functions for calculating root-mean-square deviation
#

import sys
import os
import math
from datetime import datetime
import molecule
import matrix

LastModDate = "Jun 27, 2011"




def rotate_molA_onto_molB_for_matched_atoms_with_minRMSD(molA,molB,anumA,anumB):

  ## (1) cal gravity center gA and gB 
  gA = [0.0, 0.0, 0.0] 
  gB = [0.0, 0.0, 0.0] 
  for i in range(len(anumA)):
    p = anumA[i]
    q = anumB[i]
    gA[0] += molA.atoms[p].x 
    gA[1] += molA.atoms[p].y 
    gA[2] += molA.atoms[p].z
 
    gB[0] += molB.atoms[q].x 
    gB[1] += molB.atoms[q].y 
    gB[2] += molB.atoms[q].z

  for k in range(3):
    gA[k] = gA[k]/len(anumA)
    gB[k] = gB[k]/len(anumA)

 ## (2) Cal 4x4 matrix B ##
  a = [0.0,0.0,0.0]
  b = [0.0,0.0,0.0]
  B = [[0.0 for i in range(4)] for j in range(4)] 
  U = [[0.0 for i in range(4)] for j in range(4)] 
  R = [[0.0 for i in range(3)] for j in range(3)] 
  B[0][0] = B[1][1] = B[2][2] =B[3][3] = 0.0

  for i in range(len(anumA)):
    p = anumA[i]
    q = anumB[i]

    a[0] = (molB.atoms[q].x - gB[0]) + (molA.atoms[p].x - gA[0])
    a[1] = (molB.atoms[q].y - gB[1]) + (molA.atoms[p].y - gA[1])
    a[2] = (molB.atoms[q].z - gB[2]) + (molA.atoms[p].z - gA[2])

    b[0] = (molB.atoms[q].x - gB[0]) - (molA.atoms[p].x - gA[0])
    b[1] = (molB.atoms[q].y - gB[1]) - (molA.atoms[p].y - gA[1])
    b[2] = (molB.atoms[q].z - gB[2]) - (molA.atoms[p].z - gA[2])

    B[0][0] +=  b[0]*b[0] + b[1]*b[1] + b[2]*b[2]
    B[1][1] +=  b[0]*b[0] + a[1]*a[1] + a[2]*a[2]
    B[2][2] +=  a[0]*a[0] + b[1]*b[1] + a[2]*a[2]
    B[3][3] +=  a[0]*a[0] + a[1]*a[1] + b[2]*b[2]

    B[0][1] +=  a[2]*b[1] - a[1]*b[2]
    B[0][2] += -a[2]*b[0] + a[0]*b[2]
    B[0][3] +=  a[1]*b[0] - a[0]*b[1]
    B[1][2] +=  b[0]*b[1] - a[0]*a[1]
    B[1][3] +=  b[0]*b[2] - a[0]*a[2]
    B[2][3] +=  b[1]*b[2] - a[1]*a[2]

  for i in range(4):
    for j in range(4):
      B[i][j] /= len(anumA) 

  B[1][0] = B[0][1]
  B[2][0] = B[0][2]
  B[3][0] = B[0][3]
  B[2][1] = B[1][2]
  B[3][1] = B[1][3]
  B[3][2] = B[2][3]

  #print "matrixB\n",B

  ### (3) cal eig_vector for the matrix B, and chose min_eigen_value eigen vector as the quaternion ###
#  [eval,evec] = numpy.linalg.eig(B)
  matrix.cal_eigen_vector_by_Jacobi_Wilkinson(4,B,U)

  min_eval =  -1.0
  min_i = 0
  for i in range(4):
    if (B[i][i]>=0.0): 
      if (min_eval<0.0) or (B[i][i]<min_eval):
        min_i = i
        min_eval = B[i][i]
  #print "min_i %d min_eval %f"%(min_i, min_eval)

  q = [0.0,0.0,0.0,0.0]
  for i in range(4):
    q[i] = U[i][min_i] 
 
  ### (4) calcualte rotationa matrix from the quaternion ###

  R[0][0] = 2.0*q[0]*q[0] + 2.0*q[1]*q[1] - 1.0
  R[0][1] = 2.0*q[1]*q[2] - 2.0*q[0]*q[3]
  R[0][2] = 2.0*q[1]*q[3] + 2.0*q[0]*q[2]

  R[1][0] = 2.0*q[1]*q[2] + 2.0*q[0]*q[3]
  R[1][1] = 2.0*q[0]*q[0] + 2.0*q[2]*q[2] - 1.0
  R[1][2] = 2.0*q[2]*q[3] - 2.0*q[0]*q[1]

  R[2][0] = 2.0*q[1]*q[3] - 2.0*q[0]*q[2]
  R[2][1] = 2.0*q[2]*q[3] + 2.0*q[0]*q[1]
  R[2][2] = 2.0*q[0]*q[0] + 2.0*q[3]*q[3] - 1.0

  ### (5) rotate molecule A onto molecule B ###
  pos  = [0.0, 0.0, 0.0]
  rpos = [0.0, 0.0, 0.0]
  for p in range (molA.Natom):
    pos[0] = molA.atoms[p].x - gA[0]
    pos[1] = molA.atoms[p].y - gA[1]
    pos[2] = molA.atoms[p].z - gA[2]
    matrix.prod_matrix_vec(3,rpos,R,pos)
    molA.atoms[p].x = rpos[0] + gB[0]
    molA.atoms[p].y = rpos[1] + gB[1]
    molA.atoms[p].z = rpos[2] + gB[2]
  ### (6) calculate RMS for matched molecule ##
  rmsd = 0.0 
  for i in range(len(anumA)):
    p = anumA[i]
    q = anumB[i]
    dx = molA.atoms[p].x  - molB.atoms[q].x
    dy = molA.atoms[p].y  - molB.atoms[q].y
    dz = molA.atoms[p].z  - molB.atoms[q].z
    rmsd += dx*dx + dy*dy + dz*dz
  rmsd = math.sqrt(rmsd/len(anumA))

  return(rmsd)


#############
### MAIN ####
#############

def _main():
  if (len(sys.argv)<2):
    print "rmsd.py [sdffileA] [sdffileB]"
    print " (Natom of fileA and fileB should be the same.)"
    sys.exit(1)
    pass 

  molA = molecule.Molecule()
  molA.read_in_sdf(sys.argv[1])
  
  molB = molecule.Molecule()
  molB.read_in_sdf(sys.argv[2])
  anumA = [i for i in range(molA.Natom)]
  anumB = [i for i in range(molA.Natom)]
  rmsd = rotate_molA_onto_molB_for_matched_atoms_with_minRMSD(molA,molB,anumA,anumB)
  molA.write_in_sdf("rotA.sdf")
  print "RMSD %f"%(rmsd)

if __name__ == '__main__': _main()

