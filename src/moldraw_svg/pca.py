##
## <pca.py>
##  functions for calculating PCA of the molecule and rotate it.
##  This is mainly for 'moldraw.cgi'

import sys
import os
import math
import molecule
import matrix

LastModDate = "Jun 26, 2011"

def rotate_molecule_XY_flat_plane_by_PCA(mol):
  Cmat = cal_covariance_matrix(mol)
  #print_matrix(3,Cmat,"Cmat")
  Lmat = [[0.0 for i in range(3)]for j in range(3)]
  Umat = [[0.0 for i in range(3)]for j in range(3)]
  Rmat = [[0.0 for i in range(3)]for j in range(3)]
  matrix.copy_matrix(3,Lmat,Cmat)
  matrix.cal_eigen_vector_by_Jacobi_Wilkinson(3,Lmat,Umat)
  make_rotational_matrix_XY_flat_plane(Rmat,Lmat,Umat)
  #print_matrix(3,Lmat,"Lmat")
  #print_matrix(3,Umat,"Umat")
  #print_matrix(3,Rmat,"Rmat")
  rotate_molecule(mol,Rmat,mol.gpos)

def cal_Rmat_and_gpos_molecule_XY_flat_plane_by_PCA(mol):
  Cmat = cal_covariance_matrix(mol)
  #print_matrix(3,Cmat,"Cmat")
  Lmat = [[0.0 for i in range(3)]for j in range(3)]
  Umat = [[0.0 for i in range(3)]for j in range(3)]
  Rmat = [[0.0 for i in range(3)]for j in range(3)]
  matrix.copy_matrix(3,Lmat,Cmat)
  matrix.cal_eigen_vector_by_Jacobi_Wilkinson(3,Lmat,Umat)
  make_rotational_matrix_XY_flat_plane(Rmat,Lmat,Umat)
  return(Rmat,mol.gpos)


def rotate_two_molecules_XY_flat_plane_by_PCA(molA,molB):
  gpos = [0.0, 0.0, 0.0]
  Cmat =  cal_covariance_matrix_for_two_molecules(gpos,molA,molB)
  #print_matrix(3,Cmat,"Cmat")
  Lmat = [[0.0 for i in range(3)]for j in range(3)]
  Umat = [[0.0 for i in range(3)]for j in range(3)]
  Rmat = [[0.0 for i in range(3)]for j in range(3)]
  matrix.copy_matrix(3,Lmat,Cmat)
  matrix.cal_eigen_vector_by_Jacobi_Wilkinson(3,Lmat,Umat)
  make_rotational_matrix_XY_flat_plane(Rmat,Lmat,Umat)
  #print_matrix(3,Lmat,"Lmat")
  #print_matrix(3,Umat,"Umat")
  #print_matrix(3,Rmat,"Rmat")
  rotate_molecule(molA,Rmat,gpos)
  rotate_molecule(molB,Rmat,gpos)


def cal_Rmat_and_gpos_for_two_molecules_XY_flat_plane_by_PCA(molA,molB):
  gpos = [0.0, 0.0, 0.0]
  Cmat =  cal_covariance_matrix_for_two_molecules(gpos,molA,molB)
  #print_matrix(3,Cmat,"Cmat")
  Lmat = [[0.0 for i in range(3)]for j in range(3)]
  Umat = [[0.0 for i in range(3)]for j in range(3)]
  Rmat = [[0.0 for i in range(3)]for j in range(3)]
  matrix.copy_matrix(3,Lmat,Cmat)
  matrix.cal_eigen_vector_by_Jacobi_Wilkinson(3,Lmat,Umat)
  make_rotational_matrix_XY_flat_plane(Rmat,Lmat,Umat)
  #print_matrix(3,Rmat,"Rmat")
  return(Rmat,gpos)


def cal_covariance_matrix(mol):
  mol.gpos = [0.0,0.0,0.0]
  rpos = [0.0,0.0,0.0]
  Cmat = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
  for a in range(mol.Natom):  
    mol.gpos[0] += mol.atoms[a].x
    mol.gpos[1] += mol.atoms[a].y
    mol.gpos[2] += mol.atoms[a].z
  for i in range(3):
    mol.gpos[i] /= mol.Natom 

  for a in range(mol.Natom):  
    rpos[0] = mol.atoms[a].x - mol.gpos[0] 
    rpos[1] = mol.atoms[a].y - mol.gpos[1] 
    rpos[2] = mol.atoms[a].z - mol.gpos[2] 
    for i in range(3):
      for j in range(3):
        Cmat[i][j] += rpos[i]*rpos[j]
  for i in range(3):
    for j in range(3):
      Cmat[i][j] /= mol.Natom
  return(Cmat)

def cal_covariance_matrix_for_two_molecules(gpos,molA,molB):
  rpos = [0.0,0.0,0.0]
  gpos[0] = 0.0 
  gpos[1] = 0.0 
  gpos[2] = 0.0 
  Cmat = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
  NatomAB = molA.Natom + molB.Natom
  for a in range(molA.Natom):  
    gpos[0] += molA.atoms[a].x
    gpos[1] += molA.atoms[a].y
    gpos[2] += molA.atoms[a].z
  for a in range(molB.Natom):  
    gpos[0] += molB.atoms[a].x
    gpos[1] += molB.atoms[a].y
    gpos[2] += molB.atoms[a].z

  for i in range(3):
    gpos[i] /= NatomAB

  for a in range(molA.Natom):  
    rpos[0] = molA.atoms[a].x - gpos[0] 
    rpos[1] = molA.atoms[a].y - gpos[1] 
    rpos[2] = molA.atoms[a].z - gpos[2] 
    for i in range(3):
      for j in range(3):
        Cmat[i][j] += rpos[i]*rpos[j]

  for a in range(molB.Natom):  
    rpos[0] = molB.atoms[a].x - gpos[0] 
    rpos[1] = molB.atoms[a].y - gpos[1] 
    rpos[2] = molB.atoms[a].z - gpos[2] 
    for i in range(3):
      for j in range(3):
        Cmat[i][j] += rpos[i]*rpos[j]

  for i in range(3):
    for j in range(3):
      Cmat[i][j] /= NatomAB

  return(Cmat)




def make_rotational_matrix_XY_flat_plane(Rmat,Lmat,Umat):
  # x-axis : the largest  width  (largest  eigen value)
  # y-axis : the middle   width  (middle   eigen value)
  # z-axis : the smallest width  (smallest eigen value)
  la = [0.0,0.0,0.0]
  index  = [0,1,2] 
  la[0] = Lmat[0][0]
  la[1] = Lmat[1][1]
  la[2] = Lmat[2][2]
  sindex = sorted(index,lambda x,y:cmp(Lmat[y][y],Lmat[x][x]))
  #print index
  #print sindex
  for i in range(3):
    for j in range(3):
      Rmat[i][j] = Umat[j][sindex[i]]


def make_rotational_matrix_xyz_rot_angle_degree(Rmat,th_x,th_y,th_z):
  Rx  = [[0.0 for i in range(3)]for j in range(3)]
  Ry  = [[0.0 for i in range(3)]for j in range(3)]
  Rz  = [[0.0 for i in range(3)]for j in range(3)]
  Ryx = [[0.0 for i in range(3)]for j in range(3)]
  cx = math.cos(math.pi*th_x/180.0)
  sx = math.sin(math.pi*th_x/180.0)
  
  cy = math.cos(math.pi*th_y/180.0)
  sy = math.sin(math.pi*th_y/180.0)

  cz = math.cos(math.pi*th_z/180.0)
  sz = math.sin(math.pi*th_z/180.0)

  Rx[0][0] = 1.0
  Rx[1][1] = cx 
  Rx[1][2] = -sx 
  Rx[2][1] = sx 
  Rx[2][2] = cx
 
  Ry[1][1] = 1.0
  Ry[0][0] = cy 
  Ry[0][2] = -sy
  Ry[2][0] = sy 
  Ry[2][2] = cy

  Rz[2][2] = 1.0
  Rz[0][0] = cz 
  Rz[0][1] = -sz
  Rz[1][0] = sz 
  Rz[1][1] = cz

  matrix.prod_matrix(3,Ryx,Ry,Rx)
  matrix.prod_matrix(3,Rmat,Rz,Ryx)




def cal_rot_angle_xyz_from_rot_matrix(Rmat):
  if (math.fabs(Rmat[2][0]-1.0)>0.00000001):
    th_y = 180.0*math.asin(Rmat[2][0])/math.pi
    th_z = 180.0*math.atan2(Rmat[1][0],Rmat[0][0])/math.pi
    th_x = 180.0*math.atan2(Rmat[2][1],Rmat[2][2])/math.pi
  else:
    print "woops Rmat01 %f Rmat12 %f"%(Rmat[0][1],Rmat[0][2])
    if (Rmat[2][0]>0):
      th_y =  90.0
      th_x = 180.0*math.atan2(-Rmat[0][1],-Rmat[0][2])/math.pi
    else:
      th_y = -90.0
      th_x = 180.0*math.atan2(Rmat[0][1],Rmat[0][2])/math.pi
    th_z = 0.0
  return([th_x,th_y,th_z])

def rotate_molecule(mol,Rmat,gpos):
  pos  = [0.0, 0.0, 0.0]
  rpos = [0.0, 0.0, 0.0]
  for a in range(mol.Natom):
    pos[0] = mol.atoms[a].x - gpos[0]
    pos[1] = mol.atoms[a].y - gpos[1]
    pos[2] = mol.atoms[a].z - gpos[2]
    matrix.prod_matrix_vec(3,rpos,Rmat,pos)
    mol.atoms[a].x = rpos[0] +  gpos[0]
    mol.atoms[a].y = rpos[1] +  gpos[1]
    mol.atoms[a].z = rpos[2] +  gpos[2]




#############
### MAIN ####
#############


def _main():
  if (len(sys.argv)<2):
    print "pca.py [sdf or pdb file]"
    print "pca.py [sdf or pdb file A] [sdf or pdb file B]"
    sys.exit(1)

  if (len(sys.argv)==2):
    mol = molecule.Molecule()
    if (sys.argv[1].find(".sdf")!=-1):
      mol.read_in_sdf(sys.argv[1])
    elif (sys.argv[1].find(".pdb")!=-1):
      mol.read_in_pdb(sys.argv[1])
      mol.guess_bonds_from_atom_xyz()
    elif (sys.argv[1].find(".mol2")!=-1):
      mol.read_in_mol2(sys.argv[1])
    else:
      print "#no input file"
      sys.exit(1)

    (Rmat,gpos) = cal_Rmat_and_gpos_molecule_XY_flat_plane_by_PCA(mol)
    print "Rmat0 %f %f %f"%(Rmat[0][0], Rmat[0][1], Rmat[0][2])
    print "Rmat1 %f %f %f"%(Rmat[1][0], Rmat[1][1], Rmat[1][2])
    print "Rmat2 %f %f %f"%(Rmat[2][0], Rmat[2][1], Rmat[2][2])
    print "gpos %f %f %f"%(gpos[0], gpos[1], gpos[2])
    
    rotate_molecule_XY_flat_plane_by_PCA(mol)
    mol.write_in_sdf("out.sdf")
    th_y = math.asin(Rmat[2][0])
    th_z = math.atan2(Rmat[1][0],Rmat[0][0])
    th_x = math.atan2(Rmat[2][1],Rmat[2][2])
    print("x %f y %f z %f"%(th_x,th_y,th_z)) 
    print("rotate x %.0f"%(180.0*th_x/math.pi)) 
    print("rotate y %.0f"%(180.0*th_y/math.pi)) 
    print("rotate z %.0f"%(180.0*th_z/math.pi)) 
    [th_x,th_y,th_z] = cal_rot_angle_xyz_from_rot_matrix(Rmat)
    print("rotate x %.0f"%(th_x)) 
    print("rotate y %.0f"%(th_y)) 
    print("rotate z %.0f"%(th_z)) 

    make_rotational_matrix_xyz_rot_angle_degree(Rmat,0.0,0.0,0.0)
    matrix.print_matrix(3,Rmat,"Rmat 0.0 0.0 0.0")
    [th_x,th_y,th_z] = cal_rot_angle_xyz_from_rot_matrix(Rmat)
    print("rotate x %.0f"%(th_x)) 
    print("rotate y %.0f"%(th_y)) 
    print("rotate z %.0f"%(th_z)) 

    th_x = 90.0
    th_y = -90.0
    th_z = 90.0
    make_rotational_matrix_xyz_rot_angle_degree(Rmat,th_x, th_y, th_z)
    matrix.print_matrix(3,Rmat,"Rmat %.1f %.1f %.1f"%(th_x, th_y, th_z))
    [th_x,th_y,th_z] = cal_rot_angle_xyz_from_rot_matrix(Rmat)
    print("rotate x %.0f"%(th_x)) 
    print("rotate y %.0f"%(th_y)) 
    print("rotate z %.0f"%(th_z)) 
    make_rotational_matrix_xyz_rot_angle_degree(Rmat,th_x,th_y,th_z)
    matrix.print_matrix(3,Rmat,"Rmat_regenerated")

  if (len(sys.argv)==3):
    molA = molecule.Molecule()
    if (sys.argv[1].find(".sdf")!=-1):
      molA.read_in_sdf(sys.argv[1])
    elif (sys.argv[1].find(".pdb")!=-1):
      molA.read_in_pdb(sys.argv[1])
      molA.guess_bonds_from_atom_xyz()
    else:
      print "#no input file"
      sys.exit(1)

    molB = molecule.Molecule()
    if (sys.argv[2].find(".sdf")!=-1):
      molB.read_in_sdf(sys.argv[1])
    elif (sys.argv[2].find(".pdb")!=-1):
      molB.read_in_pdb(sys.argv[2])
      molB.guess_bonds_from_atom_xyz()
    else:
      print "#no input file"
      sys.exit(1)

    (Rmat,gpos) = cal_Rmat_and_gpos_for_two_molecules_XY_flat_plane_by_PCA(molA,molB)
    print "Rmat0 %f %f %f"%(Rmat[0][0], Rmat[0][1], Rmat[0][2])
    print "Rmat1 %f %f %f"%(Rmat[1][0], Rmat[1][1], Rmat[1][2])
    print "Rmat2 %f %f %f"%(Rmat[2][0], Rmat[2][1], Rmat[2][2])
    print "gpos %f %f %f"%(gpos[0], gpos[1], gpos[2])

    rotate_two_molecules_XY_flat_plane_by_PCA(molA,molB)
    molA.write_in_sdf("outA.sdf")
    molB.write_in_sdf("outB.sdf")

#  Cmat = cal_covariance_matrix(mol)
#  Lmat = [[0.0 for i in range(3)]for j in range(3)]
#  Umat = [[0.0 for i in range(3)]for j in range(3)]
#  Rmat = [[0.0 for i in range(3)]for j in range(3)]
#  copy_matrix(3,Lmat,Cmat)
#  print_matrix(3,Lmat,"Lmat")
#  print_matrix(3,Umat,"Umat")
#  cal_eigen_vector_by_Jacobi_Wilkinson(3,Lmat,Umat)
#  print_matrix(3,Lmat,"Lmat")
#  print_matrix(3,Umat,"Umat")
#  x = [0.0, 0.0, 0.0]
#  y = [0.0, 0.0, 0.0]
#  x[0] = Umat[0][0]
#  x[1] = Umat[1][0]
#  x[2] = Umat[2][0]
#  prod_matrix_vec(3,y,Cmat,x)
#  print "x  %f %f %f"%(x[0],x[1],x[2])
#  print "y  %f %f %f"%(y[0],y[1],y[2])
#  nx =  Normalize_vec(3,x)
#  ny =  Normalize_vec(3,y)
#  lenx = length_vec(3,x) 
#  leny = length_vec(3,y) 
#  print "lenx %f leny %f"%(lenx,leny)
#  print "ny %f %f %f"%(nx[0],nx[1],nx[2])
#  print "ny %f %f %f"%(ny[0],ny[1],ny[2])
#  make_rotational_matrix_XY_flat_plane(Rmat,Lmat,Umat)
#  #set_transpose_matrix(3,Rmat,Umat)
#  print_matrix(3,Rmat,"Rmat")
#  rotate_molecule(mol,Rmat,mol.gpos)
#  mol.write_in_sdf("out.sdf")

if __name__ == '__main__': _main()

