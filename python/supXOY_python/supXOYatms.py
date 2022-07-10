#!/usr/bin/env python

##
## <supXOYatms.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import pdb
import glob

LastModDate = '2019/05/31'


def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        opt_dic[argv[i][1:]] = argv[i+1]


def pick_subset_of_pdb_Molecule(Mo,pattern):
  Msub = pdb.Molecule()
  (Atomname,Resname,Chain,Resnum) = pattern.split(':')
  Msub.Natom = 0 
  Atomname = Atomname.replace(' ','')
  Resname = Resname.replace(' ','')
  Chain = Chain.replace(' ','')
  Resnum = Resnum.replace(' ','')
  print "#pick_subset_of_pdb_Molecule(Mo:%s,pattern '%s:%s:%s:%s')"%(Mo.filename,Atomname, Resname,Chain, Resnum)
  for i in range(Mo.Natom):
    #sys.stdout.write("'%s' '%s' '%s' '%s'\n"%(Mo.atom[i].replace(' ',''), Mo.resi[i].replace(' ',''), Mo.chain[i].replace(' ',''), Mo.rnum[i].replace(' ','')))
    #if (Atomname == '') or (Mo.atom[i].replace(' ','') == Atomname): 
    #  sys.stdout.write(" atomOK")
    #if (Resname == '') or (Mo.resi[i].replace(' ','') == Resname): 
    #  sys.stdout.write(" resiOK")
    #if (Chain == '') or (Mo.chain[i].replace(' ','') == Chain): 
    #  sys.stdout.write(" chOK")
    #if (Resnum == '') or (Mo.rnum[i].replace(' ','') == Resnum): 
    #  sys.stdout.write(" rnumOK")
    #sys.stdout.write("\n") 
    if (Atomname == '') or (Mo.atom[i].replace(' ','') == Atomname): 
      if (Resname == '') or (Mo.resi[i].replace(' ','') == Resname): 
        if (Chain == '') or (Mo.chain[i].replace(' ','') == Chain): 
          if (Resnum == '') or (Mo.rnum[i].replace(' ','') == Resnum): 
            Msub.Natom += 1
            Msub.anum.append(Mo.anum[i])
            Msub.atom.append(Mo.atom[i])
            Msub.resi.append(Mo.resi[i])
            Msub.rnum.append(Mo.rnum[i])
            Msub.posX.append(Mo.posX[i])
            Msub.posY.append(Mo.posY[i])
            Msub.posZ.append(Mo.posZ[i])
            Msub.chain.append(Mo.chain[i])
            Msub.AHtype.append(Mo.AHtype[i])
            Msub.line.append(Mo.line[i])

  return(Msub)
  pass


def extract_pos3_from_mol(mol,pattern):
  if (mol.Natom !=1):
    print "#ERROR: Natom for '%s' is %d.\n"%(pattern, mol.Natom)
    sys.exit(1)
  else:
    p  = [0.0, 0.0, 0.0]
    p[0] = mol.posX[0] 
    p[1] = mol.posY[0] 
    p[2] = mol.posZ[0] 
    return(p)

def cross_vec3(a,b):   ## c = a x b ## 
  c = [0.0, 0.0, 0.0]
  c[0] = a[1]*b[2] - a[2]*b[1]
  c[1] = a[2]*b[0] - a[0]*b[2]
  c[2] = a[0]*b[1] - a[1]*b[0] 
  return(c)


def dot_prod3(a,b):
  return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])



def sub_vec3(a,b):
  c = [0.0, 0.0, 0.0]
  c[0] = a[0] - b[0]
  c[1] = a[1] - b[1]
  c[2] = a[2] - b[2] 
  return(c)

def normalize_vec3(v):
  norm = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
  norm = math.sqrt(norm)
  if (norm<=0.0):
     return(0.0);
  v[0] /= norm  
  v[1] /= norm
  v[2] /= norm;
  return(norm);

def multi_mat3x3(A,B): ## C = A B ## 
  C   = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  for i in range(3):
    for j in range(3):
      C[i][j] = 0.0;
      for k in range(3):
        C[i][j] += A[i][k] * B[k][j];
  return(C)

def make_transpose_matrix3x3(A):
  TA  = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  for i in range(3):
    for j in range(3):
       TA[i][j] = A[j][i];
  return(TA)



def multiply_mat3x3_by_vec3(R,x): ## y = R*x  ## 
  y = [0.0, 0.0, 0.0]
  for i in range(3):
    y[i] = 0.0;
    for j in range(3):
      y[i] += R[i][j]*x[j]
  return(y)



def make_rotmatrix_from_rot_axis_and_angle(axis,angle_degree):
  R   = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  S   = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  SS   = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

  a = axis[0] 
  b = axis[1]
  c = axis[2]
  S[0][0] = 0.0 
  S[0][1] = -c  
  S[0][2] = b
  S[1][0] = c
  S[1][1] = 0.0
  S[1][2] = -a
  S[2][0] = -b  
  S[2][1] =  a
  S[2][2] = 0.0;


  SS[0][0] = -c*c-b*b
  SS[0][1] = a*b
  SS[0][2] = a*c
  SS[1][0] = a*b
  SS[1][1] = -c*c-a*a
  SS[1][2] = b*c
  SS[2][0] = a*c
  SS[2][1] = b*c
  SS[2][2] = -a*a-b*b;

  angle_radian = angle_degree*math.pi/180.0 
  sint  = math.sin(angle_radian)
  ccost = 1.0 - math.cos(angle_radian)
 
  for i in range(3):
    for j in range(3):
      R[i][j] = sint*S[i][j] + ccost*SS[i][j]
      if (i==j):
        R[i][j] += 1.0

  return(R)


def rotate_y_around_x_o_axis(x,o,y,rangle_degree):
  raxis = sub_vec3(o,x)
  normalize_vec3(raxis)
  Rmat = make_rotmatrix_from_rot_axis_and_angle(raxis,rangle_degree)
  yO      = sub_vec3(y,x)
  new_y =  multiply_mat3x3_by_vec3(Rmat,yO) 
  new_y[0] = new_y[0] + x[0]
  new_y[1] = new_y[1] + x[1]
  new_y[2] = new_y[2] + x[2]
  return(new_y) 



def print_matrix3(A,title):
 print "#>%s\n"%(title);
 print "#%8.3lf %8.3lf %8.3lf"%(A[0][0],A[0][1],A[0][2])
 print "#%8.3lf %8.3lf %8.3lf"%(A[1][0],A[1][1],A[1][2])
 print "#%8.3lf %8.3lf %8.3lf"%(A[2][0],A[2][1],A[2][2])


def make_XYZ_from_xoy_atoms(x,o,y):
  X = sub_vec3(x,o)
  normalize_vec3(X)
  b = sub_vec3(y,o)
  normalize_vec3(b)
  Xdotb = dot_prod3(X,b)
  Y = [0.0, 0.0, 0.0]
  for i in range(3):
    Y[i] = b[i] - Xdotb * X[i] 
  normalize_vec3(Y)
  Z = cross_vec3(X,Y)
  return((X,Y,Z))


def make_Rmat_from_tXYZ_and_rXYZ(tXX,tYY,tZZ,rXX,rYY,rZZ):
  Ut   = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  Ur   = [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

  for i in range(3):
    Ut[i][0] = tXX[i]
    Ut[i][1] = tYY[i]
    Ut[i][2] = tZZ[i]

  for i in range(3):
    Ur[i][0] = rXX[i]
    Ur[i][1] = rYY[i]
    Ur[i][2] = rZZ[i]

  print_matrix3(Ut,"Ut")
  print_matrix3(Ur,"Ur")

  UtT  = make_transpose_matrix3x3(Ut)
  Rmat = multi_mat3x3(Ur,UtT)
  return(Rmat)

def transpose_atoms_by_Rmat_Oorig_Onew(mol,Rmat,Oorg,Onew):
  orgpos  = [0.0, 0.0, 0.0]
  newpos = [0.0, 0.0, 0.0]
 
  for i in range (mol.Natom):
    orgpos[0] = mol.posX[i]  
    orgpos[1] = mol.posY[i]  
    orgpos[2] = mol.posZ[i]  
    org_O = sub_vec3(orgpos,Oorg)
    new_O =  multiply_mat3x3_by_vec3(Rmat,org_O) 
    mol.posX[i] = new_O[0] + Onew[0] 
    mol.posY[i] = new_O[1] + Onew[1] 
    mol.posZ[i] = new_O[2] + Onew[2] 


def write_in_pdb(mol,ofname,NEWRNUM='',COMMAND='',DATE=''):
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_in_pdb()-->'%s'"%(ofname)
    of = open(ofname,'w')
#          1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM   1071  HZ2 LYS A  64     -10.286   7.104  12.073  1.00  0.00           H
#ATOM   1072  HZ3 LYS A  64     -10.003   6.378  10.583  1.00  0.00           H
  of.write("REMARK COMMAND %s\n"%(COMMAND))
  of.write("REMARK DATE    %s\n"%(DATE))
  for i in range(mol.Natom):
    #print i,mol.line[i]
    line = mol.line[i]
    head_line = mol.line[i][0:30]
    if (NEWRNUM != ''):
      head_line = "%s%4s%s"%(head_line[0:22],NEWRNUM,head_line[26:])
    tail_line = mol.line[i][54:]
    of.write("%s%8.3f%8.3f%8.3f%s\n"%(head_line,mol.posX[i],mol.posY[i],mol.posZ[i],tail_line))
  if (ofname!='-'):
    of.close()


 

################
##### MAIN #####
################

OPT = {}
OPT['iref'] = 'N-glycan_1gya_1.pdb'
OPT['itar'] = 'MAN.pdb'
OPT['tX'] = 'C1:MAN:A:1'
OPT['tO'] = 'O1:MAN:A:1'
OPT['tY'] = 'HO1:MAN:A:1'
OPT['rX'] = 'HO2:MAN:A:110'
OPT['rO'] = 'O2:MAN:A:110'
OPT['rY'] = 'C2:MAN:A:110'
OPT['opdb'] = '-'
OPT['tnewrnum'] = ''
OPT['tYrotang'] = '0.0'

if (len(sys.argv)<2):
  print "supXOYatms.py <options>"
  print "  superimpose 'tar' molecule on 'ref' molecule, using X-O-Y atom matching."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print "  -iref: input reference PDB file [%s]"%(OPT['iref']) 
  print "  -itar: input target    PDB file [%s]"%(OPT['itar']) 
  print "  -rX,-rO,-rY: reference atoms pattern (atomname):(resname):(chain):(resnum)"
  print "                rX [%s] rO [%s] rY [%s]"%(OPT['rX'],OPT['rO'],OPT['rY'])
  print "  -tX,-tO,-tY: target atoms pattern (atomname):(resname):(chain):(resnum)"
  print "                tX [%s] tO [%s] tY [%s]"%(OPT['tX'],OPT['tO'],OPT['tY'])
  print "  -tnewrnum  : target new residue number [%s]"%(OPT['tnewrnum'])
  print "  -tYrotang  : rotation angle around rX-rO axis for tY [%s]"%(OPT['tYrotang'])
  print "  -opdb  : output PDB filename [%s]"%(OPT['opdb'])
  sys.exit()

read_option(sys.argv, OPT)
if (OPT['opdb'] == '-'):
  of = sys.stdout
else:
  print "#write_to --> '%s'"%(OPT['opdb'])
  of = open(OPT['opdb'],'w')



R = pdb.Molecule()
R.read(OPT['iref'],AHtype='B')

T = pdb.Molecule()
T.read(OPT['itar'],AHtype='B')


rXatms = pick_subset_of_pdb_Molecule(R,OPT['rX'])
rX     = extract_pos3_from_mol(rXatms,OPT['rX'])

rOatms = pick_subset_of_pdb_Molecule(R,OPT['rO'])
rO     = extract_pos3_from_mol(rOatms,OPT['rO'])

rYatms = pick_subset_of_pdb_Molecule(R,OPT['rY'])
rY     = extract_pos3_from_mol(rYatms,OPT['rY'])

print "#rX",rX
print "#rO",rO
print "#rY",rY

tXatms = pick_subset_of_pdb_Molecule(T,OPT['tX'])
tX     = extract_pos3_from_mol(tXatms,OPT['tX'])

tOatms = pick_subset_of_pdb_Molecule(T,OPT['tO'])
tO     = extract_pos3_from_mol(tOatms,OPT['tO'])

tYatms = pick_subset_of_pdb_Molecule(T,OPT['tY'])
tY     = extract_pos3_from_mol(tYatms,OPT['tY'])

print "#tX",tX
print "#tO",tO
print "#tY",tY

if (OPT['tYrotang'] != '0.0'):
  new_tY = rotate_y_around_x_o_axis(tX,tO,tY,float(OPT['tYrotang']))
  print "#    tY %8.3f%8.3f%8.3f"%(tY[0],tY[1],tY[2])
  print "#new_tY %8.3f%8.3f%8.3f"%(new_tY[0],new_tY[1],new_tY[2])
  tY[0] = new_tY[0]
  tY[1] = new_tY[1]
  tY[2] = new_tY[2]
  pass


(rXX,rYY,rZZ) = make_XYZ_from_xoy_atoms(rX,rO,rY)
(tXX,tYY,tZZ) = make_XYZ_from_xoy_atoms(tX,tO,tY)

Rmat = make_Rmat_from_tXYZ_and_rXYZ(tXX,tYY,tZZ,rXX,rYY,rZZ)
print_matrix3(Rmat,"Rmat")

transpose_atoms_by_Rmat_Oorig_Onew(T,Rmat,tO,rO)

write_in_pdb(T,OPT['opdb'],NEWRNUM=OPT['tnewrnum'],COMMAND=OPT['COMMAND'],DATE=OPT['START_DATE'])

sys.exit(1)

