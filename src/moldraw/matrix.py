##
## <matrix.py>
##  functions for calculating PCA of the molecule and rotate it.
##

import sys
import os
import math

LastModDate = "June 26, 2011"


def cal_eigen_vector_by_Jacobi_Wilkinson(N,A,U):
  ## A and U are N x N matrices
  ##
  ##      [before]           [after] 
  ##  A : input matrix  --> diagonal matrix for eigen value
  ##  U : ------------  --> eigen vector matrix  
  ##    
  ##  Each column vector of U corresponds to the eigen vector for A.
  ## evec0 = |Umat[0][0]|  eval0 = A[0][0]
  ##         |Umat[1][0]|
  ##         |Umat[2][0]|
  ## evec1 = |Umat[0][1]|  eval1 = A[1][1]
  ##         |Umat[1][1]|
  ##         |Umat[2][1]|
  ## evec2 = |Umat[0][2]|  eval2 = A[2][2]
  ##         |Umat[1][2]|
  ##         |Umat[2][2]|


  R     = [[0.0 for i in range(N)] for j in range(N)] 
  TR    = [[0.0 for i in range(N)] for j in range(N)] 
  BUFF  = [[0.0 for i in range(N)] for j in range(N)] 
  sqrt2 = math.sqrt(2.0);

  set_identiry_matrix(N,U)
  max_mini = 0.00000001
  c_max  = 100  
  [mi,mj,max] =  find_max_abs_nondiagonal(N,A)
  c = 0;
  while ((max>max_mini) and (c<c_max)):
    #print "mi %d mj %d max %.10f"%(mi,mj,max)
    #print_matrix(3,A,"A")
    #print_matrix(3,U,"U")
    ## --------- SET of cos sin ----------- ##
    wa = (A[mi][mi] + A[mj][mj])/2.0
    sa = (A[mi][mi] - A[mj][mj])/2.0
    r  = math.sqrt(sa*sa + A[mi][mj]*A[mi][mj])
    if (sa>0.0):
      co =  math.sqrt(1.0+sa/r)/sqrt2
      si =  A[mi][mj]/2.0/r/co;
    else:
      co =  math.sqrt(1.0-sa/r)/sqrt2
      si = - A[mi][mj]/2.0/r/co

    ## -------- SET of Rot Matrix R and TR-----## 
    set_identiry_matrix(N,R)
    R[mi][mi] =  co
    R[mi][mj] = -si
    R[mj][mi] =  si
    R[mj][mj] =  co

    set_transpose_matrix(N,TR,R)

    prod_matrix(N,BUFF,A,R)
    prod_matrix(N,A,TR,BUFF)
    prod_matrix(N,BUFF,U,R)
    copy_matrix(N,U,BUFF)

    [mi,mj,max] = find_max_abs_nondiagonal(N,A)
    c+= 1




def find_max_abs_nondiagonal(N,A):   ## find maximum |A[mi][mj]| = max, and return [mi,mj,max] ## 
  max = math.fabs(A[0][1])
  mi = 0
  mj = 1
  for i in range(N):
    for j in range(N):
      if (i!=j):
        absA = math.fabs(A[i][j]) 
        if (absA>max):
          mi = i
          mj = j
          max = absA
  return([mi,mj,max])


def prod_matrix(N,C,A,B): ## C = A * B ## 
  for i in range(N):
    for j in range(N):
      C[i][j] = 0.0
      for k in range(N):
        C[i][j] += A[i][k] * B[k][j]


def prod_matrix_vec(N,y,A,x): ## y = A * x ## 
  for i in range(N):
    y[i] = 0.0
    for k in range(N):
      y[i] += A[i][k] * x[k]


def copy_matrix(N,A,B): ## A = B ## 
  for i in range(N):
    for j in range(N):
      A[i][j] = B[i][j]

def set_identiry_matrix(N,A):
  for i in range(N):
    for j in range(N):
      if (i==j):
        A[i][j] = 1.0
      else:
        A[i][j] = 0.0

def set_transpose_matrix(N,TA,A):
  for i in range(N):
    for j in range(N):
      TA[j][i] = A[i][j];

def print_matrix(N,A,comment):
  print ">%s"%(comment)
  for i in range(N):
    for j in range(N):
      sys.stdout.write(" %f"%(A[i][j]))
    sys.stdout.write("\n")
  sys.stdout.write("\n")

def length_vec(N,x):
  norm = 0.0
  for i in range(N):
    norm += x[i]*x[i]
  return(math.sqrt(norm))

def Normalize_vec(N,x):
  len = length_vec(N,x)
  nx = [0.0 for i in range(N)]
  for i in range(N):
    nx[i]  = x[i]/len
  return(nx)
