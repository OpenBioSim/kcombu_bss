#!/usr/bin/env python

import sys
import os
import math

def Cal_Inverse_Matrix3D_by_Cramer_Rule(InvA,A):
# 
# 
# <Cramer's Rule>
#
# Inv[A] = 1/|A| transpose[Delta_ij] 
#
# Delta_ij = (-1)**(i+j) * |Aij|
#
# Aij = matrix removing i-th row and j-th column.
#
  det =  A[0][0]*A[1][1]*A[2][2] + A[0][2]*A[1][0]*A[2][1] + A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] 

  print "def %f"%(det)

  InvA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1])
  InvA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1])
  InvA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1])

  InvA[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0])
  InvA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0])
  InvA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0])

  InvA[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
  InvA[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]);
  InvA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]);

  for i in range(3):
    for j in range(3):
      InvA[i][j] /= det;

  return(det)


### MAIN ####

a = []
b = []
c = []

a = [0.4,0.2,0.1]
b = [1,0.3,0.8]
c = [1,1,0.2]

x = [0,0,2]

A = B = C = 0.0
for i in range(3):
  A += (a[i]-x[i])*(a[i]-x[i])
  B += (b[i]-x[i])*(b[i]-x[i])
  C += (c[i]-x[i])*(c[i]-x[i])

A = math.sqrt(A)
B = math.sqrt(B)
C = math.sqrt(C)

print "a %f %f %f A %f"%(a[0],a[1],a[2],A)
print "b %f %f %f B %f"%(b[0],b[1],b[2],B)
print "c %f %f %f C %f"%(c[0],c[1],c[2],C)
print "x %f %f %f"%(x[0],x[1],x[2])

X    = [[0.0 for i in range(3)] for j in range(3)]
invX = [[0.0 for i in range(3)] for j in range(3)]


X[0][0] = b[0] - a[0]
X[0][1] = b[1] - a[1]
X[0][2] = b[2] - a[2]

X[1][0] = b[0] - c[0]
X[1][1] = b[1] - c[1]
X[1][2] = b[2] - c[2]

X[2][0] = c[0] - a[0]
X[2][1] = c[1] - a[1]
X[2][2] = c[2] - a[2]


print "X",X

Cal_Inverse_Matrix3D_by_Cramer_Rule(invX,X)

print "invX",invX
