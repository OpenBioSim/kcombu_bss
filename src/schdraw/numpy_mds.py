##
## <numpy_mds.py>
##
##  for Multi-Dimensional Scaling
##  using python "numpy" package 
##

import sys
import os
import math 
import random

import numpy
import matrix 

LastModDate = "July 4, 2011"

def metric_MDS(N,dmat,xpos,ypos,zpos,Sim2Dis='F'):
  print "#metric_MDS(N,dmat,xpos,ypos,zpos,Sim2Dis='%s')"%(Sim2Dis)

  DDmat    = [[0.0 for i in range(N)]for j in range(N)]
  #DDmat    = numpy.matrix([[0.0 for i in range(N)] for j in range(N)])
  aveDDrow = [0.0 for i in range(N)]
  #Bmat     = [[0.0 for i in range(N)]for j in range(N)]
  #Vmat     = [[0.0 for i in range(N)]for j in range(N)]
  Bmat    = numpy.matrix([[0.0 for i in range(N)] for j in range(N)])
  Vmat    = numpy.matrix([[0.0 for i in range(N)] for j in range(N)])


  ## [1] Making B[i][j] = -0.5 *(DD[i][j] - DD[i][*] - DD[*][j] + DD[*][*]) ##
  print "## [1] Making B[i][j] = -0.5 *(DD[i][j] - DD[i][*] - DD[*][j] + DD[*][*]) ##"

  for i in range(N):
    DDmat[i][i] = 0.0
    for j in range(i+1,N):
      if (Sim2Dis=='T'):
        dis = 1.0 - dmat[i][j]
      else:
        dis = dmat[i][j]
      #print "i %d j %d dmat %f"%(i,j,dmat[i][j])
      DDmat[i][j] = DDmat[j][i] = dis * dis 

  aveDDall = 0.0
  for i in range(N):
    aveDDrow[i] = 0.0
    for j in range(N):
      aveDDrow[i] += DDmat[i][j] 
      aveDDall    += DDmat[i][j] 
    aveDDrow[i] /= N 

  aveDDall /= (N*N)

  for i in range(N):
    for j in range(N):
      Bmat[i,j] = -0.5 *(DDmat[i][j] - aveDDrow[i] - aveDDrow[j] + aveDDall)

  #print "Bmat",Bmat

  ## [2] Calculate eigen vector for Bmat ##
  print "## [2] Calculate eigen vector for Bmat ##"
  #matrix.cal_eigen_vector_by_Jacobi_Wilkinson(N,Bmat,Vmat)
  (eigval,eigvec) = numpy.linalg.eig(Bmat)

  index = [i for i in range(N)] 
  sindex = sorted(index,lambda x,y:cmp(eigval[y],eigval[x]))
  for ii in range(5):
    i = sindex[ii] 
    print ">slambda%d %f"%(i,eigval[i]) 

  ## [3] Set XYZ-coordinate pos[i] = sqrt(eval) * evec[i] for the three largest eval ##
  j0 = sindex[0]
  j1 = sindex[1]
  j2 = sindex[2]
  sqrtLamX = sqrtLamY = sqrtLamZ = 0.0 
  if (eigval[j0]>0.0):
    sqrtLamX = math.sqrt(eigval[j0])
  if (eigval[j1]>0.0):
    sqrtLamY = math.sqrt(eigval[j1])
  if (eigval[j2]>0.0):
    sqrtLamZ = math.sqrt(eigval[j2])
 
  for i in range(N):
    xpos[i] = sqrtLamX*eigvec[i,j0] 
    ypos[i] = sqrtLamY*eigvec[i,j1] 
    zpos[i] = sqrtLamZ*eigvec[i,j2] 


#############
### MAIN ####
#############


def _main():
  if (len(sys.argv)<2):
    print "numpy_mds.py [Ndata] [Mdim]"
    sys.exit(1)

  N = int(sys.argv[1])
  M = int(sys.argv[2])

  pos = [[random.random() for i in range(N)]for j in range(M)]
 
  nposX = [0.0 for i in range(N)]
  nposY = [0.0 for i in range(N)]
  nposZ = [0.0 for i in range(N)]
 
  dmat   = [[0.0 for i in range(N)]for j in range(N)]
  dmat0  = [[0.0 for i in range(N)]for j in range(N)]
  for i in range (N):
    for j in range (N):
      dd = 0.0
      for k in range (M):
        dd += (pos[k][i]-pos[k][j]) *(pos[k][i]-pos[k][j])
      dmat[i][j] = dmat0[i][j] = math.sqrt(dd)
 
  metric_MDS(N,dmat,nposX,nposY,nposZ)
  print nposX 
  print nposY

  ndmat  = [[0.0 for i in range(N)]for j in range(N)]
  dmatfile = "dmat.comp"
  print "#write to '%s'"%(dmatfile)
  of = open(dmatfile,"w")
  for i in range (N):
    for j in range (i+1,N):
      dx = nposX[i] - nposX[j]
      dy = nposY[i] - nposY[j]
      ndmat[i][j] = math.sqrt(dx*dx + dy*dy) 
      of.write("%d %d %f %f\n"%(i,j,dmat0[i][j],ndmat[i][j]))
  of.close()

  xyfile = "xy.comp"
  print "#write to '%s'"%(xyfile)
  of = open(xyfile,"w")
  for i in range(N):
    of.write("%f %f %f %f\n"%(pos[0][i],pos[1][i],nposX[i],nposY[i]))
  of.close()

if __name__ == '__main__': _main()

