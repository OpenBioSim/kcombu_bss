#!/usr/bin/env python
import sys
import os
import math

LastModDate = "July 8, 2010"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def write_VRML_Sphere(of,x,y,z,radius,R,G,B,T):
  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(x,y,z))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry Sphere{\n")
  of.write("        radius %f\n"%(radius))
  of.write("      }")
  of.write("      appearance Appearance{\n")
  of.write("        material Material{\n")
  of.write("          diffuseColor %f %f %f\n"%(R,G,B))
  of.write("          transparency %s\n"%(T))
  of.write("        }\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")

def write_VRML_Wire_Octahedron(of,tx,ty,tz,radius,R,G,B,T):
  pnt = []
  pnt.append([ 0.0,  0.0,  1.0]) #0
  pnt.append([ 1.0,  0.0,  0.0]) #1
  pnt.append([ 0.0,  1.0,  0.0]) #2
  pnt.append([-1.0,  0.0,  0.0]) #3
  pnt.append([ 0.0, -1.0,  0.0]) #4
  pnt.append([ 0.0,  0.0, -1.0]) #5

  tri = []
  tri.append([0,2,1]) #0
  tri.append([0,3,2]) #1
  tri.append([0,4,3]) #2
  tri.append([0,1,4]) #3
  tri.append([1,2,5]) #4
  tri.append([2,3,5]) #5
  tri.append([3,4,5]) #6
  tri.append([4,1,5]) #7

  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(tx,ty,tz))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry IndexedLineSet{\n")
  of.write("        coord Coordinate{\n")
  of.write("          point[\n")
  for i in range(len(pnt)):
    of.write("            %.2f %.2f %.2f"%(radius*pnt[i][0],radius*pnt[i][1],radius*pnt[i][2]))
    if (i<(len(pnt)-1)):
      of.write(",")
    of.write("\n")
  of.write("          ]\n")
  of.write("        }\n")
  of.write("        coordIndex[\n")
  for i in range(len(tri)):
    of.write("           %d %d %d %d -1"%(tri[i][0],tri[i][1],tri[i][2],tri[i][0]))
    if (i<(len(tri)-1)):
      of.write(",")
    of.write("\n")
  of.write("        ]\n")
  of.write("        colorPerVertex FALSE\n")
  of.write("        color Color{\n")
  of.write("          color[%f %f %f]\n"%(R,G,B))
  of.write("        }\n")
  of.write("        colorIndex[\n")
  of.write("          ")
  for i in range(len(tri)):
    of.write(" 0")
  of.write("\n")
  of.write("        ]\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")



def write_VRML_Wire_Icosahedron(of,tx,ty,tz,radius,R,G,B,Nrecur):
  #print "orig a %f t %f l %f"%(a,t,l) 
  t = 2*math.sqrt(5)/5.0
  l = t
  pi = math.pi 

  pnt = []
  pnt.append([ 0.0,  0.0,  1.0]) #0
  for i in range(5):
    th = i * 2*pi/5.0
    x = t * math.cos(th) 
    y = t * math.sin(th)
    z = l/2.0 
    pnt.append([ x, y, z]) 
  for i in range(5):
    th = i * 2*pi/5.0 + pi/5.0
    x = t * math.cos(th) 
    y = t * math.sin(th)
    z = -l/2.0 
    pnt.append([ x, y, z]) 
  
  pnt.append([ 0.0,  0.0, -1.0]) 

  tri = []
  tri.append([0,2,1])   #0
  tri.append([0,3,2])   #1
  tri.append([0,4,3])   #2
  tri.append([0,5,4])   #3
  tri.append([0,1,5])   #4
  tri.append([1,2,6])   #5
  tri.append([2,3,7])   #6
  tri.append([3,4,8])   #7
  tri.append([4,5,9])   #8
  tri.append([5,1,10])  #9
  tri.append([2,7,6])   #10
  tri.append([3,8,7])   #11
  tri.append([4,9,8])   #12
  tri.append([5,10,9])  #13
  tri.append([1,6,10])  #14
  tri.append([6,7,11])  #15
  tri.append([7,8,11])  #16
  tri.append([8,9,11])  #17
  tri.append([9,10,11]) #18
  tri.append([10,6,11]) #19


  ## adding new points by recursive dividing triangles ##
  for r in range(Nrecur):
    trinew = []
    Npnt = len(pnt)
    Ntri = len(tri)
    np1 = [0.0,0.0,0.0]
    np2 = [0.0,0.0,0.0]
    np3 = [0.0,0.0,0.0]
    for t in range(Ntri):
      nlen1 = nlen2 = nlen3 = 0
      for k in range(3):
        np1[k] = (pnt[tri[t][0]][k] + pnt[tri[t][1]][k])/2.0
        np2[k] = (pnt[tri[t][1]][k] + pnt[tri[t][2]][k])/2.0
        np3[k] = (pnt[tri[t][2]][k] + pnt[tri[t][0]][k])/2.0
        nlen1 += np1[k]*np1[k]
        nlen2 += np2[k]*np2[k]
        nlen3 += np3[k]*np3[k]
      nlen1 = math.sqrt(nlen1)
      nlen2 = math.sqrt(nlen2)
      nlen3 = math.sqrt(nlen3)
      for k in range(3):
        np1[k] /= nlen1
        np2[k] /= nlen2
        np3[k] /= nlen3
      pnt.append([np1[0],np1[1],np1[2]])
      pnt.append([np2[0],np2[1],np2[2]])
      pnt.append([np3[0],np3[1],np3[2]])

      trinew.append([tri[t][0], Npnt,      Npnt+2])
      trinew.append([Npnt,      tri[t][1], Npnt+1])
      trinew.append([Npnt+1,    tri[t][2], Npnt+2])
      trinew.append([Npnt,Npnt+1,Npnt+2])
      Npnt += 3
    tri = trinew


  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(tx,ty,tz))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry IndexedLineSet{\n")
  of.write("        coord Coordinate{\n")
  of.write("          point[\n")
  for i in range(len(pnt)):
    of.write("            %.2f %.2f %.2f"%(radius*pnt[i][0],radius*pnt[i][1],radius*pnt[i][2]))
    if (i<(len(pnt)-1)):
      of.write(",")
    of.write("\n")
  of.write("          ]\n")
  of.write("        }\n")
  of.write("        coordIndex[\n")
  for i in range(len(tri)):
    of.write("           %d %d %d %d -1"%(tri[i][0],tri[i][1],tri[i][2],tri[i][0]))
    if (i<(len(tri)-1)):
      of.write(",")
    of.write("\n")
  of.write("        ]\n")
  of.write("        colorPerVertex FALSE\n")
  of.write("        color Color{\n")
  of.write("          color[%f %f %f]\n"%(R,G,B))
  of.write("        }\n")
  of.write("        colorIndex[\n")
  of.write("          ")
  for i in range(len(tri)):
    of.write(" 0")
  of.write("\n")
  of.write("        ]\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")







def write_VRML_Cylinder(of,x,y,z,nx,ny,nz,radius,height,R,G,B,T):
  wx  = nz 
  wy  = 0.0
  wz  = -nx
#  print "ny %f"%(ny)
  th  = math.acos(ny) 

  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(x,y,z))
  of.write("  rotation %f %f %f %f\n"%(wx,wy,wz,th))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry Cylinder{\n")
  of.write("        radius %f\n"%(radius))
  of.write("        height %f\n"%(height))
  of.write("      }")
  of.write("      appearance Appearance{\n")
  of.write("        material Material{\n")
  of.write("          diffuseColor %f %f %f\n"%(R,G,B))
  of.write("          transparency %s\n"%(T))
  of.write("        }\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")



def write_VRML_Wire_Circle(of,x,y,z,nx,ny,nz,radius,height,R,G,B,Ndiv):
  wx  = nz 
  wy  = 0.0
  wz  = -nx
#  print "ny %f"%(ny)
  th  = math.acos(ny) 

  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(x,y,z))
  of.write("  rotation %f %f %f %f\n"%(wx,wy,wz,th))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry IndexedLineSet{\n")
  of.write("        coord Coordinate{\n")
  of.write("          point[\n")
  of.write("            %.2f %.2f %.2f,\n"%(0.0,0.0,0.0))
  for i in range(Ndiv):
    phi = 2*3.14159*float(i)/float(Ndiv)
    xd = radius * math.cos(phi)
    yd = 0.0
    zd = radius * math.sin(phi)
    of.write("            %.2f %.2f %.2f"%(xd,yd,zd))
    if (i<Ndiv):
      of.write(",")
    of.write("\n")
  of.write("          ]\n")
  of.write("        }\n")
  of.write("        coordIndex[\n")
  of.write("          ")
  for i in range(Ndiv):
    of.write(" %d"%(i+1))
  of.write(" 1 -1,\n")
  for i in range(Ndiv):
    of.write("           0 %d -1"%(i+1))
    if (i<Ndiv):
      of.write(",")
    of.write("\n")
  of.write("        ]\n")
  of.write("        colorPerVertex FALSE\n")
  of.write("        color Color{\n")
  of.write("          color[%f %f %f]\n"%(R,G,B))
  of.write("        }\n")
  of.write("        colorIndex[\n")
  of.write("          ")
  for i in range(Ndiv+2):
    of.write(" 0")
  of.write("\n")
  of.write("        ]\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")


def write_VRML_Wire_Disc(of,x,y,z,nx,ny,nz,radius,height,R,G,B,Ndiv):
  wx  = nz 
  wy  = 0.0
  wz  = -nx
#  print "ny %f"%(ny)
  th  = math.acos(ny) 

  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(x,y,z))
  of.write("  rotation %f %f %f %f\n"%(wx,wy,wz,th))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry IndexedLineSet{\n")


  of.write("        coord Coordinate{\n")
  of.write("          point[\n")
  ## points on upper ring ##  
  yd = 0.5*height 
  of.write("            %.2f %.2f %.2f,\n"%(0.0,yd,0.0))
  for i in range(Ndiv):
    phi = 2*3.14159*float(i)/float(Ndiv)
    xd = radius * math.cos(phi)
    yd = 0.5*height 
    zd = radius * math.sin(phi)
    of.write("            %.2f %.2f %.2f,\n"%(xd,yd,zd))

  ## points on lower ring ##  
  yd = -0.5*height 
  of.write("            %.2f %.2f %.2f,\n"%(0.0,yd,0.0))
  for i in range(Ndiv):
    phi = 2*3.14159*float(i)/float(Ndiv)
    xd = radius * math.cos(phi)
    yd = -0.5*height 
    zd = radius * math.sin(phi)
    of.write("            %.2f %.2f %.2f"%(xd,yd,zd))
    if (i<Ndiv):
      of.write(",")
    of.write(" \n")
  of.write("          ]\n")
  of.write("        }\n")


  of.write("        coordIndex[\n")
  of.write("          ")
  ## lines on upper ring ##  
  for i in range(Ndiv):
    of.write(" %d"%(i+1))
  of.write(" 1 -1,\n")

  for i in range(Ndiv):
    of.write("           0 %d -1,\n"%(i+1))

  ## lines on lower ring ##  
  of.write("          ")
  for i in range(Ndiv):
    of.write(" %d"%(i+2+Ndiv))
  of.write(" %d -1,\n"%(Ndiv+2))

  for i in range(Ndiv):
    of.write("           %d %d -1,\n"%(Ndiv+1,Ndiv+i+2))

  ## lines connecting upper and lower rings ##  
  for i in range(Ndiv):
    of.write("           %d %d -1"%(i+1,i+2+Ndiv))
    if (i<Ndiv):
      of.write(",")
    of.write("\n")
  
  of.write("        ]\n")


  of.write("        colorPerVertex FALSE\n")
  of.write("        color Color{\n")
  of.write("          color[%f %f %f]\n"%(R,G,B))
  of.write("        }\n")
  of.write("        colorIndex[\n")
  of.write("          ")
  for i in range(3*Ndiv+2):
    of.write(" 0")
  of.write("\n")
  of.write("        ]\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")








##############
#### MAIN ####
##############

OPT = {}
OPT['tmin'] =   0.0
OPT['tmax'] = 100.0
OPT['smin'] = 0.5
OPT['smax'] = 2.0
OPT['sc']   = 'T'
OPT['ow'] = 'out.wrl'
OPT['tr'] = 0.0
OPT['rg'] = 'd'
OPT['sp'] = 'i'
OPT['sdiv'] = 'T'

if (len(sys.argv)<2):
  print "ring2wrl.py [DABLring PDB file] <options>"
  print " for converting DABLring PDB file into VRML format(*.wrl)"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -rg   : Ring type.   'C'ylinder, 'c':wired circle, 'd':wired disc[%s]"%(OPT['rg'])
  print " -sp   : Sphere type. 'S'phere,   'o'ctahedron-wire 'i'cosahedron-wire[%s]"%(OPT['sp'])
  print " -sc   : scale type 'C'onstant,'T':t-Factor[%s]"%(OPT['sc'])
  print " -sdiv : scale-dependent Ndivision control (T or F)[%s]"%(OPT['sdiv'])
  print " -tmin : tFactor min[%f]"%(OPT['tmin'])
  print " -tmax : tFactor max[%f]"%(OPT['tmax'])
  print " -smin : scale   min[%f]"%(OPT['smin'])
  print " -smax : scale   max[%f]"%(OPT['smax'])
  print " -tr   : transparency [%f]"%(OPT['tr'])
  print " -ow   : output VRML file [%s]"%(OPT['ow'])
  sys.exit(1)

ipdbfile = sys.argv[1]
read_option(sys.argv,OPT)

OPT['smin'] = float(OPT['smin'])
OPT['smax'] = float(OPT['smax'])
OPT['tmin'] = float(OPT['tmin'])
OPT['tmax'] = float(OPT['tmax'])

R = {}
G = {}
B = {}

R['DON'] = 0.0
G['DON'] = 0.0
B['DON'] = 1.0

R['ACC'] = 1.0
G['ACC'] = 0.0
B['ACC'] = 0.0

R['BOT'] = 1.0
G['BOT'] = 0.0
B['BOT'] = 1.0

R['CAL'] = 0.5
G['CAL'] = 0.5
B['CAL'] = 0.5

R['RGC'] = 0.1
G['RGC'] = 0.5
B['RGC'] = 0.1

### <<FORMAT EXAMPLE>>
#REMARK                                                  Rg          nvec_x  nvec_y  nvec_z anum atom ele
#HETATM    1  CAL ATP A 355      16.223   8.151  -1.015 0.000  0.00   0.000   0.000   0.000  2939  PG  P L
#HETATM    2  ACC ATP A 355      17.550   8.536  -0.515 0.000  0.00   0.000   0.000   0.000  2940  O1G O A
#HETATM    3  ACC ATP A 355      15.748   6.835  -0.427 0.000  0.00   0.000   0.000   0.000  2941  O2G O A
#HETATM    4  ACC ATP A 355      16.255   8.014  -2.438 0.000  0.00   0.000   0.000   0.000  2942  O3G O A
#:
#HETATM   20  BOT ATP A 355       9.942   8.562   5.764 0.000  0.00   0.000   0.000   0.000  2958  O2' O B
#HETATM   21  CAL ATP A 355       8.988   9.957   4.085 0.000  0.00   0.000   0.000   0.000  2959  C1' C L
#HETATM   22  ACC ATP A 355       6.731   9.904   1.248 0.000  0.00   0.000   0.000   0.000  2962  N7  N A
#HETATM   23  DON ATP A 355       3.616  10.001   1.128 0.000  0.00   0.000   0.000   0.000  2965  N6  N D
#HETATM   24  ACC ATP A 355       3.756  10.034   3.445 0.000  0.00   0.000   0.000   0.000  2966  N1  N A
#HETATM   25  ACC ATP A 355       5.825  10.053   4.705 0.000  0.00   0.000   0.000   0.000  2968  N3  N A
#HETATM   26  RGC ATP E 355       6.923   9.917   2.422 1.160  0.00   0.005   1.000  -0.013
#HETATM   27  RGU ATP E 355       6.929  11.077   2.407 1.160  0.00
#HETATM   28  RGL ATP E 355       6.917   8.756   2.437 1.160  0.00
#HETATM   29  RGC ATP E 355       5.112  10.009   3.454 1.368  0.00   0.043   0.999   0.026
#HETATM   30  RGU ATP E 355       5.172  11.375   3.490 1.368  0.00
#HETATM   31  RGL ATP E 355       5.053   8.642   3.418 1.368  0.00
#HETATM   29  RGC ATP E 355       5.112  10.009   3.454 1.368  0.00   0.043   0.999   0.026
#01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
#          1         2         3         4         5         6         7         8         9        
#TER


of = open(OPT['ow'],'w')
print "#output_VRML_file() -->'%s'"%(OPT['ow'])
of.write("#VRML V2.0 utf8\n")
of.write("\n")
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
f = open(ipdbfile)
for line in f:
  line = line.rstrip('\n')
  atom = line[13:16]
  if (atom=="DON")or(atom=="ACC")or(atom=="BOT")or(atom=="CAL")or(atom=="RGC"):
    of.write("#%s %s\n"%(line,atom))
    x  = float(line[30:38])
    y  = float(line[38:46])
    z  = float(line[46:54])
    rg = float(line[54:60])
    if (rg==0.0):
      rg = 1.87
    tFactor = float(line[60:66])
    tnorm = (tFactor - OPT['tmin'])/(OPT['tmax']-OPT['tmin'])
    if (OPT['sc'] == 'C'):
      scale = 1.0
    elif (OPT['sc'] == 'T'):
      scale  = OPT['smin'] + tnorm *(OPT['smax']-OPT['smin'])

    if (atom=="RGC"):
      nx  = float(line[66:74])
      ny  = float(line[74:82])
      nz  = float(line[82:90])

    if (OPT['sdiv'] == 'T'):
      if (tnorm>0.8):
        Nrec = 4
      elif (tnorm>0.6):
        Nrec = 3
      elif (tnorm>0.4):
        Nrec = 2 
      elif (tnorm>0.2):
        Nrec = 1
      else:
        Nrec = 0
      Ndiv = int(60*tnorm+12) 
    else:
      Nrec = 2
      Ndiv = 30
 
    if (atom=="RGC"):
      if (OPT['rg']=='C'):
        write_VRML_Cylinder(of,x,y,z,nx,ny,nz,rg*scale,0.2*scale,R[atom],G[atom],B[atom],OPT['tr'])
      if (OPT['rg']=='c'):
        write_VRML_Wire_Circle(of,x,y,z,nx,ny,nz,rg*scale,0.2*scale,R[atom],G[atom],B[atom],Ndiv)
      if (OPT['rg']=='d'):
        write_VRML_Wire_Disc(of,x,y,z,nx,ny,nz,rg*scale,0.2*scale,R[atom],G[atom],B[atom],Ndiv)
    else:
      if (OPT['sp']=='S'):
        write_VRML_Sphere(of,x,y,z,0.5*scale,R[atom],G[atom],B[atom],OPT['tr'])
      if (OPT['sp']=='o'):
        write_VRML_Wire_Octahedron(of,x,y,z,0.5*scale,R[atom],G[atom],B[atom])
      if (OPT['sp']=='i'):
        write_VRML_Wire_Icosahedron(of,x,y,z,0.5*scale,R[atom],G[atom],B[atom],Nrec)

of.close()
