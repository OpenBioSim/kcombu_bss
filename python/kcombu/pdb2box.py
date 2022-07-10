#!/usr/bin/python
import sys
import os
import math

LastModDate = "July 16, 2010"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



##############
#### MAIN ####
##############

OPT = {}
OPT['ow'] = 'out.wrl'
OPT['op'] = ''

if (len(sys.argv)<2):
  print "pdb2box.py [PDB file] <options>"
  print " for making covering box (*.wrl) from PDB file"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ow   : output VRML file [%s]"%(OPT['ow'])
  print " -op   : output PDB  file [%s]"%(OPT['op'])
  sys.exit(1)

ipdbfile = sys.argv[1]
read_option(sys.argv,OPT)



### <<FORMAT EXAMPLE>>
#REMARK                                                  Rg          nvec_x  nvec_y  nvec_z anum atom ele
#HETATM    1  CAL ATP A 355      16.223   8.151  -1.015 0.000  0.00   0.000   0.000   0.000  2939  PG  P L
#HETATM    2  ACC ATP A 355      17.550   8.536  -0.515 0.000  0.00   0.000   0.000   0.000  2940  O1G O A
#HETATM    3  ACC ATP A 355      15.748   6.835  -0.427 0.000  0.00   0.000   0.000   0.000  2941  O2G O A
#HETATM    4  ACC ATP A 355      16.255   8.014  -2.438 0.000  0.00   0.000   0.000   0.000  2942  O3G O A
#01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
#          1         2         3         4         5         6         7         8         9        
#TER

f = open(ipdbfile)
init = 1
for line in f:
  line = line.rstrip('\n')
  atom = line[13:16]
  if (line.startswith("ATOM"))or(line.startswith("HETATM")):
    x  = float(line[30:38])
    y  = float(line[38:46])
    z  = float(line[46:54])
    if (init==1)or(x<xmin):
      xmin = x
    if (init==1)or(x>xmax):
      xmax = x
    if (init==1)or(y<ymin):
      ymin = y
    if (init==1)or(y>ymax):
      ymax = y
    if (init==1)or(z<zmin):
      zmin = z
    if (init==1)or(z>zmax):
      zmax = z
    init = 0
f.close()

xmid = (xmin+xmax)/2.0
ymid = (ymin+ymax)/2.0
zmid = (zmin+zmax)/2.0
xw   = xmax-xmin
yw   = ymax-ymin
zw   = zmax-zmin
print "xmin %f xmax %f xmid %f xw %f"%(xmin,xmax,xmid,xw)
print "ymin %f ymax %f ymid %f yw %f"%(ymin,ymax,ymid,yw)
print "zmin %f zmax %f zmid %f zw %f"%(zmin,zmax,zmid,zw)



xpnt = []
ypnt = []
zpnt = []

xpnt.append(xmax)
ypnt.append(ymax)
zpnt.append(zmax)

xpnt.append(xmax)
ypnt.append(ymax)
zpnt.append(zmin)

xpnt.append(xmax)
ypnt.append(ymin)
zpnt.append(zmax)

xpnt.append(xmax)
ypnt.append(ymin)
zpnt.append(zmin)

xpnt.append(xmin)
ypnt.append(ymax)
zpnt.append(zmax)

xpnt.append(xmin)
ypnt.append(ymax)
zpnt.append(zmin)

xpnt.append(xmin)
ypnt.append(ymin)
zpnt.append(zmax)

xpnt.append(xmin)
ypnt.append(ymin)
zpnt.append(zmin)

### OUTPUT VRML FILE ###
if (OPT['ow'] != ''):
  of = open(OPT['ow'],'w')
  print "#output_VRML_file() -->'%s'"%(OPT['ow'])
  of.write("#VRML V2.0 utf8\n")
  of.write("\n")
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("Transform{\n")
  of.write("  translation %f %f %f\n"%(xmid,ymid,zmid))
  of.write("  children[\n")
  of.write("    Shape{\n")
  of.write("      geometry Box{\n")
  of.write("        size %f %f %f\n"%(xw,yw,zw))
  of.write("      }")
  of.write("      appearance Appearance{\n")
  of.write("        material Material{\n")
  of.write("          diffuseColor 1.0 1.0 1.0\n")
  of.write("          transparency 1.0\n")
  of.write("        }\n")
  of.write("      }\n")
  of.write("    }\n")
  of.write("  ]\n")
  of.write("}\n")
  for i in range(8):
    of.write("Transform{\n")
    of.write("  translation %f %f %f\n"%(xpnt[i],ypnt[i],zpnt[i]))
    of.write("  children[\n")
    of.write("    Shape{\n")
    of.write("      geometry Sphere{\n")
    of.write("        radius 0.05\n")
    of.write("      }")
    of.write("      appearance Appearance{\n")
    of.write("        material Material{\n")
    of.write("          diffuseColor 0.1 0.1 0.1\n")
    of.write("          transparency 1.0\n")
    of.write("        }\n")
    of.write("      }\n")
    of.write("    }\n")
    of.write("  ]\n")
    of.write("}\n")
  
  of.close()

### OUTPUT PDB FILE ###
if (OPT['op'] != ''):
  of = open(OPT['op'],'w')
  print "#output_PDB_file() -->'%s'"%(OPT['op'])
  for i in range(8):
    of.write("HETATM%5d %4s %3s %s%5d   %8.3f%8.3f%8.3f\n"%(i,"XXX","BOX",' ',1,xpnt[i],ypnt[i],zpnt[i]))

  of.close()

