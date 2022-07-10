#!/usr/bin/env python
##
## <hist2D.py>
##


import sys
import os
import re
import math

LastModDate = "June 19, 2012"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def isfloat(str):
  for i in range(len(str)):
    if (str[i].isdigit()):
      return(1)
  return(0)


def write_ix_iy_2D_data(ofname,DAT2d,ixlist,iylist,wx,wy):
  print "#write_ix_iy_2D_data()-->'%s'"%(ofname)
  of = open(ofname,"w")
  for ix in (ixlist):
    for iy in (iylist):
      val = 0
      if (DAT2d.has_key(ix)) and (DAT2d[ix].has_key(iy)):
        val = float(DAT2d[ix][iy])
      of.write("%f %f %f\n"%(wx*(ix+0.5),wy*(iy+0.5),val))
    of.write("\n")
  of.close()
  pass


def write_ix_iy_2D_tab_data(ofname,DAT2d,ixlist,iylist,wx,wy):
  print "#write_ix_iy_2D_data()-->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("\t")
  for iy in (iylist):
    of.write("%f\t"%(wy*(iy+0.5)))
  of.write("\n")
 
  for ix in (ixlist):
    of.write("%f\t"%(wx*(ix+0.5)))
    for iy in (iylist):
      val = 0.0
      if (DAT2d.has_key(ix)) and (DAT2d[ix].has_key(iy)):
        val = float(DAT2d[ix][iy])
      of.write("%f\t"%(val))
    of.write("\n")
  of.close()
  pass


def write_point_for_zero_2D_data(ofname,DAT2d,ixlist,iylist,wx,wy):
  print "#write_ix_iy_2D_data()-->'%s'"%(ofname)
  of = open(ofname,"w")
  for ix in (ixlist):
    for iy in (iylist):
      val = 0
      if (DAT2d.has_key(ix)) and (DAT2d[ix].has_key(iy)):
        val = float(DAT2d[ix][iy])
      if (val==0):
        of.write("%f %f %f\n"%(wx*(ix+0.5),wy*(iy+0.5),val))
        #of.write("%f %f 0\n"%(wx*(ix-1),wy*(iy-1)))
        #of.write("%f %f 0\n"%(wx*ix,wy*iy))
        #of.write("\n")
        #of.write("%f %f 0\n"%(wx*(ix-1),wy*iy))
        #of.write("%f %f 0\n"%(wx*ix,wy*(iy-1)))
        #of.write("\n")
        #of.write("\n")

  of.close()
  pass


###############
#### MAIN #####
###############

OPT = {}
OPT['if'] = ''
OPT['cx'] = '1'
OPT['cy'] = '2'
OPT['cz'] = '3'
OPT['wx'] = '10.0'
OPT['wy'] = '10.0'
OPT['tz'] = '2.0'
OPT['minx'] = '0.0'
OPT['maxx'] = '100.0'
OPT['miny'] = '0.0'
OPT['maxy'] = '100.0'

if (len(sys.argv)<2):
  print "hist2D.py [space-separated data file] <options>"
  print " for analyzing property data file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
#  print " -if  : input property file (space-splited)  [%s]"%(OPT['if']) 
  print " -cx : column number for X [%s]"%(OPT['cx']) 
  print " -cy : column number for Y [%s]"%(OPT['cy']) 
  print " -cz : column number for Z [%s]"%(OPT['cz']) 
  print " -minx : minimum of x [%s]"%(OPT['minx'])
  print " -maxx : maximum of x [%s]"%(OPT['maxx'])
  print " -miny : minimum of x [%s]"%(OPT['miny'])
  print " -maxy : maximum of x [%s]"%(OPT['maxy'])
  print " -wx : width for bin of X  [%s]"%(OPT['wx'])
  print " -wy : width for bin of X  [%s]"%(OPT['wy'])
  print " -tz : threshold of Z  [%s]"%(OPT['tz'])
  sys.exit(1)

read_option(sys.argv,OPT)
ifname = sys.argv[1]

cx = int(OPT['cx'])-1

cx = int(OPT['cx'])-1
cy = int(OPT['cy'])-1
cz = int(OPT['cz'])-1
wx = float(OPT['wx'])
wy = float(OPT['wy'])
tz = float(OPT['tz'])
minx = float(OPT['minx'])
maxx = float(OPT['maxx'])
miny = float(OPT['miny'])
maxy = float(OPT['maxy'])

if (os.access(ifname,os.R_OK)==0):
  print "#ERROR:Can't open '%s'"%(ifname)
  sys.exit(1)

f = open(ifname)

N2d   = {}
S2d   = {}
SS2d  = {}
POS2d = {}

### (1) Read file line one by one ###
ix_dic = {}
iy_dic = {}

for line in f:
  line = line.rstrip('\n')
  if (len(line)>10) and (line.startswith('#')==0):
    field = line.split()
    if (isfloat(field[cx])) and (isfloat(field[cy])) and (isfloat(field[cz])):
      x = float(field[cx])
      y = float(field[cy])
      z = float(field[cz])
      ix = int(x/wx)
      if (x==minx):
        ix = 0
      if (x==maxx):
        ix = int((maxx-0.5*wx)/wx) 
      iy = int(y/wy)
      if (y==miny):
        iy = 0
      if (y==maxy):
        iy = int((maxy-0.5*wy)/wy) 

      #print "x %f ix %d"%(x,ix)
      ix_dic[ix] = 1
      iy_dic[iy] = 1
      if (N2d.has_key(ix)==0):
        N2d[ix] = {}
        S2d[ix] = {}
        SS2d[ix] = {}
        POS2d[ix] = {}
      if (N2d[ix].has_key(iy)==0):
        N2d[ix][iy] = 0 
        S2d[ix][iy] = 0.0 
        SS2d[ix][iy] = 0.0 
        POS2d[ix][iy] = 0 
      N2d[ix][iy] += 1
      S2d[ix][iy] += z 
      SS2d[ix][iy] += z*z
      if (z<=tz):
        POS2d[ix][iy] += 1 
 
    #print "x %f %d y %f %d z %f"%(x,ix,y,iy,z)
f.close()

ixlist = sorted(ix_dic.keys(),lambda x,y:cmp(x,y))
for ix in (ixlist):
  print ix
iylist = sorted(iy_dic.keys(),lambda x,y:cmp(x,y))

write_ix_iy_2D_data("N.dat",N2d,ixlist,iylist,wx,wy)
write_ix_iy_2D_data("S.dat",S2d,ixlist,iylist,wx,wy)
for ix in (ixlist):
  for iy in (ixlist):
    if (N2d.has_key(ix)) and (N2d[ix].has_key(iy)) and (N2d[ix][iy]>0):
      S2d[ix][iy]   =  S2d[ix][iy]/N2d[ix][iy]
      POS2d[ix][iy] =  100.0*float(POS2d[ix][iy])/float(N2d[ix][iy])
    else:
      S2d[ix][iy]   =  -10.0 
      POS2d[ix][iy] =  0.0 

write_ix_iy_2D_data("ave.dat",S2d,ixlist,iylist,wx,wy)
write_ix_iy_2D_data("rpos.dat",POS2d,ixlist,iylist,wx,wy)
write_ix_iy_2D_tab_data("rpos_tab.txt",POS2d,ixlist,iylist,wx,wy)

write_point_for_zero_2D_data("zero.dat",N2d,ixlist,iylist,wx,wy)
