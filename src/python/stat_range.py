#!/usr/bin/env python
import sys
import os
import math

LastModDate = "Aug 27, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def read_file(ifname,colx,coly,colz,ymin,ymax,dat,width,histx,histz,histxz,avexonz):
  if not os.access(ifname,os.R_OK):
    print "#WARNING:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      X = line.split()
      x = float(X[colx]) 
      if (coly>=0):
        y = float(X[coly]) 
      #print "x '%s' %f width %f"%(X[colx],x,width)
      if ((coly<0) or ((ymin<y) and (y<=ymax))) and (X[colx] != 'Nan') and (X[colx] != 'nan'):
        intx = int(x/width)
        histx[intx] = histx.get(intx,0) + 1 
        dat['Sx']   += x 
        dat['SSx']  += x*x
        dat['Nall'] += 1
        dat['x'].append(x)
        dat['line'].append(line)
        dat['file'].append(ifname)
        if (colz>=0):
          z = float(X[colz]) 
          intz = int(z/width)
          histz[intz] = histz.get(intz,0) + 1 
          if (histxz.has_key(intx)==0):
            histxz[intx] = {}
          histxz[intx][intz] = histxz[intx].get(intz,0) + 1 
          avexonz[intz] = avexonz.get(intz,0.0) + x 
          dat['Sz']   += z 
          dat['SSz']  += z*z
          dat['Sxz']  += x*z
          dat['z'].append(z)
        #if (OPT['O']=='T'):
        #  print "#%s"%(line)
  f.close()

  return(1)


############
### MAIN ###
############

OPT = {}
OPT['x'] = 0
OPT['y'] = -1 
OPT['ymin'] = 0.0
OPT['ymax'] = 1.0
OPT['z'] = -1
OPT['w'] = 0.1
OPT['ozx'] = ''
OPT['ohx'] = ''
OPT['ohz'] = ''
OPT['ohzx'] = ''
OPT['oaxz'] = ''
OPT['oall'] = ''

if (len(sys.argv)<2):
  print "stat_range.py [file1]:[file2]:...:[fileN] <options>"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<optios>"
  print " -x   : column for 'x'(target) (1,2,...) [%d]"%(OPT['x'])
  print " -y   : column for 'y'(range restriction) (1,2,...) -1:accept everything[%d]"%(OPT['y'])
  print " -ymin: minimum value for y (>ymin)[%f]"%(OPT['ymin'])
  print " -ymax: maximum value for y (<=ymax)[%f]"%(OPT['ymax'])
  print " -z   : column for 'z'(correlated variables with 'x') (1,2,...) -1:ignore[%d]"%(OPT['z'])
  print " -w   : width for x-histogram [%f]"%(OPT['w'])
  print " -ozx : output data  z and x[%s]"%(OPT['ozx'])
  print " -oaxz: output average x on z [%s]"%(OPT['oaxz'])
  print " -ohx : output histogram file for x[%s]"%(OPT['ohx'])
  print " -ohz : output histogram file for z[%s]"%(OPT['ohz'])
  print " -ohzx: output 2D histogram file for x and z[%s]"%(OPT['ohzx'])
  print " -oall: output all the data [%s]"%(OPT['oall'])
  sys.exit(1)


read_option(sys.argv,OPT)
colx = int(OPT['x']) - 1 
coly = int(OPT['y']) - 1 
colz = int(OPT['z']) - 1 
width = float(OPT['w'])
ymin = float(OPT['ymin'])
ymax = float(OPT['ymax'])


histx   = {}
histz   = {}
histxz  = {}
avexonz = {}
dat = {}
dat['Nall'] = 0
dat['Sx']   = 0.0
dat['SSx']  = 0.0
dat['Sz']   = 0.0
dat['SSz']  = 0.0
dat['Sxz']  = 0.0
dat['x'] = []
dat['y'] = []
dat['z'] = []
dat['line'] = []
dat['file'] = []

list = sys.argv[1].split(':')
ifnamelist = []
for ifname in (list):
  if (ifname != ''): 
     ifnamelist.append(ifname)

for ifname in (ifnamelist):
  if (ifname != ''): 
    read_file(ifname,colx,coly,colz,ymin,ymax,dat,width,histx,histz,histxz,avexonz)

if (int(dat['Nall'])==0):
  sys.exit(1)

Nall = float(dat['Nall'])

aveX  = dat['Sx']/Nall
varX  = dat['SSx']/Nall-aveX*aveX
varXunbiased  = dat['SSx']/(Nall-1.0)- Nall*aveX*aveX/(Nall-1.0)

if (colz>=0):
  aveZ  = dat['Sz']/Nall
  varZ  = dat['SSz']/Nall - aveZ * aveZ
  varXZ = dat['Sxz']/Nall - aveX * aveZ
  ccXZ  = varXZ/math.sqrt(varX)/math.sqrt(varZ)

sys.stdout.write("Nfile %d %-40s Ndata %4d aveX %f (%.1f %%) varXunb %f x %s"%(len(ifnamelist),sys.argv[1],dat['Nall'],aveX,100*aveX,varXunbiased,OPT['x']))
if (colz>=0):
  sys.stdout.write(" z %s aveZ %f varZ %f ccXZ %f"%(OPT['z'],aveZ,varZ,ccXZ))
sys.stdout.write("\n")  


###############
#### OUTPUT ###
###############

if (OPT['ohx'] != ''):
  print "#Output_Histogram_for_x()-->'%s'"%(OPT['ohx'])
  of = open(OPT['ohx'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[x] [f(x)] [N(x)]\n")
  ilist = sorted(histx.keys(),lambda x,y:cmp(x,y))
  for i in (ilist): 
    of.write("%f %f %d\n"%(width*i,float(histx[i])/Nall,histx[i]))
  of.close()

if (OPT['ohz'] != ''):
  print "#Output_Histogram_for_z()-->'%s'"%(OPT['ohz'])
  of = open(OPT['ohz'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[z] [f(z)] [N(z)]\n")
  ilist = sorted(histz.keys(),lambda x,y:cmp(x,y))
  for i in (ilist): 
    of.write("%f %f %d\n"%(width*i,float(histz[i])/Nall,histz[i]))
  of.close()

if (OPT['ohzx'] != ''):
  print "#Output_2D_Histogram_for_x_and_z()-->'%s'"%(OPT['ohzx'])
  of = open(OPT['ohzx'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  ilistX = sorted(histx.keys(),lambda x,y:cmp(x,y))
  ilistZ = sorted(histz.keys(),lambda x,y:cmp(x,y))
  for i in (ilistX): 
    for j in (ilistZ): 
      of.write("%f %f %f %d\n"%(width*j,width*i,float(histxz[i].get(j,0))/Nall,histxz[i].get(j,0)))
    of.write("\n")
  of.close()


if (OPT['oaxz'] != ''):
  print "#Output_Averaged_x_on_z()-->'%s'"%(OPT['oaxz'])
  of = open(OPT['oaxz'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[z] [avex(z)]\n")
  for intz in (avexonz.keys()):
    avexonz[intz] = avexonz[intz]/histz[intz] 
  ilist = sorted(avexonz.keys(),lambda x,y:cmp(x,y))
  for i in (ilist): 
    of.write("%f %f\n"%(width*i,avexonz[i]))
  of.close()

if (OPT['ozx'] != ''):
  print "#Output_x_and_z()-->'%s'"%(OPT['ozx'])
  of = open(OPT['ozx'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[z] [x]\n")
  for i in range(len(dat['x'])):
    of.write("%f %f\n"%(dat['z'][i],dat['x'][i]))

if (OPT['oall'] != ''):
  print "#Output_all_the_data()-->'%s'"%(OPT['oall'])
  of = open(OPT['oall'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[line] [filename]\n")
  for i in range(len(dat['line'])):
    of.write("%s %s\n"%(dat['line'][i],dat['file'][i]))
