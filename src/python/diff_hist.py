#!/usr/bin/env python
import sys
import os
import math

LastModDate = "Nov 22, 2010"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

############
### MAIN ###
############

OPT = {}
OPT['x'] = 5
OPT['r'] = 4 
OPT['w'] = 1.0
OPT['oh'] = ''
OPT['O'] = 'F'

if (len(sys.argv)<2):
  print "diff_hist.py [filename] <options>"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<optios>"
  print " -x   : column for 'x'(target)    (1,2,...) [%d]"%(OPT['x'])
  print " -r   : column for 'r'(reference) (1,2,...) [%d]"%(OPT['y'])
  print " -w   : width for x-histogram [%f]"%(OPT['w'])
  print " -oh  : output histogram file [%s]"%(OPT['oh'])
  print " -O   : output line satisfid the ymin-ymax range ('T' or 'F') [%s]"%(OPT['O'])
  sys.exit(1)


read_option(sys.argv,OPT)
colx = int(OPT['x']) - 1 
colr = int(OPT['r']) - 1 
width = float(OPT['w'])

ifname = sys.argv[1]
if not os.access(ifname,os.R_OK):
  print "#ERROR:Can't open filename '%s'" %(ifname)
  sys.exit(1)

hist = {}
Nall = 0
f = open(ifname)
Sx = 0.0
SSx = 0.0
Sr = 0.0
SSr = 0.0

for line in f:
  line = line.rstrip('\n')
  if (line.startswith('#')==0):
    X = line.split()
    x = float(X[colx]) 
    r = float(X[colr]) 
    d = x - r
    #print "x %f r %f d %f"%(x,r,d)
    intd = int(d/width)
    hist[intd] = hist.get(intd,0) + 1 
    Sx  += x 
    SSx += x*x
    Sr  += r 
    SSr += r*r
    Nall += 1
    if (OPT['O']=='T'):
      print "#%s"%(line)
f.close()

ilist = sorted(hist.keys(),lambda x,y:cmp(x,y))

aveX = Sx/float(Nall)
sdX  = math.sqrt(SSx/float(Nall)-aveX*aveX)
aveR = Sr/float(Nall)
sdR  = math.sqrt(SSr/float(Nall)-aveR*aveR)

print "#COMMAND %s"%(OPT['COMMAND'])
print "#[filename] [Ndata] [aveX] [sdX] [aveR] [sdR]"
print "%s %d %f %f %f %f"%(ifname,Nall,aveX,sdX,aveR,sdR)

#### OUTPUT HISTOGRAM ###
if (OPT['oh'] != ''):
  if (OPT['oh']=='-'):
    of = sys.stdout
  else:
    print "#Output_Histogram()-->'%s'"%(OPT['oh'])
    of = open(OPT['oh'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#aveX %f sdX %f\n"%(aveX,sdX))
  of.write("#aveR %f sdR %f\n"%(aveR,sdR))
  for i in (ilist):
    of.write("%f %f %d\n"%(width*i,float(hist[i])/Nall,hist[i]))
  of.close()
