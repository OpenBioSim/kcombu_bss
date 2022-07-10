#!/usr/bin/env python
import sys
import os
import math

LastModDate = "Feb 25, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def read_table_file(ifname,dat):
  ## dat --> dat[Ndata][Nfield]:
  ## CAUTION: The last field of dat[Ndata] is the all the line string !!

  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    if (line.startswith('#')==0) and (len(line)>10):
      line = line.rstrip('\n')
      X = line.split()
      X.append(line)
      dat.append(X)
  pass

############
### MAIN ###
############

OPT = {}
OPT['x'] = 0
OPT['y'] = -1 
OPT['ymin'] = 0.0
OPT['ymax'] = 1.0
OPT['w'] = 0.1
OPT['oh'] = ''
OPT['O'] = 'F'
OPT['x']  = '1'
OPT['xo'] = '='
OPT['y']  = '-1'
OPT['ymin']  = '-1.0'
OPT['ymax']  = '+1.0'
OPT['i1']  = '-1'
OPT['i2']  = '-1'

if (len(sys.argv)<2):
  print "cmp2files.py [filenameA] [filenameB] <options>"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<optios>"
  print " -x    : focused column number (1,2,3,....) [%s]"%(OPT['x'])
  print " -xo   : operator for x '<' '>' '=' [%s]"%(OPT['xo'])
  print " -y    : ranged column number (1,2,3,....) (-1:don't care)[%s]"%(OPT['y'])
  print " -ymin : minimum value for y (>ymin)[%s]"%(OPT['ymin'])
  print " -ymax : maximum value for y (<=ymax)[%s]"%(OPT['ymax'])
  print " -O    : output line satisfid the ymin-ymax range ('T' or 'F') [%s]"%(OPT['O'])
  print " -i1   : identity column number 1 (1,2,3,....) (-1:don't care)[%s]"%(OPT['i1'])
  print " -i2   : identity column number 2 (1,2,3,....) (-1:don't care)[%s]"%(OPT['i2'])
  sys.exit(1)


read_option(sys.argv,OPT)
ifnameA = sys.argv[1]
ifnameB = sys.argv[2]
ymin = float(OPT['ymin'])
ymax = float(OPT['ymax'])
xcol = int(OPT['x'])-1
ycol = int(OPT['y'])-1
icol1 = int(OPT['i1'])-1
icol2 = int(OPT['i2'])-1

#print "#icol1 %d icol2 %d"%(icol1,icol2)
### [1] read table data ###

datA = []
datB = []
read_table_file(ifnameA,datA)
read_table_file(ifnameB,datB)



if (len(datA)!=len(datB)):
  if (OPT['i1']=='-1') and (OPT['i2']=='-1'):
    print "#ERROR:number of data is different!!"
    print "#datA('%s'):%d datB('%s'):%d"%(ifnameA,len(datA),ifnameB,len(datB))
    sys.exit(1)
  else:
    datAorig = datA
    datBorig = datB
    datA = []
    datB = []
    dicB = {}
    for j in range(len(datBorig)):
      key = datBorig[j][icol1] + ':' + datBorig[j][icol1]
      dicB[key] = j
 
    for i in range(len(datAorig)):
      key = datAorig[i][icol1] + ':' + datAorig[i][icol1]
      if (dicB.has_key(key)):
        #print "#append i %d j %d"%(i,hit_j)
        datA.append(datAorig[i])
        datB.append(datBorig[dicB[key]])
    pass


### [2] output data ###

for i in range(len(datA)):
  if (xcol<len(datA[i])) and (xcol<len(datB[i])):
    xA = datA[i][xcol]
    xB = datB[i][xcol]
    outX = 0
    if (OPT['xo']=='=') and (float(xA)==float(xB)):
      outX = 1
    if (OPT['xo']=='<') and (float(xA)< float(xB)):
      outX = 1
    if (OPT['xo']=='>') and (float(xA)> float(xB)):
      outX = 1

    outY = 1
    if (int(OPT['y'])>=1) and (ycol<len(datA[i])) and (ycol<len(datB[i])):
      yA = datA[i][ycol]
      yB = datB[i][ycol]
      if (float(yA)<=ymin) or (float(yA)>ymax):
        outY = 0 
      if (float(yB)<=ymin) or (float(yB)>ymax):
        outY = 0 

    if (outX==1) and (outY==1):
      print "%s vs %s :: %s"%(xA,xB,datA[i][-1])
      print "%s vs %s :: %s"%(xA,xB,datB[i][-1])
 
sys.exit(1)

colx = int(OPT['x']) - 1 
coly = int(OPT['y']) - 1 
width = float(OPT['w'])
ymin = float(OPT['ymin'])
ymax = float(OPT['ymax'])

if not os.access(ifname,os.R_OK):
  print "#ERROR:Can't open filename '%s'" %(ifname)
  sys.exit(1)

hist = {}
Nall = 0
f = open(ifname)
Sx = 0.0
SSx = 0.0
for line in f:
  line = line.rstrip('\n')
  if (line.startswith('#')==0):
    X = line.split()
    x = float(X[colx]) 
    if (int(OPT['y'])>=1):
      y = float(X[coly]) 
    if int(OPT['y']<=0) or ((ymin<y) and (y<=ymax)):
      intx = int(x/width)
      hist[intx] = hist.get(intx,0) + 1 
      Sx  += x 
      SSx += x*x
      Nall += 1
      if (OPT['O']=='T'):
        print "#%s"%(line)
f.close()

ilist = sorted(hist.keys(),lambda x,y:cmp(x,y))

aveX = Sx/float(Nall)
sdX  = math.sqrt(SSx/float(Nall)-aveX*aveX)

#print "#COMMAND %s"%(OPT['COMMAND'])
#print "#[filename] [Ndata] [average] [SD]"
#print "%s Ndata %d average %f SD %f"%(ifname,Nall,aveX,sdX)
print "%-33s Ndata %4d average %f"%(ifname,Nall,aveX)

#### OUTPUT HISTOGRAM ###
if (OPT['oh'] != ''):
  if (OPT['oh']=='-'):
    of = sys.stdout
  else:
    print "#Output_Histogram()-->'%s'"%(OPT['oh'])
    of = open(OPT['oh'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#Average %f\n"%(aveX))
  of.write("#SD     %f\n"%(sdX))
  for i in (ilist):
    of.write("%f %f %d\n"%(width*i,float(hist[i])/Nall,hist[i]))
  of.close()
