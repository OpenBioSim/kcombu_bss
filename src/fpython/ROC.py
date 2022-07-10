#!/usr/bin/env python
##
## <ROC.py>
##

import os
import sys 
import math
from datetime import datetime

LastModDate = "Apr 29, 2014"

def read_option(argv,opt_dic):
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



def read_tabular_file(input_file,colX,X,datalines,avoid_self_dock='F',colTar='1',colRef='2'):
  print "#read_tabular_file('%s' colX %d)"%(input_file,colX)
  if (os.access(input_file,os.R_OK)==0):
    print "#ERROR:Can't open tabular file '%s'"%(input_file)  
    sys.exit(1)

  colX0 = colX-1

  f = open(input_file)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      field = line.split()
      if (colX0<len(field)) and (colX0<len(field)):
        accept = 1
        if (avoid_self_dock=='T'):
          if (field[int(colTar)-1] == field[int(colRef)-1]):
            accept = 0 
        if (accept==1):   
          x = float(field[colX0]) 
          X.append(x)
          datalines.append(line)
    pass
  f.close()


def cal_statistics_XY(X,Y):
  sX  = 0.0
  sXX = 0.0
  sY  = 0.0
  sYY = 0.0
  sXY = 0.0
  N = len(X)
  for i in range(N):
    sX  += X[i]
    sXX += X[i]*X[i]
    sY  += Y[i]
    sYY += Y[i]*Y[i]
    sXY += X[i]*Y[i]
  avX = sX/N
  avY = sY/N
  sdX = math.sqrt(sXX/N - avX * avX)
  sdY = math.sqrt(sYY/N - avY * avY)
  varXY = sXY/N - avX * avY
  cc  = varXY/sdX/sdY
  return((avX,sdX,avY,sdY,cc))



def cal_measures_from_2x2_table(Nxy):
  # Nxy[1][1]:true  positive(TP)
  # Nxy[0][0]:true  negative(TN)
  # Nxy[1][0]:false positive(FP)
  # Nxy[0][1]:false negative(FN)

  N   = float(Nxy[0][0]+Nxy[0][1]+Nxy[1][0]+Nxy[1][1])
  Nx1 = float(Nxy[1][0]+Nxy[1][1])
  Ny1 = float(Nxy[0][1]+Nxy[1][1])
  Ny0 = float(Nxy[1][0]+Nxy[0][0])

  recall    = 0.0
  precision = 0.0
  tpr = 0.0
  fpr = 0.0
  mcc = 0.0
  fmeasure = 0.0
  if (Ny1>0):
    recall    = float(Nxy[1][1])/float(Ny1)
  if (Nx1>0):
    precision = float(Nxy[1][1])/float(Nx1)

  tpr = recall
  if (Ny0>0): 
    fpr = Nxy[1][0]/float(Ny0)
  if ((recall+precision)>0.0):
     fmeasure = 2.0*recall*precision/(recall+precision)
  mcc_denom = (Nxy[1][1]+Nxy[1][0])*(Nxy[0][0]+Nxy[0][1])*(Nxy[1][1]+Nxy[0][1])*(Nxy[0][0]+Nxy[1][0])
  if (mcc_denom >0):
    mcc = float(Nxy[1][1]*Nxy[0][0]-Nxy[0][1]*Nxy[1][0])/math.sqrt(mcc_denom)

  return((recall,precision,fpr,tpr,fmeasure,mcc))


################
##### MAIN #####
################

OPT = {}


OPT['ifX'] = ''
OPT['ifY'] = ''

OPT['cX'] = '7'
OPT['cY'] = '6'
OPT['of'] = ''
OPT['of'] = '-'
OPT['minY'] = '0.0'
OPT['maxY'] = '2.0'
OPT['alg'] = 'I'
OPT['oXY'] = ''
OPT['threX'] = '0.5'
OPT['avself'] = 'F'
OPT['cT'] = '1'
OPT['cR'] = '2'

if (len(sys.argv)<2):
  print "ROC.py  <options>"
  print " for calculating ROC curve."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ifX   : input file for variable X (prediction) [%s]"%(OPT['ifX'])
  print " -ifY   : input file for variable Y (actual)     [%s]"%(OPT['ifY'])
  print " -cX    : column number for X (1,2,..) (range_defined_variable)[%s]"%(OPT['cX'])
  print " -cY    : column number for Y (1.2,..) (target_variable)       [%s]"%(OPT['cY'])
  print " -minY  : min-value for correct range Y (minY<=Y) [%s]"%(OPT['minY'])
  print " -maxY  : max-value for correct range Y (maxY>=Y) [%s]"%(OPT['maxY'])
  print " -threX : threshold value for upper and lower X [%s]"%(OPT['threX'])
  print " -alg   : Algorithm. 'N'aive, 'I'ncremental [%s]"%(OPT['alg'])
  print " -of    : output ROC curve file [%s]"%(OPT['of'])
  print " -oXY   : output XY file [%s]"%(OPT['oXY'])
  print " -avself: avoid self docking (T or F) [%s]"%(OPT['avself']) 
  print " -cT    : column number for target    ligand [%s]"%(OPT['cT'])
  print " -cR    : column number for reference ligand [%s]"%(OPT['cR'])
  sys.exit(1)

read_option(sys.argv,OPT)


cX      = int(OPT['cX'])
cY      = int(OPT['cY'])
minY  = float(OPT['minY'])
maxY  = float(OPT['maxY'])
threX = float(OPT['threX'])


### (1) read data X,Y ###
X    = []
Y    = []
datalinesX = []
datalinesY = []
if (OPT['ifX']==''):
  print "#ERROR: the option '-ifX' is required."
  sys.exit(1)
else:
  read_tabular_file(OPT['ifX'],cX,X,datalinesX,avoid_self_dock=OPT['avself'],colTar=OPT['cT'],colRef=OPT['cR'])

if (OPT['ifY']==''):
  print "#ERROR: the option '-ifY' is required."
  sys.exit(1)
else:
  read_tabular_file(OPT['ifY'],cY,Y,datalinesY,avoid_self_dock=OPT['avself'],colTar=OPT['cT'],colRef=OPT['cR'])
N  = len(X)
Nx = len(X)
Ny = len(Y)
print "#Nx %d Ny %d"%(Nx,Ny)
if (Nx != Ny):
  print "#Nx is not equal to Ny."
  sys.exit(1)

#### (2) Simple Statistics ####
(avX,sdX,avY,sdY,cc) = cal_statistics_XY(X,Y)

Ny0 = 0
Ny1 = 0

NuppX = 0
NlowX = 0
Ny1_uppX = 0
Ny1_lowX = 0

for i in range(Nx):
  if (X[i]>=threX):
    upplowX = 'U'
  else:
    upplowX = 'L'
 
  if (minY<=Y[i]) and (Y[i]<=maxY):
    corY = 1
  else:
    corY = 0

  if (corY==1):
    Ny1 += 1
  else:
    Ny0 += 1

  if (upplowX=='U'):
    NuppX += 1
    if (corY==1):
      Ny1_uppX += 1
  else:
    NlowX += 1
    if (corY==1):
      Ny1_lowX += 1

RcorrY      = 100.0*float(Ny1)/(Ny0+Ny1)
RcorrYuppX  = 100.0*float(Ny1_uppX)/NuppX
RcorrYlowX  = 100.0*float(Ny1_lowX)/NlowX
RuppX = 100.0*float(NuppX)/(NuppX+NlowX)
RlowX = 100.0*float(NlowX)/(NuppX+NlowX)


### (3) try each threshold x ###
if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_ROCcurve_to -->'%s'"%(OPT['of'])

of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))
of.write("#[thre_x:1] [fpr:2] [tpr:3] [recall:4] [precision:5] [fmeasure:6] [mcc:7]\n")
of.write("#[Nx0y0:8] [Nx0y1:9] [Nx1y0:10] [Nx1y1:11]\n")

max_mcc        = 0.0
thre_x_max_mcc = 0.0

max_fmeasure        = 0.0
thre_x_max_fmeasure = 0.0
auc = 0.0
fpr0 = -1.0 
tpr0 = -1.0 


## 'N'aive algorithm ##
if (OPT['alg']=='N'):
  ###  make Unique X values ###
  uniqX = {}
  for x in (X):
    uniqX[x] = uniqX.get(x,0) + 1 
  uniqXlist = sorted(uniqX.keys(), lambda x,y:cmp(x,y)) 

  for thre_x in (uniqXlist):
    Nxy = [[0,0],[0,0]]
    for i in range(N):
      x = X[i] 
      y = Y[i]
      intX = 0
      intY = 0
      if (minY<=y) and (y<=maxY):
        intY = 1
      if (x>=thre_x):
        intX = 1
      Nxy[intX][intY] += 1

    (recall,precision,fpr,tpr,fmeasure,mcc) = cal_measures_from_2x2_table(Nxy)
    if (fpr0>=0.0):
      auc += math.fabs(fpr-fpr0)*(tpr+tpr0)/2.0
    of.write("%f %f %f %f %f %f %f %4d %4d %4d %4d\n"%(thre_x,fpr,tpr,recall,precision,fmeasure,mcc,Nxy[0][0],Nxy[0][1],Nxy[1][0],Nxy[1][1]))
    fpr0 = fpr
    tpr0 = tpr
    if (math.fabs(mcc) > math.fabs(max_mcc)):
      max_mcc = mcc
      thre_x_max_mcc = thre_x
    if (fmeasure > max_fmeasure):
      max_fmeasure = fmeasure 
      thre_x_max_fmeasure = thre_x


## 'I'ncremental algorithm ##
if (OPT['alg']=='I'):
  index = []
  for i in range(Nx):
    index.append(i)
 
  sindex = sorted(index, lambda x,y:cmp(X[x],X[y])) 

  Nxy = [[0,0],[0,0]]
  Nxy[0][0] = 0
  Nxy[0][1] = 0 
  Nxy[1][0] = Ny0 
  Nxy[1][1] = Ny1 

  for i in range(Nx):
    si = sindex[i]
    thre_x = X[si]
    y = Y[si]
    intY = 0
    if (minY<=y) and (y<=maxY):
      intY = 1

    if (i==0) or (X[sindex[i]] != X[sindex[i-1]]):
      (recall,precision,fpr,tpr,fmeasure,mcc) = cal_measures_from_2x2_table(Nxy)
      if (i>0):
        auc += math.fabs(fpr-fpr0)*(tpr+tpr0)/2.0
      of.write("%f %f %f %f %f %f %f %4d %4d %4d %4d\n"%(thre_x,fpr,tpr,recall,precision,fmeasure,mcc,Nxy[0][0],Nxy[0][1],Nxy[1][0],Nxy[1][1]))
      fpr0 = fpr
      tpr0 = tpr
      if (math.fabs(mcc) > math.fabs(max_mcc)):
        max_mcc = mcc
        thre_x_max_mcc = thre_x
      if (fmeasure > max_fmeasure):
        max_fmeasure = fmeasure 
        thre_x_max_fmeasure = thre_x

    Nxy[0][intY] += 1
    Nxy[1][intY] -= 1
    


## (4) Final output ##
of.write("#avX %f sdX %f avY %f sdY %f cc %f\n"%(avX,sdX,avY,sdY,cc))
of.write("#Ratio_y1 %.2f %% Ny0  %d Ny1 %d\n"%(100.0*float(Ny1)/(Ny0+Ny1),Ny0,Ny1))
of.write("#auc_for_ROCcurve %f\n"%(auc))
of.write("#max_mcc      %f thre_x_mac_mcc      %f\n"%(max_mcc,thre_x_max_mcc))
of.write("#max_fmeasure %f thre_x_mac_fmeasure %f\n"%(max_fmeasure,thre_x_max_fmeasure))
of.write("#CORRELATION  cc  %f  auc %f max_mcc %f max_fmeasure %f\n"%(cc,auc,max_mcc,max_fmeasure))


sys.stdout.write("#avX %f sdX %f avY %f sdY %f cc %f\n"%(avX,sdX,avY,sdY,cc))
sys.stdout.write("#auc_for_ROCcurve %f\n"%(auc))
sys.stdout.write("#Ratio_y1 %.2f %% Ny0  %d Ny1 %d\n"%(100.0*float(Ny1)/(Ny0+Ny1),Ny0,Ny1))
sys.stdout.write("#max_mcc      %f thre_x_mac_mcc      %f\n"%(max_mcc,thre_x_max_mcc))
sys.stdout.write("#max_fmeasure %f thre_x_mac_fmeasure %f\n"%(max_fmeasure,thre_x_max_fmeasure))
sys.stdout.write("#UPPER_LOWERX   N %d NlowX %d (%d) NuppX %d (%d)\n"%(N,NlowX,Ny1_lowX,NuppX,Ny1_uppX))
sys.stdout.write("#ACCURACY %s avY  %.3f RcrctY %5.1f %% %d RcrctYlowX %5.1f %% %d RcrctYuppX %5.1f %% %d\n"%(OPT['ifY'],avY,RcorrY,Ny1,RcorrYlowX,Ny1_lowX,RcorrYuppX,Ny1_uppX))
sys.stdout.write("#CORRELATION %s %s  cc  %.3f  auc %.3f  max_mcc %.3f  max_f %.3f\n"%(OPT['ifY'],OPT['cX'],cc,auc,max_mcc,max_fmeasure))

if (OPT['of']!='-'):
  of.close()


if (OPT['oXY'] != ''):
  of = open(OPT['oXY'],'w')
  print "#write_XY_to -->'%s'"%(OPT['oXY'])
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#[x:1] [y:2]\n")
  for i in range(N):
    x = X[i] 
    y = Y[i]
    of.write("%f %f %s\n"%(x,y,datalinesY[i]))

sys.exit(1)

