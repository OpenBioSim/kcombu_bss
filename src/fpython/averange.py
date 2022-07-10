#!/usr/bin/env python
##
## <averange.py>
##

import os
import sys 
import math
from datetime import datetime

LastModDate = "July 17, 2013"

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



def Read_Data(input_file,cX,cY,X,Y,datalines):
  #print "#Read_Data('%s' cX %d cY %d)"%(input_file,cX,cY)
  if (os.access(input_file,os.R_OK)==0):
    print "#ERROR:Can't open '%s'"%(input_file)  
    sys.exit(1)
  f = open(input_file)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      field = line.split()
      if (cX<len(field)) and (cY<len(field)):
        x = float(field[cX]) 
        y = float(field[cY]) 
        X.append(x)
        Y.append(y)
        datalines.append(line)
    pass
  f.close()


################
##### MAIN #####
################

OPT = {}

OPT['cX'] = '1'
OPT['cY'] = '5'
OPT['of'] = ''
OPT['bin'] = '1.0'
OPT['of'] = '-'
OPT['minX'] = '0.0'
OPT['maxX'] = '1.0'
OPT['minY'] = '0.0'
OPT['maxY'] = '2.0'
OPT['ineqmaxX'] = 'E'


if (len(sys.argv)<2):
  print "averange.py [input_file] <options>"
  print " for calculating average correct Y value for focused X."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -cX   : column number for X (1,2,..) (range_defined_variable)[%s]"%(OPT['cX'])
  print " -cY   : column number for Y (1.2,..) (target_variable)       [%s]"%(OPT['cY'])
  print " -minX : min-value for focused range X (minX<=X) [%s]"%(OPT['minX'])
  print " -maxX : max-value for focused range X (maxX>=X) [%s]"%(OPT['maxX'])
  print " -minY : min-value for correct range Y (minY<=Y) [%s]"%(OPT['minY'])
  print " -maxY : max-value for correct range Y (maxY>=Y) [%s]"%(OPT['maxY'])
  print " -ineqmaxX: 'E':X<=maxX, 'I':X<maxX [%s]"%(OPT['ineqmaxX'])
  sys.exit(1)

read_option(sys.argv,OPT)


input_file  = sys.argv[1]
cX      = int(OPT['cX'])-1
cY      = int(OPT['cY'])-1
bin_width = float(OPT['bin'])
minX = float(OPT['minX'])
maxX = float(OPT['maxX'])
minY = float(OPT['minY'])
maxY = float(OPT['maxY'])

X    = []
Y    = []
datalines = []
Read_Data(input_file,cX,cY,X,Y,datalines)
N = len(X)

N_ix   = {}
Sy_ix  = {}
Syy_ix = {}


if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_to -->'%s'"%(OPT['of'])

of.write("#COMMAND averange.py %s -cX %s -cY %s -minX %s -maxX %s -minY %s -maxY %s\n"%(input_file,OPT['cX'],OPT['cY'],OPT['minX'],OPT['maxX'],OPT['minY'],OPT['maxY']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))

Nfocus         = 0
Nfocus_correct = 0
for i in range(N):
  if (OPT['ineqmaxX']=='E'):
   if (minX <=X[i]) and (X[i] <= maxX):
      Nfocus += 1 
      if (minY <=Y[i]) and (Y[i] <= maxY):
        Nfocus_correct += 1 

  if (OPT['ineqmaxX']=='I'):
   if (minX <=X[i]) and (X[i] < maxX):
      Nfocus += 1 
      if (minY <=Y[i]) and (Y[i] <= maxY):
        Nfocus_correct += 1 


if (Nfocus>0):
  Rcorrect = 100.0 * Nfocus_correct/Nfocus
else:
  Rcorrect = 0.0
of.write("#[filename] [Rcorrect(%%)] [Nfocus_correct] [Nfocus]\n")
of.write("%-16s %5.1f %6d %6d\n"%(input_file,Rcorrect,Nfocus_correct,Nfocus))

if (OPT['of']!='-'):
  of.close()
