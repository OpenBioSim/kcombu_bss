#!/usr/bin/env python
##
## <avecurve.py>
##

import os
import sys 
import math
from datetime import datetime

LastModDate = "Aug 24, 2013"

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
  print "#Read_Data('%s' cX %d cY %d)"%(input_file,cX,cY)
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
OPT['cY'] = '2'
OPT['of'] = ''
OPT['bin'] = '0.1'
OPT['of'] = '-'
if (len(sys.argv)<2):
  print "avecurve.py [input_file] <options>"
  print " for drawing average value for Y curves against X."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -cX   : column number for X (1,2,..) (range_defined_variable)[%s]"%(OPT['cX'])
  print " -cY   : column number for Y (1.2,..) (target_variable)       [%s]"%(OPT['cY'])
  print " -bin  : bin width for variable X [%s]"%(OPT['bin'])

  sys.exit(1)

read_option(sys.argv,OPT)


input_file  = sys.argv[1]
cX      = int(OPT['cX'])-1
cY      = int(OPT['cY'])-1
bin_width = float(OPT['bin'])
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

of.write("#COMMAND avecurve.py %s -cX %s -cY %s -bin %s\n"%(input_file,OPT['cX'],OPT['cY'],OPT['bin']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))

for i in range(N):
  #print "%f %f %s"%(X[i],Y[i],datalines[i])
  ix = int(X[i]/bin_width)
  #print "%f %f ix %d"%(X[i],Y[i],ix)
  N_ix[ix]   = N_ix.get(ix,0) + 1
  Sy_ix[ix]  = Sy_ix.get(ix,0.0) + Y[i] 
  Syy_ix[ix] = Syy_ix.get(ix,0.0) + Y[i] * Y[i]

of.write("#[x:1] [ave_y:2] [ave_y - sd_y:3] [ave_y + sd_y:4] [sd_y:5] [Ncount:6]\n")
ix_list = sorted(N_ix.keys(), lambda x,y:cmp(x,y))
for ix in (ix_list):
  My  = Sy_ix[ix]/N_ix[ix]
  SDy = math.sqrt(Syy_ix[ix]/N_ix[ix] - My * My)
  of.write("%f %f %f %f %f %d\n"%(ix*bin_width,My,My-SDy,My+SDy,SDy,N_ix[ix]))

if (OPT['of']!='-'):
  of.close()
