#!/usr/bin/env python
##
## <close_val_row.py>
##


import sys
import os
import math
import random
from datetime import datetime

LastModDate = "Oct 15, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



###############
#### MAIN #####
###############

OPT = {}
OPT['col'] = 1
OPT['tar'] = '0.1' 

if (len(sys.argv)<3):
  print "close_val_row.py [input_table_file] <options>"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -col : focus column number [%s]"%(OPT['col']) 
  print " -tar : target value for the column '-col' [%s]"%(OPT['tar'])
  sys.exit(1)

read_option(sys.argv,OPT)

itabfile = sys.argv[1]

col_num = int(OPT['col']) - 1 
tar_val = float(OPT['tar'])

if (os.access(itabfile,os.R_OK)==0):
  print "#ERROR:can't open '%s'"%(itabfile)
  sys.exit(1)

f = open(itabfile)


# >> file format example <<
# #DATE   :2011/10/15 14:37:10
# #[thre:1] [recall(11/1*):2] [precision(11/*1):3] [fmeasure:4] [fpr(01/0*):5] [tpr(11/1*):6] [hitratio(*1/**):7]
# 0.018868 1.000000 0.000400 0.000800 0.999998 1.000000 0.999998
# 0.019240 1.000000 0.000400 0.000800 0.999996 1.000000 0.999996
# 0.020000 1.000000 0.000400 0.000800 0.999994 1.000000 0.999994
# 0.020419 1.000000 0.000400 0.000800 0.999990 1.000000 0.999990
# 0.020420 1.000000 0.000400 0.000800 0.999988 1.000000 0.999988

LINES  = []
VALUES = []
DIFFS  = []
min_diff = -1.0
i_min = 0
i = 0
for line in f:
  if (line.startswith('#')==0) and (len(line)>10): 
    line = line.rstrip('\n')
    field = line.split()
    LINES.append(line)
    VALUES.append(field)
    val = float(field[col_num])
    diff = math.sqrt((val - tar_val)*(val - tar_val))
    if (min_diff<0.0) or (diff < min_diff):
      i_min = i
      min_diff = diff  
    i += 1
f.close()

print "#%s min_diff %f"%(LINES[i_min],min_diff)

 

