#!/usr/bin/env python
##
## <mkdes_mdir.py>
##
## for making atompair descriptor file for compound files stored in the multi directories 
##



import sys
import os

LastModDate = "July 7, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
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
OPT['A'] = 'F'
OPT['dir'] = ''

if (len(sys.argv)<3):
  print "mkdes_mdir.py [sdf_file]"
  print " for making atompair descriptor file for compound files stored in the multi directories"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -dir : target directory [%s]"%(OPT['dir']) 
  print " -A   : Action (T or F) [%s]"%(OPT['A']) 

  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['dir'])
print "#Ndir %d"%(len(dirlist))

for x in (dirlist):
  print x
  comstr = "dkcombu -M S -idl %s/%s -ode %s.des"%(OPT['dir'],x,x)
  print "#command '%s'"%(comstr)
  if (OPT['A']=='T'):
    system(comstr)
