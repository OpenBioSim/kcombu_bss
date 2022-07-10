#!/usr/bin/env python
##
## <count_files.py>
##


import sys
import os

LastModDate = "Aug 31, 2011"

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
OPT['srcdir'] = 'src'
OPT['ol'] = ''

if (len(sys.argv)<3):
  print "count_files.py <options>"
  print " for counting files under '-srcdir'."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -srcdir : original    src directory [%s]"%(OPT['srcdir']) 
  print " -ol     : output file list file [%s]"%(OPT['ol']) 
  print " -A      : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['srcdir'])
#for s in (xsublist):
#  print s

filelist = []


for x in (dirlist):
  xfull = OPT['srcdir'] + '/' + x
  if (os.path.isdir(xfull)==0):
    sys.stdout.write("%s\n"%(xfull))
    filelist.append(x)
  else:
    sys.stdout.write(">%s"%(xfull))
    dirlist_under = os.listdir(xfull)
    Nfile_under_dir = 0
    for y in (dirlist_under):
      yfull = OPT['srcdir']+'/'+ x + '/' + y
      Nfile_under_dir += 1 
      #print y,yfull
      filelist.append(yfull)
    sys.stdout.write("(%d files)\n"%(Nfile_under_dir))

print "#Nfile %d"%(len(filelist))

if (OPT['ol'] != ''):
  print "#write_list_files() --> '%s'"%(OPT['ol'])
  of = open(OPT['ol'],'w') 
  for file in (filelist):
    of.write("%s\n"%(file))
  of.close() 
