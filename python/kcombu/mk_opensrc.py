#!/usr/bin/env python
##
## <mk_opensrc.py>
##


import sys
import os

LastModDate = "Nov 11, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def copy_file(origfile,newfile,xsublist):
  xsubhit = 0
  for s in (xsublist):
    if (origfile.find(s)>=0):
      xsubhit = 1 
  if (xsubhit==1):
    print "#'%s' is excluded for xsublist %s "%(origfile,xsublist)

  if (xsubhit==0):
    str = "cp %s %s"%(origfile,newfile)
    print "#%s"%(str)
    if (OPT['A']=='T'):
      os.system(str)

def mk_dir(newdir):
  if (os.path.exists(newdir)==0):
    print "#mkdir %s"%(newdir)
    if (OPT['A']=='T'):
      os.mkdir(newdir)

###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['srcdir'] = 'src'
OPT['tardir'] = 'kcombu_opensrc'
OPT['xsub']   = 'dkcombu:fkcombu:gprof:AtmPairDesc:transform:stamp_transf:enemin_transf:.o:.pyc:\.sdf|\.pdb|\.mol2|\.pdf|\.png'

if (len(sys.argv)<3):
  print "mk_opensrc.py <options>"
  print " for making open src packages for the 'kcombu' program, excluding improper files."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -srcdir : original    src directory [%s]"%(OPT['srcdir']) 
  print " -tardir : target(new) src directory [%s]"%(OPT['tardir']) 
  print " -xsub   : substring for excluding files (str1:str2)[%s]"%(OPT['xsub']) 
  print " -A      : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['srcdir'])
xsublist = OPT['xsub'].split(':')
#for s in (xsublist):
#  print s

mk_dir(OPT['tardir']+'/'+OPT['srcdir'])



for x in (dirlist):
  xfull = OPT['srcdir']+'/'+x
  print '>',x,xfull
  if (os.path.isdir(xfull)==0):
    copy_file(xfull,OPT['tardir']+'/'+xfull,xsublist)
  else:
    mk_dir(OPT['tardir']+'/'+xfull)
    dirlist_under = os.listdir(xfull)
    for y in (dirlist_under):
      yfull = xfull + '/'+y
      if (os.path.isdir(yfull)==0):
        copy_file(yfull,OPT['tardir']+'/'+yfull,xsublist)
      else:
        mk_dir(OPT['tardir']+'/'+yfull)
        dirlist_under_under = os.listdir(yfull)
        for z in (dirlist_under_under):
          zfull = yfull + '/' + z
          print '>>>',z,zfull 
          if (os.path.isdir(zfull)==0):
            copy_file(zfull,OPT['tardir']+'/'+zfull,xsublist)
