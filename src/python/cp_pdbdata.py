#!/usr/bin/env python
##
## <cp_pdbdata.py>
##


import sys
import os
import re
import math

LastModDate = "Feb 12, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def mk_dir(newdir):
  if (os.path.exists(newdir)==0):
    print "#mkdir %s"%(newdir)
    if (OPT['A']=='T'):
      os.mkdir(newdir)


def copy_file(origfile,newfile):
  str = "cp %s %s"%(origfile,newfile)
  print "#%s"%(str)
  if (OPT['A']=='T'):
    os.system(str)


###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['srcdir'] = 'PDB'
OPT['tardir'] = '/DB/PDBv3'
OPT['tar2d']  = 'T'

if (len(sys.argv)<3):
  print "cp_pdbdata.py <options>"
  print " for making open src packages for the 'kcombu' program, excluding improper files."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -srcdir : original    src directory  [%s]"%(OPT['srcdir']) 
  print " -tardir : target(new) src directory  [%s]"%(OPT['tardir']) 
  print " -tar2d  : make 2-letter target subdir (T or F)[%s]"%(OPT['tar2d'])
  print " -A      : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)

srcdirlist = os.listdir(OPT['srcdir'])

srclist          = []
srclist_fullpath = []
for s in (srcdirlist):
  if (os.path.isdir(OPT['srcdir']+'/'+s)):
    list = os.listdir(OPT['srcdir']+'/'+s)
    for p in (list):
      srclist.append(p)
      srclist_fullpath.append(OPT['srcdir']+'/'+s+'/'+p)
    pass
  else:
    srclist.append(s)
    srclist_fullpath.append(OPT['srcdir']+'/'+s)

#for i in range(len(srclist)):
#  print srclist[i],srclist_fullpath[i] 


for i in range(len(srclist)):
  srcfile = srclist[i]
  srcfile_fullpath = srclist_fullpath[i]
  ### srcfile is, for example "pdb1ubq.ent"
  field = srcfile.split('.')
  if (len(field)==2) and (field[1]=='ent'):
    id = field[0]
    if (OPT['tar2d']=='T'):
      newdir = OPT['tardir'] + '/' + id[4:6]
    else:
      newdir = OPT['tardir']
    mk_dir(newdir)
    newfile =  newdir + '/' + srcfile 
    copy_file(srcfile_fullpath,newfile)
