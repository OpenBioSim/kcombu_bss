#!/usr/bin/env python
##
## <cp_moldata.py>
##


import sys
import os
import re
import math

LastModDate = "Aug 9, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def dirname_from_ligand_id(lig_id):
#>> example of 'lig_id' <<
# '00022356-01'      --> '00020001_00025000'
# 'NS-00022175'      --> '00020001_00025000'
# 'NS-00022175.mol2' --> '00020001_00025000'
# 'D08237.mol'       --> '00005000_00010000'

  #print "#dirname_from_ligand_id(%s):"%(lig_id)
  field  = re.split('[\-\.]',lig_id)

  id = ''
  for x in (field):
    #print x
    if (re.match(r'^[A-Za-z]',x)!=''):
      x = x[1:] 
      #print x
    if (len(x)>=4) and (x.isdigit()):
      id = x
  if (id==''):
    return('null')
  num = int(id)
  sta = int(5000*(math.ceil(float(num)/5000.0)-1))+1
  end = int(5000*math.ceil(float(num)/5000.0))
  stastr = '%d'%(sta)
  stastr = '0'*(8-len(stastr)) + stastr
  endstr = '%d'%(end)
  endstr = '0'*(8-len(endstr)) + endstr
  #print "num %8d stastr %s endstr %s"%(num,stastr,endstr)
  str = "%s_%s"%(stastr,endstr)
  return(str)



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
OPT['srcdir'] = 'src'
OPT['tardir'] = 'newsrc'
OPT['xsub']   = 'dkcombu:fkcombu:gprof:AtmPairDesc:transform:.o:.pyc:\.sdf|\.pdb|\.mol2|\.pdf|\.png'
OPT['dhead'] = ''
OPT['fhead'] = ''
OPT['ftail'] = ''


if (len(sys.argv)<3):
  print "cp_moldata.py <options>"
  print " for making open src packages for the 'kcombu' program, excluding improper files."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -srcdir : original    src directory  [%s]"%(OPT['srcdir']) 
  print " -tardir : target(new) src directory  [%s]"%(OPT['tardir']) 
  print " -dhead  : head string of target dir  [%s]"%(OPT['dhead'])
  print " -fhead  : head string of target file [%s]"%(OPT['fhead'])
  print " -ftail  : tail string of target file [%s]"%(OPT['ftail'])
  print " -A      : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)

srclist = os.listdir(OPT['srcdir'])
for srcfile in (srclist):
  subdir = dirname_from_ligand_id(srcfile)
  newdir = OPT['tardir'] + '/' + OPT['dhead'] + subdir
  mk_dir(newdir)
  (id,tail) = srcfile.split('.')
  if (OPT['ftail'] != ''):
    tail = OPT['ftail'] 
  newfile =  newdir + '/' + OPT['fhead'] + id + '.' + tail 
  copy_file(OPT['srcdir'] + '/'+srcfile,newfile)
