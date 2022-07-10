
import sys
import math
import os
import re


def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def subdir_from_ligand_id(lig_id):
#>> example of 'lig_id' <<
# '00022356-01'      --> '00020001_00025000'
# 'NS-00022175'      --> '00020001_00025000'
# 'NS-00022175.mol2' --> '00020001_00025000'
# 'D08237.mol'       --> 'D00005000_00010000'

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
  if (lig_id.startswith('D')):
    str = 'D'+str
  return(str)


