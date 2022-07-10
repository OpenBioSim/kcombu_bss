#
# <ligandbox_func.py>
#

import sys
import os 
import re
import math

LastModDate = "2019/11/14"

def library_2D_full_filename_from_ligand_id(lig_id,CONF):
## >> For temporary uploaded Query Molecule <<
## if lig_id = '12345.Qmol' --> /var/www/html/ligandbox_out/TMPOUT/12345.Qmol ##

  ## If lig_id contains '..' or '/' or ' or  ", the script stops.
  if (lig_id.find('..')>=0) or (lig_id.find('/')>=0) or (lig_id.find('\'')>=0) or (lig_id.find('\"')>=0):
      print "Content-type: text/html\n";
      print "<HTML><BODY><PRE>"
      print "#ERROR:Improper ligand_id('%s')"%(lig_id)
      print "</PRE></BODY></HTML>"
      sys.exit(1)

  if (lig_id.endswith('.Qmol')):
    filename = CONF['TMPOUT_DIR'] + '/' + lig_id
    return(filename)

#>> example of 'lig_id' <<
# '00022356-01'

  subdir = subdir_from_ligand_id(lig_id)

# modified by Sumita 2014/02/09
### lig_id = 'HTS1204-03603931' --> filename = '03600000/03603931.sdf' ##
#  if (lig_id.startswith('HTS')):
#    (head,tail) = lig_id.split('-')
#    filename = CONF.get('2D_LIBRARY_DIR','') + '/' + subdir + '/' +  tail
#  else:
#    filename = CONF.get('2D_LIBRARY_DIR','') + '/' + subdir + '/' + lig_id
  filename = CONF.get('2D_LIBRARY_DIR','') + '/' + subdir + '/' + lig_id
# modified by Sumita 2014/02/09

  if (filename.endswith('.sdf')==0):
    filename = filename + '.sdf'
  return(filename)


def subdir_from_ligand_id(lig_id):
#>> example of 'lig_id' <<
# '00000000-01'       --> '00000000'
# '00000001-01'       --> '00000000'
# '00004999-01'       --> '00000000'
# '00005000-01'       --> '00005000'
# '00009999-01'       --> '00005000'
# '00010000-01'       --> '00010000'
# '00022356-01'       --> '00020000'
# 'NS-00022175'       --> '00020000'
# 'NS-00022175.mol2'  --> '00020000'
# 'D08237.mol'        -->'D00005000'
# 'C08237.mol'        -->'C00005000'
# 'PDB_2MB-01'        -->'PDB_0_9'
# 'PDB_ALA-02'        -->'PDB_A_H'
# 'ZINC55249051.mol2' --> '55245000'
# 'HTS1204-03603931'  --> '03600000'

  #print "#subdir_from_ligand_id(%s):"%(lig_id)
  field  = re.split('[\-\.]',lig_id)
  subdir  = ''
  numid = ''
  for x in (field):
    #print x
    #ZINC55
    #01234
    if (x.startswith('ZINC')):
      x = x[4:]
    elif (re.match(r'^[A-Za-z]',x)!=''):
      x = x[1:]
      #print x

    if (len(x)>=4) and (x.isdigit()):
      numid = x

  if (numid!=''):
    num = int(numid)
    sta = int(5000*(math.floor(float(num)/5000.0)))
    stastr = '%d'%(sta)
    stastr = '0'*(8-len(stastr)) + stastr
    subdir = stastr
    if (lig_id.startswith('D')):
      subdir = 'D' + subdir
    if (lig_id.startswith('C')):
      subdir = 'C' + subdir

  elif (lig_id.startswith('PDB_')):
    field  = re.split('[\-\.]',lig_id)
    a = field[0][4:5]
    if ('0'<=a) and (a <='9'):
      subdir = 'PDB_0_9'
    elif ('A'<=a) and (a <='H'):
      subdir = 'PDB_A_H'
    elif ('I'<=a) and (a <='Q'):
      subdir = 'PDB_I_Q'
    else:
      subdir = 'PDB_R_Z'
    #print "lig_id '%s' '%s' --> '%s'"%(lig_id,a,subdir)
  #print "#num %d --> %s"%(num,subdir)
  if (subdir == ''):
    return('null')
  else:
    return(subdir)


