#!/usr/bin/env python

##
## <chkChiralForIdSmi.py>
##

import sys
import os
import random
from datetime import datetime
import re
import math
import subprocess

LastModDate = "Sep 11, 2012"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        opt_dic[argv[i][1:]] = argv[i+1]

def subdir_from_ligand_id(lig_id):
#>> example of 'lig_id' <<
# '00000000-01'      --> '00000000'
# '00000001-01'      --> '00000000'
# '00004999-01'      --> '00000000'
# '00005000-01'      --> '00005000'
# '00009999-01'      --> '00005000'
# '00010000-01'      --> '00010000'
# '00022356-01'      --> '00020000'
# 'NS-00022175'      --> '00020000'
# 'NS-00022175.mol2' --> '00020000'
# 'D08237.mol'       -->'D00005000'
# 'ZINC55249051.mol2' --> '55245000'

  field  = re.split('[\-\.]',lig_id)

  id = ''
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
      id = x
  if (id==''):
    return('null')
  num = int(id)
  sta = int(5000*(math.floor(float(num)/5000.0)))
  stastr = '%d'%(sta)
  stastr = '0'*(8-len(stastr)) + stastr
  str = stastr
  if (lig_id.startswith('D')):
    str = 'D'+str
  #print "#num %d --> %s"%(num,str)
  return(str)



def read_annotation_file(ifname,dat):
#>> FILE FORMAT EXAMPLE <<
#>02454613
#PAU
#>02454944
#TMP T DT DRT 0DT
#>02454534
#DTD D1D

  print "#read_annotation_file('%s')"%(ifname)

  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING:can't open sdf file '%s'"%(ifname)
    return(0)
  f = open(ifname)
  id = ''
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('>')):
      id = line[1:]
    elif ((line.startswith('#')==0) and (len(line)>1) and (id != '')):
      dat[id] = line
  f.close()
  pass


##############
#### MAIN #####
###############

OPT = {}
OPT['ian'] = ''
OPT['idA'] = '.'
OPT['idB'] = '.'
OPT['A'] = 'F'
OPT['ste3A'] = 'F'
OPT['ste3B'] = 'F'
OPT['subdirA']  = 'F'
OPT['subdirB']  = 'F'
OPT['tailA'] = '.sdf'
OPT['tailB'] = '.sdf'
OPT['th'] = '1.0'
OPT['oan'] = 'idlist.out'
if (len(sys.argv)<3):
  print "chkChiralForIdSmi.py <options>"
  print " for checking chiralities for Identical-SMILES molecular pairs."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ian     : input ID_A-vs-equivalent ID_B annotation file [%s]"%(OPT['ian'])
  print " -idA     : input directory for molA [%s]"%(OPT['idA'])
  print " -idB     : input directory for molB [%s]"%(OPT['idB'])
  print " -subdirA : use output subdir for molA (T or F)[%s]"%(OPT['subdirA'])
  print " -subdirB : use output subdir for molB (T or F)[%s]"%(OPT['subdirB'])
  print " -ste3A   : Setup Stereo Parity from 3D for mol A (T or F)[%s]"%(OPT['ste3A'])
  print " -ste3B   : Setup Stereo Parity from 3D for mol B (T or F)[%s]"%(OPT['ste3B'])
  print " -tailA   : tail of molA ('null' is for none)[%s]"%(OPT['tailA'])
  print " -tailB   : tail of molB ('null' is for none)[%s]"%(OPT['tailB'])
  print " -th      : threshold value (max Tanimoto) [%s]"%(OPT['th'])
  print " -A       : Action (T or F)[%s]"%(OPT['A'])
  print " -oan     : output ID_A-vs-equivalent ID_B annotation file [%s]"%(OPT['oan'])
  sys.exit(1)

read_option(sys.argv,OPT)

if (OPT['tailA']=='null'):
  OPT['tailA']=''
if (OPT['tailB']=='null'):
  OPT['tailB']=''


### READ ANNOTATION FILE ###

ANNOT = {}
read_annotation_file(OPT['ian'],ANNOT)


if (OPT['A']=='T') and (OPT['oan'] != ''):
  of = open(OPT['oan'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))

for idA in (ANNOT.keys()):
  print ">%s"%(idA)
  print "%s"%(ANNOT[idA])
  if (OPT['subdirA']=='T'):
    imolfileA = "%s/%s/%s%s"%(OPT['idA'],subdir_from_ligand_id(idA),idA,OPT['tailA'])
  else:
    imolfileA = "%s/%s%s"%(OPT['idA'],idA,OPT['tailA'])
  idBs = ANNOT[idA].split(' ')
  equiv_idBs = [] 
  for idB in (idBs):
    if (OPT['subdirB']=='T'):
      imolfileB = "%s/%s/%s%s"%(OPT['idB'],subdir_from_ligand_id(idB),idB,OPT['tailB'])
    else:
      imolfileB = "%s/%s%s"%(OPT['idB'],idB,OPT['tailB'])

    command = "pkcombu -A %s -ste3A %s -B %s -ste3B %s -chste T -KV T"%(imolfileA,OPT['ste3A'],imolfileB,OPT['ste3B'])

    if (OPT['A']!='T'):
      print "#'%s'"%(command) 
    elif (OPT['A']=='T'):
      #os.system(command)
      p = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE)
      TANIMOTO = '0.0'
      while True:
        line = p.stdout.readline()
        #print line
        if (line.startswith('TANIMOTO')):
          line = line.rstrip('\n')
          (key,value) = line.split()
          TANIMOTO = value
        if not line:
          break
      print "TANIMOTO %s"%(TANIMOTO)
      if (float(TANIMOTO) >= float(OPT['th'])):
        equiv_idBs.append(idB) 
  if (OPT['A']=='T') and (OPT['oan'] != '') and (len(equiv_idBs)>0):
    of.write(">%s\n"%(idA))
    for i in range(len(equiv_idBs)):
      if (i!=0):
        of.write(" ")
      of.write("%s"%(equiv_idBs[i]))
    of.write("\n")
    pass

if (OPT['A']=='T') and (OPT['oan'] != ''):
  print "#output_equiv_moleclue() --> '%s' done"%(OPT['oan'])
  of.close()
