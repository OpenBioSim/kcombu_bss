#!/usr/bin/env python

##
## <calRMSDtriple.py>
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





def cal_RMSD_by_pkcombu(imolfile,dirfilePDB,pdbfile,dat,Action):
  command = "pkcombu -A %s -B %s/%s -fB P -KV T -rkrm T"%(imolfile,dirfilePDB,pdbfile)
  p = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE)
  dat['NHEAVYATOM_A'] = 0
  dat['NHEAVYATOM_B'] = 0
  dat['TANIMOTO'] = 0.0
  dat['RMSD']     = -1.0
  if (Action=='T'):
    while True:
      line = p.stdout.readline()
      #print line
      if (line.startswith('#')==0) and (line.startswith('[')==0) and (line.startswith(' ')==0) and (line.find(' ')>=0):
        line = line.rstrip('\n')
        (key,value) = line.split()
        if (line.startswith('TANIMOTO')):
          dat['TANIMOTO'] = float(value)
        if (line.startswith('RMSD')):
          dat['RMSD'] = float(value)
        if (line.startswith('NHEAVYATOM_A')):
          dat['NHEAVYATOM_A'] = int(value)
        if (line.startswith('NHEAVYATOM_B')):
          dat['NHEAVYATOM_B'] = int(value)
      if not line:
        break
    #print "TANIMOTO %s RMSD %s"%(TANIMOTO,RMSD)
  else:
    print "# %s"%(command)


def RMSD_statistics(simlist,thre_tanimoto):
  N = 0
  aveRMSD = 0.0
  minRMSD = 0.0
  maxRMSD = 0.0
  for x in (simlist):
    if (x['TANIMOTO'] >= thre_tanimoto):
      if (N==0):
        minRMSD = x['RMSD']
        maxRMSD = x['RMSD']
      else:
        if (x['RMSD'] < minRMSD):
          minRMSD = x['RMSD']
        if (x['RMSD'] > maxRMSD):
          maxRMSD = x['RMSD']
      aveRMSD += x['RMSD']
      N += 1 
  if (N>0):
    aveRMSD = aveRMSD/N
  return([N,minRMSD,aveRMSD,maxRMSD])


##############
#### MAIN #####
###############

OPT = {}
OPT['itri'] = ''
OPT['ian']   = ''
OPT['idA']   = '/home/WEBDB/LIGANDBOX/201110/3D_MOL2_0525'
OPT['idB']   = '/home/WEBDB/ZINC/MOL2'
OPT['idPDB'] = '/home/WEBDB/LIGAND-EXPO/lddb/cc-inst/coord-aggregates/ipdb'

OPT['A'] = 'F'
OPT['ste3A'] = 'F'
OPT['ste3B'] = 'F'
OPT['subdirA']  = 'T'
OPT['subdirB']  = 'T'
OPT['tailA'] = '-01.mol2'
OPT['tailB'] = '.mol2'
OPT['th'] = '1.0'
OPT['of'] = 'rmsd.out'


if (len(sys.argv)<3):
  print "calRMSDtriple.py <options>"
  print " for calculating RMSDs for ([idA],[idB],[PDBlig])-triples."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -itri    : input (idA,idB,PDBlig)-triple files [%s]"%(OPT['itri'])
  print " -idA     : input directory for molA [%s]"%(OPT['idA'])
  print " -idB     : input directory for molB [%s]"%(OPT['idB'])
  print " -idPDB   : input directory for PDBlig [%s]"%(OPT['idB'])
  print " -subdirA : use output subdir for molA (T or F)[%s]"%(OPT['subdirA'])
  print " -subdirB : use output subdir for molB (T or F)[%s]"%(OPT['subdirB'])
  print " -tailA   : tail of molA ('null' is for none)[%s]"%(OPT['tailA'])
  print " -tailB   : tail of molB ('null' is for none)[%s]"%(OPT['tailB'])
  print " -th      : threshold value (max Tanimoto) [%s]"%(OPT['th'])
  print " -A       : Action (T or F)[%s]"%(OPT['A'])
  print " -of      : output file [%s]"%(OPT['of'])
  sys.exit(1)

read_option(sys.argv,OPT)

if (OPT['tailA']=='null'):
  OPT['tailA']=''
if (OPT['tailB']=='null'):
  OPT['tailB']=''


### READ TRIPLE FILE ###
if (os.access(OPT['itri'],os.R_OK)==0):
  print "#ERROR:Can't open triple file '%s'"%(OPT['itri'])
  sys.exit(1)

if (OPT['A']=='T'):
  of = open(OPT['of'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#[idA:1] [idB:2] [NheavyA:3] [NheavyB:4]\n")
  of.write("#[NstrpdbA:5] [minRMSD:6]  [aveRMSD:7]  [maxRMSD:8]\n")
  of.write("#[NstrpdbB:9] [minRMSD:10] [aveRMSD:11] [maxRMSD:12] [pdbligs]\n")

fi = open(OPT['itri'])
for line in fi:
  if (line.startswith('#')==0) and (len(line)>5):
    line = line.rstrip('\n')
    fields = line.split()
    idA = fields[0]
    idB = fields[1]
    pdbligs = []
    for i in range(2,len(fields)):
      pdbligs.append(fields[i])
    if (OPT['subdirA']=='T'):
      imolfileA = "%s/%s/%s%s"%(OPT['idA'],subdir_from_ligand_id(idA),idA,OPT['tailA'])
    else:
      imolfileA = "%s/%s%s"%(OPT['idA'],idA,OPT['tailA'])
    if (OPT['subdirB']=='T'):
      imolfileB = "%s/%s/%s%s"%(OPT['idB'],subdir_from_ligand_id(idB),idB,OPT['tailB'])
    else:
      imolfileB = "%s/%s%s"%(OPT['idB'],idB,OPT['tailB'])
    
    print idA,idB,pdbligs
    print imolfileA
    print imolfileB
    simAlist = []
    simBlist = []

    for idPDB in (pdbligs):
      dirfilePDB = "%s/%s/%s"%(OPT['idPDB'],idPDB[0],idPDB)
      print dirfilePDB 
      if (os.path.isdir(dirfilePDB)): 
        pdbfiles = os.listdir(dirfilePDB)
        for pdbfile in pdbfiles:
          datA = {}
          cal_RMSD_by_pkcombu(imolfileA,dirfilePDB,pdbfile,datA,OPT['A'])
          simAlist.append(datA)
          datB = {}
          cal_RMSD_by_pkcombu(imolfileB,dirfilePDB,pdbfile,datB,OPT['A'])
          simBlist.append(datB)

    if (OPT['A'] == 'T'):
      Nheavyatom_A = 0
      Nheavyatom_B = 0
      if (len(simAlist)>0):
        Nheavyatom_A = simAlist[0]['NHEAVYATOM_A']
      if (len(simBlist)>0):
        Nheavyatom_B = simBlist[0]['NHEAVYATOM_A']

      (N_A,minRMSD_A,aveRMSD_A,maxRMSD_A) = RMSD_statistics(simAlist,float(OPT['th']))
      (N_B,minRMSD_B,aveRMSD_B,maxRMSD_B) = RMSD_statistics(simBlist,float(OPT['th']))
      if (N_A>0) and (N_B>0):
        of.write("%-12s %-12s %2d %2d %4d %5.2f %5.2f %5.2f %4d %5.2f %5.2f %5.2f %s\n"%(idA,idB,Nheavyatom_A,Nheavyatom_B,N_A,minRMSD_A,aveRMSD_A,maxRMSD_A,N_B,minRMSD_B,aveRMSD_B,maxRMSD_B,pdbligs))
 
fi.close()

if (OPT['A']=='T'):
  print "#write_RMSDs --> '%s'"%(OPT['of'])
  of.close()

