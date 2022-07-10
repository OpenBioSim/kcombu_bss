#!/usr/bin/env python

##
## <suplig_fkcombu.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import glob

LastModDate = '2019/06/02'


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


################
##### MAIN #####
################

OPT = {}
OPT['ipros'] = ''
OPT['odlig'] = 'MODEL_HighMann'
OPT['odpro_noAsn'] = 'MODEL_PDB_NO_ASN790'
OPT['ilig']  = 'HighMann_1gya_1.pdb'
OPT['A'] = 'F'
OPT['nlig'] = '10'
OPT['iam'] =  "ASN65_ASN790.am"
OPT['olog'] = '-'
OPT['nini'] = '10'
OPT['wepa'] = '0.0'
OPT['xroini'] = '180.0'
OPT['div'] = '0/1'
OPT['nligsta'] = ''
OPT['nligend'] = ''

if (len(sys.argv)<2):
  print "suplig_fkcombu.py <options>"
  print "  repeating superimposition of ligand on protein by 'fkcombu'."
  print "  LastModDate:%s"%(LastModDate)
  print " -ipros  : input protein  PDB file pattern [%s]"%(OPT['ipros']) 
  print " -ilig   : input ligand PDB file[%s]"%(OPT['ilig']) 
  print " -odlig  : output directory for ligands [%s]"%(OPT['odlig']) 
  print " -odpro_noAsn  : output directory for protein PDB file without focuesd Asns [%s]"%(OPT['odpro_noAsn']) 
  print " -nlig   : number of ligands to be generated [%s]"%(OPT['nlig'])
  print " -nligsta: start number of ligands to be generated [%s]"%(OPT['nligsta'])
  print " -nligend: end   number of ligands to be generated [%s]"%(OPT['nligend'])
  print " -iam    : input atom matching file [%s]"%(OPT['iam'])
  print " -nini   : num. of random initial conformatinon. (nini==1, only keep the input conf). [%s]"%(OPT['nini'])
  print " -xroini : initial maximum step for rotational bond (degree) [%s]"%(OPT['xroini'])
  print " -wepa   : Weight for Eprotatrct [%s]"%(OPT['wepa'])
  print " -A      : Action ('T' or 'F')[%s]"%(OPT['A'])
  print " -olog   : output filename [%s]"%(OPT['olog'])
  print " -div    : Job divisiton for 'nlig'.[bunshi]/[bunbo] [%s]"%(OPT['div'])
  sys.exit()

read_option(sys.argv, OPT)




if (OPT['olog'] == '-'):
  olog = sys.stdout
else:
  print "#write_to_log --> '%s'"%(OPT['olog'])
  olog = open(OPT['olog'],'w')

olog.write("#COMMAND %s\n"%(OPT['COMMAND']))
olog.write("#DATE    %s\n"%(OPT['START_DATE']))
for key in (OPT.keys()):
  olog.write("#KEY  %-12s  VALUE %s\n"%(key,OPT[key]))

filelist = glob.glob(OPT['ipros'])
ipro_filelist = sorted(filelist,lambda x,y:cmp(x,y))

Nligand = int(OPT['nlig'])
Nprotein = len(ipro_filelist)
if (OPT['nligsta'] != '') and (OPT['nligend'] != ''): 
  nligsta = int(OPT['nligsta'])
  nligend = int(OPT['nligend'])
elif (OPT['div'] != '0/1'): 
  [bunshi,bunbo] = OPT['div'].split('/')
  bunshi = int(bunshi)
  bunbo  = int(bunbo)
  nligsta = Nligand*bunshi/bunbo
  nligend = Nligand*(bunshi+1)/bunbo
  olog.write("#bunshi_bunbo   %d %d\n"%(bunshi,bunbo))
  print "#bunshi %d bunbo %d Nligand  %d Nsta %d Nend %d"%(bunshi,bunbo,Nligand,nligsta,nligend)
else:
  nligsta = 0
  nligend = Nligand

print "#Nligand  %d nligsta %d nligend %d"%(Nligand,nligsta,nligend)

olog.write("#Nprotein  %d\n"%(Nprotein))
olog.write("#Nligand   %d\n"%(Nligand))
olog.write("#nligsta   %d\n"%(nligsta))
olog.write("#nligend   %d\n"%(nligend))

for i in range(Nprotein):
  ipro_file = ipro_filelist[i]
  head      = ipro_file.replace('.pdb','')
  for j in range(nligsta,nligend):
#  for j in range(Nligand):
    olig_file =  OPT['odlig'] + '/' + head + '_%d.pdb'%(j+1)
    command = "fkcombu -bR F -alg F -weam 1 -wesc 1.0 -wepc 1.0 -wepa %s -srd %d -cfesc T"%(OPT['wepa'],j+1)
    command += " -iam %s -nini %s -xroini %s"%(OPT['iam'],OPT['nini'],OPT['xroini'])
    command += " -T %s"%(OPT['ilig'])
    command += " -R %s"%(ipro_file)
    command += " -P %s/%s"%(OPT['odpro_noAsn'],ipro_file)
    command += " -opdbT %s"%(olig_file)
    print "#%s"%(command)
    if (OPT['A'] == 'T'):
      print "#EXECUTE '%s'"%(command)
      olog.write("#EXECUTE %s\n"%(command))
      os.system(command)
    else:
      olog.write("#%s\n"%(command))

olog.write("#START_DATE    %s\n"%(OPT['START_DATE']))
now = datetime.now()
olog.write("#END_DATE      %s\n"%(now.strftime("%Y/%m/%d %H:%M:%S")))

if (OPT['olog'] != '-'):
  olog.close()
