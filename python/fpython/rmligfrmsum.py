#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Apr 23, 2014"

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



############
### MAIN ###
############

OPT = {}
OPT['isum']  = ''
OPT['osum']  = '-'
OPT['rmlig'] = 'D2S_4fik_1_E_1.pdb'

if (len(sys.argv)<2):
  print "rmligfrmsum.py <options>"
  print " for removing a ligand structure from summary file."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -isum   : input  summary file [%s]"%(OPT['isum'])
  print " -osum   : output summary file [%s]"%(OPT['osum'])
  print " -rmlig  : ligand name to be removed [%s]"%(OPT['rmlig'])
  sys.exit(1)


### [1] read option ###
read_option(sys.argv,OPT)

### [2] read isum (-L) ###
if (os.access(OPT['isum'],os.R_OK)==0):
  print "#ERROR:Can't open isum '%s'."%(OPT['isum'])
  sys.exit(1)

#>>> EXAMPLE OF SUMMARY FILE <<<
# #COMMAND '/home/takawaba/work/kcombu/src/fpython/sum3Dmodel.py -L CAH2.liglist -idtarp D_PRO_K -of summary/D_PRO_K.sum -A T'
# #DATE    '2013/12/28 15:07:05'
# #[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [NheavyB-NheavyA:5]
# #[rmsd:6] [volume_overlap:7] [tanimoto_volume:8]
# MFS_3hkn_1_C_1.pdb MFS_3hkn_1_C_1.pdb 47 47  0 0.478201 1685.100342 0.877452
# MFS_3hkn_1_C_1.pdb D2S_4fik_1_E_1.pdb 47 18 -29 10.682014 605.011963 0.312287
# BON_3mnu_1_E_1.pdb MFS_3hkn_1_C_1.pdb 20 47 27 3.133920 724.765503 0.387774
# TRU_1zgf_1_D_1.pdb TRU_1zgf_1_D_1.pdb 20 20  0 1.007337 655.142944 0.741655
# TRU_1zgf_1_D_1.pdb D2S_4fik_1_E_1.pdb 20 18 -2 9.038887 371.717712 0.326388
# D2S_4fik_1_E_1.pdb MFS_3hkn_1_C_1.pdb 18 47 29 9.400905 689.157227 0.372700
# D2S_4fik_1_E_1.pdb BON_3mnu_1_E_1.pdb 18 20  2 8.520524 352.766022 0.296828
# :

fi = open(OPT['isum'])
if (OPT['osum'] == '-'):
  fo = sys.stdout
else:
  fo = open(OPT['osum'],'w')

fo.write("#COMMAND_RMLIG %s"%(OPT['COMMAND']))
fo.write("#-rmlig  %s"%(OPT['rmlig']))
fo.write("#DATE_RMLIG  %s"%(OPT['START_DATE']))
remove_lines = []
for line in fi:
  if (line.startswith('#')):
    fo.write("%s"%(line))
  else:
    X = line.split()
    if (X[0]==OPT['rmlig']) or (X[1]==OPT['rmlig']):
      remove_lines.append(line)
    else: 
      fo.write("%s"%(line))

fo.write("#Nremove_line %d\n"%(len(remove_lines)))
for x in (remove_lines):
  fo.write("#DELETE %s"%(x))
fi.close()
if (OPT['osum'] != '-'):
  fo.close()
