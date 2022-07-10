#!/usr/bin/env python
##
## <self_success.py>
##

import os
import sys 
import math
from datetime import datetime

LastModDate = "May 11, 2014"

def read_option(argv,opt_dic):
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
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



OPT['crmsd'] = '6'
OPT['of'] = '-'
OPT['maxrmsd'] = '2.0'
OPT['alg'] = 'I'
OPT['oXY'] = ''
OPT['threX'] = '0.5'
OPT['cT'] = '1'
OPT['cR'] = '2'
OPT['self'] = 'T'

if (len(sys.argv)<2):
  print "self_success.py [summary_file] <options>"
  print " for calculating success rate for self-docking/cross-docking"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -crmsd   : column number for rmsd (target_variable)       [%s]"%(OPT['crmsd'])
  print " -maxrmsd : max-value for correct range Y (Y<=maxY) [%s]"%(OPT['maxrmsd'])
  print " -cT      : column number for target    ligand [%s]"%(OPT['cT'])
  print " -cR      : column number for reference ligand [%s]"%(OPT['cR'])
  print " -self    : 'T': only for target==reference, 'N':only for target!=reference(cross). 'A':accept all. [%s]"%(OPT['self'])
  sys.exit(1)

read_option(sys.argv,OPT)


input_summary_file  = sys.argv[1]
crmsd = int(OPT['crmsd'])
cT    = int(OPT['cT'])
cR    = int(OPT['cR'])
maxrmsd = float(OPT['maxrmsd'])


# #[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [NheavyB-NheavyA:5]
# #[rmsd:6] [volume_overlap:7] [tanimoto_volume:8]
# G20_2qwf_1_G_1.pdb G20_2qwf_1_G_1.pdb 24 24  0 0.271461 834.800354 0.822288
# G20_2qwf_1_G_1.pdb GNA_2qwe_1_G_1.pdb 24 23 -1 1.075027 773.756958 0.763037
# G20_2qwf_1_G_1.pdb DAN_2qwc_1_G_1.pdb 24 20 -4 1.050479 664.700317 0.661045
# G20_2qwf_1_G_1.pdb 49A_1f8e_1_H_1.pdb 24 20 -4 1.094503 655.997986 0.646492

if (os.access(input_summary_file,os.R_OK)==0):
  print "#ERROR:Can't open summary file '%s'."%(input_summary_file)
  sys.exit(1)

fi = open(input_summary_file)
Nself = 0
Nself_success = 0
for line in fi:
  if (line.startswith('#')==0) and (len(line)>10):
    line = line.rstrip('\n')
    X = line.split()
    if ((OPT['self']=='T') and (X[cT-1]==X[cR-1])) or ((OPT['self']=='N') and (X[cT-1]!=X[cR-1])) or (OPT['self']=='A'):
      Nself += 1
      rmsd = float(X[crmsd-1])
      if (rmsd <= maxrmsd):
        Nself_success += 1

success_rate = 100.0*float(Nself_success)/float(Nself)
print "#%s self %s all %d success %d success_rate %6.1f %% %5.0f %%"%(input_summary_file,OPT['self'],Nself,Nself_success,success_rate,success_rate)

