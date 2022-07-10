#!/usr/bin/env python
##
## <self_success.py>
##

import os
import sys 
import math
from datetime import datetime

LastModDate = "May 23, 2014"

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



def read_tabular_file_one_column(input_file,colX,X):
  print "#read_tabular_file_one_column('%s' colX %s)"%(input_file,colX)
  if (os.access(input_file,os.R_OK)==0):
    print "#ERROR:Can't open tabular file '%s'"%(input_file)  
    sys.exit(1)

# #[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [NheavyB-NheavyA:5]
# #[rmsd:6] [volume_overlap:7] [tanimoto_volume:8]
# G20_2qwf_1_G_1.pdb G20_2qwf_1_G_1.pdb 24 24  0 0.271461 834.800354 0.822288
# G20_2qwf_1_G_1.pdb GNA_2qwe_1_G_1.pdb 24 23 -1 1.075027 773.756958 0.763037
# G20_2qwf_1_G_1.pdb DAN_2qwc_1_G_1.pdb 24 20 -4 1.050479 664.700317 0.661045
# G20_2qwf_1_G_1.pdb 49A_1f8e_1_H_1.pdb 24 20 -4 1.094503 655.997986 0.646492

  colX0 = int(colX)-1

  f = open(input_file)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      field = line.split()
      if (colX0<len(field)):
        X.append(field[colX0])
    pass
  f.close()


def read_tabular_file_three_columns(input_file,colX,colY,colZ,X,Y,Z):
  print "#read_tabular_file_three_columns('%s' colXYZ %s %s %s)"%(input_file,colX,colY,colZ)
  if (os.access(input_file,os.R_OK)==0):
    print "#ERROR:Can't open tabular file '%s'"%(input_file)  
    sys.exit(1)

# #[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [NheavyB-NheavyA:5]
# #[rmsd:6] [volume_overlap:7] [tanimoto_volume:8]
# G20_2qwf_1_G_1.pdb G20_2qwf_1_G_1.pdb 24 24  0 0.271461 834.800354 0.822288
# G20_2qwf_1_G_1.pdb GNA_2qwe_1_G_1.pdb 24 23 -1 1.075027 773.756958 0.763037
# G20_2qwf_1_G_1.pdb DAN_2qwc_1_G_1.pdb 24 20 -4 1.050479 664.700317 0.661045
# G20_2qwf_1_G_1.pdb 49A_1f8e_1_H_1.pdb 24 20 -4 1.094503 655.997986 0.646492

  colX0 = int(colX)-1
  colY0 = int(colY)-1
  colZ0 = int(colZ)-1

  f = open(input_file)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      field = line.split()
      if (colX0<len(field)) and (colY0<len(field)) and (colZ0<len(field)):
        X.append(field[colX0])
        Y.append(field[colY0])
        Z.append(field[colZ0])
    pass
  f.close()





################
##### MAIN #####
################

OPT = {}



OPT['crmsd'] = '6'
OPT['maxrmsd'] = '2.0'
OPT['alg'] = 'I'
OPT['oXY'] = ''
OPT['threX'] = '0.5'
OPT['cT'] = '1'
OPT['cR'] = '2'
OPT['self'] = 'T'
OPT['isum'] = ''
OPT['isim'] = ''
OPT['simtype'] = 'A'
OPT['csim'] = '9'
OPT['simthre'] = '0.7'
OPT['osum'] = '-'


if (len(sys.argv)<2):
  print "sim_success.py <options>"
  print " for calculating success rate using similarities."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -isum    : input summary file of RMSD [%s]"%(OPT['isum']) 
  print " -crmsd   : column number for rmsd (target_variable)       [%s]"%(OPT['crmsd'])
  print " -maxrmsd : max-value for correct range Y (Y<=maxY) [%s]"%(OPT['maxrmsd'])
  print " -cT      : column number for target    ligand [%s]"%(OPT['cT'])
  print " -cR      : column number for reference ligand [%s]"%(OPT['cR'])
  print " -self    : 'T': only for target==reference, 'N':only for target!=reference (cross). 'A':accept all. [%s]"%(OPT['self'])
  print " -osum    : output summary file [%s]"%(OPT['osum'])
  print "<options for similarities>"
  print " -isim     : input files for similarities,corresponding summary_file [%s]"%(OPT['isim']) 
  print " -csim     : column number for similarity [%s]"%(OPT['csim'])
  print " -simthre  : threshold value for similarity  [%s]"%(OPT['simthre']) 
  print " -simtype  : 'U':(sim>=thre), 'L':(sim<thre),'A':accept all [%s]"%(OPT['simtype']) 
  sys.exit(1)

read_option(sys.argv,OPT)


crmsd = int(OPT['crmsd'])
cT    = int(OPT['cT'])
cR    = int(OPT['cR'])
maxrmsd = float(OPT['maxrmsd'])
simthre = float(OPT['simthre'])


### Read summary file ###
if (os.access(OPT['isum'],os.R_OK)==0):
  print "#ERROR:Can't open summary file '%s'."%(OPT['isum'])
  sys.exit(1)

TAR_IDs = []
REF_IDs = []
RMSDs   = []
read_tabular_file_three_columns(OPT['isum'],OPT['cT'],OPT['cR'],OPT['crmsd'],TAR_IDs,REF_IDs,RMSDs)

### Read similarity file ###
if (os.access(OPT['isim'],os.R_OK)==0):
  print "#ERROR:Can't open similarity_file '%s'."%(OPT['isim'])
  sys.exit(1)

SIMs = []
read_tabular_file_one_column(OPT['isim'],OPT['csim'],SIMs)
print "#Npair for summary file    : %d"%(len(TAR_IDs))
print "#Npair for similarity file : %d"%(len(SIMs))

### Count Nsuccess ###
if (OPT['osum'] != ''):
  of = open(OPT['osum'],'w')
  of.write("#COMMAND:%s\n"%(OPT['COMMAND']))
  of.write("#DATE:   %s\n"%(OPT['START_DATE']))
  of.write("#[similarity:1] [rmsd:2] [tar_id:3] [ref_id:4]\n")

Nall     = 0
Nsuccess = 0
for i in range(len(TAR_IDs)):
  tar_id = TAR_IDs[i]
  ref_id = REF_IDs[i]
  rmsd   = float(RMSDs[i])
  sim    = float(SIMs[i])
  if ((OPT['self']=='T') and (tar_id==ref_id)) or ((OPT['self']=='N') and (tar_id!=ref_id)) or (OPT['self']=='A'):
    if ((OPT['simtype']=='U') and (sim>=simthre)) or ((OPT['simtype']=='L') and (sim<simthre)) or (OPT['simtype']=='A'):
      Nall += 1
      if (OPT['osum'] != ''):
        of.write("%f %f %s %s\n"%(sim,rmsd,tar_id,ref_id))

      if (rmsd <= maxrmsd):
        Nsuccess += 1
       

success_rate = 100.0*float(Nsuccess)/float(Nall)
print "#%s self %s simtype %s %s all %d success %d success_rate %6.1f %%"%(OPT['isum'],OPT['self'],OPT['simtype'],OPT['simthre'],Nall,Nsuccess,success_rate)

if (OPT['osum'] != ''):
  of.close()
  print "#write_summary_file() --> '%s'"%(OPT['osum'])
