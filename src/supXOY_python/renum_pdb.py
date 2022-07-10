#!/usr/bin/env python

##
## <supXOYatms.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import pdb
import glob

LastModDate = '2019/01/31'


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




################
##### MAIN #####
################

OPT = {}
OPT['ipdb'] = ''
OPT['opdb'] = '-'
OPT['a'] = 'T'
OPT['r'] = 'T'
OPT['sta_anum'] = '1'
OPT['sta_rnum'] = '1'
OPT['ch'] = ''
OPT['new_conect'] = ''

if (len(sys.argv)<2):
  print "renum_pdb.py <options>"
  print "  renumbering atomic number and residue number."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print "  -ipdb: input PDB file [%s]"%(OPT['ipdb']) 
  print "  -opdb: output PDB filename [%s]"%(OPT['opdb'])
  print "  -ch  : new chain identifier [%s]"%(OPT['ch'])
  print "  -a: renumbering atomic  number ('T' or 'F')[%s]"%(OPT['a'])
  print "  -r: renumbering residue number ('T' or 'F')[%s]"%(OPT['r'])
  print "  -sta_anum : start atomic  number [%s]"%(OPT['sta_anum'])
  print "  -sta_rnum : start residue number [%s]"%(OPT['sta_rnum'])
  print "  -new_conect : new CONECT atom pair (anum1):(anum2) [%s]"%(OPT['new_conect'])
  sys.exit()

read_option(sys.argv, OPT)


P = pdb.Molecule()
P.read(OPT['ipdb'],AHtype='B')
#write_in_pdb(T,OPT['opdb'],NEWRNUM=OPT['tnewrnum'],COMMAND=OPT['COMMAND'],DATE=OPT['START_DATE'])

if (OPT['opdb'] == '-'):
  of = sys.stdout
else:
  print "#write_to --> '%s'"%(OPT['opdb'])
#          1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM   1071  HZ2 LYS A  64     -10.286   7.104  12.073  1.00  0.00           H
#ATOM   1072  HZ3 LYS A  64     -10.003   6.378  10.583  1.00  0.00           H

of = open(OPT['opdb'],'w')
of.write("REMARK COMMAND %s\n"%(OPT['COMMAND']))
of.write("REMARK DATE    %s\n"%(OPT['START_DATE']))

anum = int(OPT['sta_anum'])
rnum = int(OPT['sta_rnum']) - 1
rnum_str0 = ''

anum_new_frm_orig = {}

### REMARK LINES ###
f = open(OPT['ipdb'])
for line in f:
  if (line.startswith('REMARK')):
    of.write("%s"%(line))
f.close()

for i in range(P.Natom):
  #print i,mol.line[i]
  line = P.line[i]
  head_line = P.line[i][0:6]
  anum_str  = P.line[i][6:11]
  mid_line  = P.line[i][11:20]
  chain_str = P.line[i][21:22]
  rnum_str  = P.line[i][22:26]
  tail_line = P.line[i][26:]

  anum_orig = int(anum_str)
  if (OPT['a']=='T'):
    anum_str = "%5d"%(anum) 
    anum_new_frm_orig[anum_orig] = anum
  anum += 1

  if (rnum_str != rnum_str0):
    rnum += 1
  rnum_str_new = rnum_str
  if (OPT['r']=='T'):
    rnum_str_new = "%4d"%(rnum) 
   
  if (OPT['ch']!=''):
    chain_str = OPT['ch'] 
  of.write("%s%s%s %s%s%s\n"%(head_line,anum_str,mid_line,chain_str,rnum_str_new,tail_line))
  rnum_str0 = rnum_str 

#print anum_new_frm_orig

### CONECT LINES ###
if (OPT['a']=='T'):
  f = open(OPT['ipdb'])
  #          1         2         3         4
  #01234567890123456789012345678901234567890123456789012345678901234567890123456789
  #CONECT    1    2    9
  #CONECT    2    1    3    5   10
  #CONECT    3    2    4
  #CONECT    4    3
  #CONECT    5    2    6   11   12
  #CONECT    6    5    7    8
  for line in f:
    if (line.startswith('CONECT')):
      line = line.rstrip('\n')
      X = line.split()
      #print X
      first_anum_orig = int(X[1])
      if (first_anum_orig in anum_new_frm_orig): 
        for i in range(len(X)):
          if (i==0):
            of.write("CONECT")
          else:
            anum_orig = int(X[i])
            if (anum_orig in anum_new_frm_orig): 
              anum_new  = anum_new_frm_orig[anum_orig]
              of.write("%5d"%(anum_new))
        of.write('\n')

  if (OPT['new_conect'] != ''):
    (anum1,anum2) = OPT['new_conect'].split(':')
    of.write("CONECT%5s%5s\n"%(anum1,anum2))
    pass

if (OPT['opdb'] !='-'):
  print "#write_to --> '%s'"%(OPT['opdb'])
  of.close()

