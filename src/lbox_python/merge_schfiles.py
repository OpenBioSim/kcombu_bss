#!/usr/bin/env python
##
## <merge_schfiles.py>
##


import sys
import os
import random
from datetime import datetime
import kcombu

LastModDate = "Jan 10, 2012"

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



###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['fus'] = 'A'

if (len(sys.argv)<3):
  print "merge_schfiles.py [schfile0] [schfile1] [schfile2] ... <options>"
  print " for counting files under '-idir'."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -fus : fusion method. 'A'verage, 'X':max [%s]"%(OPT['fus']) 
  print " -A   : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)


schfile_list = []
i = 0
stop = 0
while (i<len(sys.argv)) and (stop==0):
  if (sys.argv[i].startswith('-')==0):
    schfile_list.append(sys.argv[i])
  if (sys.argv[i].startswith('-')):
    stop = 1
  i += 1

hitmol_dic = {}
hitmol_schfile = {}
hitmol_similarity = {}

for schfile in (schfile_list):
  print schfile 
  hitmol_dic[schfile] = []
  kcombu.read_kcombu_search_result(schfile,hitmol_dic[schfile])
  for x in (hitmol_dic[schfile]):
    #print x['molname'],x['sim_mcs']
    if (hitmol_schfile.has_key(x['molname'])==0):
      hitmol_schfile[x['molname']] = []

    hitmol_schfile[x['molname']].append(schfile)
    if (hitmol_similarity.has_key(x['molname'])==0):
      hitmol_similarity[x['molname']] = {}
    hitmol_similarity[x['molname']][schfile] = x['sim_mcs']

score_fused = {}
for m in (hitmol_schfile.keys()):
  sys.stdout.write("%s "%(m))
  score_fused[m] = 0.0 
  Nadd = 0
  for schfile in (hitmol_schfile[m]):
    #sys.stdout.write(" %s %s"%(schfile,hitmol_similarity[m][schfile]))
    if (OPT['fus']=='A'):
      score_fused[m] += float(hitmol_similarity[m][schfile])
      Nadd += 1
    if (OPT['fus']=='X'):
      if (float(hitmol_similarity[m][schfile])>score_fused[m]):
        score_fused[m] = float(hitmol_similarity[m][schfile])
  if (Nadd>0):
    score_fused[m] /= Nadd

  sys.stdout.write(" Sfus %f"%(score_fused[m]))

  for schfile in (hitmol_schfile[m]):
    sys.stdout.write(" %s %s"%(schfile,hitmol_similarity[m][schfile]))

  sys.stdout.write("\n")





