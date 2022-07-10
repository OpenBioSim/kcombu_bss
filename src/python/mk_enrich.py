#!/usr/bin/env python
import sys
import os
import random

LastModDate = "Apr 21, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


 
def read_list_file(ifname,list):
  print "#read_list_file('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     list.append(fields[0])
  f.close()



def read_all_vs_all_active_inactive_file(ifname,actlist,inactlist,dataclass,scmat):
## actlist[0,1,2,...]        : list of data (active)
## inactlist[0,1,2,...]      : list of data (inactive)
## scmat{'datanameA'}{'datanameB'} : similarity matrix.

# #>>Similarities in List format<<
# #COMMAND kcombu -M A -ill all.list -act T -alg B -con C -osl ACT_INACT/adeno.Buc.ava
# #DATE_START Apr 21,2011 7:19:35
# #DATE_END   Apr 21,2011 7:21:14
# #COMP_TIME  98.573006 seconds
# #NDATA 1332
# #DATA 0 adeno_lig/2-hexy.mol2 active
# #DATA 1 adeno_lig/2chlade.mol2 active
# #DATA 2 adeno_lig/AB-MECA.mol2 active
# #DATA 3 adeno_lig/APNEA.mol2 active
# #DATA 28 adeno_lig/sakuran.mol2 active
# #DATA 29 adeno_lig/theophy.mol2 active
# #DATA 30 adeno_lig/visnagi.mol2 active
# #DATA 31 base/0000353-01.mol2 inactive
# #DATA 32 base/0000453-01.mol2 inactive
# #DATA 33 base/0000598-01.mol2 inactive
# #DATA 34 base/0000697-01.mol2 inactive
# #DATA 1330 base/tlglu.mol2 inactive
# #DATA 1331 base/tlthr.mol2 inactive
# 0 1 0.655172
# 0 2 0.540541
# 0 3 0.473684
# :
# 30 1328 0.222222
# 30 1329 0.200000
# 30 1330 0.193548
# 30 1331 0.206897

  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open list_all_vs_all file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)

  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#NDATA')):
      N = int(field[1])
      dataname = ['' for i in range(N)]
    if (line.startswith('#DATA')):
      dataname[int(field[1])] = field[2]
      if (field[3]=='active'):
        actlist.append(field[2])
        dataclass[field[2]] = 'active'
      if (field[3]=='inactive'):
        inactlist.append(field[2])
        dataclass[field[2]] = 'inactive'
    if (line.startswith('#')==0) and (len(line)>10) and (len(field)==3):
      i = int(field[0])
      j = int(field[1])
      if (scmat.has_key(dataname[i])==0):
        scmat[dataname[i]] = {}
      if (scmat.has_key(dataname[j])==0):
        scmat[dataname[j]] = {}
      scmat[dataname[i]][dataname[j]] =  float(field[2])
      scmat[dataname[j]][dataname[i]] =  float(field[2])
  
      #print "%s %f"%(index,scmat[index])
      #print datalist[i],datalist[j],scmat[datalist[i]][datalist[j]],scmat[datalist[j]][datalist[i]]

  f.close()

def randomize_label(datalist,dataclass):
  print "#randomize_label()"
  Nlist = len(datalist)
  for r in range(Nlist):
    i = random.randint(0,Nlist-1) 
    j = random.randint(0,Nlist-1) 
    x = datalist[i]
    y = datalist[j]
    buff = dataclass[x]
    dataclass[x] = dataclass[y]
    dataclass[y] = buff

############
### MAIN ###
############

OPT = {}
OPT['of']   = 'stdout'
OPT['xred'] = 'F'
OPT['omn'] = ''
if (len(sys.argv)<2):
  print "mk_enrich.py [ava_score_list_file]<options>"
  print " for making cumulative recall (enrichment curve) for scores of acitive-inactive comparison"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -of   : output file [%s]"%(OPT['of'])
  print " -omn  : output mean cumulative file [%s]"%(OPT['omn'])
  print " -xred : exclude redundancy (T or F)[%s]"%(OPT['xred'])
  sys.exit(1)



iavafile = sys.argv[1]
read_option(sys.argv,OPT)

### [1] read similarity list file  ###
actlist = []
inactlist = []
dataclass = {}
scmat = {}
read_all_vs_all_active_inactive_file(iavafile,actlist,inactlist,dataclass,scmat)
Nact   = len(actlist)
Ninact = len(inactlist)

Nact_all = Nact - 1
Nall      = Nact - 1 + Ninact
## [2] output all the cumulative curves ##
if (OPT['of']=='stdout'):
  of = sys.stdout
else:
  print "#write_all_the_cumulative_curves() -->'%s'"%(OPT['of'])
  of = open(OPT['of'],'w')

of.write("#>>Cumurative Curves for Each Active Compound\n")
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#[ratio_data_tested] [ratio_act_found] [name_of_data] [score] [class]\n")
ave_ratio_act_found = {}

for ref in (actlist):
  of.write("#>ref_act %s\n"%(ref))
  list  = sorted(scmat[ref].keys(), lambda x,y:cmp(scmat[ref][y],scmat[ref][x]))
  Nact_hit = 0
  Ntest    = 0
  ratio_data_tested0 = 0.0
  ratio_act_found0   = 0.0
  auc = 0.0 
  for x in (list):
    Ntest += 1
    if (dataclass[x]=='active'):
      Nact_hit += 1
    ratio_data_tested  = float(Ntest)/float(Nall) 
    ratio_act_found = float(Nact_hit)/float(Nact_all)
    ave_ratio_act_found[ratio_data_tested] = ave_ratio_act_found.get(ratio_data_tested,0.0) + ratio_act_found 
    auc += (ratio_act_found + ratio_act_found0) * (ratio_data_tested - ratio_data_tested0) * 0.5
    ratio_act_found0   = ratio_act_found
    ratio_data_tested0 = ratio_data_tested
    of.write("%f %f %s %s %f %s\n"%(ratio_data_tested,ratio_act_found,ref,x,scmat[ref][x],dataclass[x]))
  of.write("#%s ref_act %s auc %f\n"%(iavafile,ref,auc))
  of.write("\n")

if (OPT['of'] != 'stdout'):
  of.close()

## [3] output average cumulative curves ##
if (OPT['omn'] != ''):
  of = open(OPT['omn'],'w')
  print "#write_mean_culmulative_curve() -->'%s'"%(OPT['omn'])
  of.write("#>>Mean Cumurative Curve\n")
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[ratio_data_tested] [averaged_ratio_act_found]\n")
  ratio_data_tested_list  = sorted(ave_ratio_act_found.keys(), lambda x,y:cmp(x,y))
  ratio_data_tested0 = 0.0
  ratio_act_found0   = 0.0
  auc = 0.0 
  for ratio_data_tested in (ratio_data_tested_list):
    ratio_act_found = ave_ratio_act_found[ratio_data_tested]/Nact
    of.write("%f %f\n"%(ratio_data_tested,ratio_act_found))
    auc += (ratio_act_found + ratio_act_found0) * (ratio_data_tested - ratio_data_tested0) * 0.5
    ratio_act_found0   = ratio_act_found
    ratio_data_tested0 = ratio_data_tested
  of.write("#%s mean_auc %f\n"%(iavafile,auc))
  of.close()

sys.exit(1)
