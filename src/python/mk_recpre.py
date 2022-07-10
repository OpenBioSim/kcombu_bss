#!/usr/bin/env python
import sys
import os
import random

LastModDate = "June 19, 2011"

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



def read_list_all_vs_all_file(ifname,datalist,dataclass,scmat):
## datalist[0,1,2,...]        : list of data
## dataclass{dataname}        : property(class) of data 
## scmat{'datanumA:datanumB'} : similarity matrix.

##>>FILE FORMAT EXAMPLE OF KCOMBU ALL-VS-ALL SIMILARITY LIST <<
##COMMAND kcombu -M A -ill smallset -fL P -osl all_vs_all/small.ava -con C -alg B
##START_DATE Apr 13,2011 13:36:59
##END_DATE   Apr 13,2011 13:36:59
##NDATA 10
##DATA 0 lig_CDK2/SupLIG/VAR_3bhvA CDK2
##DATA 1 lig_CDK2/SupLIG/DT5_2c6mA CDK2
##DATA 2 lig_CDK2/SupLIG/C85_2uzdA CDK2
##DATA 3 lig_CDK2/SupLIG/LS1_1ke5A CDK2
##DATA 4 lig_THR/SupLIG/81A_1t4uH THR
##DATA 5 lig_THR/SupLIG/MEL_1k22H THR
##DATA 6 lig_THR/SupLIG/BT3_1d3pB THR
##DATA 7 lig_THR/SupLIG/49U_2zhfH THR
##DATA 8 lig_THR/SupLIG/T42_1ai8H THR
##DATA 9 lig_THR/SupLIG/110_1g37A THR
#0 1 0.166667
#0 2 0.250000
#0 3 0.216216
#0 4 0.102041
#0 5 0.104167
#:
#5 9 0.183333
#6 7 0.129630
#6 8 0.051948
#6 9 0.082192
#7 8 0.254902
#7 9 0.319149
#8 9 0.154930
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open list_all_vs_all file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#NDATA')):
      N = int(field[1])
    if (line.startswith('#DATA')):
      datalist.append(field[2])
      dataclass[field[2]] = field[3]
    if (line.startswith('#')==0) and (len(line)>10) and (len(field)==3):
      i = int(field[0])
      j = int(field[1])
      index = "%d:%d"%(i,j)
      scmat[index] = float(field[2])
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
OPT['rand'] = 'F'
OPT['of']   = 'stdout'
OPT['xred'] = 'F'

if (len(sys.argv)<2):
  print "mk_recpre.py [ava_score_list_file]<options>"
  print " for making recall-precision plot for scores of all-vs-all comparison"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -of   : output file [%s]"%(OPT['of'])
  print " -rand : randomize label (T or F)[%s]"%(OPT['rand'])
  print " -xred : exclude redundancy (T or F)[%s]"%(OPT['xred'])
  sys.exit(1)



iavafile = sys.argv[1]
read_option(sys.argv,OPT)
datalist = []
dataclass = {}
scmat = {} 

### [1] read similarity list file  ###
read_list_all_vs_all_file(iavafile,datalist,dataclass,scmat)

if (OPT['rand']=='T'):
  randomize_label(datalist,dataclass)

Nlist  = len(datalist)
Npair = 0
Npair_idprop = 0
for i in range(Nlist):
  x = datalist[i] 
  for j in range(i+1,Nlist):
    y = datalist[j] 
    Npair += 1
    if (dataclass[x]==dataclass[y]): 
      Npair_idprop += 1


## [2] Sort index (i:j) by scmat ##

indexlist  = sorted(scmat.keys(), lambda x,y:cmp(scmat[x],scmat[y]))

## [3] Output recall and precision for (i:j) by increasing order of score ##
if (OPT['of']=='stdout'):
  of = sys.stdout
else:
  print "#write to -->'%s'"%(OPT['of'])
  of = open(OPT['of'],'w')

 
of.write("#RECALL-PRECITION TPR-FPR DATA FILE\n")
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#Nlist %d\n"%(len(datalist)))
of.write("#Npair %d Npair_idprop %d\n"%(Npair,Npair_idprop))
N = [[0 for i in range(2)] for j in range(2)] 
N[0][0] = 0
N[1][0] = 0
N[0][1] = Npair - Npair_idprop 
N[1][1] = Npair_idprop 
auc_roc = 0.0
opt_fmeasure  = 0.0
opt_thre      = 0.0
opt_recall    = 0.0
opt_precition = 0.0

fpr0 = 1.0
tpr0 = 1.0

ave_score0 = 0.0
ave_score1 = 0.0
count0 = 0
count1 = 0

of.write("#[thre:1] [recall:2] [precision:3] [F-measure:4] [false positive rate:5] [true positive rate(recall):6]\n")
#for ind in (indexlist):
for p in range(len(indexlist)):
  ind = indexlist[p]
  output = 0
  if (p==(len(indexlist)-1)):
    output = 1
  elif (OPT['xred']!='T'):
    output = 1
  elif (OPT['xred']=='T') and (scmat[indexlist[p]] != scmat[indexlist[p+1]]):
    output = 1

  thre = scmat[ind] 
  #print "%s %f"%(ind,scmat[ind])
  field = ind.split(':')
  i = int(field[0])
  j = int(field[1])
  x = datalist[i]
  y = datalist[j]
  real = 0 
  if (dataclass[x]==dataclass[y]): 
    real = 1
  N[real][0] += 1
  N[real][1] -= 1

  if (real==1):
    count1 += 1
    ave_score1 += thre
  else:
    count0 += 1
    ave_score0 += thre

  if (output==1):
    recall = 0.0
    precision = 0.0
    fmeasure = 0.0 
    fpr = 0.0 # false positive rate (N(r=0,p=1)/N(r=0))
    tpr = 0.0 # true positive rate (recall)
    if ((N[1][0]+N[1][1])>0): 
      recall    = float(N[1][1])/(float(N[1][0]+N[1][1])) 
    if ((N[0][1]+N[1][1])>0): 
      precision = float(N[1][1])/(float(N[0][1]+N[1][1])) 
  
    fpr = float(N[0][1])/float(N[0][0]+N[0][1])
    tpr = float(N[1][1])/float(N[1][0]+N[1][1])
    auc_roc += (tpr0+tpr)*(fpr0-fpr)/2.0
    fpr0 = fpr 
    tpr0 = tpr 
    if (recall>0.0) or (precision>0.0):
      fmeasure = 2*recall*precision/(recall + precision)
    #print thre,recall,precision,real,N[0][0],N[1][0],N[0][1],N[1][1]
    #of.write("%f %f %f %f %f %f\n"%(thre,recall,precision,fmeasure,fpr,tpr))
    of.write("%f %f %f %f %f %f N00 %d N01 %d N10 %d N11 %d\n"%(thre,recall,precision,fmeasure,fpr,tpr,N[0][0],N[0][1],N[1][0],N[1][1]))
    if (fmeasure > opt_fmeasure):
      opt_fmeasure  = fmeasure 
      opt_thre      = thre
      opt_recall    = recall 
      opt_precision = precision

## [4] ending procedure ##
ave_score0 /= count0 
ave_score1 /= count1 
of.write("#%s opt_fmeasure %f auc_roc %f opt_thre %f opt_recall %f opt_prection %f av_sc0 %f av_sc1 %f\n"%(iavafile,opt_fmeasure, auc_roc,opt_thre, opt_recall, opt_precision,ave_score0,ave_score1))
 
if (OPT['of'] != 'stdout'):
  of.close()
