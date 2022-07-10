#!/usr/bin/env python
import sys
import os
import random
from datetime import datetime

LastModDate = "Oct 19, 2011"


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



def read_target_predictor_valus_in_all_vs_all_file(ifname,col_tar,col_pres,unitype_pre,val_tar,val_pre):

##>>FILE FORMAT EXAMPLE OF KCOMBU ALL-VS-ALL SIMILARITY LIST <<
# #>>Similarities in List format<<
# #COMMAND /home/takawaba/work/kcombu/dkcombu -M A -ide A_txt.des -mcs T -con C -osl avaA_BuC.sl
# #DATE_START Sep 8,2011 11:31:16
# #LIBRARY_DIRECTORY /DB/LIGAND-EXPO/SDF2d
# #NMOL_IN_LIBRARY 763
# #NDATA 763
# #DATA 0 AQO
# #DATA 1 A27
# #DATA 2 ANG
# #DATA 3 AZ0
# #DATA 760 ALD
# #DATA 761 A5T
# #DATA 762 AMJ
# #[numA:1] [nameA:2] [numB:3] [nameB:4] [sim_mcs:5] [sim_pair:6] [sim_pair_norm:7] [sim_one:8] [sim_ring:9]
# 0 AQO 1 A27 0.250000 0.133333 0.365173 0.379310 0.229508
# 0 AQO 2 ANG 0.500000 0.466667 0.675835 0.600000 0.205128
# 0 AQO 3 AZ0 0.125000 0.025373 0.153464 0.148936 0.125000
# 0 AQO 4 AMZ 0.307692 0.207469 0.452996 0.478261 0.037037
# 0 AQO 5 ASY 0.031250 0.000000 0.033764 0.031250 0.046512
# 759 AC6 761 A5T 0.177778 0.057903 0.246993 0.204545 0.210526
# 759 AC6 762 AMJ 0.333333 0.192771 0.437048 0.555556 0.666667
# 760 ALD 761 A5T 0.100000 0.100101 0.279395 0.222222 0.274725
# 760 ALD 762 AMJ 0.238095 0.281513 0.529795 0.529412 0.600000
# 761 A5T 762 AMJ 0.196078 0.129683 0.347254 0.297872 0.315789
# #COMMAND    /home/takawaba/work/kcombu/dkcombu -M A -ide A_txt.des -mcs T -con C -osl avaA_BuC.sl
# #DATE_START Sep 8,2011 11:31:16
# #DATE_END   Sep 8,2011 11:41:43
# #COMP_TIME  626.676166 seconds
  print "#read_target_predictor_valus_in_all_vs_all_file('%s',col_tar:%s col_pre:%s unitype_pre:%s)"%(ifname,col_tar,col_pres, unitype_pre)
  Ncol_tar = int(col_tar) - 1
  if (col_pres != ''):
    col_pre_array = col_pres.split(':') 
    Ncol_pre_array = []
    for x in (col_pre_array):
      Ncol_pre_array.append(int(x)-1)
    #print "Ncol_pre_array",Ncol_pre_array  

  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open list_all_vs_all file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#')==0) and (len(line)>10) and (len(field)>=Ncol_tar):
      #print line
      if (Ncol_tar>=0):
        val_tar.append(float(field[Ncol_tar]))
      if (col_pres != ''):
        vals = []
        for Ncol_pre in (Ncol_pre_array):
          vals.append(float(field[Ncol_pre]))
        uni_val = unify_values(vals,unitype_pre)
      #for x in (vals):
      #  sys.stdout.write(" %f"%(x))
      #sys.stdout.write(" %f\n"%(uni_val))
        val_pre.append(uni_val)
  f.close()





def unify_values(vals,type):
  ret = 0.0

  if (type=='I'):  ## minimum ##
    ret = vals[0]
    for x in (vals):
      if (x<ret):
        ret = x
  elif (type=='X'):  ## maximum ##
    ret = vals[0]
    for x in (vals):
      if (x>ret):
        ret = x
  elif (type=='A'):  ## average ##
    ret = 0.0
    for x in (vals):
      ret += x
    ret /= len(vals)
  else:
    print "#ERROR(unify_values):improper type '%s'"%(type) 
    sys.exit(1) 
  return(ret)



############
### MAIN ###
############

OPT = {}
OPT['of'] = 'stdout'
OPT['ct'] = '5'
OPT['cp'] = '6'
OPT['tht'] = '0.6'
OPT['G'] = 'S'
OPT['up'] = 'A'
OPT['if']  = ''
OPT['ift'] = ''
OPT['rand'] = 'T'

if (len(sys.argv)<2):
  print "simcomp_recpre.py [ava_score_list_file] <options>"
  print " for comparison of similarities using all-vs-all comparison"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -if  : ava_score_list_file [%s]"%(OPT['if'])
  print " -ift : ava_score_list_file only for target  (optional) [%s]"%(OPT['ift'])
  print " -of  : output file [%s]"%(OPT['of'])
  print " -ct  : column num for the target variable (correct standard) [%s]"%(OPT['ct'])
  print " -cp  : column nums (such as 7:8) for the predictor variables [%s]"%(OPT['cp'])
  print " -up  : unifiying type for the predictors 'I':minimum 'X':maximum 'A'verage [%s]"%(OPT['up'])
  print " -tht : threshold value of target variable for the true positive class) [%s]"%(OPT['tht'])
  print " -G   : algorithm .'N'aive, 'S'orted [%s]"%(OPT['G'])
  print " -rand: random shuffle for target values [%s]"%(OPT['rand'])
  sys.exit(1)


read_option(sys.argv,OPT)
OPT['tht'] = float(OPT['tht'])

datalist = []
dataclass = {}
scmat = {} 

### [1] read similarity list file  ###
val_tar  = []
val_pre  = []
if (OPT['if']!=''):
  if (OPT['ift']==''):
    read_target_predictor_valus_in_all_vs_all_file(OPT['if'], OPT['ct'],OPT['cp'], OPT['up'],val_tar,val_pre)
  else:
    read_target_predictor_valus_in_all_vs_all_file(OPT['ift'],OPT['ct'],'',        '',       val_tar,[]      )
    read_target_predictor_valus_in_all_vs_all_file(OPT['if'], -1,       OPT['cp'], OPT['up'],[],     val_pre)


Ndata = len(val_tar)
print "#Ndata %d"%(Ndata)

if (OPT['rand']=='T'):
  random.shuffle(val_tar)

Nreal1 = 0
class_tar = []
for i in range(Ndata):
  if (val_tar[i]>=OPT['tht']):
    class_tar.append(1)
    Nreal1 += 1
  else:
    class_tar.append(0)

index = [i for i in range(Ndata)]
sindex  = sorted(index, lambda x,y:cmp(val_pre[x],val_pre[y]))





### [2] output in incresing order of val_pre[] (naive algorithm )###
if (OPT['of']=='stdout'):
  of = sys.stdout
else:
  print "#write to -->'%s'"%(OPT['of'])
  of = open(OPT['of'],'w')
of.write("#COMMAND:%s\n"%(OPT['COMMAND']))
of.write("#DATE   :%s\n"%(OPT['START_DATE']))

N = [[0 for i in range(2)] for j in range(2)] 

if (OPT['G']=='N'):
  for i in (sindex):
    thre =  val_pre[i]
    N[0][0] = N[0][1] = N[1][0] = N[1][1] = 0 
    for j in range(Ndata):
      real = class_tar[j]
      if (val_pre[j]>=thre):
        pred = 1
      else: 
        pred = 0
      N[real][pred] += 1

    if ((N[1][0]+N[1][1])>0): 
      recall    = float(N[1][1])/(float(N[1][0]+N[1][1])) 
    if ((N[0][1]+N[1][1])>0): 
      precision = float(N[1][1])/(float(N[0][1]+N[1][1])) 
    fpr = float(N[0][1])/float(N[0][0]+N[0][1])
    tpr = float(N[1][1])/float(N[1][0]+N[1][1])
    hitratio = float(N[0][1]+N[1][1])/float(N[0][0]+N[0][1]+N[1][0]+N[1][1])
    of.write("%f %f %f %f %f %f %f\n"%(val_pre[i],recall,precision,fmeasure,fpr,tpr,hitratio))
  of.close()
  sys.exit(1)

### [2] output in incresing order of val_pre[] (sorted algorithm)###
if (OPT['G']=='S'):
  of.write("#[thre:1] [recall(11/1*):2] [precision(11/*1):3] [fmeasure:4] [fpr(01/0*):5] [tpr(11/1*):6] [hitratio(*1/**):7]\n")
  N[0][0] = 0
  N[1][0] = 0
  N[0][1] = Ndata - Nreal1
  N[1][1] = Nreal1 
  thre0 = -1.0
  for i in (sindex):
    thre =  val_pre[i]
    if (class_tar[i]==1): 
      real = 1
    else:
      real = 0
    N[real][0] += 1
    N[real][1] -= 1
    #print "#N %d %d %d %d"%(N[0][0],N[0][1],N[1][0],N[1][1])
    if ((N[1][0]+N[1][1])>0): 
      recall    = float(N[1][1])/(float(N[1][0]+N[1][1])) 
    if ((N[0][1]+N[1][1])>0): 
      precision = float(N[1][1])/(float(N[0][1]+N[1][1])) 
    fpr = float(N[0][1])/float(N[0][0]+N[0][1])
    tpr = float(N[1][1])/float(N[1][0]+N[1][1])
    hitratio = float(N[0][1]+N[1][1])/float(N[0][0]+N[0][1]+N[1][0]+N[1][1])
    if (recall>0.0) or (precision>0.0):
      fmeasure = 2*recall*precision/(recall + precision)
    if (thre != thre0):
      of.write("%f %f %f %f %f %f %f\n"%(val_pre[i],recall,precision,fmeasure,fpr,tpr,hitratio))
    thre0 = thre
  of.close()
  sys.exit(1)  

