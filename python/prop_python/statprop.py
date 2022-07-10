#!/usr/bin/env python
##
## <statprop.py>
##


import sys
import os
import re
import math
from datetime import datetime

LastModDate = "June 15, 2013"

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


def isfloat(str):
  for i in range(len(str)):
    if (str[i].isdigit()):
      return(1)
  return(0)



def modify_compound_id_string(str):
  # '04270000/04274776.sdf' --> '04274776'
  modstr = str
  if (modstr.find('/')):
    (head,tail) = modstr.split('/')
    modstr = tail
  if (modstr.find('.')):
    (head,tail) = modstr.split('.')
    modstr = head 
  return(modstr)


def read_list_file_as_dic(ifname,dic):
  print "#read_list_file_as_dic('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     #list.append(fields[0])
     dic[fields[0]] = 1
  f.close()


def get_statistics_value_list(list):
  init = 1
  S = 0.0
  SS = 0.0
  N = len(list) 
  for x in (list):
    f = float(x)
    S += f
    SS += f * f
    if (init==1):
      MIN = f 
      MAX = f 
    else:
      if (f<=MIN):
        MIN = f
      if (f>=MAX):
        MAX = f
    init = 0
  AVE = S/N
  VAR = SS/N - AVE*AVE
  return((AVE,VAR,MIN,MAX))
 

###############
#### MAIN #####
###############

OPT = {}
OPT['qn0'] = ''
OPT['qv0'] = ''
OPT['qo0'] = ''
OPT['c']   = ''
OPT['of'] = ''
OPT['oh'] = ''
OPT['og'] = 'out.gnu'
OPT['opng'] = 'out.png'
OPT['qvl'] = ''
OPT['qcid'] = 'F'
OPT['type'] = 'C'
OPT['bin'] = 'N';
OPT['Nbin'] = '50';
OPT['bw'] = '1.0'
OPT['per'] = '0.1'
OPT['title'] = ''
OPT['xlabel'] = ''
OPT['all'] = 'F'
OPT['gnu'] = 'T'
OPT['bwint'] = 'F'

if (len(sys.argv)<2):
  #x = '04270000/04274776.sdf'
  #print "'%s' '%s'"%(x,modify_compound_id_string(x))
  print "statprop.py [property_file] <options>"
  print " for analyzing property data file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print "<focusing column>" 
  print " -c     : column number for focused property (1,2,...)[%s]"%(OPT['c'])
  print " -type  : type of focusing value 'C'ontinuous, 'D'iscrete [%s]"%(OPT['type'])   
  print "<histogram>" 
  print " -oh    : output histogram file (for continuous data)[%s]"%(OPT['oh'])
  print " -bin   : policy for 'bin width' decision. 'N' from 'Nbin'. 'B' from just '-bw' [%s]"%(OPT['bin'])
  print " -Nbin  : number of bins  [%s]"%(OPT['Nbin'])
  print " -bw    : bin width for histogram [%s]"%(OPT['bw'])
  print " -bwint : enforce bin width into integer (T or F)[%s]"%(OPT['bwint']) 
  print " -per   : small/large triming percentile for focusing range (0..100) [%s]"%(OPT['per'])
  print " -title : title for gnuplot [%s]"%(OPT['title']) 
  print " -xlabel: xlabel for gnuplot [%s]"%(OPT['xlabel']) 
  print " -all   : use the plot for the all as control (T or F)[%s]"%(OPT['all']) 
  print "<gnuplot>" 
  print " -og    : output gnuplot script file (for continuous data)[%s]"%(OPT['oh'])
  print " -opng  : output graph image file (*.png) using gnuplot [%s]"%(OPT['opng'])
  print " -gnu   : execute gnuplot (T or F)[%s]"%(OPT['gnu'])
  print "<frequncies of discrete state>" 
  print " -of   : output frequency file (for discrete data)  [%s]"%(OPT['of'])
  print "<filtering condition>" 
  print " -qn0  : property id number (1,2,...) [%s]"%(OPT['qn0']) 
  print " -qv0  : value/string for i-th query [%s]"%(OPT['qv0']) 
  print " -qo0  : operator for i-th query. 'U'pper, 'L'ower, 'S'ubstring [%s]"%(OPT['qo0']) 
  print " -qvl  : input list file for the substring [%s]"%(OPT['qvl'])
  print " -qcid : compound_id modification '04270000/04274776.sdf'->'04274776'(T or F) [%s]"%(OPT['qcid'])
  sys.exit(1)

read_option(sys.argv,OPT)
ipropfile = sys.argv[1]
print "#QUERY:'%s' '%s' '%s'"%(OPT['qn0'],OPT['qo0'],OPT['qv0'])
if (OPT['qn0'] != ''):
  qn0  = int(OPT['qn0']) - 1
cn   = int(OPT['c']) - 1
OPT['bw'] = float(OPT['bw'])
ifname = sys.argv[1]
f = open(ipropfile)
Ntotal = 0
Nhit   = 0
if (OPT['qvl'] != ''):
  querydic = {} 
  hit_querydic = {} 
  read_list_file_as_dic(OPT['qvl'],querydic)
  Nquerydic  = len(querydic.keys())
  print "#Nquerydic %d"%(Nquerydic) 
  #sys.exit(1)
  #for x in (querydic.keys()):
  #  print "#querydic_key '%s'"%(x)
  #sys.exit(1)


##############################################################
### [1] Read file line one by one and check query patterns ###
##############################################################
FOCUS_FIELD_LIST     = []
FOCUS_FIELD_LIST_ALL = []

for line in f:
  line = line.rstrip('\n')
  if (len(line)>10) and (line.startswith('#')==0):
    Ntotal += 1
    field = line.split('\t')
    for i in range(len(field)):
      if (field[i].find(':')>=0):
        items = field[i].split(':')
        field[i] = items 
    #print field
    #print line
    ### (1) Check Query Pattern ###
    hit = 1
    #print "qn0 '%s' len %d\n"%(OPT['qn0'],len(field))
    if (OPT['qn0'] !='') and (len(field)>qn0):
      if (isfloat(field[qn0])==1):
        val = field[qn0]
        if (OPT['qo0'] == 'U') and (float(val)>float(OPT['qv0'])):
          hit = 0 
        if (OPT['qo0'] == 'L') and (float(val)<float(OPT['qv0'])):
          hit = 0 
      if (OPT['qo0'] == 'S'):
        if (OPT['qcid']=='T'):
          field[qn0] = modify_compound_id_string(field[qn0])
        if (OPT['qvl'] == ''):
          if (isinstance(field[qn0],str))  and (field[qn0].find(OPT['qv0'])==-1): 
            hit = 0
          if (isinstance(field[qn0],list)):
            for x in (field[qn0]):
              if (x.find(OPT['qv0'])==-1):
                hit = 0 
        elif (OPT['qvl'] != ''):
          if (isinstance(field[qn0],str)):
            if (querydic.has_key(field[qn0])==0):
              hit = 0
            #print "%s hit %d"%(field[qn0],hit)
            else:
              hit_querydic[field[qn0]] = 1
          if (isinstance(field[qn0],list)):
            hitlist = 0
            for x in (field[qn0]):
              if (querydic.has_key(x)==1):
                hitlist = 1
                hit_querydic[x] = 1
            if (hitlist==0):
              hit = 0

    ### (2) Column for focusing ###
    if (hit==1):
      Nhit += 1
      valc = field[cn]
      if (OPT['type']=='C'): 
        if (isfloat(valc)==1):
          FOCUS_FIELD_LIST.append(field[cn])
      if (OPT['type']=='D'): 
        #print field[cn]
        if (isinstance(field[cn],list)):
          for y in (field[cn]):
            FOCUS_FIELD_LIST.append(y)
            #print "#%s"%(y)
        else:
          FOCUS_FIELD_LIST.append(field[cn])

    ### (2) Column for all ###
    if (OPT['all'] == 'T'):
      valc = field[cn]
      if (OPT['type']=='C'): 
        if (isfloat(valc)==1):
          FOCUS_FIELD_LIST_ALL.append(field[cn])
      if (OPT['type']=='D'): 
        if (isinstance(field[cn],list)):
          for y in (field[cn]):
            FOCUS_FIELD_LIST_ALL.append(y)
        else:
          FOCUS_FIELD_LIST_ALL.append(field[cn])

f.close()

sys.stdout.write("#Nhit %d Ntotal %d\n"%(Nhit,Ntotal))
### FOR CONTINUOUS VALUES ###
if (OPT['type'] == 'C'):
  (AVEq,VARq,MINq,MAXq) = get_statistics_value_list(FOCUS_FIELD_LIST)
  FOCUS_FIELD_LIST = sorted(FOCUS_FIELD_LIST, lambda x,y:cmp(float(x),float(y)))
  MIN_RANGEq = float(FOCUS_FIELD_LIST[int(len(FOCUS_FIELD_LIST) * float(float(OPT['per'])*0.01))])
  MAX_RANGEq = float(FOCUS_FIELD_LIST[int(len(FOCUS_FIELD_LIST) * float(1.0 - float(OPT['per'])*0.01))-1])
  print "#FOCUS MIN %f MAX %f AVE %f VAR %f MIN_RANGE %f MAX_RANGE %f"%(MINq,MAXq,AVEq,VARq,MIN_RANGEq,MAX_RANGEq)
  MIN_RANGE = MIN_RANGEq
  MAX_RANGE = MAX_RANGEq

  if (OPT['all']=='T'):
    (AVEa,VARa,MINa,MAXa) = get_statistics_value_list(FOCUS_FIELD_LIST_ALL)
    FOCUS_FIELD_LIST_ALL = sorted(FOCUS_FIELD_LIST_ALL, lambda x,y:cmp(float(x),float(y)))
    MIN_RANGEa = float(FOCUS_FIELD_LIST_ALL[int(len(FOCUS_FIELD_LIST_ALL) * float(float(OPT['per'])*0.01))])
    MAX_RANGEa = float(FOCUS_FIELD_LIST_ALL[int(len(FOCUS_FIELD_LIST_ALL) * float(1.0 - float(OPT['per'])*0.01))-1])
    print "#ALL MIN %f MAX %f AVE %f VAR %f MIN_RANGE %f MAX_RANGE %f"%(MINa,MAXa,AVEa,VARa,MIN_RANGEa,MAX_RANGEa)
    MIN_RANGE = MIN_RANGEa
    MAX_RANGE = MAX_RANGEa

  if (OPT['bin']=='N'):
    OPT['bw'] = (MAX_RANGE - MIN_RANGE)/(int(OPT['Nbin']))
    if (OPT['bwint']=='T'):
      if (OPT['bw']<1.0):
        OPT['bw'] = 1.0
      else:
        OPT['bw'] =int(OPT['bw']) 

  print "#RANGE %f %f bw %f"%(MIN_RANGE,MAX_RANGE,OPT['bw'])

#######################################
### [2] Check FOCUS_FIELD_LIST again ### 
#######################################
Ncount = {}
Nhist_count = {}
Ncount_all = {}
Nhist_count_all = {}

if (OPT['all']=='T'):
  for x in (FOCUS_FIELD_LIST_ALL):
    if (OPT['type']=='C'): 
      fx = float(x) 
      if (fx==0.0):
        print "'%s' %f"%(x,fx)
      if (OPT['oh'] != ''):
        index = int(fx/OPT['bw'])       
        Nhist_count_all[index] = Nhist_count_all.get(index,0) + 1
    elif (OPT['type']=='D'): 
      if (OPT['of'] != ''):
        Ncount_all[x] = Ncount_all.get(x,0) + 1
  FOCUS_FIELD_LIST_ALL = []

for x in (FOCUS_FIELD_LIST):
  if (OPT['type']=='C'): 
    fx = float(x) 
    if (OPT['oh'] != ''):
      index = int(fx/OPT['bw'])       
      Nhist_count[index] = Nhist_count.get(index,0) + 1
  elif (OPT['type']=='D'): 
    if (OPT['of'] != ''):
      Ncount[x] = Ncount.get(x,0) + 1

FOCUS_FIELD_LIST = []




######################################
### [3] Output Histogram/Frequency ### 
######################################

if (OPT['oh']!=''):
  print "#output_histogram-->'%s'"%(OPT['oh'])
  of = open(OPT['oh'],'w')

  if (OPT['all']=='T'):
    index_list = sorted(Nhist_count_all.keys(), lambda x,y:cmp(x,y))
  else:
    index_list = sorted(Nhist_count.keys(), lambda x,y:cmp(x,y))

  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#Ntotal %d\n"%(Ntotal))
  of.write("#Nhit   %d\n"%(Nhit))
  of.write("#[state] [Frequency] [Count]\n")
  for x in (index_list):
    of.write("%f %f %d"%(x*OPT['bw'],float(Nhist_count.get(x,0))/Nhit,Nhist_count.get(x,0)))
    if (OPT['all'] == 'T'): 
      of.write(" %f %d\n"%(float(Nhist_count_all[x])/Ntotal,Nhist_count_all[x]))
    else:
      of.write("\n")
  of.close()

if (OPT['og']!='') and (OPT['type']=='C'):
  print "#output_gnuplot_script -->'%s'"%(OPT['og'])
  of = open(OPT['og'],'w')
  of.write("#\"%s\"\n"%(OPT['COMMAND']))
#  of.write("set title \"%s\"\n"%(OPT['title']))
  of.write("set title \"RED:ALL, BLUE:%s\"\n"%(OPT['qv0']))
  of.write("set xrange [%f:%f]\n"%(MIN_RANGE,MAX_RANGE))
  of.write("set xlabel \"%s\"\n"%(OPT['xlabel']))
  of.write("set ylabel \"Frequency\"\n")
  of.write("unset key\n")
  if (OPT['opng'] != ''):
    of.write("set term png large\n")
    of.write("set output \"%s\"\n"%(OPT['opng']))
  if (OPT['all']=='T'):
    of.write("plot \"%s\" u 1:2 w lp 3,\"%s\" u 1:4 w lp 1\n"%(OPT['oh'],OPT['oh']))
  else: 
    of.write("plot \"%s\" w lp 1\n"%(OPT['oh']))
  if (OPT['opng'] != ''):
    of.write("quit\n")
  of.close()

  if (OPT['gnu']=='T'):
    str = "gnuplot %s"%(OPT['og'])
    print "#%s"%(str)
    os.system(str)
    if (OPT['opng'] != ''):
      print "#gnuplot -->'%s'"%(OPT['opng'])

for x in (querydic.keys()):
  if (hit_querydic.has_key(x)==0):
    print "NO_HIT",x


if (OPT['of']!=''):
  of = open(OPT['of'],'w')
  print "#output_frequency-->'%s'"%(OPT['of'])
  if (OPT['all']=='T'):
    index_list = sorted(Ncount_all.keys(), lambda x,y:cmp((Ncount.get(y,0),Ncount_all[y]),(Ncount.get(x,0),Ncount_all[x])))
  else:
    index_list = sorted(Ncount.keys(), lambda x,y:cmp(Ncount[y],Ncount[x]))
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#Ntotal %d\n"%(Ntotal))
  of.write("#Nhit   %d\n"%(Nhit))
  if (OPT['all']!='T'):
    of.write("#[state] [Frequency_hit(%%)] [Count_hit]\n")
  elif (OPT['all']=='T'):
    of.write("#[state] [Count_hit] [Count_all] [Hit_ratio(%%)] [Frequency_hit(%%)] [Frequency_all(%%)]\n")

  for x in (index_list):
    if (OPT['all']!='T'):
      of.write("%-20s %5.2f %7d\n"%(x,float(100.0*Ncount.get(x,0))/Nhit,Ncount.get(x,0)))
    if (OPT['all']=='T'):
      of.write("%-20s %7d %7d  %5.2f %5.2f %5.2f\n"%(x,Ncount.get(x,0),Ncount_all.get(x,0),float(100.0*Ncount.get(x,0))/Ncount_all.get(x,0),float(100.0*Ncount.get(x,0))/Nhit,float(100.0*Ncount_all.get(x,0))/Ntotal))
  of.close()

