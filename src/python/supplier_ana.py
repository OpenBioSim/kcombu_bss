#!/usr/bin/env python
##
## <supplier_ana.py>
##


import sys
import os
import re
import math

LastModDate = "Mar 8, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
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


###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['if'] = ''
OPT['c'] = '7'
OPT['sc'] = '16'
OPT['qs'] = ''
OPT['of'] = 'out'

if (len(sys.argv)<2):
  #print "'%s' %d"%('abc',isfloat('abc'))
  #print "'%s' %d"%('012',isfloat('012'))
  #print "'%s' %d"%('1.12',isfloat('1.12'))
  #print "'%s' %d"%('.12',isfloat('.12'))
  print "supplier_ana.py <options>"
  print " for analyzing property data file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -if  : input property file (tab-splited)  [%s]"%(OPT['if']) 
  print " -A   : Action (T or F) [%s]"%(OPT['A']) 
  print " -c   : column number for focused property (1,2,...)[%s]"%(OPT['c'])
  print " -sc  : column number for suppliers (1,2,...)[%s]"%(OPT['sc'])
  print " -qs  : query string of property (1,2,...)[%s]"%(OPT['qs'])
  print " -of  : output file [%s]"%(OPT['of'])
  sys.exit(1)

read_option(sys.argv,OPT)
OPT['sc'] = int(OPT['sc']) - 1
OPT['c']  = int(OPT['c'])  - 1

S = SS = 0.0
N = 0
ifname = sys.argv[1]
f = open(OPT['if'])
Ssup  = {}
SSsup = {}
Nsup  = {}
MINsup = {}
MAXsup = {}

Stotal = SStotal = 0.0
Ntotal = 0
MINtotal = 0
MAXtotal = 0
### (1) Read file line one by one ###
for line in f:
  line = line.rstrip('\n')
  if (len(line)>10) and (line.startswith('#')==0):
    field = line.split('\t')
    for i in range(len(field)):
      if (field[i].find(':')>=0):
        items = field[i].split(':')
        field[i] = items 
    #print field
    #print line
    supplier_list = []
    if (isinstance(field[OPT['sc']],str)):
      supplier_list.append(field[OPT['sc']])
    if (isinstance(field[OPT['sc']],list)):
      for x in (field[OPT['sc']]):
        supplier_list.append(x)
    ### (2) Check focused propert ##
    if (isfloat(field[OPT['c']])==1):
      val = float(field[OPT['c']])
      if (Ntotal==0):
        MAXtotal = val
        MINtotal = val
      else:
        if (val > MAXtotal):
          MAXtotal = val
        if (val < MINtotal):
          MINtotal = val
      Ntotal += 1
      Stotal += val
      SStotal += val * val
     
      for x in (supplier_list):
        Nsup[x]  = Nsup.get(x,0) + 1
        Ssup[x]  = Ssup.get(x,0.0) + val 
        SSsup[x] = SSsup.get(x,0.0) + val*val
        if (MINsup.has_key(x)==0) or (val<MINsup[x]):
          MINsup[x] = val
        if (MAXsup.has_key(x)==0) or (val>MAXsup[x]):
          MAXsup[x] = val
f.close()


### OUTPUT ###
Mtotal  = Stotal/Ntotal
SDtotal = math.sqrt(SStotal/Ntotal - Mtotal * Mtotal)

if (OPT['of']!=''):
  print "#write_query_string_histogram -->'%s'"%(OPT['of'])
  of = open(OPT['of'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  supplier_list = sorted(Nsup.keys(),lambda y,x:cmp(Nsup[x],Nsup[y]))
  #of.write("#supplier\tNcount\tMean\tSD\tmin\tmax\n")
  #of.write("%s\t%s\t%f\t%f\t%f\t%f\n"%('TOTAL',Ntotal,Mtotal,SDtotal,MINtotal,MAXtotal))
  of.write("#supplier Ncount Mean SD min max\n")
  of.write("%-20s %6d %8.3f %8.3f %8.3f %8.3f\n"%('TOTAL',Ntotal,Mtotal,SDtotal,MINtotal,MAXtotal))
  for x in (supplier_list):
    M = Ssup[x]/Nsup[x]
    SD = math.sqrt(SSsup[x]/Nsup[x]- M*M)
    xx = x
    xx = xx.replace(' ','_')
    of.write("%-20s %8d %8.3f %8.3f %8.3f %8.3f\n"%(xx,Nsup[x],M,SD,MINsup[x],MAXsup[x]))
  of.close()

