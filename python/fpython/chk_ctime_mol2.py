#!/usr/bin/env python
import sys
import os
from datetime import datetime
import time
from stat import *
LastModDate = "Oct 4, 2013"


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


def pickup_focused_line(ifname,pattern_head):
  if (os.access(ifname,os.R_OK)==0):
    return('')
  f = open(ifname)
  for line in f:
    if (line.startswith(pattern_head)):
      line = line.rstrip('\n')
      return(line)
  f.close()
  #print "#WARNING:no line for pattern_head '%s' in the file '%s'"%(pattern_head,ifname)
  return('')


def pickup_and_count_focused_line(ifname,pickup_head,count_head):
  Ncount = 0
  pickup_line = '' 
  if (os.access(ifname,os.R_OK)==0):
    return(Ncount)
  f = open(ifname)
  for line in f:
    if (line.startswith(count_head)):
      Ncount += 1
    if (line.startswith(pickup_head)):
      line = line.rstrip('\n')
      pickup_line = line 
  f.close()
  return((pickup_line,Ncount))




############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'

OPT['idmol2']  = 'tmpout'
OPT['itail'] = '.mol2'
OPT['ohist'] = 'tmp.hist'
OPT['bin'] = '0.1'
OPT['O'] = 'F'

if (len(sys.argv)<2):
  print "chk_ctime_mol2.py <options>"
  print " for making summary of computational time in mol2 file,written such as:"
  print "@<TRIPOS>COMMENT"
  print "  COMP_TIME = 8.373514 seconds"
  print "  coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -idmol2 : input directories for MOL2 files 'idir1:idir2:idir3' [%s]"%(OPT['idmol2']) 
  print " -itail  : tail of input files [%s]"%(OPT['itail']) 
  print " -ohist  : output histogram file [%s]"%(OPT['ohist']) 
  print " -bin    : bin of histogram [%s]"%(OPT['bin']) 
  print " -O      : output COMP_TIME and Nconf for each file (T or F)[%s]"%(OPT['O'])
  sys.exit(1)


### [1] read option ###
read_option(sys.argv,OPT)
OPT['bin'] = float(OPT['bin'])

idmol2_list = OPT['idmol2'].split(':')


sec_last_mod = {}
minT = 0.0
maxT = 0.0
aveT = 0.0
N = 0
Ncount = {}
minNconf = 0
maxNconf = 0
aveNconf = 0

for idmol2 in (idmol2_list):
  
  filelist = os.listdir(idmol2)
  
  for file in (filelist):
    if (file.endswith(OPT['itail'])):
      file_full = idmol2 + '/' + file
#      ctime_line = pickup_focused_line(file_full,"  COMP_TIME = ")
      (ctime_line,Nconf) = pickup_and_count_focused_line(file_full,"  COMP_TIME = ","@<TRIPOS>MOLECULE")
      fields = ctime_line.split()
      if (ctime_line != '') and (len(fields)>=2):
        ctime = float(fields[2])
        if (OPT['O']=='T'):
          print file_full,ctime_line,ctime,"Nconf:",Nconf
        aveT += ctime
        if (N==0) or (ctime > maxT): 
          maxT = ctime
        if (N==0) or (ctime < minT): 
          minT = ctime
    
        aveNconf += Nconf
        if (N==0) or (Nconf > maxNconf): 
          maxNconf = Nconf 
        if (N==0) or (Nconf < minNconf): 
          minNconf = Nconf 

        N += 1
        itime = int(ctime/OPT['bin'])
        Ncount[itime] = Ncount.get(itime,0) + 1

aveT     = aveT/N
aveNconf = float(aveNconf)/N
print "#count %d ave_time %f min_time %f max_time %f aveNconf %f minNconf %d maxNconf %d"%(N,aveT,minT,maxT,aveNconf,minNconf,maxNconf)

#### [4] output histogram ###
if (OPT['ohist'] != ''):
  of = open(OPT['ohist'],'w')
  print "#write_histogram() --> '%s'"%(OPT['ohist'])
  of.write("#HISTOGRAM OF COMPUTATIONAL TIME\n")
  of.write("#COMMAND '%s'\n"%(OPT['COMMAND']))
  of.write("#DATE    '%s'\n"%(OPT['START_DATE']))
  of.write("#count %d ave_time %f min_time %f max_time %f\n"%(N,aveT,minT,maxT))
  minI = int(minT/OPT['bin'])
  maxI = int(maxT/OPT['bin'])
  for i in range(minI,maxI):
    of.write("%f %f %d\n"%(i*OPT['bin'],Ncount.get(i,0.0)/float(N),Ncount.get(i,0)))
  of.close()  


