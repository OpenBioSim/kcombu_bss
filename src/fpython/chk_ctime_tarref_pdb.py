#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Nov 6, 2013"


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

############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'

OPT['idtarp']  = 'tmpout'
OPT['ohist'] = 'tmp.hist'
OPT['bin'] = '0.1'

if (len(sys.argv)<2):
  print "chk_ctime_tarref_pdb.py <options>"
  print " for making summary of computational time for target_reference PDB files."
  print " written such as 'REMARK COMP_TIME  1.393790 seconds'."
  print "  coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -L      : lists of library molecules 'list1:list2:list3' [%s]"%(OPT['L'])
  print " -idtarp : input dirs for predicted target 3D model in PDB 'dir1:dir2:dir3' [%s]"%(OPT['idtarp'])
  print " -ohist  : output histogram file [%s]"%(OPT['ohist']) 
  print " -bin    : bin of histogram [%s]"%(OPT['bin']) 
  sys.exit(1)

### [1] read option ###
read_option(sys.argv,OPT)
OPT['bin'] = float(OPT['bin'])

L_lists      = OPT['L'].split(':')
idtarp_lists = OPT['idtarp'].split(':')

NL_lists      = len(L_lists)
Nidtarp_lists = len(idtarp_lists)

if (NL_lists != Nidtarp_lists):
  print "#ERROR: length of L(%d) and idtarp(%d) are not equal."%(NL_lists,Nidtarp_lists)
  sys.exit(1)

minT = 0.0
maxT = 0.0
sumT = 0.0
N = 0
Ncount = {}

for n in range(NL_lists):
  listfile = L_lists[n]
  idtarp   = idtarp_lists[n]

  ### [2] read liglist (-L) ###
  liglist = []
  read_list_file(listfile,liglist)

  ### [3] read models #####
  for ligA in (liglist):
    fieldA = ligA.split('_')
    ligname3A = fieldA[0]
    for ligB in (liglist):
      fieldB = ligB.split('_')
      ligname3B = fieldB[0]
      modelfile = "%s/%s_%s"%(idtarp,ligname3A,ligB)
      if (os.access(modelfile,os.R_OK)==0):
        modelfile = "%s/%s/%s_%s"%(idtarp,ligname3A,ligname3A,ligB)
        if (os.access(modelfile,os.R_OK)==0):
          print "#ERROR:Can't find modelfile '%s'"%(modelfile)
          sys.exit(1)

      timeline = pickup_focused_line(modelfile,"REMARK COMP_TIME")
      ## exmple of timeline:  'REMARK COMP_TIME  16.512306 seconds'
      fields = timeline.split()
      time = float(fields[2])
      sumT += time
      if (N==0) or (time<minT):
        minT = time
      if (N==0) or (time>maxT):
        maxT = time
      N += 1
      itime = int(time/OPT['bin'])
      Ncount[itime] = Ncount.get(itime,0) + 1
      #print "modelfile '%s' %s %f"%(modelfile,timeline,time)

aveT = sumT/N

print "#count %d ave_time %f min_time %f max_time %f"%(N,aveT,minT,maxT)
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


