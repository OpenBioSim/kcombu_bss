#!/usr/bin/env python
import sys
import os
import math

LastModDate = "Mar 24, 2010"

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
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     list.append(fields[0])
  f.close()


##############
#### MAIN ####
##############

OPT = {}
OPT['A'] = 'F'
OPT['L'] = ''
OPT['QL'] = ''
OPT['idl'] = '.'
OPT['idq'] = '.'
OPT['ods'] = '.'
OPT['tam'] = '0.9' 
OPT['div'] = '0/1'
OPT['fQ'] = '-'
OPT['fL'] = '-'
OPT['con'] = 'P'

if (len(sys.argv)<2):
  print "scan_kcombu.py <options>"
  print " for query-library search using kcombu program"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L  : list of library molecules[%s]"%(OPT['L'])
  print " -QL : list of query molecules[%s]"%(OPT['QL'])
  print " -idl: input directory for library molecules[%s]"%(OPT['idl'])
  print " -idq: input directory for query molecules[%s]"%(OPT['idq'])
  print " -con: Type of connectivity of graph. 'C'onnected, 'U'nconnected, 'P'rogressive_unconnected [%s]"%(OPT['con'])
  print " -tam: threshold for Tanimoto_molecula_formula filtering [%s]"%(OPT['tam'])
  print " -ods: output directory for search results[%s]"%(OPT['ods'])
  print " -fQ : file format of query   'S'df or 'P'db [%s]"%(OPT['fQ'])
  print " -fL : file format of library 'S'df or 'P'db [%s]"%(OPT['fL'])
  print " -div: division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -A  : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)


### (1) read option ###
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)
print "bunshi %d bunbo %s"%(bunshi,bunbo)

### (2) read list of query compounds ###

query_list = []
read_list_file(OPT['QL'],query_list)

Nlist = len(query_list);
Nstart = bunshi*int(Nlist/bunbo);
Nend   = (bunshi+1)*int(Nlist/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Nlist
print "#Nlist %d bunshi/bunbo %d/%d start %d end %d"%(Nlist,bunshi,bunbo,Nstart,Nend)


### (3) do kcombu ###
for i in range(Nstart,Nend):
  q = query_list[i]
  print q
  command = "kcombu -M L -L %s -idl %s -con %s -tam %s -Q %s/%s -fQ %s -fL %s -osc %s/%s"%(OPT['L'],OPT['idl'],OPT['con'],OPT['tam'],OPT['idq'],q,OPT['fQ'],OPT['fL'],OPT['ods'],q)
  if (OPT['A'] == 'T'):
    os.system(command)
  else:
    print "#%s"%(command)

print "#Nlist %d bunshi/bunbo %d/%d start %d end %d"%(Nlist,bunshi,bunbo,Nstart,Nend)

