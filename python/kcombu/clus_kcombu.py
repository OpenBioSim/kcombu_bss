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

def mark_neighbor(x,cluster_num,cluster_member):
  if (mark[x]==-1):
    mark[x]  = cluster_num
    cluster_member.append(x)
    for y in (NEIGHBOR[x]):
      mark_neighbor(y,cluster_num,cluster_member)



##############
#### MAIN ####
##############

OPT = {}
OPT['L'] = '.'
OPT['iav'] = '.'
OPT['om'] = ''
OPT['or'] = 'rep.out'

if (len(sys.argv)<2):
  print "clus_kcombu.py [list_of_query_compound_file] <options>"
  print " for single linkage clustering from kcombu all-vs-all results."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print "-L   : list of library molecules[%s]"%(OPT['L'])
  print "-iav : input file for similar pairs in all-vs-all comparison. [%s]"%(OPT['iav'])
  print "-om  : output file for members of clusters [%s]"%(OPT['om'])
  print "-or  : output file for representative members of clusters [%s]"%(OPT['or'])
  sys.exit(1)


### (1) read option ###
read_option(sys.argv,OPT)
iqlistfile = sys.argv[1]

### (2) read list of compounds ###
complist = []
read_list_file(OPT['L'],complist)

Nlist = len(complist);

NEIGHBOR = {}
mark = {}
for x in (complist):
  NEIGHBOR[x] = [] 
  mark[x] = -1

### (3) read similar pairs file ###
of = open(OPT['iav'])
##[molA] [molB] [tanimoto]
#208264 208262 1.000000
#102529 13904 1.000000
for line in of:
  if (line.startswith('#')==0) and (len(line)>5):
    line = line.rstrip('\n')
    (A,B,tanimoto) = line.split()
    #print "'%s' '%s' %f"%(A,B,float(tanimoto))
    NEIGHBOR[A].append(B)
    NEIGHBOR[B].append(A)

#### (4) Do single-linkage clustering ####
cluster_num = 0
CLUSTER_MEMBER = [] 
for x in (complist):
  if (mark[x]==-1):
    member = []
    mark_neighbor(x,cluster_num,member)
    CLUSTER_MEMBER.append(member)
    cluster_num += 1

#for x in (complist):
#  print "%s mark %d"%(x,mark[x])

#### (5) Output ####
if (OPT['om'] != ''):
  print "#write_members_of_clusters -->'%s'"%(OPT['om'])
  of = open(OPT['om'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  for c in range (len(CLUSTER_MEMBER)):
    of.write(">%d Nmember %d\n"%(c,len(CLUSTER_MEMBER[c])))
    for x in (CLUSTER_MEMBER[c]):
      of.write("%s\n"%(x))
  of.close()


if (OPT['or'] != ''):
  print "#write_representatives_of_clusters -->'%s'"%(OPT['or'])
  of = open(OPT['or'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#[representative] [cluster_number] [Number of members]\n")
  for c in range (len(CLUSTER_MEMBER)):
    of.write("%s %d %d\n"%(CLUSTER_MEMBER[c][0],c,len(CLUSTER_MEMBER[c])))
  of.close()
