#!/usr/bin/env python
import sys
import os


LastModDate = "June 18,2012"


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
     list.append(fields[1])
  f.close()


if (len(sys.argv)<2):
  print "get_uniseq.py [uniid_list_file] [sequence file]"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  sys.exit(1)

ilistfile = sys.argv[1]
iseqfile = sys.argv[2]

idlist = []
read_list_file(ilistfile,idlist)
iddic = {}
for id in (idlist):
  iddic[id] = iddic.get(id,0) + 1
  #print id,iddic[id]


f = open(iseqfile)
out = 0
for line in f:
  if (line.startswith('>')):
    field = line.split()
    entry = field[0][1:]
    (head,tail) = entry.split('_')
    if (entry.endswith('HUMAN1')):
      id = head + '_1'
    elif (entry.endswith('HUMAN2')):
      id = head + '_2'
    else:
      id = head 
    #print entry,id
    if (iddic.has_key(id)):
      out = 1
      line = '>' + id + ' ' + line[1:]
    else:
      out = 0
  if (out==1):
    sys.stdout.write("%s"%(line))

