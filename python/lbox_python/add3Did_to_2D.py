#!/usr/bin/env python
#
# <add3Did_to_2D.py>
#

import sys
import math
import os
import re
import lbox_func
LastModDate = "Aug 11, 2011"


def Append_IDs_at_the_tail_of_file(ifname,idlist,Action):
  filename = ifname 
  if (os.access(filename,os.R_OK)==0):
    filename = ifname + '.sdf'
    if (os.access(filename,os.R_OK)==0):
      filename = ifname + '.mol'
      if (os.access(filename,os.R_OK)==0):
        filename = ifname + '.mol2'
        if (os.access(filename,os.R_OK)==0):
          print "#WARNING:Can't find '%s' (.sdf|.mol|.mol2)"%(ifname)
          return(0) 
  if (Action=='T'):
    print "Append() -->'%s'"%(filename)
    of = open(filename,'a')
    of.write("\n")
    of.write("> <LIGANDBOX_ID>\n")
    for i in range(len(idlist)):
      id = idlist[i]
      if (i>0):
        of.write(" ")
      of.write("%s"%(id))
    of.write("\n")
    of.write("\n")
    of.close() 
  else:
    print "-->'%s'"%(filename)
    of = sys.stdout
    of.write("> <LIGANDBOX_ID>\n")
    for i in range(len(idlist)):
      id = idlist[i]
      if (i>0):
        of.write(" ")
      of.write("%s"%(id))
    of.write("\n")
  return(1)


###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['D'] = 'O'
OPT['d2D'] = 'tmpout'

if (len(sys.argv)<2):
  print "add3Did_to_2D.py [2Did - 3Dids file]"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -D   :output directory style. 'O'ne directory,'S'plited (such as '00020001_00025000/') [%s]"%(OPT['D'])
  print " -d2D :directory of 2D file [%s]"%(OPT['d2D'])
  print " -A   :Action. 'T' or 'F' [%s]"%(OPT['A']);
  sys.exit(1)

lbox_func.read_option(sys.argv,OPT)

if2Dto3Dfile = sys.argv[1]

f = open(if2Dto3Dfile)
for line in f:
 # print line
  line = line.rstrip('\n')
  field = line.split()
  id2D = field[0]
  id3Dlist = []
  for i in range(1,len(field)):
    id3Dlist.append(field[i])
  subdir = lbox_func.subdir_from_ligand_id(id2D)
  print "%s subdir %s"%(id2D,subdir)
  filename = OPT['d2D'] + '/' + subdir + '/' + id2D 
  Append_IDs_at_the_tail_of_file(filename,id3Dlist,OPT['A'])
f.close()


