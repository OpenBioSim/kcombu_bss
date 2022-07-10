#!/usr/bin/env python
#
# <keggdrug_add3Did_to_2D.py>
#

import sys
import math
import os
import re
import lbox_func
LastModDate = "Oct 12, 2011"


def read_IDNUMBER_from_mol2_file(ifname):
  f = open(ifname)
  comment_start = 0;
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('@<TRIPOS>COMMENT')):
      comment_start = 1
    elif (line.startswith('@<TRIPOS>')):
      comment_start = 0

    if (comment_start==1):
      if (line.startswith('    IDNUMBER = ')):
        field = line.split()
        f.close()
        return(field[2])
  f.close()

def read_contents_of_file(ifname):
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
  f = open(filename)
  contents = ''
  for line in f:
    contents = contents + line
  f.close() 
  return(contents)


def write_content2D_and_id3Ds(ofname,contents,id3Dlist):
  print "#write_content2D_and_id3Ds()-->'%s'"%(ofname)
  of = open(ofname,'w')
  of.write("%s\n"%(contents))
  of.write("\n")
  of.write("> <SUPPLIERNAME_1> (KEGG_DRUG)\n")
  of.write("KEGG_DRUG\n")
  of.write("\n")
  of.write("> <LIGANDBOX_ID>\n")
  for i in range(len(id3Dlist)):
    id = id3Dlist[i]
    if (i>0):
      of.write(" ")
    of.write("%s"%(id))
  of.write("\n")
  of.write("\n")
  of.close() 




###############
#### MAIN #####
###############

OPT = {}
OPT['id3D'] = '/DB/LigandBox/KEGG_DRUG_3D'
OPT['id2D'] = '/DB/kegg/medicus/drug/mol'
OPT['od2D'] = 'tmpout'
OPT['cp2D'] = 'T'
OPT['D'] = 'S'
OPT['A'] = 'F'

if (len(sys.argv)<2):
  print "keggdrug_add3Did_to_2D.py <options>"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -id3D :input directory of 3D file [%s]"%(OPT['id3D'])
  print " -id2D :input directory of 2D file [%s]"%(OPT['id2D'])
  print " -od2D :ouput directory of 2D file [%s]"%(OPT['od2D'])
  print " -cp2D :copy all the file in 'id2D' into the 'od2D'.  'T' or 'F' [%s]"%(OPT['cp2D']);
  print " -D   :output directory style. 'O'ne directory,'S'plited (such as '00020001_00025000/') [%s]"%(OPT['D'])
  print " -A   : Action (T or F) [%s]"%(OPT['A'])
  print "<procedure>"
  print " [1] read all the files in the directory 'id3D'" 
  print " [2] copy all the file in the directory 'id2D' into 'od2D'."
  print " [3] rewrite new 2D file in the directory 'od2D' with the list of id3D"
  sys.exit(1)

lbox_func.read_option(sys.argv,OPT)

### [1] read all the files in the directory 'id3D' ##
flist3D = os.listdir(OPT['id3D'])
id3D_from_id2D = {}
for id3Dfile in (flist3D):
  id3D = id3Dfile.replace('.mol2','')
  field = id3D.split('-')
  id2D = field[1]
  #id2D = read_IDNUMBER_from_mol2_file(OPT['id3D']+'/'+id3Dfile)
  
  print "'%s'-->'%s'"%(id3D,id2D)
  if (id3D_from_id2D.has_key(id2D)==0):
    id3D_from_id2D[id2D] = []
  id3D_from_id2D[id2D].append(id3D)


### [2] copy all the file in the directory 'id2D' into 'od2D'.
if (OPT['cp2D']=='T'):
  flist2D = os.listdir(OPT['id2D'])
  for id2Dfile in (flist2D):
    id2D = id2Dfile.replace('.mol','')

    if (OPT['D']=='S'):
      odir = OPT['od2D'] + '/' + lbox_func.subdir_from_ligand_id(id2D)
      if (os.access(odir,os.R_OK)==0):
         com = "mkdir %s"%(odir)
         print "#%s"%(com)
         if (OPT['A']=='T'):
           os.system(com)
    else:
      odir = OPT['od2D']

    new_id2Dfile = odir  + '/' + id2D + '.sdf'
    str = "cp %s/%s %s"%(OPT['id2D'],id2Dfile,new_id2Dfile) 
    print "#%s"%(str)
    if (OPT['A']=='T'):
      os.system(str)

### [3] write new 2D file in the directory 'od2D' with id3Dlist

for id2D in (id3D_from_id2D.keys()):
  print ">%s"%(id2D)
  content2D = read_contents_of_file(OPT['id2D']+'/'+id2D)
  #print content2D
  if (OPT['D']=='S'):
    odir = OPT['od2D'] + '/' + lbox_func.subdir_from_ligand_id(id2D)
  else:
    odir = OPT['od2D']
  ofname2D = odir  + '/' + id2D + '.sdf'
  id3Dlist = sorted(id3D_from_id2D[id2D],lambda x,y:cmp(x,y))
  for id3D in (id3Dlist):
    sys.stdout.write(" %s"%(id3D))
  sys.stdout.write('\n')
  if (OPT['A']=='T'):
    write_content2D_and_id3Ds(ofname2D,content2D,id3Dlist)

if (OPT['cp2D']=='T'):
  print "#Num_of_files_in id2D '%s' : %d"%(OPT['id2D'],len(flist2D))
print "#Num_of_files_with_id3D in od2D '%s' : %d"%(OPT['od2D'],len(id3D_from_id2D.keys()))

