#!/usr/bin/env python
#
# <split_mol2.py>
#

import sys
import math
import os
import re
import lbox_func

LastModDate = "Aug 11, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



def write_contents_of_molecule(OPT,contents,prop):
  if (OPT['id']=='N'):
    lig_name = prop['NAMIKI_ID']
  elif (OPT['id']=='L'):
    lig_name = prop['LIGANDBOX_ID']
  else: 
    lig_name = prop['LIGANDBOX_ID']
  if (OPT['D']=='O'):
    odir = OPT['od'] 
  if (OPT['D']=='S'):
    odir = OPT['od'] + '/' + lbox_func.subdir_from_ligand_id(lig_name)
  ofname = odir + '/' + lig_name
  if (OPT['add']=='T') and (lig_name.endswith('.mol2')==0):
    ofname = ofname + '.mol2'
  print("'%s' ofname '%s' len %d %s %s"%(lig_name,ofname,len(contents),prop.get('LIGANDBOX_ID',''),prop.get('NAMIKI_ID','')))
  if (os.access(odir,os.R_OK)==0):
    com = "mkdir %s"%(odir)
    print "#%s"%(com)
    if (OPT['A']=='T'):
      os.system(com)   
  if (OPT['A']=='T'):   
    fo = open(ofname,"w")
    fo.write("%s"%(contents))   
    fo.close()
  pass

def write_id3D_at_the_bottom_of_2D_data(OPT,id2D,id3Dlist):
  fname2D = OPT['d2D'] + '/' + lbox_func.subdir_from_ligand_id(id2D) + '/' + id2D
  filename = fname2D
  if (os.access(filename,os.R_OK)==0):
    filename = filename + '.sdf'
    if (os.access(filename,os.R_OK)==0):
      filename = filename + '.mol'
      if (os.access(fname2D,os.R_OK)==0):
        print "#WARNING:Can't open %s + '.sdf' or '.mol'"%(fname2D)
        return(0)

  print "# write_id3D_at_the_bottom_of_2D_data()--> '%s'"%(filename)
  ## (1) read contents of file ##
  fi = open(filename) 
  lines = []
  for line in fi:
    lines.append(line) 
  fi.close()

  ## (2) check the line '>  <LIGANDBOX_ID> ##
#> <LIGANDBOX_ID>
#D00003210-01 D00003210-02 D00003210-03 D00003210-04
  Nline_lbox = -1
  for i in range(len(lines)):
    line = lines[i]
    if (line.startswith('> <LIGANDBOX_ID>')):
      Nline_lbox = i 

  ## (3) Add id3D into the next line '> <LIGANDBOX_ID>'
  if (Nline_lbox>=0):
    id_dic = {} 
    dataline = lines[Nline_lbox+1]
    dataline = dataline.rstrip('\n') 
    dataline = dataline.rstrip('\r') 
    fields = dataline.split(' ')
    for field in (fields):
       id_dic[field] = 1
    for id3D in (id3Dlist):
       id_dic[id3D] = 1
    lines[Nline_lbox+1] = ''
    new_id3Dlist = sorted(id_dic.keys(),lambda x,y:cmp(x,y))
    for i in range(len(new_id3Dlist)):
      if (i>0):
        lines[Nline_lbox+1] += " "
      lines[Nline_lbox+1] += "%s"%(new_id3Dlist[i]) 
    lines[Nline_lbox+1] += "\n"
  else:
    lines.append('\n')   
    lines.append('> <LIGANDBOX_ID>\n')
    str = ''
    new_id3Dlist = sorted(id3Dlist,lambda x,y:cmp(x,y))
    for i in range(len(new_id3Dlist)):
      if (i>0):
        str += ' '
      str += new_id3Dlist[i]
    lines.append(str+'\n')   
    lines.append('\n')   

  ## (4) Write again ##
  fo = open(filename,'w') 
  for line in (lines):
    fo.write("%s"%(line))
  fo.close()
  #for line in (lines):
  #  sys.stdout.write("%s"%(line))
  pass






###############
#### MAIN #####
###############

OPT = {}
OPT['od'] = 'tmpout'
OPT['id'] = 'L'
OPT['A'] = 'F'
OPT['D'] = 'S'
OPT['d2D'] = ''
OPT['add'] = 'T'
OPT['ol'] = ''

if (len(sys.argv)<3):
  print "split_mol2.py [mol2_file]"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -D   :output directory style. 'O'ne directory,'S'plited (such as '00020001_00025000/') [%s]"%(OPT['D'])
  print " -od  :output directory [%s]"%(OPT['od'])
  print " -id  :location of ID.'N':NAMIKI_ID,'L':LIGANDBOX_ID [%s]"%(OPT['id'])
  print " -d2D :directory of 2D molecule (namiki2D) [%s]"%(OPT['d2D'])  
  print " -add :add '*.mol2' into the output file (T or F) [%s]"%(OPT['add'])
  print " -ol  :output list of id2D vs id3D [%s]"%(OPT['ol']) 
  print " -A   :Action. 'T' or 'F' [%s]"%(OPT['A']);
  sys.exit(1)

read_option(sys.argv,OPT)

fname = sys.argv[1]

fi = open(fname)


contents = ''
Nlig_hash = {}
Nline = 0
lig_name = ''
status = ''
ID2D_to_3D = {} 
for line in fi:
  dline = line.rstrip('\n')
  dline = dline.rstrip('\r')

  if (line.startswith("@")):
    Nline = 0
  else:
    Nline += 1

  if (line.startswith("@<TRIPOS>COMMENT")):
    ### output contents ###
    if (contents != ''):
      write_contents_of_molecule(OPT,contents,prop)
      id2D = prop.get('NAMIKI_ID','')
      id3D = prop.get('LIGANDBOX_ID','')
      if (id3D=='00000010-02'):
         print contents 
      if (id2D != '') and (id3D != ''):
        if (ID2D_to_3D.has_key(id2D)==0):
          ID2D_to_3D[id2D] = []
        ID2D_to_3D[id2D].append(id3D)
    contents = line 
    lig_name = ''
    prop = {}

  else:
    contents = contents + line 

  if (status=='MOLECULE') and (line.startswith('@')==0):
    if (Nline==1):
      lig_name = dline.replace(' ','')
  if (status=='COMMENT') and (line.startswith('@')==0):
    line = line.rstrip('\n')
    line = line.rstrip(' ')
    line = line.lstrip(' ')
    line = line.lstrip('\n')
    field = line.split(' = ')
    if (len(field)==2):
      head  = field[0]
      value = field[1]
      #print field,head,value
      prop[head] = value
  if (line.startswith("@<TRIPOS>MOLECULE")):
    status = 'MOLECULE'
  elif (line.startswith("@<TRIPOS>COMMENT")):
    status = 'COMMENT'
  elif (line.startswith("@")):
    status = ''

fi.close()

if (len(contents)>10):
  ### output contents ###
  write_contents_of_molecule(OPT,contents,prop)
  id2D = prop.get('NAMIKI_ID','')
  id3D = prop.get('LIGANDBOX_ID','')
  if (id2D != '') and (id3D != ''):
    if (ID2D_to_3D.has_key(id2D)==0):
      ID2D_to_3D[id2D] = []
    ID2D_to_3D[id2D].append(id3D)
    print "end add %s"%(id3D)

if (OPT['ol'] != ''):
  of = open(OPT['ol'],'w')

for id2D in (ID2D_to_3D.keys()):

  print id2D
  for id3D in (ID2D_to_3D[id2D]):
    sys.stdout.write("%s "%(id3D)) 
  sys.stdout.write("\n")

  if (OPT['ol'] != ''):
    of.write("%s"%(id2D))
    for id3D in (ID2D_to_3D[id2D]):
      of.write(" %s"%(id3D)) 
    of.write('\n')

  write_id3D_at_the_bottom_of_2D_data(OPT,id2D,ID2D_to_3D[id2D])


if (OPT['ol'] != ''):
  of.close()
