#!/usr/bin/env python
#
# <split_mol2.py>
#

import sys
import math
import os
import re

LastModDate = "July 14, 2011"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def dirname_from_ligand_id(lig_id):
#>> example of 'lig_id' <<
# '00022356-01'      --> '00020001_00025000'
# 'NS-00022175'      --> '00020001_00025000'
# 'NS-00022175.mol2' --> '00020001_00025000'

  field  = re.split('[\-\.]',lig_id)

  id = ''
  for x in (field):
    #print x
    if (len(x)>=8) and (x.isdigit()):
      id = x
  if (id==''):
    return('null')
  num = int(id)
  sta = int(5000*(math.ceil(float(num)/5000.0)-1))+1
  end = int(5000*math.ceil(float(num)/5000.0))
  stastr = '%d'%(sta)
  stastr = '0'*(8-len(stastr)) + stastr
  endstr = '%d'%(end)
  endstr = '0'*(8-len(endstr)) + endstr
  #print "num %8d stastr %s endstr %s"%(num,stastr,endstr)
  str = "%s_%s"%(stastr,endstr)
  return(str) 


def write_contents_of_molecule(OPT,contents,prop,fol):
  if (OPT['id']=='N'):
    lig_name = prop['NAMIKI_ID']
  elif (OPT['id']=='L'):
    lig_name = prop['LIGANDBOX_ID']
  else: 
    lig_name = prop['LIGANDBOX_ID']
  if (OPT['D']=='O'):
    odir = OPT['od'] 
  if (OPT['D']=='S'):
    odir = OPT['od'] + '/' + dirname_from_ligand_id(lig_name)
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
    if (OPT['ol'] != ''):
      fol.write("%s %s %s\n"%(lig_name,OPT['cl'],prop.get('NAMIKI_ID','')))
  pass



###############
#### MAIN #####
###############

OPT = {}
OPT['od'] = 'tmpout'
OPT['ol'] = 'outlist'
OPT['cl'] = ''
OPT['id'] = 'L'
OPT['add'] = 'F'
OPT['A'] = 'F'
OPT['D'] = 'O'

if (len(sys.argv)<3):
  print "split_mol2.py [mol2_file]"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -D  :output directory style. 'O'ne directory,'S'plited (such as '00020001_00025000/') [%s]"%(OPT['D'])
  print " -od :output directory [%s]"%(OPT['od'])
  print " -ol :output list file [%s]"%(OPT['ol'])
  print " -id :location of ID.'N':NAMIKI_ID,'L':LIGANDBOX_ID [%s]"%(OPT['id'])
  print " -cl :name of class for '-ol' list file [%s]"%(OPT['cl'])
  print " -add:add '*.mol2' into the output file (T or F)[%s]"%(OPT['add'])
  print " -A  :Action. 'T' or 'F' [%s]"%(OPT['A']);
  sys.exit(1)

read_option(sys.argv,OPT)

fname = sys.argv[1]

fi = open(fname)
if (OPT['ol'] != ''):
  fol = open(OPT['ol'],'w')
  print "#write_compound_list --> '%s'"%(OPT['ol'])


contents = ''
Nlig_hash = {}
Nline = 0
lig_name = ''
status = ''
for line in fi:
  dline = line.rstrip('\n')
  dline = dline.rstrip('\r')

  if (line.startswith("@")):
    Nline = 0
  else:
    Nline += 1

  if (line.startswith("@<TRIPOS>COMMENT")):
    if (contents != ''):
      write_contents_of_molecule(OPT,contents,prop,fol)
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

if (contents != ''):
  write_contents_of_molecule(OPT,contents,prop,fol)

if (OPT['ol'] != ''):
  fol.close()


