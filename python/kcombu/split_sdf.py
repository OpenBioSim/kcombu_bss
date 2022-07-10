#!/usr/bin/env python
##
## <split_sdf.py>
##
## for splitting a multi-compound SDF file into many one-compound SDF files
##
## <FileFormat>
##
## (1)Separator line between compounds
##
##  $$$$
##
## (2) Compound name  
##   i) For most of the case, the compound name is the first line of the contents.
##      such as: 
## $$$$
## 53225002                 <-- this is the compound name
##   -OEChem-07041100532D
## 
##  28 29  0     0  0  0  0  0  0999 V2000
##
##   ii) For the ZINC database, the first line is "blank'. The example is as follows:
##    
## M  END
## > <id>
## ZINC03814459
## 
## > <Cluster>
## 2
## 
## $$$$
## 
##   SciTegic07110717182D
##



import sys
import re
import math
import os

LastModDate = "July 13, 2011"

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



###############
#### MAIN #####
###############

OPT = {}
OPT['od'] = 'tmpout'
OPT['ol'] = 'outlist'
OPT['cl'] = ''
OPT['id'] = 'F'
OPT['add'] = 'F'
OPT['A'] = 'F'
OPT['D'] = 'O'

if (len(sys.argv)<3):
  print "split_sdf.py [sdf_file]"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -D  :output directory style. 'O'ne directory,'S'plited (such as '00020001_00025000/') [%s]"%(OPT['D'])
  print " -od :output directory [%s]"%(OPT['od']) 
  print " -ol :output list file [%s]"%(OPT['ol'])
  print " -id :location of ID. 'F'irst_line,'I'> <id> tag(ZINC),'N'> <NAMIKI_ID> [%s]"%(OPT['id']) 
  print " -cl :name of class for '-ol' list file [%s]"%(OPT['cl']) 
  print " -add:add '*.sdf' into the output file (T or F)[%s]"%(OPT['add']) 
  print " -A  :Action. 'T' or 'F' [%s]"%(OPT['A']);  
  sys.exit(1)

read_option(sys.argv,OPT)

ifname = sys.argv[1]
if not os.access(ifname,os.R_OK):
  print "#ERROR:Can't open sdffile '%s'"%(ifname)
  sys.exit(1)

if (ifname.endswith('.gz')):
  command = "zcat %s"%(ifname)
  fi = os.popen(command)
  print "#os.popen('%s')"%(command)
else:
  fi = open(ifname)

if (OPT['ol'] != ''):
  fol = open(OPT['ol'],'w')
  print "#write_compound_list --> '%s'"%(OPT['ol'])


contents = ''
lig_hash = {}
idread = 0
lig_name = ''
Nout = 0
for line in fi:
  #sys.stdout.write("%s"%(line))
  dline = line.rstrip('\n')
  dline = dline.rstrip('\r')
  if (line.startswith("$$$$")):
    if (contents != ''):
      if (lig_name == ''):
        print "#ERROR:Compound name is not defined."
 #       print contents 
        sys.exit(1)
      if (lig_hash.has_key(lig_name)):
        print "#ERROR:Compound name '%s' is already used."%(lig_name)
        sys.exit(1)
      if (OPT['D']=='O'):
        odir = OPT['od']
      if (OPT['D']=='S'):
        odir = OPT['od'] + '/' + dirname_from_ligand_id(lig_name)
      ofname = odir + '/' + lig_name

      if (OPT['add']=='T') and (lig_name.endswith('.sdf')==0):
        ofname = ofname + '.sdf'
      print("'%s' len %d -->'%s'"%(lig_name,len(contents),ofname))

      if (os.access(odir,os.R_OK)==0):
        com = "mkdir %s"%(odir)
        print "#%s"%(com)
        if (OPT['A']=='T'):
          os.system(com)

      if (OPT['A']=='T'):
        fo = open(ofname,"w")
        fo.write("%s"%(contents))   
        fo.close()
      Nout += 1
      if (OPT['ol'] != ''):
        fol.write("%-20s %s\n"%(lig_name,OPT['cl']))
      lig_hash[lig_name] = 1
      lig_name = ''
    contents = ''

  else:
    if (OPT['id']=='F'):
      if (contents == '') and (len(dline)>0):
        lig_name = dline 

    if (OPT['id']=='I'):
      if (idread==1):
        lig_name = dline
      if (line.startswith('> <id>')):
        idread = 1
      else:
        idread = 0 

    if (OPT['id']=='N'):
      if (idread==1):
        lig_name = dline
      if (line.startswith('> <NAMIKI_ID>')):
        idread = 1
      else:
        idread = 0 


    contents = contents + line 

fi.close()
if (OPT['ol'] != ''):
  fol.close()

print "#Nout %d"%(Nout)


