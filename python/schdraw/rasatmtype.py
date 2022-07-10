#!/usr/bin/env python

import sys
import os
from datetime import datetime


LastModDate = "Jan 26, 2012"


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



def read_pdb_file(ifname,atomnum_list,element_list, atomtype_list):
  print "#read_pdb_file('%s')"%(ifname)
  if (os.access(ifname,os.R_OK)==0):  
    print "#ERROR:Can't open pdbfile '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
#          1         2         3         4         5         6         7      
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#MODEL        1
#HETATM    1  C   MOL A   1      -6.968   0.820   0.000  0.00  0.00           C   R
#HETATM    2  C   MOL A   1      -8.400   1.268   0.000  0.00  0.00           C   L
#MODEL        1
#HETATM    1  C   MOL A   1      -6.968   0.820   0.000  0.00  0.00           C  C@
#HETATM    2  C   MOL A   1      -8.400   1.268   0.000  0.00  0.00           C  C@
#

  for line in f:
    if (line.startswith('HETATM') or line.startswith('ATOM')):
      line = line.rstrip('\n')
      anumstr = line[6:11]
      anum = int(line[6:11])
      tailstr = line[66:] 
      field = tailstr.split()
      ele = ' '
      type = ' ' 
      if (len(field)>0):
        ele  = field[0]
      if (len(field)>1):
        type  = field[1]
 
      atomnum_list.append(anum)
      element_list.append(ele)
      atomtype_list.append(type)
      pass 
  pass


def set_atomtype_color(type,color):
  if (type=='D'):
    color['D'] = 'blue'
    color['A'] = 'red'
    color['B'] = 'purple'
    color['R'] = '[0,128,0]'
    color['L'] = 'gray'



##########
## MAIN ##
##########

OPT = {}
OPT['of'] = '-'
OPT['at'] = 'D'
if (len(sys.argv)<2):
  print "rasatmtype.py [pdbfile] <options>"
  print " for making rasmol script file for atomtype coloring a PDB file."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>" 
  print " -of : output file [%s]"%(OPT['of'])
  print " -at : atom type coling 'D'abrl [%s]"%(OPT['at'])
  sys.exit(1)

read_option(sys.argv,OPT)

ipdbfile = sys.argv[1]

atomtype_color = {}
set_atomtype_color(OPT['at'],atomtype_color)

## [1] read PDB file ##
atomnum_list  = []
element_list  = []
atomtype_list = []

read_pdb_file(ipdbfile,atomnum_list,element_list, atomtype_list)


## [2] Output Rasmol script ##

if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_rasmol_script() --> '%s'"%(OPT['of'])

of.write("#COMMAND %s\n"%(OPT['COMMAND']))

of.write("echo \">>DABRL coloring <<\"\n")
of.write("echo \"'elemmen(CNO)+DABRL style'\"\n")
of.write("echo \"Donor    =>blue,  Accepror=>red,       Both=>green\"\n")
of.write("echo \"aLyphatic=>gray, aRomatic=>dark green.\"\n")

for i in range(len(atomnum_list)):
  color = atomtype_color.get(atomtype_list[i],'white')
  of.write("#%d '%s' '%s'\n"%(atomnum_list[i],element_list[i],atomtype_list[i]))
  of.write("select atomno==%d\n"%(atomnum_list[i]))
  of.write("color %s\n"%(color))

of.write("connect true\n")
of.write("set background white\n")
of.write("select all\n")
of.write("label %a\n")
of.write("wireframe 30\n")
of.write("select !hydrogen\n")
of.write("spacefill 100\n")
if (OPT['of']!='-'):
  of.close()
