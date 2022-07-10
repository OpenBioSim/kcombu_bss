#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Jan 12, 2012"

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

def Action(comstr,action):
  if (action=='T'):
    print "#'%s'"%(comstr)
    os.system(comstr)  
  else:
    print "#'%s'"%(comstr)



def write_renumbered_pdb_file(ifname,ofname):
  print "#write_renumbered_pdb_file('%s') --> '%s'"%(ifname,ofname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open input pdbfile '%s'" %(ifname)
    return(0)
  f  = open(ifname)
  of = open(ofname,"w")
 
#          1         2         3         4 
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#HETATM 8845  PG  ATP A1294     -16.660 200.221 120.089  1.00  0.00           P
#HETATM 8846  O1G ATP A1294     -18.104 200.605 119.813  1.00  0.00           O
#HETATM 8847  O2G ATP A1294     -15.715 201.401 120.222  1.00  0.00           O
#HETATM 8848  O3G ATP A1294     -16.512 199.227 121.224  1.00  0.00           O
#HETATM 8849  PB  ATP A1294     -14.888 199.775 117.852  1.00  0.00           P
#HETATM 8850  O1B ATP A1294     -14.756 198.706 116.794  1.00  0.00           O
#HETATM 8851  O2B ATP A1294     -13.705 200.003 118.772  1.00  0.00           O
#HETATM 8852  O3B ATP A1294     -16.178 199.392 118.770  1.00  0.00           O
 
  Natom = 0
  for line in f:
    if (line.startswith('HETATM') or line.startswith('ATOM ')):
      Natom += 1
      head = line[0:6]
      tail = line[11:]
      newline = "%s%5d%s"%(head,Natom,tail)
      of.write("%s"%(newline))
  f.close()
  of.close()


############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['idpdb']  = 'SupLIG'
OPT['odmol2'] = 'tmpout'
OPT['gen3d'] = 'F'
OPT['short'] = 'F'
if (len(sys.argv)<2):
  print "mkMOL2frmPDB.py <options>"
  print " for making 3D mol2 file from 3D PDB file using OpenBabel"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L        : list of library molecules[%s]"%(OPT['L'])
  print " -idpdb    : input directory for 3D PDB molecules[%s]"%(OPT['idpdb'])
  print " -odmol2   : output directory for 3D MOL2 molecules[%s]"%(OPT['odmol2'])
  print " -short    : If 'T':shorten mol2 name for output, ('SIA_1mweA'-->'SIA')[%s]"%(OPT['short'])
  print " -gen3d    : add '--gen3d' option (T or F)[%s]"%(OPT['gen3d'])
  print " -div      : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -A        : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)


### [1] read option ###
PID = os.getpid()
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

### [2] read liglist (-L) ###
liglist = []
read_list_file(OPT['L'],liglist)
print "#PID %d"%(PID)

Nlig = len(liglist)
Nsta = Nlig*bunshi/bunbo
Nend = len(liglist)*(bunshi+1)/bunbo
print "#bunshi %d bunbo %d Npair %d Nsta %d Nend %d"%(bunshi,bunbo,Nlig,Nsta,Nend)
if (Nend>Nlig):
  Nend = Nlig

for i in range(Nsta,Nend):
  lig = liglist[i]
  (ligname3,pdbch) = lig.split('_')
  print "%s '%s' '%s'"%(lig,ligname3,pdbch)
  ipdbfile   = "%s/%s"%(OPT['idpdb'],lig)
  if (OPT['short']=='T'):
    omol2file  = "%s/%s.mol2"%(OPT['odmol2'],ligname3)
  else:
    omol2file  = "%s/%s.mol2"%(OPT['odmol2'],lig)
  tmppdbfile = lig 

  if (OPT['A']=='T'):
    write_renumbered_pdb_file(ipdbfile,tmppdbfile)
  else:
    print "#write_renumbered_pdb_file('%s') --> '%s'"%(ipdbfile,tmppdbfile)

  comstr = "babel -h "
  if (OPT['gen3d']=='T'):
    comstr += " --gen3d "

  comstr += " -ipdb %s -omol2 %s"%(tmppdbfile,omol2file)
  Action(comstr,OPT['A'])
  
  comstr = "rm %s"%(tmppdbfile)
  Action(comstr,OPT['A'])
