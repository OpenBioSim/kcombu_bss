#!/usr/bin/env python
##
## <sele_rand_files.py>
##


import sys
import os
import random
from datetime import datetime

LastModDate = "Dec 12, 2011"

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


def read_Natom_in_sdf_file(ifname):
#>> EXAMPLE 1 << (/DB/LIGAND-EXPO/SDF2d)
#/DB/LIGAND-EXPO/SDF/ATP
#     RDKit          2D
#
# 31 33  0  0  0  0  0  0  0  0999 V2000
#   -7.0413    1.9118    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
#   -8.3230    2.6910    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#:
#>> EXAMPLE 2 << (/DB/kegg/medicus/drug/mol)
# 
# 
# 
# 19 22  0  0  0  0  0  0  0  0999 V2000
#   42.6766  -20.8611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   41.8327  -21.9863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#:
#>> EXAMPLE 3 << (/DB/LigandBox/namiki2010_2D)
#
#  Symyx   10231017462D 1   1.00000     0.00000     0
#
# 16 17  0     0  0            999 V2000
#   -1.1917   -0.3875    0.0000 O   0  0  0  0  0  0           0  0  0
#   -2.1250   -0.0625    0.0000 C   0  0  0  0  0  0           0  0  0
#:

  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING:Can't open '%s'"%(ifname)
    return(0)
  f = open(ifname) 
  n = 0
  Natom = 0
  for line in f:
    n += 1
    if ((n==4) and (len(line)>10)):
      str1 = line[0:3]
      Natom = int(str1)
      break 
  f.close()
  return(Natom)




def read_Natom_in_mol2_file(ifname):
#>>EXAMPLE OF MOL2 FILE FORMAT <<
#@<TRIPOS>MOLECULE
#5950
# 13 12 0 0 0
#SMALL
#GASTEIGER
#Energy = 0
#
#@<TRIPOS>ATOM
#      1 O           5.1350   -0.2500    0.0000 O.3     1  ALA1       -0.4795
#      2 OXT         4.2690    1.2500    0.0000 O.2     1  ALA1       -0.2493
#      3 N           2.5369    0.2500    0.0000 N.3     1  ALA1       -0.3186

  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING:Can't open '%s'"%(ifname)
    return(0)
  f = open(ifname) 
  Nline = 0
  Natom = 0
  status = ' ' 
  for line in f:
    #print line
    Nline += 1
    if (status=='M') and (Nline==2):
        field = line.split()
        Natom = int(field[0])
    if (line.startswith("@<TRIPOS>MOLECULE")):
      status = 'M'
      Nline = 0
    elif (line.startswith("@<TRIPOS>ATOM")):
      status = 'A'
    elif (line.startswith("@<TRIPOS>BOND")):
      status = 'B'

  f.close()
  print "filename '%s' Natom %d"%(ifname,Natom)
  return(Natom)







###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/DB/LIGAND-EXPO/SDF2d'
OPT['ol']   = 'out.list'
OPT['nout'] = -1
OPT['seed'] = 0
OPT['tail'] = ''
OPT['mina'] = 0 
OPT['fL'] = 'S'
OPT['odir'] = '' 
OPT['omm']  = '' 


if (len(sys.argv)<3):
  print "sele_rand_files.py <options>"
  print " for counting files under '-idir'."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir : input src directory [%s]"%(OPT['idir']) 
  print " -tail : tail string required for file [%s]"%(OPT['tail'])
  print " -nout : number of selecting molecule [%d]"%(OPT['nout'])
  print " -seed : seed number for random [%d]"%(OPT['seed'])
  print " -mina : minimum atom number [%d]"%(OPT['mina'])  
  print " -fL   : filetype 'S'df, '2':mol2. [%s]"%(OPT['fL'])
  print " -ol   : output list file for selected molecules [%s]"%(OPT['ol']) 
  print " -odir : output dir list  for selected molecules [%s]"%(OPT['odir']) 
  print " -omm  : output multi molecular file  for selected molecules [%s]"%(OPT['omm']) 
  sys.exit(1)

read_option(sys.argv,OPT)
random.seed(int(OPT['seed']))

OPT['nout'] = int(OPT['nout'])
OPT['mina'] = int(OPT['mina'])

dirlist = os.listdir(OPT['idir'])
#for s in (xsublist):
#  print s

FILELIST = []

### (1) Getting FILELIST[] in the director 'idir' ###
Ndir = 0
for x in (dirlist):
  xfull = OPT['idir'] + '/' + x
  if (os.path.isdir(xfull)==0):
    #sys.stdout.write("%s\n"%(xfull))
    if (OPT['tail']=='') or (x.endswith(OPT['tail'])):
      FILELIST.append(x)
  else:
    Ndir += 1
    #sys.stdout.write(">%s"%(xfull))
    dirlist_under = os.listdir(xfull)
    Nfile_under_dir = 0
    for y in (dirlist_under):
      #yfull = OPT['idir']+'/'+ x + '/' + y
      yfull = x + '/' + y
      #print y,yfull
      if (OPT['tail']=='') or (y.endswith(OPT['tail'])):
        FILELIST.append(yfull)
        Nfile_under_dir += 1 
    sys.stdout.write("%s (%d files)\n"%(xfull,Nfile_under_dir))

Nfile = len(FILELIST)
print "#Ndir %d Nfile %d"%(Ndir,Nfile)

if (Nfile==0):
  print "#ERROR:no file is found."
  sys.exit(1)

### (2) Choose output FILELIST and SELECT_MARK_FILE[i] := 1 ###
if (OPT['nout'] < 0):
  SELECT_MARK_FILE = [1 for i in range(len(FILELIST))] 
else:
  SELECT_MARK_FILE = [0 for i in range(len(FILELIST))] 
  for i in range(OPT['nout']):
    ok = 0
    Nround = 0
    while (ok==0): 
      Nround += 1
      k = random.randint(0,len(FILELIST)) 
      ok = 0
      if (SELECT_MARK_FILE[k]==0):
        ok = 1
        if (OPT['mina']>0):
          if (OPT['fL']=='2'):
            Natom = read_Natom_in_mol2_file(OPT['idir'] + '/' + FILELIST[k])
          else:
            Natom = read_Natom_in_sdf_file(OPT['idir'] + '/' + FILELIST[k])
          if (Natom<OPT['mina']):
            ok = 0 
     
      if (ok==1):   
        SELECT_MARK_FILE[k] = 1
      if (Nround >= len(FILELIST)):
        print "#ERROR:too much Nround %d"%(Nround) 
        sys.exit(1)

### (3) write output FILELIST in OPT['ol'] ##
if (OPT['ol'] != ''):
  print "#write_list_of_selected_molecules(Nout:%d) --> '%s'"%(OPT['nout'],OPT['ol'])
  of = open(OPT['ol'],'w') 
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#NFILE %d NFILEOUT %d\n"%(Nfile,OPT['nout']))
  for i in range(len(FILELIST)):
    if (SELECT_MARK_FILE[i]==1):
      of.write("%s\n"%(FILELIST[i]))
  of.close() 

if (OPT['odir'] != ''):
  print "#write_selected_molecules(Nout:%d) --> in the directory '%s/' "%(OPT['nout'],OPT['odir'])
  for i in range(len(FILELIST)):
    if (SELECT_MARK_FILE[i]==1):
      ifilename = OPT['idir'] + '/' + FILELIST[i]
      ofilename = OPT['odir'] + '/' + FILELIST[i]
      command = "cp %s %s"%(ifilename,ofilename)     
      print "#'%s'"%(command) 
      os.system(command)


if (OPT['omm'] != ''):
  print "#write_multi_molecular_files_for_selected_molecules(Nout:%d) --> '%s' "%(OPT['nout'],OPT['omm'])
  ofm = open(OPT['omm'],'w') 
  for i in range(len(FILELIST)):
    if (SELECT_MARK_FILE[i]==1):
      ifilename = OPT['idir'] + '/' + FILELIST[i]
      if (os.access(ifilename,os.R_OK)==0):
        print "#WARNING:Can't open '%s'"%(ifilename)
      else:
        fi = open(ifilename)   
        for line in fi:
          ofm.write("%s"%(line))
        fi.close()
        if (OPT['fL']=='S'):
          ofm.write("$$$\n")
  ofm.close() 


