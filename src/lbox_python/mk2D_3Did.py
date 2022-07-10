#!/usr/bin/env python
##
## <mk2D_3Did.py>
##


import sys
import os
import random
from datetime import datetime

LastModDate = "Oct 4, 2011"

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



def read_comment_in_mol2_file(fname,dat):
  #print "#read_comment_in_mol2_file('%s'):"%(fname)

#@<TRIPOS>COMMENT
#    LIGANDBOX_ID = 00071284-02
#    SOURCE = Namiki1010
#    NAMIKI_ID = NS-00108162
#    SUPPLIER = ButtPark
#    IDNUMBER = 107\40-54
#    MOLECULAR FORMULA = C20H23N7FS
#    MOLECULAR WEIGHT = 412.516
#    MOLECULAR CHARGE = 1
#    NUM OF DONOR = 2
#    NUM OF ACCEPTOR = 4
#    HOMO = -10.47
#    LUMO = -4.027
#    NUM OF CHIRAL ATOMS = 0
#    DELTAG = 3.4427
#    ATOMIC DELTAG = 0.0662
#    NOTE = ButtPark_107\40-54;SPECS_AO-476/42871082
#
#@<TRIPOS>MOLECULE
# NS-00108162.mol2
# 52 55 0 0 0
#SMALL
#GASTEIGER
  if (os.access(fname,os.R_OK)==0):
    print "#WARNING:can't open '%s'"%(fname) 
    return(0)
  f = open(fname)
  status = ' '
  for line in f:
    if (line.startswith('@<TRIPOS>COMMENT')):
      status = 'C'   
    elif (line.startswith('@<TRIPOS>')):
      if (status=='C'):
        return(0)
      status = ' '

    if (status == 'C') and (line.startswith('@')==0):
      line = line.rstrip('\n')
      field = line.split('=')
      if (len(field)>=2):
        head = field[0]
        value = field[1]
        head = head.lstrip(' ')
        head = head.rstrip(' ')
        value = value.lstrip(' ')
        value = value.rstrip(' ')
        #print "'%s' = '%s'"%(head,value) 
        dat[head] = value
  f.close()
  return(1)


###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/DB/LigandBox/namiki2010_3D'
OPT['ol']   = 'out.list'
OPT['tail'] = ''

if (len(sys.argv)<3):
  print "mk2D_3Did.py <options>"
  print " for counting files under '-idir'."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir : input src directory [%s]"%(OPT['idir']) 
  print " -tail : tail string required for file [%s]"%(OPT['tail'])
  print " -ol   : output file list file [%s]"%(OPT['ol']) 
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['idir'])
#for s in (xsublist):
#  print s

filelist = []

### (1) Getting filelist[] in the director 'idir' ###
Ndir = 0
for x in (dirlist):
  xfull = OPT['idir'] + '/' + x
  if (os.path.isdir(xfull)==0):
    #sys.stdout.write("%s\n"%(xfull))
    if (OPT['tail']=='') or (x.endswith(OPT['tail'])):
      filelist.append(x)
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
        filelist.append(yfull)
        Nfile_under_dir += 1 
    sys.stdout.write("%s (%d files)\n"%(xfull,Nfile_under_dir))

Nfile = len(filelist)
print "#Ndir %d Nfile %d"%(Ndir,Nfile)

if (Nfile==0):
  print "#ERROR:no file is found."
  sys.exit(1)

### (2) Making 3Didlist for each 2D id ###
id3Dfrm2D = {}

for x in (filelist):
  #print x
  dat = {}
  read_comment_in_mol2_file(OPT['idir']+'/'+x, dat)
  #print "%s %s"%(dat.get('LIGANDBOX_ID',''),dat.get('NAMIKI_ID',''))
  id3D = dat.get('LIGANDBOX_ID','')
  id2D = dat.get('NAMIKI_ID','')
  if (id2D != '') and (id3D != ''):
    if (id3Dfrm2D.has_key(id2D)==0):
      id3Dfrm2D[id2D] = []
    id3Dfrm2D[id2D].append(id3D)
 
### (3) output id3Dfrm2D[][] ###
if (OPT['ol'] != ''):
  of = open(OPT['ol'],'w')
  print "#write_id2D--id3Ds -->'%s'"%(OPT['ol'])
else:
  of = sys.stdout
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))
for id2D in id3Dfrm2D.keys():
  of.write("%s"%(id2D))
  for id3D in (id3Dfrm2D[id2D]):
    of.write(" %s"%(id3D))
  of.write("\n")
of.close()
