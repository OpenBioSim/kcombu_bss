#!/usr/bin/env python
##
## <mkPROPERTYlist.py>
##


import sys
import os
import random
from datetime import datetime

LastModDate = "Nov 14, 2011"

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

def read_comment_in_sdf_file(fname,dat):
#  18 17  2  0     0  0
#  19 16  1  0     0  0
#  20 18  1  0     0  0
#  13  9  2  0     0  0
#  19 20  2  0     0  0
# M  END
# > <NAMIKI_ID> (NS-00000475)
# NS-00000475
# 
# > <SUPPLIERNAME_1> (Bionet)
# Bionet
# 
# > <SUPPLIERID_1> (5W-0887)
# 5W-0887
# 
# > <SUPPLIERNAME_2> (Pharmeks(Natural))
# Pharmeks(Natural)
# 
# > <SUPPLIERID_2> (P2000N-48472)
# P2000N-48472


  if (os.access(fname,os.R_OK)==0):
    print "#WARNING:can't open sdf file '%s'"%(fname) 
    return(0)
  f = open(fname)
  status = ' '
  item = ''
  for line in f:
    line = line.rstrip('\n')
    line = line.rstrip('\r')
    if line.startswith('>'):
      field = line.split(' ')
      if (len(field)>1):
        item  = field[1]
        item  = item.lstrip('<')
        item  = item.rstrip('>')
        status = 'C'
    elif (status=='C'):
      if (line.startswith('>')==0) and (len(line)>1) and (item != ''):
        value = line
        # 'Pharmeks(Natural)' --> 'Pharmeks'
        #if value.endswith('(Natural)'):
        #  value = value.replace('(Natural)','')
        dat[item] = value
      status = ' '
      item = ''
 
  f.close()
  return(1)





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
OPT['of']   = ''
OPT['ol']   = ''
OPT['odl'] = ''
OPT['tail'] = '.mol2'
OPT['A'] = 'F'

if (len(sys.argv)<3):
  print "mkPROPERTYlist.py <options>"
  print " for making compound list satisfying given property conditions in mol2 file"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir : input mol2 file directory [%s]"%(OPT['idir']) 
  print " -tail : tail string required for molecular file [%s]"%(OPT['tail'])
  print " -of   : output file for molecules and suppliers [%s]"%(OPT['of']) 
  print " -ol   : output file for supplier list file [%s]"%(OPT['ol']) 
  print " -odl  : output dir for each supplier molecular list  [%s]"%(OPT['odl']) 
  print " -A    : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)




### (1) Getting filelist[] in the director 'idir' ###
dirlist = os.listdir(OPT['idir'])
#for s in (xsublist):
#  print s
filelist = []

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

### (2) Read a comment in each SDF file ###
Nsuppl = {}
MOL_SUPPL = {}
index_suppl = {}
for x in (filelist):
  #print x
  dat = {}
  #read_comment_in_sdf_file(OPT['idir']+'/'+x, dat)
  read_comment_in_mol2_file(OPT['idir']+'/'+x, dat)
  hit = 0
  if (int(dat['NUM OF DONOR'])<=10) and (int(dat['NUM OF ACCEPTOR'])<=10) and (float(dat['MOLECULAR WEIGHT'])<=500.0) and (float(dat['DELTAG'])<5.0):
    hit = 1 
  print ">%s HIT %d"%(x, hit)
  #for ind in (dat.keys()):
  #  print " '%s' '%s'"%(ind,dat[ind])

