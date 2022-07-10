#!/usr/bin/env python
##
## <gen2Dfrm3D.py>
##


import sys
import os
import random
from datetime import datetime

LastModDate = "Apr 3, 2012"

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



###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/DB/LigandBox/namiki2010_2D'
OPT['odir'] = '.'
OPT['of']   = ''
OPT['ol']   = ''
OPT['odl'] = ''
OPT['tail'] = '01.mol2'
OPT['newtail'] = ''
OPT['A'] = 'F'
OPT['rmconf'] = 'T'
OPT['rpname'] = 'T'
OPT['mkodir'] = 'T'
OPT['div'] = '0/1'

if (len(sys.argv)<3):
  print "gen2Dfrm3D.py <options>"
  print " for making compound list for each supplier."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir    : input molecular file directory [%s]"%(OPT['idir']) 
  print " -tail    : tail string required for molecular file [%s]"%(OPT['tail'])
  print " -newtail : new tail string replaced for the tail [%s]"%(OPT['newtail'])
  print " -odir    : output molecular file directory [%s]"%(OPT['odir']) 
  print " -rmconf  : remove conformer number (such as '-01','02')[%s]"%(OPT['rmconf'])
  print " -rpname  : replace molecular name by input file head string (T or F)[%s]"%(OPT['rpname'])
  print " -mkodir  : make output directory, if it does not exist. (T or F)[%s]"%(OPT['mkodir'])
  print " -A       : Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)

[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)



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


Nsta = Nfile*bunshi/bunbo
Nend = Nfile*(bunshi+1)/bunbo
print "#bunshi %d bunbo %d Nmodel %d Nsta %d Nend %d"%(bunshi,bunbo,Nfile,Nsta,Nend)

### (2) make and write 2D-generated SDF file using 'ckcombu' ###
for i in range(Nsta,Nend):
#for imol in (filelist):
  imol = filelist[i] 
  print imol 
  if (imol.find('/')>=0):
    (idir,ifile) = imol.split('/')
  else:
    ifile = imol
    idir = ''
  (ihead,itail) = ifile.split('.')

  print "idir '%s' ifile '%s' ihead '%s' itail '%s'"%(idir,ifile,ihead,itail)
  ifile_full = OPT['idir'] + '/' + imol


  if (OPT['rmconf'] == 'T'):
    if (ihead.endswith('-01')):
      ihead = ihead.replace('-01','') 

  odir_full  = OPT['odir'] + '/' + idir
  ofile_full = odir_full + '/' +  ihead + '.sdf' 

  comstr = ''
  if (OPT['mkodir']=='T') and (os.path.exists(odir_full)==0):
    comstr = "mkdir %s"%(odir_full)
    print "#%s"%(comstr) 
    if (OPT['A'] == 'T'):
      os.system(comstr)

  comstr = "ckcombu -A %s -2dA T -prp T -ansm T -osA %s"%(ifile_full,ofile_full)
  if (OPT['rpname'] == 'T'):
    comstr += " -nameA %s"%(ihead)

  print "#%s"%(comstr) 
  if (OPT['A'] == 'T'):
    os.system(comstr)



  pass

sys.exit(1)

