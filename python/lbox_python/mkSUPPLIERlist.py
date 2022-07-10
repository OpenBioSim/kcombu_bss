#!/usr/bin/env python
##
## <mkSUPPLIERlist.py>
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



###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/DB/LigandBox/namiki2010_2D'
OPT['of']   = ''
OPT['ol']   = ''
OPT['odl'] = ''
OPT['tail'] = ''
OPT['A'] = 'F'

if (len(sys.argv)<3):
  print "mkSUPPLIERlist.py <options>"
  print " for making compound list for each supplier."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir : input SDF file directory [%s]"%(OPT['idir']) 
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
if (OPT['of']=='-'):
  of = sys.stdout.out
elif (OPT['of']!=''):
  of = open(OPT['of'],'w')
  print "#write_molecule_and_suppliers() --> '%s'"%(OPT['of'])

Nsuppl = {}
MOL_SUPPL = {}
index_suppl = {}
for x in (filelist):
  #print x
  dat = {}
  read_comment_in_sdf_file(OPT['idir']+'/'+x, dat)
  #print ">%s"%(x)
  MOL_SUPPL[x] = [] 
  suppl_dic = {}
  for ind in (dat.keys()):
    if (ind.startswith('SUPPLIERNAME')):
      #print "'%s' '%s'"%(ind,dat[ind])
      suppl = dat[ind]
      if (index_suppl.has_key(suppl)==0):
        index_suppl[suppl] = len(index_suppl.keys())  + 1
      suppl_dic[suppl] = 1

  for suppl in (suppl_dic.keys()):
    MOL_SUPPL[x].append(index_suppl[suppl])
    Nsuppl[suppl] = Nsuppl.get(suppl,0) + 1
  if (OPT['of'] != ''):
    of.write("%s"%(x))      
    for suppl in (suppl_dic.keys()):
      of.write(" %d"%(index_suppl[suppl])) 
    of.write(" \n") 



### (3) output supllier list file ###
if (OPT['ol']=='-'):
  ol = sys.stdout.out
elif (OPT['ol']!= ''):
  ol = open(OPT['ol'],'w')
  print "#write_supplier_list() --> '%s'"%(OPT['ol'])
else:
  ol = ''


Nmolecule  = len(filelist)

slist =  sorted(Nsuppl.keys(), lambda x,y : cmp(Nsuppl[y],Nsuppl[x]))
#slist =  sorted(Nsuppl.keys(), lambda x,y : cmp(index_suppl[x],index_suppl[y]))
Ncomp_supp = 0
for y in (slist):
  Ncomp_supp += Nsuppl[y]

if (OPT['ol'] != ''):
  ol.write("#COMMAND  %s\n"%(OPT['COMMAND']))
  ol.write("#DATE     %s\n"%(OPT['START_DATE']))
  ol.write("#Nmolecule %d\n"%(Nmolecule))
  for y in (slist):
    ol.write("#[%3d] %8d %6.1f %% %s\n"%(index_suppl[y],Nsuppl[y],100.0*float(Nsuppl[y])/Nmolecule,y))
  ol.write("#%7d Nchem_total\n"%(Ncomp_supp))
  ol.close()

if (OPT['of'] != ''):
  of.write("#Nmolecule %d\n"%(Nmolecule))
  for y in (slist):
    of.write("#[%3d] %8d %6.1f %% %s\n"%(index_suppl[y],Nsuppl[y],100.0*float(Nsuppl[y])/Nmolecule,y))
  of.write("#%7d Nchem_total\n"%(Ncomp_supp))
  of.close()


### (4) Output each supplier molecular list ###
if (OPT['odl'] != ''):
  for suppl in (slist):
    s = suppl
    s = s.replace('/','_')
    s = s.replace('.','_')
    s = s.replace('\\','_')
    s = s.replace('+','_')
    s = s.replace('*','_')
    s = s.replace(' ','_')
    s = s.replace('(','_')
    s = s.replace(')','_')
    ofname = OPT['odl'] + '/' + s + '.list'
    print "'%s'-->'%s'"%(suppl,ofname)
    if (OPT['A']=='T'):
      f = open(ofname,'w')
      f.write("#COMMAND  %s\n"%(OPT['COMMAND']))
      f.write("#DATE     %s\n"%(OPT['START_DATE']))
      f.write("#SUPPLIER %s\n"%(suppl))
      f.write("#SUPPLIER_STRING %s\n"%(s))
      f.write("#N_MOLECULE_FOR_SUPPLIER %d\n"%(Nsuppl[suppl]))
      list = []
      for x in (filelist):
        if (index_suppl[suppl] in MOL_SUPPL[x]):
          list.append(x)

      for x in (sorted(list,lambda x,y:cmp(x,y))):
        f.write("%s\n"%(x))
      f.close()
