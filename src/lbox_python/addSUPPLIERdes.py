#!/usr/bin/env python
##
## <mkSUPPLIERdes.py>
##


import sys
import os
import random
from datetime import datetime

LastModDate = "Oct 7, 2011"

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


def make_new_headline_adding_supplier_number_array(headline,num_supplier,Ncount_supplier):
## >> head line example ##
#>00000001_00005000/NS-00002705.sdf 23 C19_O3_CL C@ 12 C 7 O 2 O1 1 X 1
#>00000001_00005000/NS-00001572.sdf 17 C14_N3 C@ 12 N@ 3 C 2 

  line = headline.rstrip('\n')
  line = line.lstrip('>')
  field = line.split()
  isdffile = OPT['idl'] + '/' + field[0]
  #print "isdffile '%s'"%(isdffile)
  dat = {}
  supp_array =  [] 
  chem_supp = {}
  read_comment_in_sdf_file(isdffile,dat)
  for key in (dat.keys()):
    if (key.startswith('SUPPLIERNAME')):
      supp = dat[key]
      if (num_supplier.has_key(supp)==0):
        num_supplier[supp] = len(num_supplier.keys())
      #print "%s %d"%(supp,num_supplier[supp])
      supp_array.append(num_supplier[supp])
      chem_supp[num_supplier[supp]] = 1
      Ncount_supplier[supp] = Ncount_supplier.get(supp,0) + 1

  supp_array = chem_supp.keys()
  supp_array.sort()
  supp_str = ''
  for x in (supp_array):
    if (supp_str != ''):
      supp_str = supp_str + ':'
    supp_str = supp_str + "%d"%(x)
  new_headline = ">%s %%0 %s\n"%(line,supp_str)
  return(new_headline)


###############
#### MAIN #####
###############

OPT = {}
OPT['ides']   = ''
OPT['odes']   = 'out.des'
OPT['okey']   = 'out.key'
OPT['idl'] = ''


if (len(sys.argv)<2):
  print "mkSUPPLIERdes.py <options>"
  print " for counting files under '-idir'."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ides  : input descriptor file [%s]"%(OPT['ides']) 
  print " -idl   : input library directory (in SDF) [%s]"%(OPT['idl']) 
  print " -odes  : output descriptor file [%s]"%(OPT['odes']) 
  print " -okey  : output key number file [%s]"%(OPT['okey']) 
  sys.exit(1)

read_option(sys.argv,OPT)

if (os.access(OPT['ides'],os.R_OK)==0):
  print "#WARNING:can't open sdf file '%s'"%(idesfile) 
  sys.exit(1)
fi = open(OPT['ides'],"rb")
fo = open(OPT['odes'],"wb")
print "#write_descriptor_file() -->'%s'"%(OPT['odes'])

b0 = '\n' 
b = ' '
status = ' '
Nsupplier = 0
num_supplier = {}
Ncount_supplier = {}

while (b != ''):
  b = fi.read(1)
 
  if (status=='D') and (ord(b0)==255) and (b=='\n'): ## for binary mode ##
    status = ' '
  
  if (status=='D') and (b0=='/') and (b=='\n'):     ## for text mode ##
    status = ' '

  if (status == ' ') and (b0=='\n') and (b=='>'):
    status = '>'
    headline = ''

  if (status=='>'):
    headline = headline + b 
  else:
    fo.write(b)   

  if (status=='>') and (b=='\n'):
    ## MODIFY headline ###
    new_headline =  make_new_headline_adding_supplier_number_array(headline,num_supplier,Ncount_supplier)
    #sys.stdout.write("%s"%(new_headline))
    fo.write("%s"%(new_headline))
    status = 'D'

  b0 = b
fi.close()


### OUTPUT KEY NUMBER DESCRIPTION ###
fk = open(OPT['okey'],'w')

fo.write("#KEY %0 DESCRIPTION\n")


fk.write("#COMMAND %s\n"%(OPT['COMMAND']))
fk.write("#DATE    %s\n"%(OPT['START_DATE']))
fk.write("#KEY %0 DESCRIPTION\n")
for supp in (sorted (num_supplier.keys(), lambda x,y:cmp(num_supplier[x],num_supplier[y]))):
  fo.write("#%d %s\n"%(num_supplier[supp],supp))
  fk.write("%d %s\n"%(num_supplier[supp],supp))
  sys.stdout.write("#%d %s\n"%(num_supplier[supp],supp))


fo.close()
