#!/usr/bin/env python
import sys
import os
from datetime import datetime


LastModDate = "Feb 18, 2011"

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


def read_atommatch_file_and_make_inverse_str(ifname):
  print "#read_atommatch_file_and_make_inverse_str(%s)"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  invstr = ''
  status = ''

  for line in f:
    if (line.startswith('#COMMAND')):
      line += "#COMMANDinv '%s'\n"%(OPT['COMMAND'])
    if (line.startswith('#DATE')):
      line += "#DATEinv '%s'\n"%(OPT['START_DATE'])
    if (line.startswith('#MoleculeA')):
      line = line.replace('#MoleculeA','#MoleculeB')
    elif (line.startswith('#MoleculeB')):
      line = line.replace('#MoleculeB','#MoleculeA')
    elif (line.startswith('#NatomA')):
      line = line.replace('#NatomA','#NatomB')
    elif (line.startswith('#NatomB')):
      line = line.replace('#NatomB','#NatomA')
    elif (line.startswith('#NheavyatomA')):
      line = line.replace('#NheavyatomA','#NheavyatomB')
    elif (line.startswith('#NheavyatomB')):
      line = line.replace('#NheavyatomB','#NheavyatomA')
    elif (line.startswith('#NpermuA')):
      line = line.replace('#NpermuA','#NpermuB')
    elif (line.startswith('#NpermuB')):
      line = line.replace('#NpermuB','#NpermuA')

##[numA] [num_in_fileA] [atomnameA] --- [numB] [num_in_fileB] [atomnameB] [numpair] [atomtype] [ECA] [ECB] [ECdiff] [Nnei_diff]\n");
#>1
#1     2024   N1  --- 17    2078   N17  1 N1    6  6  0  0
#2     2025   S1  --- 1     2062   S1   2 S    21 20  1  0
#3     2026   O1  --- 19    2080   O19  3 O1    6  6  0  0
#4     2027   O2  --- 18    2079   O18  4 O1    6  6  0  0
#//
    if (line.startswith('/')==1):    
      status = ''
    if (status=='a') and (line.startswith('#')==0):
      line = line.rstrip('\n')
      #print "%s\n"%(line)
      X = line.split(' --- ')
      Y = X[1][16:].split()
      line = "%s --- %s %2s %-4s %2s %2s %2s %2s\n"%(X[1][0:16],X[0],Y[0],Y[1],Y[3],Y[2],Y[4],Y[5])
      #print "%s\n"%(line)

    if (line.startswith('>')==1):    
      status = 'a'
    invstr += line
  f.close()
  return(invstr)







############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['idam'] = ''
OPT['odam'] = ''
OPT['A'] = ''

if (len(sys.argv)<2):
  print "cpinv_atommatch.py <options>"
  print " for making atom match file for all-vs-all comparison for a given molecule list."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L    : list of library molecules[%s]"%(OPT['L'])
  print " -idam : input  directory for atom_maching file[%s]"%(OPT['idam'])
  print " -odam : output directory for atom_maching file[%s]"%(OPT['odam'])
  print "<options for target matching>"
  
  sys.exit(1)


### [1] read option ###

read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
liglist = []
read_list_file(OPT['L'],liglist)


### [3] Copy pairs ##

Ncopy = 0
Ninv  = 0
Nnot  = 0
for i in range(len(liglist)):
  ligA = liglist[i]
  for j in range(i+1,len(liglist)):
    ligB = liglist[j]
    pair = {}
    pair['A'] = ligA
    pair['B'] = ligB
    infile      = "%s/%s_%s"%(OPT['idam'],ligA,ligB)
    infile_inv  = "%s/%s_%s"%(OPT['idam'],ligB,ligA)
    outfile     = "%s/%s_%s"%(OPT['odam'],ligA,ligB)
    if (os.path.isfile(infile)==1):
      command = "cp %s %s"%(infile,outfile)
      print "#%s"%(command)
      if (OPT['A']=='T'):
        os.system(command)
      Ncopy += 1
    elif (os.path.isfile(infile_inv)==1):
      print "#inverse exists (%s --> %s)."%(infile_inv,outfile)
      invstr = read_atommatch_file_and_make_inverse_str(infile_inv)
      print "%s"%(invstr)
      if (OPT['A']=='T'):
        of = open(outfile,'w')
        of.write("%s"%(invstr))
        of.close()
      Ninv += 1
    else:
      print "#not exists (%s and %s)."%(infile,infile_inv)
      Nnot += 1

print "#Ncopy %d Ninv %d Nnot %d"%(Ncopy,Ninv,Nnot)

