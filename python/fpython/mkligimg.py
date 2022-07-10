#!/usr/bin/env python
import sys
import os
from datetime import datetime
import time

LastModDate = "July 18, 2013"


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


def read_file_string(ifname):
  print "#read_file_string('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  content = ''
  for line in f:
    content += line
  f.close()
  return(content)

def write_file_string(ofname,content):
  print "#write_file_string('%s')"%(ofname)
  f = open(ofname,'w')
  f.write("%s"%(content))
  f.close()

def write_file_string_with_MODEL(ofname,fmode,content,model_number):
  print "#write_file_string_with_MODEL('%s' model_num %d)"%(ofname,model_number)
  if (fmode == 'w'):
    f = open(ofname,'w')
  if (fmode == 'a'):
    f = open(ofname,'a')
  f.write("MODEL        %d\n"%(model_number))
  f.write("%s"%(content))
  f.write("ENDMDL\n")
  f.close()

############
### MAIN ###
############

OPT = {}
OPT['L']   = 'ligand.list'
OPT['idlig'] = 'SupLIG'
OPT['odimg'] = 'outdir'
OPT['A'] = 'F'
OPT['iras'] = ''
OPT['iref'] = ''

if (len(sys.argv)<2):
  print "mkligimg.py <options>"
  print " for making images of ligand structure using RasMol."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L     :list of library molecules[%s]"%(OPT['L'])
  print " -idlig  :input  dir for molecules[%s]"%(OPT['idlig'])
  print " -odimg :output dir for images of ligands [%s]"%(OPT['odimg'])
  print " -iras  :input common script for RasMol[%s]"%(OPT['iras'])
  print " -iref  :input reference PDB molecule file[%s]"%(OPT['iref'])
  print " -A     :Action (T or F) [%s]"%(OPT['A'])
  sys.exit(1)


### [1] read option ###
read_option(sys.argv,OPT)

PID = os.getpid()

tmprasfile = "%d.ras"%(PID)
tmpppmfile = "%d.ppm"%(PID)
tmppdbfile = "%d.pdb"%(PID)


print "#tmprasfile '%s'"%(tmprasfile)  
print "#tmpppmfile '%s'"%(tmpppmfile)  

common_rasmol_script = ''
if (OPT['iras']!=''):
  common_rasmol_script = read_file_string(OPT['iras'])

ref_pdb_content = ''
if (OPT['iref']!=''):
  ref_pdb_content = read_file_string(OPT['iref'])


#write_file_string(tmprasfile,common_rasmol_script)
#sys.exit(1)
### [2] read liglist (-L) ###
liglist = []
read_list_file(OPT['L'],liglist)

Nlig  = len(liglist)

for lig in (liglist):  
  (head,tail) = lig.split('.')
  iligfile = OPT['idlig'] + '/' + lig
  oimgfile = OPT['odimg'] + '/' + head + '.png' 
  print lig,iligfile,oimgfile


  lig_pdb_content = read_file_string(iligfile)


  ## make rasmol_script ##  
  rasmol_script = 'load "%s"\n'%(tmppdbfile)
  rasmol_script += common_rasmol_script
  rasmol_script += 'select */2\n'
  rasmol_script += 'spacefill false\n'
  rasmol_script += 'wireframe false\n'
  rasmol_script += 'select */2 and !hydrogen\n'
  rasmol_script += 'wireframe\n'
  rasmol_script += 'color green\n'
  #rasmol_script += 'color [128,128,128]\n'
  rasmol_script += 'write ppm "%s"\n'%(tmpppmfile)
  rasmol_script += 'quit\n'
 
  rasmol_command = "rasmol -nodisplay -insecure -script %s"%(tmprasfile)
  convert_command = "convert %s %s"%(tmpppmfile,oimgfile)

  if (OPT['A']=='T'):
    write_file_string(tmprasfile,rasmol_script)
    write_file_string_with_MODEL(tmppdbfile,'w',lig_pdb_content,1)
    write_file_string_with_MODEL(tmppdbfile,'a',ref_pdb_content,2)


    os.system(rasmol_command)
    os.system(convert_command)
  else:
    print "%s"%(rasmol_script)
    print "#%s"%(rasmol_command)
    print "#%s"%(convert_command)


if (OPT['A']=='T'):
  os.system("rm %s"%(tmprasfile))
  os.system("rm %s"%(tmppdbfile))
  os.system("rm %s"%(tmpppmfile))
