#!/usr/bin/env python
import sys
import os
from datetime import datetime
import time


LastModDate = "Sep 2, 2013"


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
    os.system(comstr)  
  else:
    print "#%s"%(comstr)


def add_comp_time_to_sdf(ifname,comp_time):
  print "#add_comp_time_to_sdf('%s',comp_time)"%(ifname)
  lines = []
  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING(add_comp_time_to_sdf):Can't open '%s'."%(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    lines.append(line)
  f.close()

  of = open(ifname,'w')
  #REMARK COMP_TIME  3.891029 seconds
  for i in range(len(lines)):
    line = lines[i]
    of.write("%s"%(line))
  of.write(">  <comp_time>\n")
  of.write("%f seconds\n"%(comp_time))
  of.write("\n")
  of.write("$$$$\n")
  of.close()
  return(1)



def add_comp_time_to_mol2(ifname,comp_time):
  print "#add_comp_time_to_mol2('%s',comp_time)"%(ifname)
  lines = []
  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING(add_comp_time_to_mol2):Can't open '%s'."%(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    lines.append(line)
  f.close()

  of = open(ifname,'w')
  for i in range(len(lines)-1):
    line = lines[i]
    of.write("%s"%(line))
  of.write("@<TRIPOS>COMMENT\n")
  of.write("  COMP_TIME = %f seconds\n"%(comp_time))
  of.write("\n")
  of.close()
  return(1)




############
### MAIN ###
############

OPT = {}
OPT['L'] = 'ligand.list'
OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['idsdf']  = ''
OPT['idpdb']  = ''
OPT['idsmi']  = ''
OPT['odmol2'] = ''
OPT['odsdf']  = ''

OPT['c'] = '20'
OPT['i'] = '20'
OPT['nGenerations'] = '200'
OPT['randomSeed'] = '1'
OPT['singleconf'] = 'F'
OPT['onlycharge'] = 'F'
OPT['H'] = 'F'
OPT['f'] = '/home/takawaba/tool/balloon_x86_64/MMFF94.mff'
OPT['iligname3'] = 'T'
OPT['oligname3'] = 'T'

if (len(sys.argv)<2):
  print "mk3Dfrm2D_balloon.py <options>"
  print " for making 3D mol2 file from 2D sdf file using balloon"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<basic usage>"
  print "mk3Dfrm2D_ballon.py -L [ligand_list] -odmol2 [MOL2_3D_BALLOON]"
  print "mk3Dfrm2D_ballon.py -L [ligand_list] -odsdf  [SDF_3D_BALLOON]"
  print "<options for input/output>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print "** '-isdf' or '-idsmi' should be assigned." 
  print " -idsdf  : input directory for SDF molecules[%s]"%(OPT['idsdf'])
  print " -idpdb  : input directory for PDB molecules[%s]"%(OPT['idpdb'])
  print " -idsmi  : input directory for SMILES for molecules[%s]"%(OPT['idsmi'])
  print "** '-odmol2' or '-iodsdf' should be assigned." 
  print " -odmol2 : output directory for 3D MOL2 molecules[%s]"%(OPT['odmol2'])
  print " -odsdf  : output directory for 3D SDF molecules[%s]"%(OPT['odsdf'])
  print "<other options>"
  print " -iligname3: input  file name is just ligname3 or not 'T' or 'F' [%s]"%(OPT['iligname3'])
  print " -oligname3: output file name is just ligname3 or not 'T' or 'F' [%s]"%(OPT['oligname3'])
  print " -div  : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -nGenerations : max number of generation [%s]"%(OPT['nGenerations'])
  print " -c    : Number of conformers to generate [%s]"%(OPT['c'])
  print " -i    : Maximum iterations for conjugate gradient [%s]"%(OPT['i'])
  print " -randomSeed : seed of random number generator [%s]"%(OPT['randomSeed'])
  print " -H    :Add hydrogens to the structures ('T' or 'F') [%s]"%(OPT['H']) 
  print " -f    :location of energy parameter file [%s]"%(OPT['f']) 
  print " -singleconf :Output only the lowest-energy conformation (T or F) [%s]"%(OPT['singleconf'])
  print " -onlycharge :Only assign partial atomic charges to the input structures(T or F)[%s]"%(OPT['onlycharge'])
  print " -A    : Action (T or F) [%s]"%(OPT['A'])
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
  field = lig.split('_')
  ligname3 = field[0]
  (ligfilecore,tail) = lig.split('.')
  #(ligname3,pdbch) = lig.split('_')
  #print "%s '%s' '%s'"%(lig,ligname3,pdbch)
 
  if (OPT['odmol2'] != ''): 
    outfile = "%s/%s.mol2"%(OPT['odmol2'],ligname3)
  if (OPT['odsdf'] != ''): 
    outfile = "%s/%s.sdf"%(OPT['odsdf'],ligname3)

  if (OPT['oligname3']=='T'): 
    ofilecore = ligname3
  else:
    ofilecore = ligfilecore

  if (OPT['odmol2'] != ''):
    outfile = "%s/%s.mol2"%(OPT['odmol2'],ofilecore) 
  if (OPT['odsdf'] != ''):
    outfile = "%s/%s.sdf"%(OPT['odsdf'],ofilecore) 
  
  if (os.path.exists(outfile)):
    print "#WARNING:outfile '%s' already exists."%(outfile)
    comstr = "rm %s"%(outfile) 
    Action(comstr,OPT['A'])



  if (OPT['iligname3']=='T'): 
    ifilecore = ligname3
  else:
    ifilecore = ligfilecore

  if (OPT['idsdf'] != ''):
    infile = "%s/%s.sdf"%(OPT['idsdf'],ifilecore) 
  elif (OPT['idpdb'] != ''):
    infile = "%s/%s.pdb"%(OPT['idpdb'],ifilecore) 
  elif (OPT['idsmi'] != ''):
    infile = "%s/%s.smi"%(OPT['idsmi'],ifilecore) 

  comstr = "balloon --input-file %s --output-file %s"%(infile,outfile)

  comstr += " --nGeneration %s -c %s -i %s"%(OPT['nGenerations'],OPT['c'],OPT['i']) 

  if (OPT['randomSeed'] != ''):
    comstr += " --randomSeed %s"%(OPT['randomSeed'])

  if (OPT['singleconf'] == 'T'):
    comstr += " --singleconf"

  if (OPT['onlycharge'] == 'T'):
    comstr += " --onlycharge"

  if (OPT['f'] != ''):
    comstr += " -f %s"%(OPT['f'])
  
  if (OPT['H'] == 'T'):
    comstr += " -H"

  sta_sec = time.time()
  Action(comstr,OPT['A'])
  end_sec = time.time()
  delta_sec = end_sec - sta_sec
  if (OPT['A']=='T') and (OPT['odsdf']!=''):
    add_comp_time_to_sdf(outfile,delta_sec)
  if (OPT['A']=='T') and (OPT['odmol2']!=''):
    add_comp_time_to_mol2(outfile,delta_sec)
