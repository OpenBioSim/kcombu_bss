#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Apr 16, 2014"

sys.path.append('/home/takawaba/work/kcombu/src/moldraw/')
import molecule

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



def read_list_file(ifname,list,prop):
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
     prop[fields[0]] = fields[1]
     list.append(fields[0])
  f.close()




############
### MAIN ###
############

OPT = {}

OPT['A'] = 'F'
OPT['of'] = '-'
OPT['isdf'] = ''
OPT['osdf'] = 'out.sdf'
OPT['olog'] = 'out.log'

if (len(sys.argv)<2):
  print "delHsulfonamide.py <options>"
  print "-SO2NH2 ==> -SO2NH"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " (1) Second hydrogen in SO2NH2 will be removed."
  print " (2) If two (SO2NH2) groups are contined, the first one will be dehydrogenated."
  print "<options>"
  print " -isdf : input SDF file (with SO2NH2) [%s]"%(OPT['isdf']) 
  print " -osdf : outpu SDF file (with SO2NH ) [%s]"%(OPT['osdf']) 
  print " -olog : outpu LOG file [%s]"%(OPT['olog']) 
  sys.exit(1)


read_option(sys.argv,OPT)
olog = open(OPT['olog'],'a')
olog.write("#COMMAND %s\n"%(OPT['COMMAND']))
olog.write("#DATE    %s\n"%(OPT['START_DATE']))

### [1] read SDF molecule (OPT['isdf']) ###
m = molecule.Molecule()
m.read_in_sdf(OPT['isdf'])
m.write_in_sdf("temp.sdf")

## [2] Find -SO2NH2 
##       O  H1
##       |  |
##      -S--N--H2
##       | 
##       O
##
## * H2 will be deleted.


Nsulfonamide = 0
numHdelete = -1 
for i in range(m.Natom):
  if (m.atoms[i].element=='S'):
    numS = i
    print "%d %d %s"%(i,m.atoms[i].num_in_file,m.atoms[i].atomname)
    ## check neighbors of 'S' ##
    NconS = {}
    numN = -1
    for j in range(m.Natom):
      if (j!=i) and (m.contable[i][j] != '0'):
        NconS[m.atoms[j].element] = NconS.get(m.atoms[j].element,0) + 1 
        if (m.atoms[j].element=='N'):
          numN = j

    if (NconS.get('O',0)>=2) and (NconS.get('N',0)==1):
      print "NconS",NconS
      ## check neighbors of 'N' ##
      numH1 = -1
      numH2 = -1
      NconN = {}
      for k in range(m.Natom):
        if (k!=numS) and (k!=numN) and (m.contable[numN][k] != '0'):
          NconN[m.atoms[k].element] = NconN.get(m.atoms[k].element,0) + 1 
          if (m.atoms[k].element=='H'):
            if (numH1==-1):
              numH1 = k
            else:
              numH2 = k
      print "NconN",NconN 
      if (NconN.get('H',0)==2):
        ### FIND Nsulfonamide ###
        print "#SULFONAMIDE S:%d N:%d H1:%d H2:%d"%(m.atoms[numS].num_in_file, m.atoms[numN].num_in_file, m.atoms[numH1].num_in_file, m.atoms[numH2].num_in_file)
        Nsulfonamide += 1
        if (numHdelete==-1):
          numHdelete = numH2

print "#%s Nsulfonamide %d"%(OPT['isdf'],Nsulfonamide)
olog.write("%s Nsulfonamide %d. "%(OPT['isdf'],Nsulfonamide))
if (numHdelete >=0):
  olog.write(" Hydrogen %d is removed."%(m.atoms[numHdelete].num_in_file))
olog.write("\n")

## [3]  Delete the hydrogen (H2) ### 
if (numHdelete >=0):
  newm = molecule.Molecule()
  newm.Natom      = m.Natom - 1
  newm.Nheavyatom = m.Nheavyatom
  newm.filename = m.filename
  natom = 0
  newatomnum_from_oldatomnum = {}
  for a in range(m.Natom):
    if (a != numHdelete):
      atom = molecule.Atom()
      atom.x = m.atoms[a].x
      atom.y = m.atoms[a].y
      atom.z = m.atoms[a].z
      atom.num = natom  
      atom.num_in_file = natom  + 1 
      atom.element = m.atoms[a].element 
      newatomnum_from_oldatomnum[a] = atom.num
      newm.atoms.append(atom)
      natom += 1
  newm.Nbond = 0
  for b in range(m.Nbond):
    if (m.bonds[b].natom1!=numHdelete) and (m.bonds[b].natom2!=numHdelete):
      bond = molecule.Bond()
      bond.num = newm.Nbond
      bond.natom1     = newatomnum_from_oldatomnum[m.bonds[b].natom1] 
      bond.natom2     = newatomnum_from_oldatomnum[m.bonds[b].natom2] 
      bond.bondtype   = m.bonds[b].bondtype
      bond.bondstereo = m.bonds[b].bondstereo
      newm.bonds.append(bond)
      newm.Nbond += 1
else:
  newm = molecule.Molecule()
  newm.read_in_sdf(OPT['isdf'])


### [4] write newm to OPT['osdf'] ###
if (numHdelete >=0):
  newm.write_in_sdf(OPT['osdf'])
olog.close()
