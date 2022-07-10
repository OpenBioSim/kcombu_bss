#!/usr/bin/env python
import sys
import os


LastModDate = "June 18,2012"


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
     list.append(fields[1])
  f.close()

def read_ligand_pdb_sws_list(ifname,ligpdblist,prop):
  print "#read_ligand_pdb_sws_list('%s')"%(ifname)
  # #[ligand]_[pdb][chain] [SWS_ID] [residue number for ligand] [chaindID for ligand]
  # LW4_3iw4A  KPCA 901 A
  # 07U_3txoA  KPCL 1 A
  # ANP_3g51A  KS6A3_1 480 A
  # ANP_3kn5A  KS6A5_2 400 A
  # STU_3a60A  KS6B1 400 A
  # ATP_2y4iB  KSR2 1932 B
  # STU_3s95A  LIMK1 1 A


  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open ligand_pdb_sws_list '%s'" %(ifname)
    sys.exit(1)
  f = open(ifname)

  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0) and (len(line)>10):
      field = line.split()
      ligpdblist.append(field[0])
      (lig,pdbch) = field[0].split('_')
      sws = rnum = chain = ''
      if (len(field)>1):
        sws = field[1]
      if (len(field)>2):
        rnum = field[2]
      if (len(field)>3):
        chain = field[3]
      prop[field[0]] = {}
      prop[field[0]]['LINE']   = line 
      prop[field[0]]['LIG']   = lig 
      prop[field[0]]['PDBCH'] = pdbch 
      prop[field[0]]['SWS']   = sws
      prop[field[0]]['RNUM']  = rnum
      prop[field[0]]['CHAIN'] = chain
      #print "'%s' '%s' '%s' '%s' '%s'"%(line,field[0],sws,rnum,chain)
  print "#len(ligpdblist):%d"%(len(ligpdblist))



if (len(sys.argv)<2):
  print "get_multilig.py [uniid_list_file] [sequence file]"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  sys.exit(1)

ilistfile = sys.argv[1]
ligpdblist = []
ligpdbdat = {}
read_ligand_pdb_sws_list(ilistfile,ligpdblist,ligpdbdat)

Nlig = {}
for x in (ligpdblist):
  lig = ligpdbdat[x]['LIG']
  Nlig[lig] = Nlig.get(lig,0) + 1

for x in (ligpdblist):
  lig = ligpdbdat[x]['LIG']
  if (Nlig[lig]>1):
    print ligpdbdat[x]['LINE'],Nlig[lig]


