#!/usr/bin/env python

import sys
import os
import math
import re
import pgdb
from datetime import datetime
import time
import glob

LastModDate = "May 7, 2014"

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
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open '%s'."%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      fields = line.split()
      list.append(fields[0])
  f.close()


##### MAIN ####
OPT = {}
OPT['il'] = ''
OPT['dbname'] = 'pdb_mmcif_dc'
OPT['user']   = 'takawaba'
OPT['ol']    = 'out.list'
OPT['oerr'] = 'out.err'
OPT['maxheavy'] = '170'

if (len(sys.argv)<2):
  print "mkligprolist.py <options>"
  print "  for getting ligand-protein list from PDB_ID list."
  print "  coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -il       : input file for pdb_id [%s]"%(OPT['il'])
  print " -ol       : output file for ligand-receptor list [%s]"%(OPT['ol'])
  print " -oerr     : output file for error PDB_ID list [%s]"%(OPT['oerr'])
  print " -dbname   : database of PDB_MMCIF RDB [%s]"%(OPT['dbname']) 
  print " -user     : username of PDB_MMCIF RDB [%s]"%(OPT['user']) 
  print " -maxheavy : maximum number of heavy atom for ligand [%s]"%(OPT['maxheavy']) 
  sys.exit(1)

read_option(sys.argv,OPT)
pdb_id_list = []
read_list_file(OPT['il'],pdb_id_list)
print "#Npdb_id %d"%(len(pdb_id_list))
of = open(OPT['ol'],'w')
print "write_to_file() -->'%s'"%(OPT['ol'])

of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))

db = pgdb.connect(database=OPT['dbname'],user=OPT['user'])
cur = db.cursor()
error_pdbids_dic = {} 


for pdb_id in (pdb_id_list):
  prot_asym_ids = []
  prot_dic = {} 
  lig_asym_ids = []
  lig_dic = {}

  ## get asym_id for 'polypeptide'
  command = "SELECT asym_id,pdbx_description FROM unitmol WHERE pdb_id='%s' and poly_type='polypeptide(L)' ORDER BY nresidue DESC;"%(pdb_id)
  print "#command %s"%(command) 
  cur.execute(command)
  db.commit()
  rows = cur.fetchall()
  for row in (rows):
    prot_asym_ids.append(row[0])
    prot_dic[row[0]] = {}
    prot_dic[row[0]]['description'] = row[1]

  if (len(prot_asym_ids)==0):
    error_pdbids_dic[pdb_id] = error_pdbids_dic.get(pdb_id,'') + ' LOST_PDB_ID'

  ## get asym_id for 'ligand'
  command = "SELECT asym_id,comp_id,nheavyatom FROM unitmol WHERE pdb_id='%s' and nheavyatom>1 and nheavyatom<= %s ORDER BY nheavyatom DESC;"%(pdb_id,OPT['maxheavy'])
  print "#command %s"%(command) 
  cur.execute(command)
  db.commit()
  rows = cur.fetchall()
  for row in (rows):
    asym_id = row[0]
    lig_asym_ids.append(asym_id)
    lig_dic[asym_id] = {}
    lig_dic[asym_id]['comp_id'] = row[1]
    if (row[1]==''):
      lig_dic[asym_id]['comp_id'] = pdb_id+asym_id 


    lig_dic[asym_id]['nheavyatom'] = row[2]


  if (len(prot_asym_ids)>0) and (len(lig_asym_ids)==0):
    error_pdbids_dic[pdb_id] = error_pdbids_dic.get(pdb_id,'') + ' NO_LIGAND'
  of.write("#PDB_ID %s Nligand %d Nprotein %d\n"%(pdb_id,len(lig_asym_ids),len(prot_asym_ids)))

  for lig_asym_id in (lig_asym_ids):
    of.write("#LIGAND %s_%s_1_%s_1.pdb nheavy %s\n"%(lig_dic[lig_asym_id]['comp_id'],pdb_id,lig_asym_id,lig_dic[lig_asym_id]['nheavyatom']))     
  for prot_asym_id in (prot_asym_ids):
    of.write("#PROTEIN %s_1_%s_1.pdb %s\n"%(pdb_id,prot_asym_id,prot_dic[prot_asym_id]['description']))     

  ## CHECK CONTACT ### 
  Nconpair = 0
  for lig_asym_id in (lig_asym_ids):
    for prot_asym_id in (prot_asym_ids):
      command = "SELECT * FROM contact WHERE pdb_id='%s' and assembly_id='1' and asym_id_a='%s' and asym_id_b='%s';"%(pdb_id,lig_asym_id,prot_asym_id)
      cur.execute(command)
      db.commit()
      rowsA = cur.fetchall()
      command = "SELECT * FROM contact WHERE pdb_id='%s' and assembly_id='1' and asym_id_b='%s' and asym_id_a='%s';"%(pdb_id,lig_asym_id,prot_asym_id)
      cur.execute(command)
      db.commit()
      rowsB = cur.fetchall()
      if (len(rowsA)>0) or (len(rowsB)>0):
        of.write("%s_%s_1_%s_1.pdb "%(lig_dic[lig_asym_id]['comp_id'],pdb_id,lig_asym_id))
        of.write(" %s_1_%s_1.pdb"%(pdb_id,prot_asym_id))     
        of.write("\n")
        Nconpair += 1

db.close()
of.close()
print "write_to_file() -->'%s'"%(OPT['ol'])

print "#Nerror_pdbds %d"%(len(error_pdbids_dic.keys()))

if (OPT['oerr']!=''):
  print "write_error_pdbids() -->'%s'"%(OPT['oerr'])
  of = open(OPT['oerr'],'w')
  of.write("#LOST PDBIDs\n")
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#Nerror_pdbds %d\n"%(len(error_pdbids_dic.keys())))
  for pdbid in (error_pdbids_dic.keys()):
    of.write("%s %s\n"%(pdbid,error_pdbids_dic[pdbid]))
  of.close()
