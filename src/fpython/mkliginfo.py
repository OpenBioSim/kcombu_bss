#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Apr 14, 2014"

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
  #print "#read_list_file('%s')"%(ifname)
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


def read_ligand_info_in_PDB(ipdbfile,DAT):
  f = open(ipdbfile)
#REMARK     [1]:auth_asym_id [2]:label_seq_id
#REMARK     [3]:label_asym_id [4]:auth_seq_id [5]:type_symbol
#REMARK              [1] [2]                                       [3]   [4][5]
#          1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#HETATM 2474  C17 Y8L A   1       5.548  23.033   6.437  0.50  0.00  B  1299  C
#HETATM 2476  C18 Y8L A   1       5.706  21.502   6.484  0.50  0.00  B  1299  C
#HETATM 2478  N19 Y8L A   1       4.841  20.873   7.494  0.50  0.00  B  1299  N
#HETATM 2480  C20 Y8L A   1       5.005  21.589   8.751  0.50  0.00  B  1299  C
#:
#HETATM 2248  C22 BYP A   1       4.673  22.552   7.714  0.50  0.00  B   299  C
#HETATM 2249  N1  BYP A   1       1.850  30.748   6.638  0.50  0.00  B   299  N
#HETATM 2250  C2  BYP A   1       0.745  31.369   7.121  0.50  0.00  B   299  C
#HETATM 2251  C6  BYP A   1       2.074  29.484   7.035  0.50  0.00  B   299  C
#HETATM 2252  C4  BYP A   1       0.194  29.368   8.379  0.50  0.00  B   299  C
#HETATM 2253  C3  BYP A   1      -0.145  30.702   8.026  0.50  0.00  B   299  C
#HETATM 2254  N5  BYP A   1       1.294  28.803   7.878  0.50  0.00  B   299  N
#:
  for line in f:
    if (line.startswith('HETATM') or line.startswith('ATOM')):
      DAT['lig3name']      = line[17:20]
      DAT['auth_asym_id']  = line[21:22]
      DAT['label_seq_id']  = line[25:26]
      DAT['label_asym_id'] = line[68:69]
      DAT['auth_seq_id']   = line[71:75]
  f.close()





def perform_kcombu_and_get_key_values(command,valkey):
#>>Example of -KV output of pkcombu<<
# NHEAVYATOM_A 28
# NHEAVYATOM_B 26
# NATOM_MATCH  15
# TANIMOTO     0.384615
# RMSD         0.194683
  f = os.popen(command)
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#')==0) and (len(field)>=2):
      key   = field[0]
      value = field[1]
      valkey[key] = value
    pass
  f.close()


############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'
OPT['idlig']  = 'SupLIG'
OPT['idpro']  = 'SupPRO'

OPT['A'] = 'F'
OPT['of'] = '-'
OPT['proname'] = 'NEU'

if (len(sys.argv)<2):
  print "mkliginfo.py <options>"
  print " for making ligand information summary text file."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -L       : list of library molecules[%s]"%(OPT['L'])
  print " -idlig   : input dir for ligands  in PDB [%s]"%(OPT['idlig'])
  print " -idpro   : input dir for proteins in PDB [%s]"%(OPT['idpro'])
  print " -A       : Action (T or F) [%s]"%(OPT['A'])
  print " -of      : output summary file [%s]"%(OPT['of']) 
  print " -proname : name of receptor protein [%s]"%(OPT['proname']) 
  sys.exit(1)


### [1] read option ###
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_summary_file() --> '%s'"%(OPT['of'])

of.write("## DATA LIST FOR SUPERPOSITION OF MOLECULES USED IN THE PAPER OF FKCOMBU (2014).##\n")
of.write("## provided by Takeshi Kawabata (kawabata@protein.osaka-u.ac.jp)                ##\n")
of.write("#COMMAND '%s'\n"%(OPT['COMMAND']))
of.write("#DATE    '%s'\n"%(OPT['START_DATE']))
of.write('#[1:receptor_protein_name]\n')
of.write('#[2:ligand_filename] [3:PDB_ID] [4:ligand 3-letter name] [5:ligand author_asym_id(chain_id)]\n')
of.write('#[6:ligand author_seq_id] [7:ligand label_asym_id]\n')
of.write('#[8:protein_filename] [9:protein author_asym_id(chain_id)] [10:protein label_asym_id]\n')
of.write('#[11:Nheavyatom] [12:Nrotbond]\n')
liglist = []
receptor = {}
read_list_file(OPT['L'],liglist,receptor)


for lig in (liglist):
  iligfile = "%s/%s"%(OPT['idlig'],lig)
  iprofile = "%s/%s"%(OPT['idpro'],receptor[lig])
  field = lig.split('_')
  lig3name = field[0] 
  pdbid    = field[1] 
  LIG = {}
  read_ligand_info_in_PDB(iligfile,LIG)
  PRO = {}
  read_ligand_info_in_PDB(iprofile,PRO)
 # print iligfile,LIG
  valkey = {}
  command = "fkcombu -T %s -R %s -S N -KV T"%(iligfile,iligfile)
  perform_kcombu_and_get_key_values(command,valkey)
  of.write("%s "%(OPT['proname']))
  of.write("%s %4s %s %s %s %s %s %s %s"%(lig,pdbid,lig3name,LIG['auth_asym_id'],LIG['auth_seq_id'],LIG['label_asym_id'],receptor[lig],PRO['auth_asym_id'],PRO['label_asym_id']))
  of.write(" %s %s"%(valkey['NHEAVYATOM_T'],valkey['NROTBOND_T']))
  of.write("\n")
if (OPT['of'] != '-'):
  of.close()
sys.exit(1)
