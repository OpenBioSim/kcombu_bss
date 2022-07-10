#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "Apr 5, 2014"


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


def pickup_focused_line(ifname,pattern_head):
  if (os.access(ifname,os.R_OK)==0):
    return('')
  f = open(ifname)
  for line in f:
    if (line.startswith(pattern_head)):
      line = line.rstrip('\n')
      return(line)
  f.close()
  #print "#WARNING:no line for pattern_head '%s' in the file '%s'"%(pattern_head,ifname)
  return('')

def count_Nheavyatom_in_PDBfile(ipdbfile):
  f = open(ipdbfile)
#          1         2         3
#0123456789012345678901234567890123456789
#HETATM   20  O   MOL A   1      27.169  14.191  62.446  0.00  0.00  O1  O x - A 1 0 0 0 0 0   1   3   4   0   0
#HETATM   21  H   MOL A   1      22.638  17.272  63.699  0.00  0.00      H x - H 1 0 0 0 0 0   0   0   0   0   0
#HETATM 4008  C10 EQP A 500      27.445  14.640  61.896  0.95  0.00
#HETATM 4010  N5  EQP A 500      26.556  15.627  61.837  0.77  0.00
#HETATM 4013  O3P EQP A 500      22.655  20.364  65.637  1.00  0.00
#HETATM 4019  O10 EQP A 500      28.023  14.307  62.949  0.83  0.00
#HETATM 4020  H2  EQP A 500      24.274  19.039  63.113  1.00  0.00
#HETATM 4021  H31 EQP A 500      23.256  17.174  64.524  1.00  0.00
#HETATM 4022  H32 EQP A 500      24.901  16.973  65.199  1.00  0.00
#HETATM 4030 H111 EQP A 500      28.343  14.492  59.955  1.00  0.00
#HETATM 4031 H112 EQP A 500      26.751  13.647  60.120  1.00  0.00
#HETATM 4032 H113 EQP A 500      28.237  12.939  60.857  1.00  0.00
  Nheavyatom = 0
  for line in f:
    if (line.startswith('HETATM') or line.startswith('ATOM')):
      atomstr = line[12:16].replace(' ','')
      if (atomstr[0]!='H'):
        Nheavyatom += 1
      pass
  f.close()
  return(Nheavyatom)




def perform_kcombu_and_get_key_values(command,valkey):
#>>Example of -KV output of pkcombu<<
# NHEAVYATOM_A 28
# NHEAVYATOM_B 26
# NATOM_MATCH  15
# TANIMOTO     0.384615
# RMSD         0.194683
#
#>>Example of -KV output of fkcombu<<
# NHEAVYATOM_T 28
# NHEAVYATOM_R 28
# NATOM_MATCH  28
# TANIMOTO     1.000000
# Etotal     1135.069824
# Eatommatch 1133.180298
# Eselfcrash 0.188955
# Ercptcrash 0.000000
# Evolmovlap -418.469696
# RMSD_MATCH 8.953962
# TANIMOTO_VOL 0.219865
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

def read_REMARK_COMMENT_in_pdb(ifname,valkey):
#          1         2         3
#0123456789012345678901234567890123456789
#REMARK COMMENT: [rank 1] ninit 3 mcs T1K Nmatch 9 tani 0.1837 Etotal +9.92 Eam
#REMARK COMMENT: +9.60 Esc +0.31 Epc +0.00 Evo -524.02 rmsd_match 1.745 tani_vol
#REMARK COMMENT: 0.3014
  comment = ''
  if (os.access(ifname,os.R_OK)==0):
    return('')
  f = open(ifname)
  for line in f:
    if (line.startswith('REMARK COMMENT')):
      line = line.rstrip('\n')
      line = line[16:]
      if (comment != ''):
        comment += ' ' 
      comment += line
  f.close()

  ## translating 'REMARK COMMENT' line ##
  fields = comment.split()
  n = 2
  while (n < len(fields)):
    key   = fields[n]
    value = fields[n+1]
    valkey[key] = value
    n += 2
 
  #print valkey 

  return(comment)




############
### MAIN ###
############

OPT = {}
OPT['L']      = 'ligand.list'
OPT['idref']  = 'SupLIG'

OPT['A'] = 'F'
OPT['idtarp']  = 'tmpout'
OPT['self'] = 'F'
OPT['permsd'] = 'T'
OPT['fT'] = 'P'
OPT['tailm'] = ''
OPT['of'] = '-'


if (len(sys.argv)<2):
  print "sum3Dmodel.py <options>"
  print " for making summary of accuracry of predicted simimilariy-based 3D model."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ** The summary is shown in 'stdout'" 
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idref  : input dir for mols with 'correct' 3D conf in PDB [%s]"%(OPT['idref'])
  print " -idtarp : input dir for predicted target 3D model [%s]"%(OPT['idtarp'])
  print " -fT     : file type for 3D model for target. 'P':pdb,'2':mol2 [%s]"%(OPT['fT'])
  print " -tailm  : tail of file for 3D model. [%s]"%(OPT['tailm'])
  print " -permsd : Permutations for getting minimum RMSD. 'T':permute only target,'B'oth_permute, 'F'alse [%s]"%(OPT['permsd'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print " -self   : only self modeling (T or F)[%s]"%(OPT['self'])
  print " -of     : output summary file [%s]"%(OPT['of']) 
  sys.exit(1)


### [1] read option ###
PID = os.getpid()
read_option(sys.argv,OPT)

### [2] read liglist (-L) ###
if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write_summary_file() --> '%s'"%(OPT['of'])

of.write("#COMMAND '%s'\n"%(OPT['COMMAND']))
of.write("#DATE    '%s'\n"%(OPT['START_DATE']))
of.write("#[ligA:1] [ligB:2] [NheavyA:3] [NheavyB:4] [NheavyB-NheavyA:5]\n")
of.write("#[rmsd:6] [volume_overlap:7] [tanimoto_volume:8]\n")
of.write("#[rmsd_match:9]\n")
liglist = []
read_list_file(OPT['L'],liglist)
print "#PID %d"%(PID)

for ligA in (liglist):
  fieldA = ligA.split('_')
  ligname3A = fieldA[0]
  for ligB in (liglist):
    (ligheadB,tail) = ligB.split('.')
    fieldB = ligheadB.split('_')
    ligname3B = fieldB[0]
    print "ligheadB '%s' tail '%s'"%(ligheadB,tail)
    if ((OPT['self']!='T') or (ligA == ligB)):
      if (OPT['tailm']==''):
        modelfile = "%s/%s_%s"%(OPT['idtarp'],ligname3A,ligB)
      else:
        #modelfile = "%s/%s_%s.%s"%(OPT['idtarp'],ligname3A,ligB,OPT['tailm'])
        modelfile = "%s/%s_%s.%s"%(OPT['idtarp'],ligname3A,ligheadB,OPT['tailm'])


      if (os.access(modelfile,os.R_OK)==0):
        if (OPT['tailm']==''):
          modelfile = "%s/%s/%s_%s"%(OPT['idtarp'],ligname3A,ligname3A,ligB)
        else:
          modelfile = "%s/%s/%s_%s.%s"%(OPT['idtarp'],ligname3A,ligname3A,ligheadB,OPT['tailm'])

        if (os.access(modelfile,os.R_OK)==0):
          print "#ERROR:Can't find modelfile '%s'"%(modelfile)
          sys.exit(1)


      command_rmsd = "fkcombu -T %s -fT %s -R %s/%s  -fR P -con C -permsd %s -S N -KV T -at E"%(modelfile,OPT['fT'],OPT['idref'],ligA,OPT['permsd'])
      command_vol  = "fkcombu -T %s -fT %s -R %s/%s  -fR P -con C -mcs F -S N -KV T"%(modelfile,OPT['fT'],OPT['idref'],ligB)

      print "#%s"%(command_rmsd)
      print "#%s"%(command_vol)
      mdat = {}
      commentline = read_REMARK_COMMENT_in_pdb(modelfile,mdat)

      print "#%s"%(commentline)
      if (OPT['A']=='T'):
        fdat = {}
        perform_kcombu_and_get_key_values(command_rmsd,fdat)

        vdat = {}
        perform_kcombu_and_get_key_values(command_vol,vdat)
        diffBA = int(vdat.get('NHEAVYATOM_R',0))-int(vdat.get('NHEAVYATOM_T',0))

        of.write("%s %s %s %s %2d"%(ligA,ligB,vdat.get('NHEAVYATOM_T',0),vdat.get('NHEAVYATOM_R',0),diffBA))
        of.write(" %s"%(fdat.get('RMSD_MATCH',0.0)))
        of.write(" %f %s"%(-float(vdat.get('Evolmovlap',0.0)),vdat.get('TANIMOTO_VOL',0.0)))
        of.write(" %s"%(mdat.get('rmsd_match',-1.0)))
        of.write("\n")


if (OPT['of'] != '-'):
  of.close()


