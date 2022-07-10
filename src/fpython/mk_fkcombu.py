#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "May 9, 2014"


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

def read_ligand_protein_list_file(ifname,list):
  print "#read_ligand_protein_list_file('%s')"%(ifname)
#>> FILE EXAMPLE 
# #[ligand_file] [protein_file]
# G20_2qwf_1_G_1.pdb 2qwf_1_A_1.pdb
# GNA_2qwe_1_G_1.pdb 2qwe_1_A_1.pdb
# DAN_2qwc_1_G_1.pdb 2qwc_1_A_1.pdb
# 49A_1f8e_1_H_1.pdb 1f8e_1_A_1.pdb
# 9AM_1f8d_1_H_1.pdb 1f8d_1_A_1.pdb
# G39_2qwh_1_E_1.pdb 2qwh_1_A_1.pdb
# SIA_2c4a_1_B_1.pdb 2c4a_1_A_1.pdb

  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR(read_ligand_protein_list_file): Can't open '%s'."%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0) and (len(line)>1):
      fields = line.split()
      ligand  = fields[0]
      protein = fields[1]
      dat = {}
      dat['lig'] = ligand
      dat['pro'] = protein
      list.append(dat)
  f.close()




############
### MAIN ###
############

OPT = {}
OPT['L']   = 'ligand.list'
OPT['idref'] = 'SupLIG'

OPT['A'] = 'F'
OPT['E'] = 'A'
OPT['S'] = 'F'
OPT['div'] = '0/1'
OPT['idsdf']  = 'SDF_3D'
OPT['idmol2'] = 'MOL2_3D'
OPT['rms'] = 'T'
OPT['stp'] = 'T'
OPT['chsp3'] = 'F'
OPT['stprng'] = 'T'
OPT['gra'] = 'T'
OPT['mcs'] = 'T'


OPT['con'] = 'C'
OPT['mtd'] = '-1'
OPT['odtar'] = 'tmpout'
OPT['avex'] = 'F'
OPT['self'] = 'F'
OPT['idpro'] = ''

OPT['nini'] = '10'
OPT['nout'] = '1'
OPT['weam'] = 1.0
OPT['wesc'] = 1.0
OPT['wepc'] = 1.0
OPT['wevo'] = 0.0
OPT['avo'] = 'I' 

OPT['fix'] = 'F'
OPT['rgd'] = 'F'
OPT['pca'] = 'F'
OPT['ich'] = 'F'

OPT['at'] = 'K'
OPT['protype'] = 'R';
OPT['amste3D'] = 'T'
OPT['damrl'] = 'H'
OPT['fR'] = 'P'
OPT['nk'] = '100'
OPT['mxam'] = '1'
OPT['sort'] = 'E'
OPT['mcs0'] = ''
OPT['mcs1'] = ''
OPT['mcs2'] = ''
OPT['mcs3'] = ''
OPT['mcs4'] = ''
OPT['mcs5'] = ''
OPT['mcs6'] = ''
OPT['mcs7'] = ''
OPT['mcs8'] = ''
OPT['mcs9'] = ''
OPT['w']   = 'D'
OPT['wtd'] = '0.0'
OPT['subdir'] = 'F'
OPT['fixstp'] = 'R'



if (len(sys.argv)<2):
  print "mk_fkcombu.py <options>"
  print " for making similariy-based 3D model using 'fkcombu' program."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -idref  : input directory for template reference molecules[%s]"%(OPT['idref'])
  print " -fR     : type of reference molecule 'P'db, '2':mol2, 'S'df [%s]"%(OPT['fR']) 
  print " -idsdf  : input directory for target 3D SDF  molecules[%s]"%(OPT['idsdf'])
  print " -idmol2 : input directory for target 3D MOL2 molecules[%s]"%(OPT['idmol2'])
  print "           (One of '-idsdf' and '-isdmol2' should be assigned.)"
  print " -odtar  : output directory for predicted (transformed) 3D model in PDB format [%s]"%(OPT['odtar'])
  print " -subdir : use target ligand name as output subdirectory. (T or F)[%s]"%(OPT['subdir'])
  print " -idpro  : input directory for receptor protein PDB file (optional) [%s]"%(OPT['idpro'])
  print " -protype: type of receptor protein. use 'R'eference_protein, 'T'arget_protein [%s]"%(OPT['protype'])
  print " -div    : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print " -avex   : avoid existing file pairs in '-odtar' (T or F)[%s]"%(OPT['avex'])
  print " -self   : only self modeling (T or F)[%s]"%(OPT['self'])
  print "<options for general stargety of transformation>"
  print " -E    : Energy-strategy.'A'tom-match, 'V'olume-overlap, 'B'oth, 'X':by detailed options.[%s]"%(OPT['E'])
  print " -S    : Search-strategy.'F'lexible-refine,'R'igid-refine, 'f'lexible-fast,'r'idid-fast 'N':do nothing 'X':by detailed options [%s]"%(OPT['S'])
  print "<options for fkcombu>"
  print " -rms    : rmsd-min rigid-body fitting of molA onto molB (T or F) [%s]"%(OPT['rms'])
  print " -stprng : Stamp 5- or 6-ring conformations (T or F) [%s]"%(OPT['stprng'])
  print " -chsp3 : Chilarity fit of molA into molB (T or F) [%s]"%(OPT['chsp3'])
  print " -stp   : Stamping rotatable bond into molA from molB (T or F) [%s]"%(OPT['stp'])
  print " -gra   : Gradient-based fitting (T or F) [%s]"%(OPT['gra'])
  print " -mcs   : Do MCS. ('T' or 'F'). If 'F', do not use any atom match [%s]"%(OPT['mcs'])
  print " -pca   : PCA-Evolmovlap fitting (T or F) [%s]"%(OPT['pca'])
  print " -ich   : initial chirarity change ('T' or 'F') [%s]"%(OPT['ich'])
  print ""
  print " -amste3D: adjust atom match to agree 3D chirality by permutations (T or F)[%s]"%(OPT['amste3D'])
  print " -sort  : sort value for multiple target 3D models. 'E':Etotal, 'V':tanimoto_volume [%s]"%(OPT['sort'])
  print " -mcs[x]: for [x]-th type for multiple MCS matches (x=0,1,2,..). [con]:[mtd]:[at]:"
  print " -nini  : number of random initial conformation [%s]"%(OPT['nini'])
  print " -nout  : number of output conformations [%s]"%(OPT['nout'])
  print " -weam  : Weight for Eatommatch [%s]"%(OPT['weam'])
  print " -wesc  : Weight for Eselfcrash [%s]"%(OPT['wesc'])
  print " -wepc  : Weight for Eprotcrash [%s]"%(OPT['wepc'])
  print " -wevo  : Weight for Evolmovlap [%s]"%(OPT['wevo'])
  print " -mxam  : maximum number of atom matches for the 3D-modelling [%s]"%(OPT['mxam'])
  print " -avo   : atom match type for Evolmovlap. 'I'gnore,'C'are [%s]"%(OPT['avo'])
  print " -fixstp: fixing rule for stamped bond for steepest descent initials. 'F':fix. 'R'otatable. [%s]"%(OPT['fixstp'])
  print "<options for MCS>"
  print " -con   : Connectivity.'D'isconnected,'C'onnected, 'S'ubstructure(A<=B), 'I'somorphic (A=B),"
  print "        :        'T'opo_constrained_disconnected  't'opo_constrained_connected)[%s]"%(OPT['con'])
  print " -mtd   : max difference of topological distance(num of bonds in the shortest path).<0:don't care. [%s]"%(OPT['mtd'])
 
  print " -at    : atomtype classification. 'X':don't care 'E'lement(C,N,O) ele+'R'ing (C,C@,N,N@)"
  print "        :   ele+'B'ond_num (C1,C2,C3) 'T':ele+bond+ring (C2,C2@,C3,C3@)"
  print "        :   'K':combu-recommend(C,C@,O,O@,O1,N,N@,N1,S,S@,P,X) 'D'abrl (D,A,B,R,L) 'k'cf [%s]"%(OPT['at'])
  print " -damrl :Assignment of DAMRL type.  '3':by 3D conf, 'H':by hydrogens, 'X':pdb->'3',other_format->'H' [%s]"%(OPT['damrl'])
  print " -nk    : Nkeep for build up [%s]"%(OPT['nk'])
  print " -w     : weighting scheme. 'D'efault weighting. 'S'pecified by user  [%s]"%(OPT['w'])
  print " -wtd   : weight_select_dis for topological distance       [%s]"%(OPT['wtd'])
  sys.exit(1)


### [1] read option ###
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

### [2] read liglist (-L) ###
ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)

### [3] count Npair ###
Npair = 0
for T in (ligpro_list):   ## target ##
  ligT = T['lig']
  fieldT = ligT.split('_')
  ligname3T = fieldT[0]
  for R in (ligpro_list): ## reference ##
    ligR = R['lig']
    fieldR = ligR.split('_')
    ligname3R = fieldR[0]

    if (OPT['subdir']=='T'):
      odir_target = "%s/%s"%(OPT['odtar'],ligname3T)
      outpdbfile = "%s/%s_%s"%(odir_target,ligname3T,ligR)
      if (os.path.isdir(odir_target)==0):
        print "subdir '%s' does not exist."%(odir_target)
        mkdir_command = "mkdir %s"%(odir_target)
        if (OPT['A']=='T'):
          os.system(mkdir_command)
    else:
      outpdbfile = "%s/%s_%s"%(OPT['odtar'],ligname3T,ligR)

    #outpdbfile = "%s/%s_%s"%(OPT['odtar'],ligname3T,ligR)
    if ((OPT['self']!='T') or (ligT == ligR)) and ((OPT['avex']=='F') or (os.access(outpdbfile,os.R_OK)==0)):
      Npair += 1

#Npair = Nlig * Nlig

Nsta = Npair*bunshi/bunbo
Nend = (Npair*(bunshi+1))/bunbo
print "#Nlig %d Npair %d bunshi %d bunbo %d Nsta %d Nend %d"%(Nlig,Npair,bunshi,bunbo,Nsta,Nend)
npair = 0

### [4] do fkcombu ###
  
for T in (ligpro_list):   ## target ##
  ligT = T['lig']
  fieldT = ligT.split('_')
  ligname3T = fieldT[0]
  for R in (ligpro_list): ## reference ##
    ligR = R['lig']
    fieldR = ligR.split('_')
    ligname3R = fieldR[0]
    (ligheadR,ligtailR) = ligR.split('.') 

    if (OPT['subdir']=='T'):
      odir_target = "%s/%s"%(OPT['odtar'],ligname3T)
      outpdbfile = "%s/%s_%s"%(odir_target,ligname3T,ligR)
      if (os.path.isdir(odir_target)==0):
        print "subdir '%s' does not exist."%(odir_target)
        mkdir_command = "mkdir %s"%(odir_target)
        if (OPT['A']=='T'):
          os.system(mkdir_command)  
    else:
      outpdbfile = "%s/%s_%s"%(OPT['odtar'],ligname3T,ligR)


    if ((OPT['self']!='T') or (ligT == ligR)) and ((OPT['avex']=='F') or (os.access(outpdbfile,os.R_OK)==0)):
      if (Nsta<=npair) and (npair<Nend):
        if (OPT['idmol2'] != ''):
          tarmolstr = "%s/%s.mol2"%(OPT['idmol2'],ligname3T)
        if (OPT['idsdf'] != ''):
          tarmolstr = "%s/%s.sdf"%(OPT['idsdf'],ligname3T)
        if (OPT['fR'] == 'P'):
          refmolstr = "%s/%s.pdb"%(OPT['idref'],ligheadR)
        if (OPT['fR'] == '2'):
          refmolstr = "%s/%s.mol2"%(OPT['idref'],ligheadR)
        if (OPT['fR'] == 'S'):
          refmolstr = "%s/%s.sdf"%(OPT['idref'],ligheadR)


        comstr = "fkcombu -T %s -R %s -fR %s -con %s -mtd %s -nk %s -E %s -S %s -at %s -rms %s -stprng %s -chsp3 %s -stp %s -gra %s -mcs %s -pca %s -nini %s -nout %s -weam %s -wesc %s -wepc %s -wevo %s -avo %s -ich %s -amste3D %s -opdbT %s"%(\
          tarmolstr,refmolstr,OPT['fR'],OPT['con'],OPT['mtd'],OPT['nk'],OPT['E'],OPT['S'],OPT['at'],\
          OPT['rms'],OPT['stprng'],OPT['chsp3'],OPT['stp'],OPT['gra'],OPT['mcs'],OPT['pca'],\
          OPT['nini'],OPT['nout'],\
          OPT['weam'],OPT['wesc'],OPT['wepc'],OPT['wevo'],OPT['avo'],OPT['ich'],OPT['amste3D'],outpdbfile)
        
        comstr = comstr + " -sort %s"%(OPT['sort']) 

        if (OPT['mxam'] != '1'):
          comstr = comstr + " -mxam %s"%(OPT['mxam']) 

        if (OPT['idpro'] != ''):
          if (OPT['protype']=='R'):
            comstr = comstr + " -P %s/%s -fP P"%(OPT['idpro'],R['pro']) 
          elif (OPT['protype']=='T'):
            comstr = comstr + " -P %s/%s -fP P"%(OPT['idpro'],T['pro']) 
          else:
            print "#ERROR: improper protype:'%s'"%(OPT['protype'])

        if (OPT['at']=='D'):
          comstr = comstr + " -damrl %s"%(OPT['damrl']) 
        if (OPT['w']=='S'):
          comstr = comstr + " -w S -wtd %s"%(OPT['wtd']) 
        if (OPT['fixstp']!='F') and (OPT['fixstp']!=''):
          comstr = comstr + " -fixstp %s"%(OPT['fixstp']) 

        if (OPT['mcs0']!=''):
          comstr = comstr + " -mcs0 %s"%(OPT['mcs0']) 
        if (OPT['mcs1']!=''):
          comstr = comstr + " -mcs1 %s"%(OPT['mcs1']) 
        if (OPT['mcs2']!=''):
          comstr = comstr + " -mcs2 %s"%(OPT['mcs2']) 
        if (OPT['mcs3']!=''):
          comstr = comstr + " -mcs3 %s"%(OPT['mcs3']) 
        if (OPT['mcs4']!=''):
          comstr = comstr + " -mcs4 %s"%(OPT['mcs4']) 
        if (OPT['mcs5']!=''):
          comstr = comstr + " -mcs5 %s"%(OPT['mcs5']) 
        if (OPT['mcs6']!=''):
          comstr = comstr + " -mcs6 %s"%(OPT['mcs6']) 
        if (OPT['mcs7']!=''):
          comstr = comstr + " -mcs7 %s"%(OPT['mcs7']) 
        if (OPT['mcs8']!=''):
          comstr = comstr + " -mcs8 %s"%(OPT['mcs8']) 
        if (OPT['mcs9']!=''):
          comstr = comstr + " -mcs9 %s"%(OPT['mcs9']) 

        print "#'%s'"%(comstr)

        if (OPT['A']=='T'):
          os.system(comstr)  
      npair += 1 
