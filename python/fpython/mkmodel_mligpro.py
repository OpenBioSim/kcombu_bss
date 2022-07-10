#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "June 18, 2012"


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
      sws = rnum = chain = ''
      if (len(field)>1):
        sws = field[1]
      if (len(field)>2):
        rnum = field[2]
      if (len(field)>3):
        chain = field[3]
      prop[field[0]] = {} 
      prop[field[0]]['SWS']   = sws
      prop[field[0]]['RNUM']  = rnum
      prop[field[0]]['CHAIN'] = chain 
      #print "'%s' '%s' '%s' '%s' '%s'"%(line,field[0],sws,rnum,chain)
  print "#len(ligpdblist):%d"%(len(ligpdblist))
  f.close() 


def mk_pdbfile(pdbch):
  pdb   = pdbch[0:len(pdbch)-1]
  dir   = pdbch[1:3]
  ch    = pdbch[len(pdbch)-1:len(pdbch)]
  pdbfile = OPT['idpdb'] + "/"+dir+"/pdb"+pdb+".ent"
  return(pdbfile)


def extract_specific_ligand_from_PDB_file(ipdbfile,opdbfile,ligname,rnumlig,chainlig):
  print "#extract_specific_ligand_from_PDB_file('%s',name '%s' rnum '%s' chain '%s')-->'%s'"%(ipdbfile,ligname,rnumlig,chainlig,opdbfile)
  if not os.access(ipdbfile,os.R_OK):
    print "#ERROR:Can't open ipdbfile '%s'" %(ipdbfile)
    sys.exit(1)
#          1         2
#0123456789012345678901234567890123456789
#HETATM 8689  C2  N76 A1298      14.298   4.174  -6.086  1.00  0.00           C   pdb1oiu A
#HETATM 8690  C8  N76 A1298      14.626  -0.177  -5.242  1.00  0.00           C   pdb1oiu A
#HETATM 8691  C10 N76 A1298      14.003   4.616  -2.053  1.00  0.00           C   pdb1oiu A
#HETATM 8692  C11 N76 A1298      13.519   4.455  -0.627  1.00  0.00           C   pdb1oiu A
#HETATM 8693  C12 N76 A1298      13.532   5.831   0.036  1.00  0.00           C   pdb1oiu A
  extract_lines = []
  header_lines  = []
  f = open(ipdbfile)
  for line in f:
    if (line.startswith('ATOM  ')) or (line.startswith('HETATM')):
      resname = line[17:20]
      chain   = line[21]
      rnum    = line[22:26].replace(' ','')
      if (resname == ligname) and ((chain == chainlig)or (chainlig=='')) and ((rnum==rnumlig)or(rnumlig=='')):
        extract_lines.append(line)
    if (line.startswith('REMARK    COMMENT ')) or (line.startswith('REMARK    DATE ')) or (line.startswith('REMARK    COMMAND ')):
      header_lines.append(line)

  f.close()

  if (len(extract_lines)>0):
    of = open(opdbfile,"w")
    for line in (header_lines):
      of.write("%s"%(line))
    for line in (extract_lines):
      of.write("%s"%(line))
    of.close()
  return(len(extract_lines))
 



def extract_line_with_specific_header(ifname,header_str):
  print "#extract_line_with_specific_header('%s','%s'):"%(ifname,header_str)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open ipdbfile '%s'" %(ifname)
  f = open(ifname)
  for line in f:
    if (line.startswith(header_str)):
      f.close()
      return(line.rstrip('\n'))
  f.close()
  return("")
  pass



############
### MAIN ###
############

OPT = {}

OPT['div'] = '0/1'
OPT['rms'] = 'T'
OPT['stp'] = 'T'
OPT['chr'] = 'T'
OPT['rng'] = 'T'
OPT['gra'] = 'T'
OPT['mcs'] = 'T'


OPT['con'] = 'C'
OPT['mtd'] = '-1'

OPT['nini'] = '10'
OPT['nout'] = '1'
OPT['weam'] = 1.0
OPT['wesc'] = 1.0
OPT['werc'] = 1.0
OPT['wevo'] = 0.0
OPT['aevo'] = 'I' 

OPT['fix'] = 'F'
OPT['rgd'] = 'F'
OPT['pca'] = 'F'
OPT['ich'] = 'F'

OPT['at'] = 'K'


OPT['ilig']     = ''
OPT['iapo']     = ''
OPT['imol2']    = 'LIG_MOL2'
OPT['idpdb']    = '/DB/PDBv3'
OPT['otmp']     = 'tmpout'
OPT['olQonP0']  = 'LIG_TMP'
OPT['olPmdl']   = 'LIG_TARmdl'
OPT['olPcor']   = 'LIG_TARcor'
OPT['of']       = 'summary.out'
OPT['A']   = 'F'
OPT['EVA'] = 'F'
OPT['avex'] = 'F'


if (len(sys.argv)<=2):
  print "mkmodel_mprolig.py <options>"
  print " for making ligand-protein complex model using 'MATRAS' and 'fkcombu' program."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " P:target, Q:template"
  print "<options>"
  print " -ilig    : input list of lig_pdb swiss [%s]"%(OPT['ilig'])
  print " -iapo    : input list of apo pdb swiss [%s]"%(OPT['iapo'])
  print " -idpdb   : input dir for PDB [%s]"%(OPT['idpdb'])
  print " -imol2   : input directory for 3D MOL2 molecules[%s]"%(OPT['imol2'])
  print " -olQonP0 : outdir. for superimposed ligand of templateQ on apo targetP[%s]"%(OPT['olQonP0'])
  print " -olPmdl  : outdir. for modeled target ligandP on apo targetP[%s]"%(OPT['olPmdl'])
  print " -olPcor  : outdir. for correct target ligandP (superimposed on apo target)[%s]"%(OPT['olPcor'])
  print " -otmp    : outdir. for temporary files [%s]"%(OPT['otmp'])
  print " -div     : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -A       : Action (T or F) [%s]"%(OPT['A'])
  print " -EVA     : Action only for evaluation (T or F) [%s]"%(OPT['EVA'])
  print " -avex    : avoid writing, if file already exists. (T or F)[%s]"%(OPT['avex']) 
  print " -of      : output summary file [%s]"%(OPT['of'])
  print " -h       : detailed help"
 
  if (len(sys.argv)==2) and (sys.argv[1]=='-h'): 
    print "<options for fkcombu>"
    print " -rms : rmsd-min rigid-body fitting of molA onto molB (T or F) [%s]"%(OPT['rms'])
    print " -rng : Stamping ring structure (T or F) [%s]"%(OPT['rng'])
    print " -chr : Chilarity fit of molA into molB (T or F) [%s]"%(OPT['chr'])
    print " -stp : Stamping rotatable bond into molA from molB (T or F) [%s]"%(OPT['stp'])
    print " -gra : Gradient-based fitting (T or F) [%s]"%(OPT['gra'])
    print " -fix  : fixing the conformation, not do any transformations (T or F)[%s]"%(OPT['fix'])
    print " -rgd  : rigid transform.in other words, just '-rms T'. (T or F)[%s]"%(OPT['rgd'])
    print " -mcs  : Do MCS. ('T' or 'F'). If 'F', do not use any atom match [%s]"%(OPT['mcs'])
    print " -pca  : PCA-Evolmovlap fitting (T or F) [%s]"%(OPT['pca'])
    print " -ich  : initial chirarity change ('T' or 'F') [%s]"%(OPT['ich'])
   
    print ""
    print " -nini : number of random initial conformation [%s]"%(OPT['nini'])
    print " -nout : number of output conformations [%s]"%(OPT['nout'])
    print " -weam : Weight for Eatommatch [%s]"%(OPT['weam'])
    print " -wesc : Weight for Eselfcrash [%s]"%(OPT['wesc'])
    print " -werc : Weight for Erecpcrash [%s]"%(OPT['werc'])
    print " -wevo : Weight for Evolmovlap [%s]"%(OPT['wevo'])
    print " -aevo : atom match type for Evolmovlap. 'I'gnore,'C'are [%s]"%(OPT['aevo'])
    print " -con : Connectivity.'D'isconnected,'C'onnected, 'S'ubstructure(A<=B), 'I'somorphic (A=B),"
    print "      :        'T'opo_constrained_disconnected  't'opo_constrained_connected)[%s]"%(OPT['con'])
    print " -mtd : max difference of topological distance(num of bonds in the shortest path).<0:don't care. [%s]"%(OPT['mtd'])
   
    print "-at  : atomtype classification. 'X':don't care 'E'lement(C,N,O) ele+'R'ing (C,C@,N,N@)"
    print "     :   ele+'B'ond_num (C1,C2,C3) 'T':ele+bond+ring (C2,C2@,C3,C3@)"
    print "     :   'K':combu-recommend(C,C@,O,O@,O1,N,N@,N1,S,S@,P,X) 'D'abrl (D,A,B,R,L) 'k'cf [%s]"%(OPT['at'])
  sys.exit(1)


### [0] read option ###
PID = os.getpid()
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

if (OPT['fix']=='T'):
  OPT['rms'] = OPT['rng'] = OPT['chr'] = OPT['stp'] = OPT['gra'] = 'F'
if (OPT['rgd']=='T'):
  OPT['rms'] = 'T'
  OPT['rng'] = OPT['chr'] = OPT['stp'] = OPT['gra'] = 'F'



### [1] read lists of ligand_pdb_sws files & apo_pdb_sws file ###
ligpdblist = []
ligpdb_prop = {}
read_ligand_pdb_sws_list(OPT['ilig'],ligpdblist,ligpdb_prop)

apopdblist = []
apopdb_prop = {}
read_ligand_pdb_sws_list(OPT['iapo'],apopdblist,apopdb_prop)

apopdb_frm_sws = {}
for x in (apopdblist):
  apopdb_frm_sws[apopdb_prop[x]['SWS']] = x 


#### [2] MAKING LIST OF TARGET-TEMPLATE PAIR ###
TAR_TEMP_PAIR = []

## P:target, Q:template
 
Nmodel = 0
for ligpdbP in (ligpdblist):
  (ligP,pdbP) = ligpdbP.split('_')
  swsP    = ligpdb_prop[ligpdbP]['SWS']
#  print ligpdbP,swsP
  if (apopdb_frm_sws.has_key(swsP)):
    apoP = apopdb_frm_sws[swsP]
    for ligpdbQ in (ligpdblist):
      swsQ    = ligpdb_prop[ligpdbQ]['SWS']
      (ligQ,pdbQ) = ligpdbQ.split('_')
      pair = {}
      pair['TAR_LIG']     = ligP
      pair['TAR_PDB']     = pdbP
      pair['TAR_LIG_PDB'] = ligpdbP
      pair['TAR_APO_PDB'] = apoP
      pair['TAR_SWS']     = swsP
      pair['TMP_LIG']     = ligQ 
      pair['TMP_LIG_PDB'] = ligpdbQ 
      pair['TMP_PDB']     = pdbQ 
      pair['TMP_SWS']     = swsQ 
      oligPmodel = "%s/%s_on_%s_by_%s_%s"%(OPT['olPmdl'],ligP,apoP,ligQ,pdbQ)
      if (OPT['avex']!='T') or (os.path.isfile(oligPmodel)==0):
        TAR_TEMP_PAIR.append(pair) 
        Nmodel += 1

print "#Nmodel %d %d"%(Nmodel, len(TAR_TEMP_PAIR))
Nsta = Nmodel*bunshi/bunbo
Nend = Nmodel*(bunshi+1)/bunbo
print "#bunshi %d bunbo %d Nmodel %d Nsta %d Nend %d"%(bunshi,bunbo,Nmodel,Nsta,Nend)


#### [2] Do superimposition one by one ###
if (OPT['A']=='T') or (OPT['EVA']=='T'):
  of = open(OPT['of'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#bunshi %d bunbo %d Nmodel %d Nsta %d Nend %d\n"%(bunshi,bunbo,Nmodel,Nsta,Nend))
  of.write("#[TANIMOTOligP_ligQ:1] [SqIDproP_proQapo:2] [RMSDlig_TAR_TARcor:3]\n")
  of.write("#[ligTAR:4] [pdbTAR:5] [apoTAR:6] [sws:7] [ligTMP:8] [pdbTMP:9]\n")
  of.write("#[RMSDproP_proQapo:10] [RDISproP_proQapo:11] [SqIDproP_proPapo:12] [RMSDproP_proPapo:13]\n")

for i in range(Nsta,Nend):
  pair = TAR_TEMP_PAIR[i] 
  ## 'P':TAR:target, 'Q':TMP:template ##
  ligP     = pair['TAR_LIG']
  pdbP     = pair['TAR_PDB']
  ligpdbP  = pair['TAR_LIG_PDB']
  apoP     = pair['TAR_APO_PDB']
  swsP     = pair['TAR_SWS']
  ligQ     = pair['TMP_LIG']
  pdbQ     = pair['TMP_PDB']
  ligpdbQ  = pair['TMP_LIG_PDB']
  swsQ     = pair['TMP_SWS']

  print ">>>[%d/%d] target lig %s pdb %s apo %s sws %s template lig %s pdb %s sws %s"%(i,Nmodel,ligP,pdbP,apoP,swsP,ligQ,pdbQ,swsQ)

  ## [(1): matras for prediction] Superimpose pdbQ on apoP, and extract superimposed ligQ on apoP 
  otmppdbfile = "%s/%d.pdb"%(OPT['otmp'],PID)
  comstr = "Matras P -A %s -Ac %s -B %s -Bc %s -opdbA %s -OA T"%(mk_pdbfile(pdbQ),pdbQ[4],mk_pdbfile(apoP),apoP[4],otmppdbfile)
  oligQonPapo = "%s/%s_%s_on_%s"%(OPT['olQonP0'],ligQ,pdbQ,apoP)
  print "#%s"%(comstr)
  print "#extract_specific_ligand_from_PDB_file('%s' -->'%s','%s','%s','%s')"%(otmppdbfile,oligQonPapo,ligQ,ligpdb_prop[ligpdbQ]['RNUM'],ligpdb_prop[ligpdbQ]['CHAIN'])

  if (OPT['A'] == 'T') and ((OPT['avex']!='T') or (os.path.isfile(oligQonPapo)==0)):
    os.system(comstr)
    okstep1 = extract_specific_ligand_from_PDB_file(otmppdbfile,oligQonPapo,ligQ,ligpdb_prop[ligpdbQ]['RNUM'],ligpdb_prop[ligpdbQ]['CHAIN'])
    if (okstep1==0):
       print "#EXTRACT_ERROR!! ('%s' -->'%s','%s','%s','%s')"%(otmppdbfile,oligQonPapo,ligQ,ligpdb_prop[ligpdbQ]['RNUM'],ligpdb_prop[ligpdbQ]['CHAIN'])
  if (OPT['A'] == 'T') or (OPT['EVA']=='T'):
    line = extract_line_with_specific_header(oligQonPapo,"REMARK    COMMENT ")
    #REMARK    COMMENT proA pdb3a4o X proB pdb2zv7 A Ncomp 255 SeqID 95.69 % CRMS 1.387 Rdis 85.68 %
    field = line.split()
    SqIDproP_proQapo = field[11]
    RMSDproP_proQapo = field[14]
    RDISproP_proQapo = field[16]

  ### [(2) fkcombu for prediction] Transform ligP (in mol2) onto the superimposed ligQ to generate ligP model 
  oligPmodel = "%s/%s_on_%s_by_%s_%s"%(OPT['olPmdl'],ligP,apoP,ligQ,pdbQ)
  otmpkcombufile = "%s/%d.kmb"%(OPT['otmp'],PID)
  comstr = "fkcombu -A %s/%s.mol2 -B %s -fB P -con %s -mtd %s -at %s -rms %s -rng %s -chr %s -stp %s -gra %s -mcs %s -pca %s -nini %s -nout %s -weam %s -wesc %s -werc %s -wevo %s -aevo %s -ich %s -opA %s > %s"%(\
          OPT['imol2'],ligP,oligQonPapo,OPT['con'],OPT['mtd'],OPT['at'],\
          OPT['rms'],OPT['rng'],OPT['chr'],OPT['stp'],OPT['gra'],OPT['mcs'],OPT['pca'],\
          OPT['nini'],OPT['nout'],\
          OPT['weam'],OPT['wesc'],OPT['werc'],OPT['wevo'],OPT['aevo'],OPT['ich'],oligPmodel,otmpkcombufile)
  print "#'%s'"%(comstr)
  if (OPT['A'] == 'T') and ((OPT['avex']!='T') or (os.path.isfile(oligPmodel)==0)):
    os.system(comstr)
  if (OPT['A'] == 'T') or (OPT['EVA']=='T'):
    line = extract_line_with_specific_header(oligPmodel,"REMARK COMMENT NhvatmA")
    #REMARK COMMENT NhvatmA 22 NhvatmB 22 Natmmatch 22 tanimoto 1.000000 Etotal +73.44 Eatommatch +73.37 Eselfcrash +0.07 Ercptcrash +0.00 rmsd_match 2.58
    if (line != ''):
      field = line.split()
      TANIMOTOligP_ligQ = "%.2f"%(100.0*float(field[9]))
      okstep2 = 1
    else:
      okstep2 = 0

  ## [(3): Matras for evaluation] MATRAS apoP and pdbP, and extract correct superimposed ligP
  otmppdbfile = "%s/%d.pdb"%(OPT['otmp'],PID)
  comstr = "Matras P -A %s -Ac %s -B %s -Bc %s -opdbA %s -OA T"%(mk_pdbfile(pdbP),pdbP[4],mk_pdbfile(apoP),apoP[4],otmppdbfile)
  oligPonPapo = "%s/%s_on_%s_by_%s_%s"%(OPT['olPcor'],ligP,apoP,ligP,pdbP)
  print "#%s"%(comstr)
  print "#extract_specific_ligand_from_PDB_file('%s' -->'%s','%s','%s','%s')"%(otmppdbfile,oligPonPapo,ligP,ligpdb_prop[ligpdbP]['RNUM'],ligpdb_prop[ligpdbP]['CHAIN'])
  if (OPT['A'] == 'T') and ((OPT['avex']!='T') or (os.path.isfile(oligPonPapo)==0)):
    os.system(comstr)
    extract_specific_ligand_from_PDB_file(otmppdbfile,oligPonPapo,ligP,ligpdb_prop[ligpdbP]['RNUM'],ligpdb_prop[ligpdbP]['CHAIN'])
  if (OPT['A'] == 'T') or (OPT['EVA']=='T'):
    line = extract_line_with_specific_header(oligPonPapo,"REMARK    COMMENT ")
    #REMARK    COMMENT proA pdb3a4o X proB pdb2zv7 A Ncomp 255 SeqID 95.69 % CRMS 1.387 Rdis 85.68 %
    #REMARK    COMMENT proA pdb2cgu A proB pdb1ia8 A Ncomp 263 SeqID 100.00 % CRMS 0.998 Rdis 92.11 %
    field = line.split()
    SqIDproP_proPapo = field[11]
    RMSDproP_proPapo = field[14]

  ### [(4): fkcombu for evaluation] compare correct superimposed ligP and ligP model 
  comstr = "fkcombu -A %s -B %s -fA P -fB P -fix T > %s"%(oligPmodel,oligPonPapo,otmpkcombufile)
  print "#'%s'"%(comstr)
  if (OPT['A'] == 'T') or (OPT['EVA']=='T'):
    os.system(comstr)
    line = extract_line_with_specific_header(otmpkcombufile,"COMMENT NhvatmA")
    #COMMENT NhvatmA 22 NhvatmB 22 Natmmatch 22 tanimoto 1.000000 Etotal +73.44 Eatommatch +73.37 Eselfcrash +0.07 Ercptcrash +0.00 rmsd_match 2.58
    if (line != ''):
      field = line.split()
      RMSDlig_TAR_TARcor = field[18]
      if (OPT['A']=='T') or (OPT['EVA']=='T'):
        of.write("%s %s %s %s %s %s %s %s %s %s %s %s %s\n"%(TANIMOTOligP_ligQ,SqIDproP_proQapo,RMSDlig_TAR_TARcor,ligP,pdbP,apoP,swsP,ligQ,pdbQ,RMSDproP_proQapo,RDISproP_proQapo,SqIDproP_proPapo,RMSDproP_proPapo))


if (OPT['A']=='T') or (OPT['EVA']=='T'):
  now = datetime.now()
  OPT['END_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  of.write("#bunshi %d bunbo %d Nmodel %d Nsta %d Nend %d\n"%(bunshi,bunbo,Nmodel,Nsta,Nend))
  of.write("#START_DATE    %s\n"%(OPT['START_DATE']))
  of.write("#END_DATE      %s\n"%(OPT['END_DATE']))
  print "#output_result() --> '%s'"%(OPT['of'])
  of.close()
