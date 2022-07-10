#!/usr/bin/env python
import sys
import os
import kcombu

LastModDate = "June 14, 2011"

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


############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['idl'] = 'SupLIG'
OPT['alg'] = 'B'
OPT['nk']  = 40
OPT['at']  = 'K'
OPT['con']  = 'C'
OPT['A'] = 'F'
OPT['del'] = 'T'
OPT['fA'] = 'P'
OPT['fB'] = 'P'
OPT['doAm'] = 'AtomMatch'
OPT['div'] = '0/1'
OPT['mco']  = '-1'
OPT['mtd'] = '-1'
OPT['per'] = 'F'
OPT['ctd'] = 'F'
OPT['wna'] = '1.0'
OPT['wec'] = '1.0'
OPT['wtd'] = '1.0'
OPT['wco'] = '0.0'
OPT['xsec'] = '-1.0'
OPT['avex'] = 'F'
OPT['dosm'] = 'F'
OPT['w'] = 'S'
OPT['sbi'] = '100'
OPT['nrs'] = 'T'
OPT['xtd'] = '4'
OPT['ec'] = '2'

if (len(sys.argv)<2):
  print "mk_atommatch.py <options>"
  print " for making atom match file for all-vs-all comparison for a given molecule list."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -L    : list of library molecules[%s]"%(OPT['L'])
  print " -idl  : input directory for library molecules[%s]"%(OPT['idl'])
  print " -doAm : output directory for atom_maching filelibrary molecules[%s]"%(OPT['doAm'])
  print " -div  : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -A    : Action (T or F) [%s]"%(OPT['A'])
  print " -avex : avoid existing files in -doAm (T or F)[%s]"%(OPT['avex'])
  print " -dosm : do only existing atom-matching file with small (0 or 1) atom pair in -doAm (T or F)[%s]"%(OPT['dosm'])
  print "<options for target matching>"
  print " -alg : Algorithm. 'B'uild-up, 'P'ca-3D-fit '3'D(raw_fit)"
  print "       : 'X':exact match(Bron-Kerbosch-pivot) 'x':exact match(Bron-Kerbosch-ver1) [%s]"%(OPT['alg'])
  print " -nk   : Nkeep for build up [%d]"%(OPT['nk'])
  print " -at   : atomtype_classification. 'X':don't care (X) 'E'lement(C,N,O) ele+'R'ing (C,CR,N,NR)"
  print "       :           ele+'B'ond_num (C1,C2,C3) 'T':ele+bond+ring (C2,C2R,C3,C3R) 'D'abrl (D,A,B,R,L) [%s]"%(OPT['at'])
  print " -con   : Connectivity.'D'isconnected,'C'onnected 'S'ubstructure(A<=B)"
  print "        :              'T'opo_constrained_disconnected  't'opo_constrained_connected)[%s]"%(OPT['con'])
  print " -mco   : maximum number of connected component. <=0:don't care [%s]"%(OPT['mco'])
  print " -mtd   : max difference of topological distance(num of bonds in the shortest path).<0:don't care. [%s]"%(OPT['mtd'])
  print " -ctd   : calculate topological distance (T or F) [%s]"%(OPT['ctd'])
  print " -w     : weighting. 'D'efault weighting. 'S'pecified by user  [%s]"%(OPT['w'])
  print " -wna   : weight_select_dis for neighoboring atom type     [%s]"%(OPT['wna'])
  print " -wec   : weight_select_dis for extended connectivity      [%s]"%(OPT['wec'])
  print " -wtd   : weight_select_dis for topological distance             [%s]"%(OPT['wtd'])
  print " -wco   : weight_select_dis for num of connected component [%s]"%(OPT['wco'])
  print " -fA,-fB: file formats. 'P'db,'S'df [%s%s]"%(OPT['fA'],OPT['fB'])
  print " -nrs : Nonredundancy check using selection distance score ('T' or 'F')[%s]"%(OPT['nrs'])
  print " -xtd : maximum topological distance for considering Dtopodis [%s]"%(OPT['xtd'])
  print "  -ec  : level of EC for scoring ('0','1','2','3','4') [%s]"%(OPT['ec'])
  print " -per   : Equivalent Atom Permutation. 'F':false, 'A':make all, 'N'on-redundant 'B'uild-up_filter [%s]"%(OPT['per'])
  print " -sbi   : sub integer [%s]"%(OPT['sbi'])
  print " -xsec  : maximum computational time (in second) [%s]"%(OPT['xsec'])
  
  sys.exit(1)


### [1] read option ###

kcombu.read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

### [2] read liglist (-L) ###
liglist = []
read_list_file(OPT['L'],liglist)



### [3] make all the pairs ##
pairlist = []

for i in range(len(liglist)):
  ligA = liglist[i]
  for j in range(i+1,len(liglist)):
    ligB = liglist[j]
    pair = {}
    pair['A'] = ligA
    pair['B'] = ligB
    if (OPT['avex']=='T'):
      outfile = "%s/%s_%s"%(OPT['doAm'],ligA,ligB)
      if (os.path.isfile(outfile)==0):
        pairlist.append(pair)
    elif (OPT['dosm']=='T'):
      mlist = []
      mdat  = {}
      outfile   = "%s/%s_%s"%(OPT['doAm'],ligA,ligB)
      outfilegz = "%s/%s_%s.gz"%(OPT['doAm'],ligA,ligB)
      ok = 0
      if ((os.access(outfile,os.R_OK)) and (os.path.getsize(outfile)<5000)) or ((os.access(outfilegz,os.R_OK)) and (os.path.getsize(outfilegz)<600)):
        ok = kcombu.read_atom_match_file(outfile,mlist,mdat)
      #print "%s %s ok %d len_mlist %d len_mlist[0]['A'] %d"%(ligA,ligB,ok,len(mlist),len(mlist[0]['A']))
      if (ok!=0) and ((len(mlist)==0) or (len(mlist[0]['A'])<=1)):
        pairlist.append(pair)
    else: 
      pairlist.append(pair)



Npair = len(pairlist)

Nsta = Npair*bunshi/bunbo
Nend = Npair*(bunshi+1)/bunbo
print "#bunshi %d bunbo %d Npair %d Nsta %d Nend %d"%(bunshi,bunbo,Npair,Nsta,Nend)

### [4] All vs all comparison in two different methods ###
command_common = "kcombu -fA %s -fB %s -alg %s -con %s -at %s -nk %s -w %s -wec %s -wna %s -wtd %s -wco %s -mco %s -mtd %s -xtd %s -nrs %s -ec %s"%(OPT['fA'],OPT['fB'],OPT['alg'],OPT['con'],OPT['at'],OPT['nk'],OPT['w'],OPT['wec'],OPT['wna'],OPT['wtd'],OPT['wco'],OPT['mco'],OPT['mtd'],OPT['xtd'],OPT['nrs'], OPT['ec'])

if (OPT['ctd']=='T'):
  command_common = command_common + ' -ctd T'
if (OPT['per']!='F'):
  command_common = command_common + ' -per ' + OPT['per']
if (float(OPT['xsec'])>0.0):
  command_common = command_common + ' -xsec ' + OPT['xsec'] 


for i in range(Nsta,Nend):
#for p in (pairlist):
  ligA = pairlist[i]['A']
  ligB = pairlist[i]['B']
  oAmfile = ligA + '_' + ligB
  print ">>>> %s vs %s <<<<"%(ligA,ligB)
  ## (3-1) comparison for the correct reference ##
  command = "%s -A %s/%s -B %s/%s -oAm %s/%s"%(command_common,OPT['idl'],ligA,OPT['idl'],ligB,OPT['doAm'],oAmfile)
  print "#%s"%(command) 
  if (OPT['A']=='T'):
    os.system(command) 
