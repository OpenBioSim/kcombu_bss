#!/usr/bin/env python
##
## <sumCmpTwoGrp.py>
##


import sys
import os

LastModDate = "Feb 7, 2012"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_comparison_bwn_two_groups(ifname,dataQlist,dataLlist,sim):
  print "#read_comparison_bwn_two_groups('%s')"%(ifname)
#<<FILE FORMAT EXAMPLE>>
# #NMOL_IN_LIBRARY 15493
# #DATAQ 0 3OX
# #DATAQ 1 IPM
# #DATAQ 2 OOK
# #DATAQ 3 DDO
# #DATAQ 4 IC8
# #DATAQ 5 VAW
# #DATAQ 6 331
# #DATAQ 13580 GRD
# #DATAQ 13581 EPQ
# #DATAL 0 C04146.mol
# #DATAL 1 C10027.mol
# #DATAL 2 C17987.mol
# #DATAL 3 C18878.mol
# #DATAL 15492 C16504.mol
# 1 /DB/LIGAND-EXPO/SDF/IPM 19 C06144.mol 0.500000
# 1 /DB/LIGAND-EXPO/SDF/IPM 37 C05682.mol 0.538462
# 1 /DB/LIGAND-EXPO/SDF/IPM 62 C00767.mol 0.529412
# 1 /DB/LIGAND-EXPO/SDF/IPM 80 C16436.mol 0.538462
# 1 /DB/LIGAND-EXPO/SDF/IPM 88 C00631.mol 0.533333
# 13581 /DB/LIGAND-EXPO/SDF/EPQ 15206 C01259.mol 0.500000
# 13581 /DB/LIGAND-EXPO/SDF/EPQ 15241 C01404.mol 0.500000
# 13581 /DB/LIGAND-EXPO/SDF/EPQ 15320 C05708.mol 0.533333
  Nvalue = 0
  f = open(ifname)  
  for line in f:
    line = line.rstrip('\n')
    #print line
    field = line.split()
    if (line.startswith('#DATAQ')):
      dataQlist.append(field[2])
    elif (line.startswith('#DATAL')):
      dataLlist.append(field[2])
    elif (line.startswith('#')==0) and (len(field)>=3):
      if (len(field)==3):
        numQ = int(field[0])
        numL = int(field[1])
        val  = float(field[2])
      elif (len(field)==5):
        numQ = int(field[0])
        numL = int(field[2])
        val  = float(field[4])
      if (sim.has_key(numQ)==0):
        sim[numQ] = {}

      sim[numQ][numL] = val 
      Nvalue += 1
  print "#Nsimilarity %d"%(Nvalue)
  f.close()

###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['th'] = 0.9
OPT['osi'] = 'sim.out'
OPT['ods'] = 'dis.out'

if (len(sys.argv)<2):
  print "sumCmpTwoGrp.py [comp_output_file] <options>"
  print " for making open src packages for the 'kcombu' program, excluding improper files."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -th  : similarity threshold [%s]"%(OPT['th']) 
  print " -osi : output file for similar molecules [%s]"%(OPT['osi']) 
  print " -ods : output file for disimilar molecules [%s]"%(OPT['ods']) 
  print " -A   : Action (T or F) [%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)
OPT['th'] = float(OPT['th'])
icmpfile = sys.argv[1]
dataQlist = []
dataLlist = []
sim_matrix = {}
read_comparison_bwn_two_groups(icmpfile,dataQlist,dataLlist,sim_matrix)
print "#NdataQ %d NdataL %d"%(len(dataQlist),len(dataLlist))

fos = open(OPT['osi'],'w')
fod = open(OPT['ods'],'w')

fos.write("#QUERY_COMPOUND_SIMILAR_TO_ONE_OF_LIBRARY_MOLECULE\n")
fos.write("#COMMAND '%s'\n"%(OPT['COMMAND']))
fod.write("#QUERY_COMPOUND_DISIMILAR_TO_ALL_OF_THE_LIBRARY_MOLECULES\n")
fod.write("#COMMAND '%s'\n"%(OPT['COMMAND']))

print "#write_similar_molecule   -->'%s'"%(OPT['osi'])
print "#write_disimilar_molecule -->'%s'"%(OPT['ods'])
Nout_sim = Nout_dis = 0

for q in range(len(dataQlist)):
  hit = 0
  max_sim = 0.0
  lib_max_sim = ''
  if (sim_matrix.has_key(q)):
    for l in (sim_matrix[q].keys()):
      sim = sim_matrix[q][l] 
      if (sim> max_sim):
        max_sim = sim 
        lib_max_sim = dataLlist[l]
  if ((q%100)==0):
    print "%s max_sim %f id %s"%(dataQlist[q],max_sim, lib_max_sim) 
  if (max_sim >= OPT['th']):
    fos.write("%s %s %f\n"%(dataQlist[q],lib_max_sim,max_sim))
    Nout_sim += 1
  else:
    fod.write("%s\n"%(dataQlist[q]))
    Nout_dis += 1
print "#Nout_sim %d Nout_dis %d"%(Nout_sim,Nout_dis)
print "#write_similar_molecule   -->'%s'"%(OPT['osi'])
print "#write_disimilar_molecule -->'%s'"%(OPT['ods'])
fos.close()
fod.close()
