#!/usr/bin/env python
import sys
import os
import kcombu

LastModDate = "Jun 15, 2011"

 
def read_list_file(ifname,list):
  print "#read_list_file('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#WARNING:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     list.append(fields[0])
  f.close()


def average_Ncomponent(matchlist):
  aveNcomponent = 0.0
  for m in (matchlist):
    aveNcomponent += m.get('Ncomponent',0)
  if (len(matchlist)>0):
    aveNcomponent /= len(matchlist)
  return(aveNcomponent)

def average_Maxdiff_topodis(matchlist):
  aveMaxdiff_topodis = 0.0
  for m in (matchlist):
    aveMaxdiff_topodis += m.get('Maxdiff_topodis',0)
  if (len(matchlist)>0):
    aveMaxdiff_topodis /= len(matchlist)
  return(aveMaxdiff_topodis)


def recall_precision_atom_match(test,refe):
  NpairT = len(test['A'])
  NpairR = len(refe['A'])
  NTandR = 0
  for t in range(NpairT):
    hit = 0
    for r in range(NpairR):
      if (test['A'][t]==refe['A'][r])and (test['B'][t]==refe['B'][r]):
        hit = 1
    if (hit==1):
      NTandR += 1
  recall = precision = 0.0
  if (NpairR>0):
    recall       = float(NTandR)/float(NpairR)
  if (NpairT>0):
    precision    = float(NTandR)/float(NpairT)

  agreement = 0.0
  if ((float(NpairT)+float(NpairR)-float(NTandR))>0.0):
    agreement = float(NTandR)/(float(NpairT)+float(NpairR)-float(NTandR))
  return([NTandR,agreement,recall,precision,NpairT,NpairR]) 

 

def recall_precision_atom_match_top_vs_ave(test,refelist):
  NTandR_ave = 0.0
  NTandR_max = -1 
  NTandR_min = -1 
  NpairT_ave = 0.0
  NpairR_ave = 0.0
  Nall = 0
  maxrankT = 0
  maxrankR = 0
 
  for r in range(len(refelist)):
    refe = refelist[r] 
    [NTandR,agreement,recall,precision,NpairT,NpairR] = recall_precision_atom_match(test,refe)
    NTandR_ave     += NTandR
    NpairT_ave    += NpairT
    NpairR_ave    += NpairR
    Nall += 1
    if (NTandR_max==-1):
      NTandR_max = NTandR_min = NTandR
      maxrankR = r 
    else:
      if (NTandR_max<NTandR):
        NTandR_max = NTandR
        maxrankR = r 
      if (NTandR_min>NTandR):
        NTandR_min = NTandR
 
  NTandR_ave /= Nall
  NpairT_ave /= Nall
  NpairR_ave /= Nall

  return([NTandR_ave,NTandR_min,NTandR_max,maxrankR,maxrankT,NpairT_ave,NpairR_ave]) 

 
def recall_precision_atom_match_ave_vs_top(testlist,refe):
  NTandR_ave = 0.0
  NTandR_max = -1 
  NTandR_min = -1 
  NpairT_ave = 0.0
  NpairR_ave = 0.0
  Nall = 0
  maxrankT = 0
  maxrankR = 0
  for tt in range(len(testlist)):
    test = testlist[t]
    [NTandR,agreement,recall,precision,NpairT,NpairR] = recall_precision_atom_match(test,refe)
    NTandR_ave     += NTandR
    NpairT_ave    += NpairT
    NpairR_ave    += NpairR
    Nall += 1
    if (NTandR_max==-1):
      NTandR_max = NTandR_min = NTandR
    else:
      if (NTandR_max<NTandR):
        NTandR_max = NTandR
        maxrankT = t
      if (NTandR_min>NTandR):
        NTandR_min = NTandR
 
  NTandR_ave /= Nall
  NpairT_ave /= Nall
  NpairR_ave /= Nall
 
  return([NTandR_ave,NTandR_min, NTandR_max, maxrankT, maxrankR, NpairT_ave,NpairR_ave]) 



def recall_precision_atom_match_ave_vs_ave(testlist,refelist):
  NTandR_ave = 0.0
  NTandR_max = -1 
  NTandR_min = -1 
  NpairT_ave = 0.0
  NpairR_ave = 0.0
  Nall = 0
  maxrankT = 0
  maxrankR = 0
  for t in range(len(testlist)):
    test = testlist[t]
    for r in range(len(refelist)):
      refe = refelist[r]
      [NTandR,agreement,recall,precision,NpairT,NpairR] = recall_precision_atom_match(test,refe)
      NTandR_ave     += NTandR
      NpairT_ave    += NpairT
      NpairR_ave    += NpairR
      Nall += 1
      if (NTandR_max==-1):
        NTandR_max = NTandR_min = NTandR
      else:
        if (NTandR_max<NTandR):
          NTandR_max = NTandR
          maxrankT = t
          maxrankR = r 
        if (NTandR_min>NTandR):
          NTandR_min = NTandR
 
  NTandR_ave /= Nall
  NpairT_ave /= Nall
  NpairR_ave /= Nall
 
  return([NTandR_ave,NTandR_min, NTandR_max, maxrankT, maxrankR, NpairT_ave,NpairR_ave]) 



def recall_precision_atom_match_ave_vs_max(testlist,refelist):
  NTandR_ave = 0.0
  NTandR_max = -1 
  NTandR_min = -1 
  maxrankT = 0
  maxrankR = 0

  NpairT_ave = 0.0
  for test in (testlist):
    NpairT_ave += len(test['A'])
  NpairT_ave /= len(testlist)

  NpairR_ave = 0.0
  for refe in (refelist):
    NpairR_ave += len(refe['A'])
  NpairR_ave /= len(refelist)


  for t in range(len(testlist)):
    test = testlist[t]
    ntandr_max = -1 
    for r in range(len(refelist)):
      refe = refelist[r]
      [NTandR,agreement,recall,precision,NpairT,NpairR] = recall_precision_atom_match(test,refe)
      #print "t %d r %d NTandR %d NpairT %d NpairR %d"%(t,r,NTandR,NpairT,NpairR)
      if (NTandR>ntandr_max):
        ntandr_max = NTandR

    #print "ntandr_max %d"%(ntandr_max)
    NTandR_ave    += ntandr_max 
    if (NTandR_max<0) or (ntandr_max > NTandR_max):
      NTandR_max  = ntandr_max 
    if (NTandR_min<0) or (ntandr_max < NTandR_min):
      NTandR_min  = ntandr_max 

  NTandR_ave /= len(testlist) 
  #print "NTandR_ave %f min %f max %f NpairT_ave %f NpairR_ave %f test %d refe %d"%(NTandR_ave,NTandR_min,NTandR_max,NpairT_ave,NpairR_ave,len(testlist),len(refelist)) 
  return([NTandR_ave,NTandR_min, NTandR_max, maxrankT, maxrankR, NpairT_ave,NpairR_ave]) 







############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['idl'] = 'SupLIG'
OPT['A'] = 'F'
OPT['idR'] = 'AtmMch';
OPT['idT'] = 'AtmMch';
OPT['idS'] = '';
OPT['cmp'] = 'aa'

if (len(sys.argv)<2):
  print "cmp_atommatch.py [matchfileA] [matchfileB] <options>"
  print " for comparing making atom match file for all-vs-all comparison for a given molecule list."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>" 
  print " -cmp  :comparison(test_vs_ref) 'tt':top-vs-top,'aa':average-vs-average" 
  print "       :  'at':average-vs-top,'ta':top-vs-average 'am':ave-vs-max [%s]"%(OPT['cmp']) 
  
  print "    #note: If you compare 'build-up' and 'exact', use 'ta', and evaluate by |8:agreement_max|"
  print "    #      If you compare 'build-up' and '3D',    use 'ta', and evaluate by |8:agreement_max|"
  print "    #      If you compare 'exact'    and '3D',    use 'am', and evaluate by |7:agreement_ave|"
  print "<options for 'P':pairwise comparison>"
  print " -iaT  :input atom match file (Test)      [%s]"%(OPT['idT'])
  print " -iaR  :input atom match file (Reference) [%s]"%(OPT['idR'])
  print " -iaS  :input atom match file (Reference2) (optional) [%s]"%(OPT['idS'])
  print "<options for 'L':library comparison>"
  print " -L    :list of library molecules[%s]"%(OPT['L'])
  print " -idl  :input directory for library molecules[%s]"%(OPT['idl'])
  print "  -A   :Action (T or F) [%s]"%(OPT['A'])
  print " -idT  :input directory (Test)      for atom_maching file [%s]"%(OPT['idT'])
  print " -idR  :input directory (Reference) for atom_maching file [%s]"%(OPT['idR'])
  print " -idS  :input directory (Reference2)for atom_maching file (optional) [%s]"%(OPT['idS'])
  sys.exit(1)


### [1] read option ###

kcombu.read_option(sys.argv,OPT)


oAmfileRlist = []
oAmfileTlist = []
oAmfileSlist = []
ligAlist = []
ligBlist = []
### [2L] make oAmfile_lists (library mode) ###
if (OPT['L'] != ''):
  liglist = []
  read_list_file(OPT['L'],liglist)
  for i in range(len(liglist)):
    ligA = liglist[i]
    for j in range(i+1,len(liglist)):
      ligB = liglist[j]
      ligAlist.append(ligA) 
      ligBlist.append(ligB) 
      oAmfileRlist.append(OPT['idR']+'/'+ligA + '_' + ligB)
      oAmfileTlist.append(OPT['idT']+'/'+ligA + '_' + ligB)
      if (OPT['idS'] != ''):
        oAmfileSlist.append(OPT['idS']+'/'+ligA + '_' + ligB)
      else: 
        oAmfileSlist.append(OPT['idR']+'/'+ligA + '_' + ligB)

### [2P] make oAmfile_lists (pairwise mode) ###
elif (OPT['iaT']!='') and (OPT['iaR'] != ''):
  oAmfileRlist.append(OPT['iaR'])
  oAmfileTlist.append(OPT['iaT'])
  ligAlist.append(OPT['iaR'])
  ligBlist.append(OPT['iaT'])
  if (OPT['idS'] != ''):
    oAmfileSlist.append(OPT['iaS'])
  else:
    oAmfileSlist.append(OPT['iaR'])
else:
  print "#ERROR: option (-L) or (-iaT, -iaR) are necessary."
  sys.exit(1)


### [3] All vs all comparison in two different methods ###
Msame_NpairR_NpairT = 0    
Ndiff_NpairR_NpairT = 0    
Nligpair = 0;
AVEagreement_ave = 0.0
AVEagreement_max = 0.0

print "#COMMAND %s"%(OPT['COMMAND'])
print "#DATE    %s"%(OPT['START_DATE'])
print "#1:taniS|2:taniR|3:taniT|4:NpairR|5:NpairT|6:NRandT|7:agreement_ave|8:agreement_max|9:maxrk_R|10:maxrk_T|"
print "#11:ligA|12:ligB|13:NheavyatomA|14:NheavyatomB|15:Natompair|16:COMP_TIME|"
print "#17:Natom_matchR|18:Natom_matchT|19:aveNcompoR|20:aveNcompoT|21:aveMaxdifGrhdisR|22:aveMaxdifGrhdisT"
for i in range(len(oAmfileTlist)):
  oAmfileT = oAmfileTlist[i]
  oAmfileR = oAmfileRlist[i]
  oAmfileS = oAmfileSlist[i]
  ligA     = ligAlist[i]
  ligB     = ligBlist[i]

  mlistR = []
  mlistT = []
  mlistS = []
  mdatR = {}
  mdatT = {}
  mdatS = {}
 
  okR = kcombu.read_atom_match_file(oAmfileR,mlistR,mdatR)
  if (okR==0):
    print "#WARNING:oAmfileR '%s' has no content."%(oAmfileR)   
    #print "len_matchlistR %d"%(len(mlistR))
    #for i in range (len(mlistR)):
    #  print "[%d] Npair %d tanimoto %f Npair %d"%(i,len(mlistR[i]['A']),mlistR[i]['tanimoto'],mlistR[i]['Npair_atom'])
    #  for j in range (len(mlistR[i]['A'])):
    #    print " %d-%d"%(mlistR[i]['A'][j],mlistR[i]['B'][j])
  okT = kcombu.read_atom_match_file(oAmfileT,mlistT,mdatT)
  if (okT==0):
    print "#WARNING:oAmfileT '%s' has no content."%(oAmfileT)
 
  if (OPT['idS'] != ''):
    oAmfileS = OPT['idS']+'/'+ligA + '_' + ligB
    okS  = kcombu.read_atom_match_file(oAmfileS,mlistS,mdatS)
    tanimotoS = mlistS[0]['tanimoto'] 
  elif (okR>0):
    okS = 1
    if (len(mlistR)>0):
      tanimotoS = mlistR[0].get('tanimoto',0.0) 
    else:
      tanimotoS = 0.0
  else:
      okS = 0

  if (okR>0) and (okT>0) and (okS>0):
    if (len(mlistT)==0) or (len(mlistR)==0):
      NTandR = 0
      agreement = 0
      recall = 0.0
      precision = 0.0
      NpairT = 0
      NpairR = 0
      NTandRave = 0
      NTandRmin = 0
      NTandRmax = 0
      agreement_ave = 0.0
      agreement_min = 0.0
      agreement_max = 0.0
      aveNcomponentR = 0.0
      aveNcomponentS = 0.0
      aveNcomponentT = 0.0
#      print "%f %f %f %2d %2d %2d %f %f %f %s %s %2d %2d %3d %f %4d %4d %.2f %.2f"%(tanimotoS,0.0,0.0,NpairR,NpairT,NTandRave,agreement_ave,agreement_min, agreement_max,ligA,ligB,mdatT['NheavyatomA'],mdatT['NheavyatomB'],mdatT['TotalNatompair'],mdatT['COMP_TIME'],len(mlistR),len(mlistT),aveNcomponentR,aveNcomponentT)
      print "%f %f %f %2d %2d %2d %f %f %4d %4d %s %s %2d %2d %3d %f %4d %4d %.2f %.2f %.2f %.2f"%(tanimotoS,0.0, 0.0,NpairR,NpairT,NTandRave,agreement_ave,agreement_max, maxrankR+1, maxrankT+1, ligA,ligB,mdatT['NheavyatomA'],mdatT['NheavyatomB'],mdatT['TotalNatompair'],mdatT['COMP_TIME'],len(mlistR),len(mlistT),aveNcomponentR,aveNcomponentT,aveMaxdiff_topodisR,aveMaxdiff_topodisT)
    else:
      if (OPT['cmp']=='tt'):
        [NTandR,agreement,recall,precision,NpairT,NpairR] = recall_precision_atom_match(mlistT[0],mlistR[0])
        NTandRave = NTandR
        NTandRmin = NTandR
        NTandRmax = NTandR
        maxrankR = 0  
        maxrankT = 0  

      if (OPT['cmp']=='ta'):
        [NTandRave,NTandRmin,NTandRmax,maxrankR, maxrankT,NpairT,NpairR] = recall_precision_atom_match_top_vs_ave(mlistT[0],mlistR)
      if (OPT['cmp']=='at'):
        [NTandRave,NTandRmin,NTandRmax,maxrankR, maxrankT,NpairT,NpairR] = recall_precision_atom_match_ave_vs_top(mlistT,mlistR[0])
      if (OPT['cmp']=='aa'):
        [NTandRave,NTandRmin,NTandRmax,maxrankR, maxrankT,NpairT,NpairR] = recall_precision_atom_match_ave_vs_ave(mlistT,mlistR)
      if (OPT['cmp']=='am'):
        [NTandRave,NTandRmin,NTandRmax,maxrankR, maxrankT,NpairT,NpairR] = recall_precision_atom_match_ave_vs_max(mlistT,mlistR)
  
      
      if (float(NpairT + NpairR - NTandRave)>0.0):
        agreement_ave = (NTandRave)/float(NpairT + NpairR - NTandRave)
        agreement_min = (NTandRmin)/float(NpairT + NpairR - NTandRmin)
        agreement_max = (NTandRmax)/float(NpairT + NpairR - NTandRmax)
      else:
        agreement_ave = 0.0
        agreement_min = 0.0
        agreement_max = 0.0
      AVEagreement_ave += agreement_ave
      AVEagreement_max += agreement_max
 
      aveNcomponentR = average_Ncomponent(mlistR)
      aveNcomponentS = average_Ncomponent(mlistS)
      aveNcomponentT = average_Ncomponent(mlistT)

      aveMaxdiff_topodisR = average_Maxdiff_topodis(mlistR)
      aveMaxdiff_topodisS = average_Maxdiff_topodis(mlistS)
      aveMaxdiff_topodisT = average_Maxdiff_topodis(mlistT)
      if (OPT['cmp']=='ta') or (OPT['cmp']=='tt'):
        aveMaxdiff_topodisT = mlistT[0].get('Maxdiff_topodis',0)

      print "%f %f %f %2d %2d %2d %f %f %4d %4d %s %s %2d %2d %3d %f %4d %4d %.2f %.2f %.2f %.2f"%(tanimotoS,mlistR[0]['tanimoto'], mlistT[0]['tanimoto'],NpairR,NpairT,NTandRave,agreement_ave,agreement_max, maxrankR+1, maxrankT+1, ligA,ligB,mdatT['NheavyatomA'],mdatT['NheavyatomB'],mdatT['TotalNatompair'],mdatT['COMP_TIME'],len(mlistR),len(mlistT),aveNcomponentR,aveNcomponentT,aveMaxdiff_topodisR,aveMaxdiff_topodisT)
     
    Nligpair += 1;
    Ndiff_NpairR_NpairT += (NpairR - NpairT);    
    if (NpairR==NpairT):
      Msame_NpairR_NpairT += 1    

### ENDING PROCEDURE ###
if (Nligpair>0):
  AVEdiff_NpairR_NpairT = float(Ndiff_NpairR_NpairT)/Nligpair
  Rsame_NpairR_NpairT   = float(Msame_NpairR_NpairT)/Nligpair
  AVEagreement_ave    /= Nligpair 
  AVEagreement_max    /= Nligpair 
else:
  AVEdiff_NpairR_NpairT = 0.0
  Rsame_NpairR_NpairT   = 0.0
  AVEagreement_ave    = 0.0
  AVEagreement_max    = 0.0

print "#%s %s Nligpair %d Msame_NpairR_NpairT %d Rsame_NpairR_NpairT %f AVEdiff_NpairR_NpairT %f AVEagreement_ave %f max %f"%(OPT['idT'],OPT['idR'],Nligpair,Msame_NpairR_NpairT,Rsame_NpairR_NpairT,AVEdiff_NpairR_NpairT,AVEagreement_ave,AVEagreement_max)
