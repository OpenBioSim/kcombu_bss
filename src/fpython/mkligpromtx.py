#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "June 19, 2012"


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
     list.append(fields[0])
  f.close()

def read_lig_pro_file(ifname,list):
  #print "#read_list_file('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
# 3RA_3rawA  CLK3 1 A 
# V25_2wu7A  CLK3 1483 A
# DKI_2wu6A  CLK3 1484 A
# ANK_3an0A  MP2K6 341 A6
# STU_3fmeA  MP2K6 1 A M
# Q9G_2xyuA  EPHA4 1898 A
# 1N1_2y6oA  EPHA4 1892 A
# 03Q_3pp0A  ERBB2 1 A
# 03P_3rcdA  ERBB2 9001 A


  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     ligpro = {}
     field = line.split()
     (lig,pdbch) = field[0].split('_')
     uniprot = field[1]
     ligpro['lig'] = lig
     ligpro['pdbch'] = pdbch 
     ligpro['uniprot'] = uniprot 
     list.append(ligpro)
  f.close()


#################
##### MAIN ######
#################

OPT = {}
OPT['iligpro'] = ''
OPT['ipro'] = ''
OPT['ilig'] = ''
OPT['otab'] = ''
OPT['minNliglabel'] = '1'
OPT['minNligpnt']   = '0'
OPT['ognu'] = ''
OPT['olabel'] = ''
OPT['opnt'] = ''

if (len(sys.argv)<2):
  print "mkligpromtx.py <options>"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -iligpro : input lig pdb file [%s]"%(OPT['iligpro'])
  print " -ipro    : input protein file [%s]"%(OPT['ipro'])
  print " -ilig    : input protein file [%s]"%(OPT['ilig'])
  print " -otab    : output tab-splited file [%s]"%(OPT['otab'])
  print " -ognu    : output gnuplot file [%s]"%(OPT['ognu'])
  print " -opnt    : output gnuplot point file [%s]"%(OPT['opnt'])
  print " -olabel  : output gnuplot label file [%s]"%(OPT['olabel'])
  print " -minNliglabel : minimum Nlig for labeling [%s]"%(OPT['minNliglabel'])
  print " -minNligpnt  : minimum Nlig for opnt [%s]"%(OPT['minNligpnt'])
  sys.exit(1)


read_option(sys.argv,OPT)

ligprolist = []
read_lig_pro_file(OPT['iligpro'],ligprolist)

prolist = []
read_list_file(OPT['ipro'],prolist)

liglist = []
read_list_file(OPT['ilig'],liglist)

#### cal statistic of ligpro_mtx ###

ligpro_mtx = {}
for lig in (liglist):
  ligpro_mtx[lig] = {}
  for pro in (prolist):
    ligpro_mtx[lig][pro] = '' 

for x in (ligprolist):
  lig   = x['lig']
  pro   = x['uniprot']        
  pdbch = x['pdbch']
  if (ligpro_mtx.has_key(lig)==0):
    ligpro_mtx[lig] = {}
  ligpro_mtx[lig][pro] = pdbch


Nlig_count = {}
for lig in (liglist):
  Nlig_count[lig] = 0 
  for pro in (prolist):
    if (ligpro_mtx[lig][pro] != ''):
      Nlig_count[lig] += 1 


Npro_count = {}
for pro in (prolist):
  Npro_count[pro] = 0 
  for lig in (liglist):
    if (ligpro_mtx[lig][pro] != ''):
      Npro_count[pro] += 1 




### OUTPUT TAB-SPLITED TABLE ###
if (OPT['otab'] != ''):
  print "#write_tab_splited_table() -->'%s'"%(OPT['otab'])
  of = open(OPT['otab'],'w')
  of.write("\t")
  for lig in (liglist):
    of.write("%s\t"%(lig))
  of.write("\n") 

  for pro in (prolist):
    of.write("%s\t"%(pro))
    for lig in (liglist):
      if (ligpro_mtx[lig][pro] != ''):
        out = '1'
      else: 
        out = ''
      of.write("%s\t"%(out))
    of.write("\n") 
  of.close()

offsetX = 50
offsetY = 10

if (OPT['ognu'] != ''):
  print "#write_gnuplot_file() -->'%s'"%(OPT['ognu'])
  of = open(OPT['ognu'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  for i in range(len(liglist)+1):
    lig  = ''
    lig0 = ''
    if (i<len(liglist)):
      lig = liglist[i]
    if (i>0):
      lig0 = liglist[i-1]
    for j in range(len(prolist)+1):
      pro = ''
      pro0 = ''
      if (j<len(prolist)):
        pro = prolist[j]
      if (j>0):
        pro0 = prolist[j-1]

      if (lig != '') and (pro != '') and (ligpro_mtx[lig][pro] != ''):
        out = '1'
      elif (lig0 != '') and (pro != '') and (ligpro_mtx[lig0][pro] != ''):
        out = '1'
      elif (lig  != '') and (pro != '') and (ligpro_mtx[lig][pro] != ''):
        out = '1'
      elif (lig0 != '') and (pro0 != '') and (ligpro_mtx[lig0][pro0] != ''):
        out = '1'
      else: 
        out = '0'
      if (out=='1'):
        print "'%s' '%s' '%s' %s' out %s"%(lig,lig0,pro,pro0,out)
      if (lig != '') and (pro != ''):
        print "'%s' '%s' '%s'"%(lig,pro,ligpro_mtx[lig][pro])
      of.write("%d %d %s\n"%(i+offsetX,j+offsetY,out))
    of.write("\n") 
  of.close()


if (OPT['olabel'] != ''):
  ofname = OPT['olabel'] + '.xlabel'
  print "#write_gnuplot_label_file() -->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("set xrange [%d:%d]\n"%(0,len(liglist)+offsetX))
  of.write("unset xtics\n");
  for i in range(len(liglist)):
    lig = liglist[i]
    if (Nlig_count[lig]>=5):
      of.write("set label \"%s (%d)\" at first %d,first %.0f rotate by 90\n"%(lig,Nlig_count[lig],i+offsetX,0.0*offsetY))
    elif (Nlig_count[lig]>=int(OPT['minNliglabel'])):
      of.write("set label \"%s (%d)\" at first %d,first %.0f rotate by 90\n"%(lig,Nlig_count[lig],i+offsetX,0.4*offsetY))
    elif (lig=='IRE'):
      of.write("set label \"%s (%d)\" at first %d,first %.0f rotate by 90\n"%(lig,Nlig_count[lig],i+offsetX,0.4*offsetY))
  of.close()  

  ofname = OPT['olabel'] + '.ylabel'
  print "#write_gnuplot_label_file() -->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("set yrange [%d:%d]\n"%(0,len(prolist)+offsetY))
  of.write("unset ytics\n");
  for j in range(len(prolist)):
    pro = prolist[j]
    of.write("set label \"%s (%d)\" at first %d,first %d\n"%(pro,Npro_count[pro],0,j+offsetY))
  of.close()  


if (OPT['opnt'] != ''):
  print "#write_gnuplot_file() -->'%s'"%(OPT['opnt'])
  of = open(OPT['opnt'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  for i in range(len(liglist)):
    lig = liglist[i]
    for j in range(len(prolist)):
      pro = prolist[j]
      if (ligpro_mtx[lig][pro] != '') and (Nlig_count[lig]>=int(OPT['minNligpnt'])):
        of.write("%d %d\n"%(i+offsetX,j+offsetY))
    of.write("\n") 
  of.close()
