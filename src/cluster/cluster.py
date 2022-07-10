#!/usr/bin/env python
import sys
import os
import math
from datetime import datetime

LastModDate = "July 19, 2011"

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


def set_end_date(opt_dic):
  now = datetime.now()
  opt_dic['END_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

class Entry:
  def __init__(self):
    self.name     = ''           ## name of the entry
    self.num      = 0            ## entry number.
    self.neighbor_nums    = []   ## a list of the number of neighbor entries 
    self.similarities     = []   ## a list of similaritiest of neighbor entries 
    self.cluster_num  = -1       ## cluster number. (>=0). if (<0), it is not assined yet.
    self.Nneighbor = 0           ## number of neighbors
    self.worst_sim =  1.1        ## worst similarity for neighbors_nums[]
    self.worst_num = -1          ## neighbors_num with the worst similarity
    self.Nneighbor_in_cluster=0  ## number of neighbors in the clusters
   
class Cluster:
  def __init__(self):
    self.num      = 0   ## cluster number 
    self.Nmember  = 0   ## number of memmbers
    self.members  = []  ## a list of number of member entries
    self.rep_entry_num  = -1 ## number of representative entry


def read_similarities_in_list_format_threshold(ifname,entry_list,sim_thre):
#>>FILE FORMAT EXAMPLE <<
# #>>Similarities in List format<<
# #COMMAND fkcombu -M A -ide temp.des -osl temp.osl
# #DATE_START Jun 22,2011 13:55:48
# #DATE_END   Jun 22,2011 13:55:48
# #COMP_TIME  0.035351 seconds
# #NDATA 131
# #DATA 0 3OX
# #DATA 1 IPM
# #DATA 2 OOK
# #DATA 3 DDO
# :
# #DATA 129 GVO
# #DATA 130 AHO
#0 1 0.045455
#0 2 0.219466
#0 3 0.030380
#0 4 0.164859
#0 5 0.344063
#0 6 0.385759
#0 7 0.022613
#:
#128 130 0.089189
#129 130 0.010526
  print "#read_similarities_in_list_format_threshold('%s',entry_list,sim_thre)"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#WARNING:Can't open similarity-list file '%s'" %(ifname)
    return(0)
  f = open(ifname)

  Ndata = 0
  ratio0 = 0

  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#NDATA')):
      (head,value) = line.split()
      Ndata = int(value)
      Npair = Ndata*(Ndata+1)/2
      print "#Ndata %d"%(Ndata)
      for i in range(Ndata):
        entry = Entry()
        entry.num = i
        entry.neighbor_nums = [] 
        entry_list.append(entry)
    elif (line.startswith('#DATA')):
      field = line.split()
      num = int(field[1])
      if (num>=0) and (num<Ndata):   
        entry_list[num].name = field[2] 
    elif (line.startswith('#')==0) and (len(line)>4):
      field = line.split() 
      i = int(field[0])
      j = int(field[1])
      sim = float(field[2])
      #print i,j,sim
      if (i!=j) and (sim>=sim_thre):
        npair = (Ndata-i)*(Ndata-i+1)/2
        ratio = int(100.0*float(Npair - npair)/float(Npair));
        if ((ratio!=ratio0) and (ratio%5)==0):
          print "#read %d %%"%(ratio)
        entry_list[i].neighbor_nums.append(j)
        entry_list[j].neighbor_nums.append(i)
        ratio0 = ratio
  f.close()

  aveNneighbor = 0
  maxNneighbor = -1 
  minNneighbor = -1 
  for x in (entry_list):
    x.Nneighbor = len(x.neighbor_nums)
    aveNneighbor += x.Nneighbor
    if (minNneighbor<0) or (x.Nneighbor < minNneighbor):
      minNneighbor = x.Nneighbor 
    if (maxNneighbor<0) or (x.Nneighbor > maxNneighbor):
      maxNneighbor = x.Nneighbor 
  aveNneighbor = float(aveNneighbor)/len(entry_list)
  print "#aveNeighbor %f min %d max %d"%(aveNneighbor, minNneighbor, maxNneighbor) 



def read_similarities_in_list_format_K_neighbors(ifname,entry_list,K):
#>>FILE FORMAT EXAMPLE <<
# #>>Similarities in List format<<
# #COMMAND fkcombu -M A -ide temp.des -osl temp.osl
# #DATE_START Jun 22,2011 13:55:48
# #DATE_END   Jun 22,2011 13:55:48
# #COMP_TIME  0.035351 seconds
# #NDATA 131
# #DATA 0 3OX
# #DATA 1 IPM
# #DATA 2 OOK
# #DATA 3 DDO
# :
# #DATA 129 GVO
# #DATA 130 AHO
#0 1 0.045455
#0 2 0.219466
#0 3 0.030380
#0 4 0.164859
#0 5 0.344063
#0 6 0.385759
#0 7 0.022613
#:
#128 130 0.089189
#129 130 0.010526
  print "#read_similarities_in_list_format_K_neighbors('%s',entry_list,sim_thre)"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#WARNING:Can't open similarity-list file '%s'" %(ifname)
    return(0)
  f = open(ifname)

  Ndata = 0
  ratio0 = 0

  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#NDATA')):
      (head,value) = line.split()
      Ndata = int(value)
      Npair = Ndata*(Ndata+1)/2
      print "#Ndata %d"%(Ndata)
      for i in range(Ndata):
        entry = Entry()
        entry.num = i
        entry.neighbor_nums = [] 
        entry_list.append(entry)
    elif (line.startswith('#DATA')):
      (head,num,name) = line.split()
      num = int(num)
      if (num>=0) and (num<Ndata):   
        entry_list[num].name = name
    elif (line.startswith('#')==0) and (len(line)>4):
      (i,j,sim) = line.split()
      i = int(i)
      j = int(j)
      sim = float(sim)
      if (i!=j):
        accept_i = 0
        accept_j = 0

        if (len(entry_list[i].neighbor_nums)<K):
          entry_list[i].neighbor_nums.append(j)
          entry_list[i].similarities.append(sim)
          accept_i = 1
        elif (sim > entry_list[i].worst_sim):
          entry_list[i].worst_sim = sim
          entry_list[i].neighbor_nums[entry_list[i].worst_num] = j
          entry_list[i].similarities[entry_list[i].worst_num] = sim
          accept_i = 1

        if (accept_i==1):
          npair = (Ndata-i)*(Ndata-i+1)/2
          ratio = int(100.0*float(Npair - npair)/float(Npair));
          if ((ratio!=ratio0) and (ratio%5)==0):
            print "#read %d %%"%(ratio)
          ratio0 = ratio
          entry_list[i].worst_sim = 1.1 
          worst_k = -1 
          for k in range(len(entry_list[i].neighbor_nums)):
            if (entry_list[i].similarities[k]<entry_list[i].worst_sim):
              entry_list[i].worst_sim = entry_list[i].similarities[k]
              worst_k = k 
          entry_list[i].worst_num = worst_k 
 
        if (len(entry_list[j].neighbor_nums)<K):
          entry_list[j].neighbor_nums.append(i)
          entry_list[j].similarities.append(sim)
          accept_j = 1
        elif (sim > entry_list[j].worst_sim):
          entry_list[j].worst_sim = sim
          entry_list[j].neighbor_nums[entry_list[j].worst_num] = i
          entry_list[j].similarities[entry_list[j].worst_num] = sim
          accept_j = 1

        if (accept_j==1):
          entry_list[j].worst_sim = 1.1
          worst_k = -1 
          for k in range(len(entry_list[j].neighbor_nums)):
            if (entry_list[j].similarities[k] < entry_list[j].worst_sim):
              entry_list[j].worst_sim = entry_list[j].similarities[k]
              worst_k = k 
          entry_list[j].worst_num = worst_k 
 

  f.close()

  aveNneighbor = 0
  maxNneighbor = -1 
  minNneighbor = -1 
  for x in (entry_list):
    x.Nneighbor = len(x.neighbor_nums)
    aveNneighbor += x.Nneighbor
    if (minNneighbor<0) or (x.Nneighbor < minNneighbor):
      minNneighbor = x.Nneighbor 
    if (maxNneighbor<0) or (x.Nneighbor > maxNneighbor):
      maxNneighbor = x.Nneighbor 
  aveNneighbor = float(aveNneighbor)/len(entry_list)
  print "#aveNeighbor %f min %d max %d"%(aveNneighbor, minNneighbor, maxNneighbor) 






def mark_neighbor_nums(entry_list,focus_i,cluster_num,cluster):
  #print "#mark_neighbor_nums(focus_i %d cluster_num %d)"%(focus_i,cluster_num)
  entry_list[focus_i].cluster_num = cluster_num
  cluster.members.append(focus_i)
  for j in (entry_list[focus_i].neighbor_nums):
    if (entry_list[j].cluster_num==-1):
      mark_neighbor_nums(entry_list,j,cluster_num,cluster)


def set_Nneighbor_in_cluster(entry_list, cluster_list):
  for x in (entry_list):
    x.Nneighbor_in_cluster = 0
    for j in (x.neighbor_nums):
      if (x.cluster_num == entry_list[j].cluster_num):
        x.Nneighbor_in_cluster += 1
  i = 0
  for c in (cluster_list):
    c.num = i 
    #of.write(">[%d] %d\n"%(c.num,c.Nmember)) 
    maxj = c.members[0]
    Nmax = entry_list[maxj].Nneighbor_in_cluster 
    for j in (c.members):
      if (entry_list[j].Nneighbor_in_cluster>Nmax):
        Nmax = entry_list[j].Nneighbor_in_cluster
        maxj = j
    c.rep_entry_num = maxj
    #print "cluster %d maxj %d"%(c.num,maxj)

def write_clusters(ofname,cluster_list,entry_list,opt_dic):
  print "#write_clusters()-->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("#COMMAND %s\n"%(opt_dic['COMMAND']))
  of.write("#START_DATE %s\n"%(opt_dic['START_DATE']))
  of.write("#END_DATE   %s\n"%(opt_dic['END_DATE']))
  of.write("#Nentry     %d\n"%(len(entry_list)))
  of.write("#Ncluster   %d\n"%(len(cluster_list)))
  i = 0
  for c in (cluster_list):
    c.num = i 
    of.write(">CLUSTER %d Nmember %d rep_entry %s\n"%(c.num,c.Nmember,entry_list[c.rep_entry_num].name)) 
    i += 1
    for j in (c.members):
      of.write("%s %d"%(entry_list[j].name,entry_list[j].Nneighbor_in_cluster))
      #for k in (entry_list[j].neighbor_nums):
      #  of.write(" %s"%(entry_list[k].name))
      of.write("\n")
  of.close()

def write_representatives(ofname,cluster_list,entry_list,opt_dic):
  print "#write_representatives()-->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("#COMMAND %s\n"%(opt_dic['COMMAND']))
  of.write("#START_DATE %s\n"%(opt_dic['START_DATE']))
  of.write("#END_DATE   %s\n"%(opt_dic['END_DATE']))
  of.write("#Nentry     %d\n"%(len(entry_list)))
  of.write("#Ncluster   %d\n"%(len(cluster_list)))
  i = 0
  for c in (cluster_list):
    of.write("%-16s cluster %5d Nmember %5d\n"%(entry_list[c.rep_entry_num].name,c.num,c.Nmember))
  of.close()




def mark_neighbor_nums(entry_list,focus_i,cluster_num,cluster):
  #print "#mark_neighbor_nums(focus_i %d cluster_num %d)"%(focus_i,cluster_num)
  entry_list[focus_i].cluster_num = cluster_num
  cluster.members.append(focus_i)
  for j in (entry_list[focus_i].neighbor_nums):
    if (entry_list[j].cluster_num==-1):
      mark_neighbor_nums(entry_list,j,cluster_num,cluster)


def set_Nneighbor_in_cluster(entry_list, cluster_list):
  for x in (entry_list):
    x.Nneighbor_in_cluster = 0
    for j in (x.neighbor_nums):
      if (x.cluster_num == entry_list[j].cluster_num):
        x.Nneighbor_in_cluster += 1
  i = 0
  for c in (cluster_list):
    c.num = i 
    #of.write(">[%d] %d\n"%(c.num,c.Nmember)) 
    maxj = c.members[0]
    Nmax = entry_list[maxj].Nneighbor_in_cluster 
    for j in (c.members):
      if (entry_list[j].Nneighbor_in_cluster>Nmax):
        Nmax = entry_list[j].Nneighbor_in_cluster
        maxj = j
    c.rep_entry_num = maxj
    #print "cluster %d maxj %d"%(c.num,maxj)

def write_clusters(ofname,cluster_list,entry_list,opt_dic):
  print "#write_clusters()-->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("#COMMAND %s\n"%(opt_dic['COMMAND']))
  of.write("#START_DATE %s\n"%(opt_dic['START_DATE']))
  of.write("#END_DATE   %s\n"%(opt_dic['END_DATE']))
  of.write("#Nentry     %d\n"%(len(entry_list)))
  of.write("#Ncluster   %d\n"%(len(cluster_list)))
  i = 0
  for c in (cluster_list):
    c.num = i 
    of.write(">CLUSTER %d Nmember %d rep_entry %s\n"%(c.num,c.Nmember,entry_list[c.rep_entry_num].name)) 
    i += 1
    for j in (c.members):
      of.write("%s %d"%(entry_list[j].name,entry_list[j].Nneighbor_in_cluster))
      #for k in (entry_list[j].neighbor_nums):
      #  of.write(" %s"%(entry_list[k].name))
      of.write("\n")
  of.close()

def write_representatives(ofname,cluster_list,entry_list,opt_dic):
  print "#write_representatives()-->'%s'"%(ofname)
  of = open(ofname,"w")
  of.write("#COMMAND %s\n"%(opt_dic['COMMAND']))
  of.write("#START_DATE %s\n"%(opt_dic['START_DATE']))
  of.write("#END_DATE   %s\n"%(opt_dic['END_DATE']))
  of.write("#Nentry     %d\n"%(len(entry_list)))
  of.write("#Ncluster   %d\n"%(len(cluster_list)))
  i = 0
  for c in (cluster_list):
    of.write("%-16s cluster %5d Nmember %5d\n"%(entry_list[c.rep_entry_num].name,c.num,c.Nmember))
  of.close()



##############
#### MAIN ####
##############

OPT = {}
OPT['A'] = 'F'
OPT['isl'] = ''
OPT['thre'] = 0.9
OPT['ocl'] = ''
OPT['ore'] = ''
OPT['alg'] = 'S'
OPT['k']    = 10
OPT['kmin'] = 5 

if (len(sys.argv)<2):
  print "cluster.py <options>"
  print " for clustering using Single Linkage Clustering or Jarvis-Patrick methods"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -isl   : input similarities in list format[%s]"%(OPT['isl'])
  print " -alg   : 'S'ingle_linkage, 'R'eccurent single linkage, 'J'arvis-Patrick [%s]"%(OPT['alg']) 
  print " -ocl   : output file for clusters [%s]"%(OPT['ocl'])
  print " -ore   : output file for representatives [%s]"%(OPT['ore'])
  print "<parameters for 'S' and 'R'>"
  print " -thre  : threshold similarity for SLC [%f]"%(OPT['thre'])
  print "<parameters for 'J'>"
  print " -k     : number of neighbors [%d]"%(OPT['k'])
  print " -kmin  : minimum number of neighbors for merging into one cluster [%d]"%(OPT['kmin'])
  sys.exit(1)

read_option(sys.argv,OPT)

### [1] Read similarities ###
entry_list = [] 
if (OPT['alg']=='S') or (OPT['alg']=='R'):
  read_similarities_in_list_format_threshold(OPT['isl'],entry_list,float(OPT['thre']))
if (OPT['alg']=='J'):
  read_similarities_in_list_format_K_neighbors(OPT['isl'],entry_list,int(OPT['k']))
print "#Nlist %d"%(len(entry_list))

#for x in (entry_list):
#  print ">%s"%(x.name)
#  for j in range(len(x.neighbor_nums)):
#    print "  %s %f"%(entry_list[x.neighbor_nums[j]].name,x.similarities[j]) 
#
#sys.exit(1)

### 'R':Single Linkage Clustering using recursive calls ### 
if (OPT['alg']=='R'):
  Ncluster = -1
  cluster_list = []
  for i in range(len(entry_list)):
    if (entry_list[i].cluster_num==-1):
      Ncluster += 1 
      clus = Cluster()
      clus.num = Ncluster
      mark_neighbor_nums(entry_list,i,Ncluster,clus)
      clus.Nmember = len(clus.members)
      cluster_list.append(clus)


### 'S': Single Linkage Clustering without recursive call ### 
if (OPT['alg']=='S'):
   
  Ncluster = -1
  Nround = 0
  while (Nround==0) or (Nchange>0):
    Nchange = 0
    for x in (entry_list):
      if (x.cluster_num == -1):
        Ncluster += 1 
        x.cluster_num = Ncluster
      for i in (x.neighbor_nums):
        y = entry_list[i]
        if (y.cluster_num==-1) or (y.cluster_num > x.cluster_num):
          y.cluster_num = x.cluster_num
          Nchange += 1
    Nround += 1  
    print "#Nround %d Nchange %d Ncluster %d"%(Nround, Nchange,Ncluster)
  
  new_clus_num = {}
  cluster_list = []
  Ncluster = -1
  for x in (entry_list):
    if (new_clus_num.has_key(x.cluster_num)==0):
      Ncluster += 1 
      new_clus_num[x.cluster_num] = Ncluster
      clus = Cluster()
      clus.num = Ncluster
      cluster_list.append(clus)
  print "Ncluster %d"%(Ncluster)

  for x in (entry_list):
    x.cluster_num =  new_clus_num[x.cluster_num]
    cluster_list[x.cluster_num].members.append(x.num)
    cluster_list[x.cluster_num].Nmember += 1


### 'J': Jarvis-Patrick algorithm ### 
if (OPT['alg']=='J'):
  Ncluster = -1
  Nround = 0
  Nchange = 0
  while (Nround==0) or (Nchange>0):
    Nchange = 0
    for x in (entry_list):
      if (x.cluster_num==-1):
        Ncluster += 1
        x.cluster_num = Ncluster
      for y in (entry_list):
        if (x.num != y.num) and ((y.cluster_num==-1) or (y.cluster_num > x.cluster_num)):
          x_is_neiy = 0
          y_is_neix = 0
          Ncommon = 0
          for m in (x.neighbor_nums):
            if (m==y.num):
              y_is_neix  = 1

          if (y_is_neix == 1):
            for n in (y.neighbor_nums):
              if (n==x.num):
                x_is_neiy  = 1

          if (x_is_neiy==1) and (y_is_neix==1):
            for m in (x.neighbor_nums):
              for n in (y.neighbor_nums):
                if (m==n):
                  Ncommon += 1

          if (x_is_neiy==1) and (y_is_neix==1) and (Ncommon >= int(OPT['kmin']) ):
            Nchange += 1
            #print "ycluster_num %d -> %d"%(y.cluster_num,x.cluster_num)
            y.cluster_num = x.cluster_num
    Nround += 1  
    print "#Nround %d Nchange %d Ncluster %d"%(Nround, Nchange,Ncluster)

  new_clus_num = {}
  cluster_list = []
  Ncluster = -1
  for x in (entry_list):
    if (new_clus_num.has_key(x.cluster_num)==0):
      Ncluster += 1 
      new_clus_num[x.cluster_num] = Ncluster
      clus = Cluster()
      clus.num = Ncluster
      cluster_list.append(clus)
  print "Ncluster %d"%(Ncluster)

  for x in (entry_list):
    x.cluster_num =  new_clus_num[x.cluster_num]
    cluster_list[x.cluster_num].members.append(x.num)
    cluster_list[x.cluster_num].Nmember += 1



### [3] Output ###
#for i in range(len(entry_list)):
#  print ">%d [%d] %s Nnei %d"%(entry_list[i].num,entry_list[i].cluster_num,entry_list[i].name,len(entry_list[i].neighbor_nums))
#  for j in range(len(entry_list[i].neighbor_nums)):
#    k = entry_list[i].neighbor_nums[j] 
#    print "  %d [%d] %s %f"%(entry_list[k].num,entry_list[k].cluster_num,entry_list[k].name,entry_list[i].similarities[j])

set_Nneighbor_in_cluster(entry_list,cluster_list)

scluster_list = sorted(cluster_list,lambda x,y:cmp(y.Nmember,x.Nmember))

set_end_date(OPT)

if (OPT['ocl']!=''):
  write_clusters(OPT['ocl'],scluster_list,entry_list,OPT)
if (OPT['ore']!=''):
  write_representatives(OPT['ore'],scluster_list,entry_list,OPT)
