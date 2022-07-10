##
## <slc.py>
##
##  for single linkage clustering
##

import sys
import os

import kcombu_func


LastModDate = "July 5, 2011"

# class LibraryMoleculeAVA:
#   def __init__(self):
#     self.mol  = ''  ## class Molecule (malloc later)
#     self.num       = 0
#     self.molname   = ''
#     self.molform   =  ''
#     self.class_str = ''
# def read_list_all_vs_all_file(ifname,libmolava_list,scmat,avadat):




class Cluster:
  def __init__(self):
    self.num      = 0   ## cluster number 
    self.Nmember  = 0   ## number of memmbers
    self.members  = []  ## a list of number of member entries
    self.rep_entry_num  = -1 ## number of representative entry


def set_neighbor_nums(libmol_list, scmat, scthre):
  for a in (libmol_list): 
    a.neighbor_nums   = []
    a.cluster_num = -1 
    for b in (libmol_list): 
      if (scmat[a.num][b.num]>=scthre):
        a.neighbor_nums.append(b.num) 


def mark_neighbor_nums(libmol_list,focus_i,cluster_num,cluster):
  #print "#mark_neighbor_nums(focus_i %d cluster_num %d)"%(focus_i,cluster_num)
  libmol_list[focus_i].cluster_num = cluster_num
  cluster.members.append(focus_i)
  for j in (libmol_list[focus_i].neighbor_nums):
    if (libmol_list[j].cluster_num==-1):
      mark_neighbor_nums(libmol_list,j,cluster_num,cluster)


def single_linkage_clustering(libmol_list,cluster_list):
  Ncluster = -1
  for i in range(len(libmol_list)):
    if (libmol_list[i].cluster_num==-1):
      Ncluster += 1
      clus = Cluster()
      clus.num = Ncluster
      mark_neighbor_nums(libmol_list,i,Ncluster,clus)
      clus.Nmember = len(clus.members)
      cluster_list.append(clus)


def write_clusters(ofname,cluster_list,libmol_list):
  print "#write_clusters()-->'%s'"%(ofname)
  of = open(ofname,"w")
  i = 0
  for c in (cluster_list):
    c.num = i
    of.write(">CLUSTER %d Nmember %d rep_entry %s\n"%(c.num,c.Nmember,libmol_list[c.rep_entry_num].molname))
    i += 1
    for j in (c.members):
      of.write("%s"%(libmol_list[j].molname))
      of.write("\n")
  of.close()



#############
### MAIN ####
#############


def _main():
  if (len(sys.argv)<2):
    print "slc.py  [ava similarity list file] [sc_thre(0.0--1.0)]"
    sys.exit(1)
 
  iavafile = sys.argv[1]
  sc_thre  = float(sys.argv[2])
  libmol_list = []

  scmat =  []
  avadat = {}
  kcombu_func.read_list_all_vs_all_file(iavafile,libmol_list,scmat,avadat)
  print "#Nmol_in_lib %d"%(len(libmol_list))
  set_neighbor_nums(libmol_list, scmat, sc_thre)
  cluster_list = []
  single_linkage_clustering(libmol_list,cluster_list)
  print "Ncluster %d"%(len(cluster_list))
  write_clusters("cluster.out",cluster_list,libmol_list)
if __name__ == '__main__': _main()
