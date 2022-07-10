#!/usr/bin/python2.5

##
## <avadraw.py>
##

import sys
import os 

from pyx import *

import math
import random

import molecule 
import molfunc
import pca 
import rmsd 
import kcombu_func
import mds
import numpy_mds

import image_mol
import pyx_mol

import slc

LastModDate = "July 9, 2011"

def draw_text_box(cvs,text_list,Xbottom_left,Ybottom_left,box_width,box_height):
# The box is defined as follows:
# 
# XXXX
# XXXX
# XXXX
# 'X':region for text draing
# (Xbottom_left, Ybottom_left) is located at the bottom left point of the box.
  Nline = len(text_list)
  text_height = box_height/Nline
  for i in range(len(text_list)):
    textstr = text_list[i]
    textstr = text.replace('_','\_') 
    cvs.text(Xbottom_left,Ybottom_left + box_height - text_height*i ,textstr,[color.rgb(0,0,0)])


def draw_molecular_box(cvs,mol,textstr,Xbottom_left,Ybottom_left,box_width,box_height,molmargin,textheight,ave_scale,lwidth,red=1.0,green=1.0,blue=1.0):

# The box is defined as follows:
#
#  ttttt
#  ooooo
#  oXXXo   
#  oXXXo
#  ooooo
#
# 'X':region for molecular drawing, 'o':molmargin, 't':textheight,
#
# (Xbottom_left, Ybottom_left) is located at the bottom left point of the box.
#  
#  (gminX,gminY) is defined as the bottom left point in XXX.
#                                                       XXX
#
  gminX = Xbottom_left + molmargin
  gminY = Ybottom_left + molmargin
  mol_scale = mol.cal_scale_for_gwidth_gheight(box_width-2*molmargin,box_height-2*molmargin-textheight)
  if (mol_scale > ave_scale):
    scale = ave_scale
  else: 
    scale = mol_scale
  cvs.fill(path.rect(Xbottom_left,Ybottom_left,box_width,box_height),[color.rgb(red,green,blue)])
  pyx_mol.draw_bonds(OPT,mol,cvs,scale,gminX,gminY,lwidth)
  pyx_mol.draw_atoms(OPT,mol,cvs,scale,gminX,gminY)
  textstr = textstr.replace('_','\_') 
  #cvs.text(Xbottom_left,Ybottom_left + box_height - textheight,textstr,[color.rgb(0,0,0)])
  cvs.text(Xbottom_left,Ybottom_left + box_height - textheight,textstr,[text.size(sizename=OPT['fontsize'])])
  #cvs.text(2, 2, "Hello, world!",[text.size(sizename="Huge")])
  pass



def draw_frame_of_box(cvs,Xbottom_left,Ybottom_left,box_width,box_height,p,x,y,mnum_by_pos,scmat,simthre1,simthre2):
# The box is defined as follows:
#
#  ttttt
#  ooooo
#  oXXXo   
#  oXXXo
#  ooooo
#
# 'X':region for molecular drawing, 'o':molmargin, 't':textheight,
#
# (Xbottom_left, Ybottom_left) is located at the bottom left point of the box.
#  
#  (gminX,gminY) is defined as the bottom left point in XXX.
#                                                       XXX
#
  up = bottom = left = right = 1
  mup = mbottom = mleft = mright = -1 
  
  m = mnum_by_pos[p][x][y] 
  xn = x + 1
  xb = x - 1
  yn = y + 1
  yb = y - 1
  #lw = 0.02
  bxL = Xbottom_left 
  bxR = Xbottom_left + box_width
  byU = Ybottom_left + box_height
  byB = Ybottom_left

  lw_bold  = 0.05
  lw_thin = lw_bold * 0.1

  lw = lw_bold 
  if (yb>=0):
    mup     = mnum_by_pos[p][x][yb]
    #print "m %d mup %d"%(m,mup) 
    if (mup>=0) and (scmat[m][mup]>=simthre1):
      up = 0
    if (mup>=0) and (scmat[m][mup]>=simthre2):
      lw = lw_thin

  if (up==1):  
    cvs.stroke(path.line(bxL,byU,bxR,byU),[style.linewidth(lw),color.rgb(0.0,0.0,0.0)])


  lw = lw_bold 
  if (yn<Nrow):
    mbottom = mnum_by_pos[p][x][yn]
    if (mbottom>=0) and (scmat[m][mbottom]>=simthre1):
      bottom = 0
    if (mbottom>=0) and (scmat[m][mbottom]>=simthre2):
      lw = lw_thin

  if (bottom==1):  
    cvs.stroke(path.line(bxL,byB,bxR,byB),[style.linewidth(lw),color.rgb(0.0,0.0,0.0)])
  
  lw = lw_bold 
  if (xb>=0):
    mleft   = mnum_by_pos[p][xb][y]
    if (mleft>=0) and (scmat[m][mleft]>=simthre1):
      left = 0
    if (mleft>=0) and (scmat[m][mleft]>=simthre2):
      lw = lw_thin

  if (left==1):  
    cvs.stroke(path.line(bxL,byU,bxL,byB),[style.linewidth(lw),color.rgb(0.0,0.0,0.0)])

  lw = lw_bold 
  if (xn<Ncol):
    mright  = mnum_by_pos[p][xn][y]
    if (mright>=0) and (scmat[m][mright]>=simthre1):
      right = 0
    if (mright>=0) and (scmat[m][mright]>=simthre2):
      lw = lw_thin
 
 
  if (right==1):  
    cvs.stroke(path.line(bxR,byU,bxR,byB),[style.linewidth(lw),color.rgb(0.0,0.0,0.0)])
  pass




def total_neighbor_similarities(mnum_by_pos,Npage,Ncol, Nrow,scmat):

  S = 0.0
  for p in range(Npage):
    for x in range(Ncol):
      for y in range(Nrow):
        m  = mnum_by_pos[p][x][y]
        x1 = x + 1
        y1 = y + 1
        if (x1<Ncol) and (mnum_by_pos[p][x1][y]>=0):
          mx1 = mnum_by_pos[p][x1][y]
          S += scmat[m][mx1]
        if (y1<Nrow) and (mnum_by_pos[p][x][y1]>=0):
          my1 = mnum_by_pos[p][x][y1]
          S += scmat[m][my1]
  return(S)


def monte_carlo_optimization(mnum_by_pos,Npage,Ncol, Nrow,scmat, Nrepeat):
  print "#monte_carlo_optimization(mnum_by_pos,Npage,Ncol, Nrow,scmat, Nrepeat:%d)"%(Nrepeat)
  mnum_by_pos_new = [[[-1 for i in range(Nrow)] for j in range(Ncol)] for k in range(Npage)]


  for p in range(Npage):
    for x in range(Ncol):
      for y in range(Nrow):
        mnum_by_pos_new[p][x][y] = mnum_by_pos[p][x][y]

  S = total_neighbor_similarities(mnum_by_pos,Npage,Ncol,Nrow,scmat)

  for r in range (Nrepeat):
    p  = random.randint(0,Npage-1) 
    xa = random.randint(0,Ncol-1) 
    ya = random.randint(0,Nrow-1) 
    xb = random.randint(0,Ncol-1) 
    yb = random.randint(0,Nrow-1) 
    #xb = xa + random.randint(0,2)-1  
    #yb = ya + random.randint(0,2)-1 
    if (xb>=0) and (yb>=0) and (xb<Ncol) and (yb<Nrow) and ((xa!=xb) or (ya!=yb)):
      ma = mnum_by_pos[p][xa][ya]
      mb = mnum_by_pos[p][xb][yb]
      mnum_by_pos_new[p][xa][ya] = mb 
      mnum_by_pos_new[p][xb][yb] = ma 
      Snew = total_neighbor_similarities(mnum_by_pos_new,Npage,Ncol,Nrow,scmat)
      print "[%d] p %d a %d (%d %d) b %d (%d %d) Snew %f S %f"%(r,p,ma,xa,ya,mb,xb,yb,Snew,S)
      if (Snew>S):
        print "Accept !!"
        mnum_by_pos[p][xa][ya] = mb 
        mnum_by_pos[p][xb][yb] = ma 
        S = Snew 



def set_pos_greedy_from_xy_mdsrank(libmol_list,mnum_by_pos,Npage,Ncol,Nrow):

  print "#set_pos_greedy_from_xy_mdsrank(libmol_list,mnum_by_pos,Npage,Ncol,Nrow)"
  gridpnt_list = []
  Ngrid = 0
  for i in range(Ncol): 
    for j in range(Nrow): 
      gridpnt = {}
      gridpnt['num'] =  Ngrid 
      gridpnt['x'] = i 
      gridpnt['y'] = j 
      gridpnt_list.append(gridpnt)
      Ngrid += 1

  pair_list = [] 
  for libmol in (libmol_list):
    for gridpnt in (gridpnt_list):
      dx = libmol.x_mdsrank - gridpnt['x']
      dy = libmol.y_mdsrank - gridpnt['y']
      pair = {}
      pair['gnum'] = gridpnt['num'] 
      pair['mnum'] = libmol.num
      pair['dis']  = math.sqrt(dx*dx + dy*dy)
      pair['used'] = 0
      pair_list.append(pair)

  spair_list = sorted(pair_list,lambda x,y:cmp(x['dis'],y['dis']))

  for pair in (spair_list):
    libmol = libmol_list[pair['mnum']]
    gridpnt = gridpnt_list[pair['gnum']]
    #print "%d %d (%f %f) (%f %f) dis %f"%(pair['mnum'],pair['gnum'],libmol.x_mdsrank, libmol.y_mdsrank,gridpnt['x'],gridpnt['y'],pair['dis']) 

  libmol_use = [0 for i in range(Nlib)]

  for r in range (Ncol*Nrow):
    i = 0
    find = 0
    while (i<len(spair_list)) and (find==0):
      libmol = libmol_list[spair_list[i]['mnum']]
      gridpnt = gridpnt_list[spair_list[i]['gnum']]
      if (spair_list[i]['used']==0) and (libmol_use[libmol.num]==0) and (mnum_by_pos[0][gridpnt['x']][gridpnt['y']]==-1):
        mnum_by_pos[0][gridpnt['x']][gridpnt['y']]= libmol.num
        libmol_use[libmol.num]= 1
        find = 1 
        spair_list[i]['used'] =1  
        print "greedy:mnum_by_pos[%d][%d]-->'%s' dis %f"%(gridpnt['x'],gridpnt['y'],libmol.molname,spair_list[i]['dis'])
      i += 1


def set_positions_by_clustering_and_1D_MDS(mnum_by_pos,libmol_list,cluster_list,Ncol,Nrow):
  print "#set_positions_by_clustering_and_1D_MDS():"
  ## [1] set averaged MDS x-value into each cluster
  for clus in (cluster_list):
    clus.aveX = 0.0
    for num in (clus.members): 
      clus.aveX += libmol_list[num].posX
      #print "  %d %f"%(num,libmol_list[num].posX)
    clus.aveX /= len(clus.members)

  cindex = [i for i in range(len(cluster_list))]
  scindex = sorted(cindex,lambda x,y:cmp(cluster_list[x].aveX,cluster_list[y].aveX))
 
  Nrank = 0   
  for i in (scindex):
    clus = cluster_list[i]
    print ">cluster [%d] Nmember %d xpos %f Nrank %d"%(i,len(cluster_list[i].members),cluster_list[i].aveX,Nrank)
    for num in (clus.members): 
      libmol_list[num].rank = Nrank
      Nrank += 1 
 
  for mol in (libmol_list):
    p         = (mol.rank)/(Ncol*Nrow)
    page_pos  = (mol.rank)%(Ncol*Nrow)
    x  = page_pos % Ncol 
    y  = page_pos / Ncol
    ## for zig-zag location ##
    if ((y%2)==1):
      x = Ncol - x -1
    #print "x %d Ncol %d"%(x,Ncol) 
    mnum_by_pos[p][x][y] = mol.num 
  pass



def rgb_from_cluster_num(c):
  if ((lib.cluster_num%8)==0):
    r = 1.0
    g = 0.9
    b = 0.9
  elif ((lib.cluster_num%8)==1):
    r = 0.9
    g = 1.0
    b = 0.9
  elif ((lib.cluster_num%8)==2):
    r = 0.9
    g = 0.9
    b = 1.0
  elif ((lib.cluster_num%8)==3):
    r = 0.9
    g = 1.0
    b = 1.0
  elif ((lib.cluster_num%8)==4):
    r = 1.0
    g = 0.9
    b = 1.0
  elif ((lib.cluster_num%8)==5):
    r = 1.0
    g = 1.0
    b = 0.9
  elif ((lib.cluster_num%8)==6):
    r = 0.9
    g = 0.9
    b = 0.9
  elif ((lib.cluster_num%8)==7):
    r = 1.0
    g = 0.8
    b = 0.8
  return([r,g,b])


#############
#### MAIN ###
#############

OPT = {}
 

## read options ##
OPT['isdf']  = '' 
OPT['ipdb']  = '' 
OPT['scp']  = 1.0
OPT['mgn'] = 2.0
OPT['romgn'] = 0.1 
OPT['thgt'] = 0.5
OPT['mk'] = ''
OPT['op'] = 'out.pdf'
OPT['cir'] = 'M'
OPT['txt'] = 'F'
OPT['col'] = 'C'
OPT['tcol'] = 'C'
OPT['bcol'] = 'C'
OPT['iam'] = ''
OPT['ttl'] = 'T'
OPT['pca'] = 'F'
OPT['rms'] = 'F'
OPT['btype'] = 'N';
OPT['lw'] = '0.01'
OPT['rcir'] = '0.3'
OPT['fontsize'] = 0.3 
OPT['isc'] = ''
OPT['fQ'] = ''
OPT['fL'] = ''
OPT['lay'] = 'L'
OPT['ncol'] = 2
OPT['nrow'] = 4
OPT['iava'] = ''
OPT['nrep'] = 100
OPT['sthL'] = 0.9
OPT['sthS'] = 0.5
OPT['init'] = 'F'
OPT['eig'] = 'N'
OPT['slc'] = 'F'
OPT['thslc'] = 0.8
OPT['fontsize']= 'normalsize'
OPT['pfontsize']= 0.3 

if (len(sys.argv)<2):
  print "avadraw.py <options>"
  print " make image file from molecules"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " -iava  : input kcombu search all-vs-all comparison  file [%s]"%(OPT['iava']) 
  print " -lay   : lay out style. 'L'inear, 'M'ds [%s]"%(OPT['lay']) 
  print " -nrep  : number of repeat for monte carlo [%d]"%(OPT['nrep'])
  print " -sthL  : similarity threshold for drawing the line [%f]"%(OPT['sthL'])
  print " -sthS  : similarity threshold for drawing the line [%f]"%(OPT['sthS'])
  print " -init  : initial configuration 'F'ile_order, 'M'ds [%s]"%(OPT['init'])
  print " -pca   : xy-flat orientation by PCA ('T' or 'F') [%s]"%(OPT['pca'])
  print " -rms   : rmsd-rotation for molA onto molB ('T' or 'F') [%s]"%(OPT['rms'])
  print " -eig   : eigen value calculation 'N'umpy, 'P'ython [%s]"%(OPT['eig']) 
  print " -slc   : do single linkage clustering (T or F) [%s]"%(OPT['slc']) 
  print " -thslc : threshold for single linkage clustering [%f]"%(OPT['thslc']) 
  print "<options for display>"
  print " -mgn   : margin (cm) [%f]"%(OPT['mgn']) 
  print " -romgn : ratio against BoxWidth of margin for each molecule  [%f]"%(OPT['romgn']) 
  print " -thgt  : height for molecular title string (cm) [%f]"%(OPT['thgt']) 
  print " -cir   : draw circle on the atom. 'H':heavyatom, 'M'arked_atom [%s]"%(OPT['cir'])
  print " -rcir  : radius of circle (pixel/angstrom) [%s]"%(OPT['rcir'])
  print " -txt   : draw text on the atom. 'A'll_atom_element, 'N':element only non-carbon, 'M'atch_num,'n'um_in_file 'F'alse, [%s]"%(OPT['txt'])
  print " -col   : coloring: '-':black_all 'C'pk_color, 'G'ray,'M'arked color(red),'R'ainbow_marked_color [%s]"%(OPT['col'])
  print " -bcol  : bond coloring: 'b':black_all 'A'tom_color [%s]"%(OPT['bcol'])
  print " -tcol  : text coloring: 'b':black_all 'A'tom_color [%s]"%(OPT['tcol'])
  print " -btype : bond type(single,double,...) display 'I'gnore [%s]"%(OPT['btype'])
  print " -lw    : line width [%s]"%(OPT['lw'])
  print " -ttl   : show title (T or F)[%s]"%(OPT['ttl']) 
  print "<options only for pdf/eps (G='P')>"
  print " -op    : output pdf/eps file) [%s]"%(OPT['op'])
  print " -scp   : scale (pixel/angstrom) for PyX   [%f]"%(OPT['scp']) 
  print " -fontsize  : fontsize for pyx [%s]"%(OPT['fontsize'])
  print "   (normalsize, large, Large, LARGE, huge, Huge, tiny, scriptsize, footnotesize, small)" 
  print " -fQ    : file type for moleculeQ. 'S'df,'P'db 'K'cf,'2':mol2[%s]"%(OPT['fQ'])
  print " -fL    : file type for moleculeL. 'S'df,'P'db 'K'cf,'2':mol2[%s]"%(OPT['fL'])
  print " -ncol  : number of column for lay='M' [%d]"%(OPT['ncol']) 
  print " -nrow  : number of row    for lay='M' [%d]"%(OPT['nrow']) 
  sys.exit(1)

molfunc.read_option(sys.argv,OPT)
  
OPT['scp']  = float(OPT['scp'])
OPT['mgn'] = float(OPT['mgn'])
OPT['romgn'] = float(OPT['romgn'])
OPT['thgt'] = float(OPT['thgt'])
### [1] Read kcomu Search Result file (OPT['isc']) ##

if (OPT['iava'] != ''):
  libmol_list = [] 
  scmat =  [] 
  avadat = {}
  kcombu_func.read_list_all_vs_all_file(OPT['iava'],libmol_list,scmat,avadat)
  print "#Nmol_in_lib %d"%(len(libmol_list))
else:
  print "#ERROR:input search result file (-isc) is obligatory."
  sys.exit(1)

Nlib = len(libmol_list)
print "#Nlib %d"%(Nlib)
print "LIB_DIR  %s"%(avadat['LIBRARY_DIRECTORY']) 
print "LIB_FILE %s"%(avadat['LIBRARY_FILE']) 


ave_WX =  0.0
ave_WY =  0.0

Ncol = int(OPT['ncol'])
Nrow = int(OPT['nrow']) 
Npage = int(math.ceil(float(Nlib)/float(Ncol*Nrow)))

### [3] Read library molecules ##
for lib in (libmol_list):
  lib.mol = molecule.Molecule()
  fname = "%s/%s"%(avadat['LIBRARY_DIRECTORY'],lib.molname)
  if (molfunc.read_molecule_in_various_formats(lib.mol,fname,OPT['fL'])==0):
    print "#ERROR:file type for '%s' cannot be determined."%(fname)
    sys.exit(1)
  print "#%s Natom %d"%(lib.molname, lib.mol.Natom)
  lib.mol.rotate_180_degree_around_Xaxis()
  lib.mol.marked_atoms = ['' for i in range(lib.mol.Natom)]
  molfunc.set_rgb_to_atoms(lib.mol,OPT['col'],1.0)
  lib.mol.set_min_max_XY_molecule()
  ave_WX += lib.mol.maxX - lib.mol.minX
  ave_WY += lib.mol.maxY - lib.mol.minY

if (len(libmol_list)>0):
 ave_WX /= (len(libmol_list))
 ave_WY /= (len(libmol_list))

mnum_by_pos = [[[-1 for i in range(Nrow)] for j in range(Ncol)] for k in range(Npage)]

### Do SLC ###
if (OPT['slc']=='T'):
  slc.set_neighbor_nums(libmol_list, scmat, float(OPT['thslc']))
  cluster_list = []
  slc.single_linkage_clustering(libmol_list,cluster_list)
  print "Ncluster %d"%(len(cluster_list))
  slc.write_clusters("cluster.out",cluster_list,libmol_list)

#### Do MDS ###
xpos = [0.0 for i in range(Nlib)]
ypos = [0.0 for i in range(Nlib)]
zpos = [0.0 for i in range(Nlib)]
if (OPT['eig']=='P'):
  mds.metric_MDS(Nlib,scmat,xpos,ypos,zpos,Sim2Dis='T')
else:
  numpy_mds.metric_MDS(Nlib,scmat,xpos,ypos,zpos,Sim2Dis='T')
olabfile = "label.gnu"
oxyzfile = "xyz.dat"
print "output_xyz()-->'%s'"%(oxyzfile)
print "output_label()-->'%s'"%(olabfile)
for i in range(len(libmol_list)):
  libmol_list[i].posX = xpos[i] 
  libmol_list[i].posY = ypos[i] 
  libmol_list[i].posZ = zpos[i] 

### set lattice position by MDS ###
index = [ i for i in range(Nlib)]
xindex = sorted(index,lambda x,y:cmp(xpos[x],xpos[y]))
index = [ i for i in range(Nlib)]
yindex = sorted(index,lambda x,y:cmp(ypos[x],ypos[y]))

for i in range(Nlib):
  ix = xindex[i]
  iy = yindex[i]
  libmol_list[ix].x_mdsrank = float(i)/float(Nrow);
  libmol_list[iy].y_mdsrank = float(i)/float(Ncol);

for i in range(Nlib):
  print "%f %f %f %f"%(xpos[i],ypos[i],libmol_list[i].x_mdsrank, libmol_list[i].y_mdsrank)

oxyz = open(oxyzfile,"w")
olb = open(olabfile,"w")
for i in range(len(libmol_list)):
  #oxyz.write("%f %f %f %s\n"%(xpos[i],ypos[i],zpos[i],libmol_list[i].molname))
  #olb.write("set label '%s' at %f,%f,%f\n"%(libmol_list[i].molname,xpos[i],ypos[i],zpos[i]))
  oxyz.write("%f %f %f %f %f %s\n"%(xpos[i],ypos[i],zpos[i],libmol_list[i].x_mdsrank, libmol_list[i].y_mdsrank,libmol_list[i].molname))
  olb.write("set label '%s' at %f,%f\n"%(libmol_list[i].molname,libmol_list[i].x_mdsrank, libmol_list[i].y_mdsrank))
olb.close()
oxyz.close()


#### [4] Set up postion of molecules
print "Nlib %d Ncol %d Nrow %d Npage %d"%(Nlib,Ncol,Nrow,Npage)

if (OPT['init']=='F'):
  print "#initialize_file_order"
  for lib in (libmol_list):
    p         = (lib.num)/(Ncol*Nrow)
    page_pos  = (lib.num)%(Ncol*Nrow)
    x  = page_pos % Ncol 
    y  = page_pos / Ncol
    print "mnum %d p %d/%d x %d/%d y %d/%d"%(lib.num,p,Npage,x,Ncol,y,Nrow)   
    print "p x y %d %d %d mnum %d"%(0,2,4,mnum_by_pos[0][2][4])
    mnum_by_pos[p][x][y] = lib.num   
if (OPT['init']=='M'):
  set_pos_greedy_from_xy_mdsrank(libmol_list,mnum_by_pos,Npage,Ncol,Nrow)
    
#for p in range(Npage): 
#  for x in range(Ncol): 
#    for y in range(Nrow): 
#      print "p x y %d %d %d mnum %d"%(p,x,y,mnum_by_pos[p][x][y])

S = total_neighbor_similarities(mnum_by_pos,Npage,Ncol, Nrow,scmat)
print "#InitScore %f"%(S)
if (int(OPT['nrep'])>0):
  monte_carlo_optimization(mnum_by_pos,Npage,Ncol, Nrow,scmat, int(OPT['nrep']))

if (OPT['slc']=='T'):
  set_positions_by_clustering_and_1D_MDS(mnum_by_pos,libmol_list,cluster_list,Ncol,Nrow)

#### [5] Making PDF file ####

pagelist = [] 
PaperWidth  = 21.00 ## cm (for A4) ##
PaperHeight = 29.70 ## cm (for A4) ##



#
# >> Layout Scheme <<
# 'm':OPT['mgn']  (page margin) 
# 't':OPT['thgt'] (text height)
# 'o':OPT['omgn'] (molecular margin)
# 'X':Region for molecular Drawing
#
# Ncol=2, Nrow = 3
#
# mmmmmmmmmmmmm
# mtttttttttttm
# mooooooooooom
# moXXXooXXXXom
# moXXXooXXXXom
# mooooooooooom
# mtttttttttttm
# mooooooooooom
# moXXXooXXXXom
# moXXXooXXXXom
# mooooooooooom
# mtttttttttttm
# mooooooooooom
# moXXXooXXXXom
# moXXXooXXXXom
# mooooooooooom
# mmmmmmmmmmmmm
#
# This is rewritten as:
#
# mmmm 
# mBBm
# mBBm
# mBBm
# mmmm
#
# The box 'B' is defined as follows:
#  ttttt
#  ooooo
#  oXXXo
#  oXXXo
#  ooooo


pagelist = [] 

PaperWidth  = 21.00 ## cm (for A4) ##
PaperHeight = 29.70 ## cm (for A4) ##


if (OPT['lay']=='L'):
  BoxWidth  = (PaperWidth  - 2*OPT['mgn'])/Ncol
  BoxHeight = (PaperHeight - 2*OPT['mgn'])/Nrow
  MolMargin = float(BoxWidth) * float(OPT['romgn'])
  MolWidth   = BoxWidth  - 2*MolMargin
  MolHeight  = BoxHeight - 2*MolMargin - OPT['thgt']
  ave_scX = MolWidth/ave_WX
  ave_scY = MolHeight/ave_WY
  if (ave_scX<ave_scY):
    ave_sc = ave_scX
  else:
    ave_sc = ave_scY

  for p in range(Npage): 
    cvs = canvas.canvas()
    for x in range(Ncol): 
      for y in range(Nrow): 
        m = mnum_by_pos[p][x][y]
        if (m>=0):
          lib = libmol_list[m] 
          #text = "%s"%(lib.molname) 
          textstr = "%s[%d]"%(lib.molname,lib.cluster_num) 
          Xbottom_left = OPT['mgn'] + BoxWidth * x 
          Ybottom_left = OPT['mgn'] + BoxHeight * (Nrow - 1 - y)
          (r,g,b) = rgb_from_cluster_num(lib.cluster_num)
          if (cluster_list[lib.cluster_num].Nmember<2):
            r = 1.0
            g = 1.0
            b = 1.0

          draw_molecular_box(cvs,lib.mol,textstr,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,MolMargin,OPT['thgt'],ave_sc,float(OPT['lw']),red=r,green=g,blue=b)
          draw_frame_of_box(cvs,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,p,x,y,mnum_by_pos,scmat,float(OPT['sthL']),float(OPT['sthS']))
      
    textstr = "-%d-"%(len(pagelist)+1) 
    cvs.text(PaperWidth/2,OPT['mgn']/2,textstr,[color.rgb(0,0,0)])
    page = document.page(cvs,paperformat=document.paperformat.A4,centered=0)
    pagelist.append(page)
    print "#add page"
    cvs = ''


if (OPT['lay']=='M'):
  BoxWidth  = (PaperWidth  - 2*OPT['mgn'])/Ncol
  BoxHeight = (PaperHeight - 2*OPT['mgn'])/Nrow
  GBoxWidth  = BoxWidth * 0.7
  GBoxHeight = BoxHeight * 0.7
  MolMargin = float(GBoxWidth) * float(OPT['romgn'])
  MolWidth   = GBoxWidth  - 2*MolMargin
  MolHeight  = GBoxHeight - 2*MolMargin - OPT['thgt']
  ave_scX = MolWidth/ave_WX
  ave_scY = MolHeight/ave_WY
  if (ave_scX<ave_scY):
    ave_sc = ave_scX
  else:
    ave_sc = ave_scY
  cvs = canvas.canvas()
  for lib in (libmol_list):
    #textstr = "%s"%(lib.molname) 
    textstr = "%s[%d]"%(lib.molname,lib.cluster_num) 
    Xbottom_left = OPT['mgn'] + BoxWidth * lib.x_mdsrank 
    Ybottom_left = OPT['mgn'] + BoxHeight * (Nrow - 1 - lib.y_mdsrank)
    #print "ave_sc %f"%(ave_sc) 
    #ave_sc = 0.046
    cvs.stroke(path.rect(Xbottom_left,Ybottom_left,GBoxWidth,GBoxHeight),[color.rgb(0,0,0)])
    cvs.fill(path.rect(Xbottom_left,Ybottom_left,GBoxWidth,GBoxHeight),[color.rgb(1,1,1)])
    (r,g,b) = rgb_from_cluster_num(lib.cluster_num)
    if (cluster_list[lib.cluster_num].Nmember<2):
      r = 1.0
      g = 1.0
      b = 1.0

    draw_molecular_box(cvs,lib.mol,textstr,Xbottom_left,Ybottom_left,GBoxWidth,GBoxHeight,MolMargin,OPT['thgt'],ave_sc,float(OPT['lw']),red=r,green=g,blue=b)
    #draw_frame_of_box(cvs,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,p,x,y,mnum_by_pos,scmat,float(OPT['sthL']),float(OPT['sthS']))
      
  page = document.page(cvs,paperformat=document.paperformat.A4,centered=0)
  pagelist.append(page)
  print "#add page"

if (OPT['lay']=='S'):
  BoxWidth  = (PaperWidth  - 2*OPT['mgn'])/Ncol
  BoxHeight = (PaperHeight - 2*OPT['mgn'])/Nrow
  MolMargin = float(BoxWidth) * float(OPT['romgn'])
  MolWidth   = BoxWidth  - 2*MolMargin
  MolHeight  = BoxHeight - 2*MolMargin - OPT['thgt']
  ave_scX = MolWidth/ave_WX
  ave_scY = MolHeight/ave_WY
  if (ave_scX<ave_scY):
    ave_sc = ave_scX
  else:
    ave_sc = ave_scY

  for p in range(Npage): 
    cvs = canvas.canvas()
    for x in range(Ncol): 
      for y in range(Nrow): 
        m = mnum_by_pos[p][x][y]
        if (m>=0):
          lib = libmol_list[m] 
          #textstr = "%s"%(lib.molname) 
          textstr = "%s[%d]"%(lib.molname,lib.cluster_num) 
          Xbottom_left = OPT['mgn'] + BoxWidth * x 
          Ybottom_left = OPT['mgn'] + BoxHeight * (Nrow - 1 - y)
          draw_molecular_box(cvs,lib.mol,textstr,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,MolMargin,OPT['thgt'],ave_sc,float(OPT['lw']))
          draw_frame_of_box(cvs,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,p,x,y,mnum_by_pos,scmat,float(OPT['sthL']),float(OPT['sthS']))
      
    textstr = "-%d-"%(len(pagelist)+1) 
    cvs.text(PaperWidth/2,OPT['mgn']/2,textstr,[color.rgb(0,0,0)])
    page = document.page(cvs,paperformat=document.paperformat.A4,centered=0)
    pagelist.append(page)
    print "#add page"
    cvs = ''


#### [5] Making document and output as the PDF file ####

print "#Npage %d"%(len(pagelist))
d = document.document(pages=pagelist)
if (OPT['op'] != ''):
  print "#-->'%s'"%(OPT['op'])
  d.writePDFfile(OPT['op'])
