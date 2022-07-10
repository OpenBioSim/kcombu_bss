#!/usr/bin/env python

##
## <schdraw.py>
##

import sys
import os 

from pyx import *

import math

import molecule 
import molfunc
import pca 
import rmsd 
import kcombu_func

import image_mol
import pyx_mol

import extDBMkegg

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
    text = text_list[i]
    text = text.replace('_','\_') 
    cvs.text(Xbottom_left,Ybottom_left + box_height - text_height*i ,text,[color.rgb(0,0,0)])
   



def draw_molecular_box(cvs,mol,text,Xbottom_left,Ybottom_left,box_width,box_height,molmargin,textheight,ave_scale,lwidth):

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
  pyx_mol.draw_bonds(OPT,mol,cvs,scale,gminX,gminY,lwidth)
  pyx_mol.draw_atoms(OPT,mol,cvs,scale,gminX,gminY)
  text = text.replace('_','\_') 
  cvs.text(Xbottom_left,Ybottom_left + box_height - textheight,text,[color.rgb(0,0,0)])
  bxL = Xbottom_left 
  bxR = Xbottom_left + box_width
  byU = Ybottom_left
  byB = Ybottom_left + box_height
  cvs.stroke(path.line(bxL,byU,bxR,byU),[style.linewidth(0.02),color.rgb(0.0,0.0,0.0)])
  cvs.stroke(path.line(bxL,byB,bxR,byB),[style.linewidth(0.02),color.rgb(0.0,0.0,0.0)])
  cvs.stroke(path.line(bxL,byU,bxL,byB),[style.linewidth(0.02),color.rgb(0.0,0.0,0.0)])
  cvs.stroke(path.line(bxR,byU,bxR,byB),[style.linewidth(0.02),color.rgb(0.0,0.0,0.0)])
  pass

#############
#### MAIN ###
#############

OPT = {}
 

## read options ##
OPT['isdf']  = '' 
OPT['ipdb']  = '' 
OPT['scp']  = 1.0
OPT['mgn'] = 2.0
OPT['omgn'] = 0.5 
OPT['thgt'] = 0.5
OPT['mk'] = ''
OPT['opdf'] = 'out.pdf'
OPT['oeps'] = ''
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
OPT['lw'] = '0.05'
OPT['rcir'] = '0.3'
OPT['ifontsize'] = 1.0
OPT['pfontsize'] = 0.3 
OPT['isc'] = ''
OPT['fQ'] = ''
OPT['fL'] = ''
OPT['lay'] = 'Q'
OPT['ncol'] = 2
OPT['nrow'] = 4
OPT['kegg'] = 'F' 

if (len(sys.argv)<2):
  print "schdraw.py <options>"
  print " make image file from molecules"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " -isc   : input kcombu search result file [%s]"%(OPT['isc']) 
  print " -lay   : lay out style. 'O'ne column, 'M'ulti-column 'Q'uery-lib-two column[%s]"%(OPT['lay']) 
  print " -pca   : xy-flat orientation by PCA ('T' or 'F') [%s]"%(OPT['pca'])
  print " -rms   : rmsd-rotation for molA onto molB ('T' or 'F') [%s]"%(OPT['rms'])
  print "<options for display>"
  print " -mgn   : margin (cm) [%f]"%(OPT['mgn']) 
  print " -omgn  : margin for each molecule (cm) [%f]"%(OPT['omgn']) 
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
  print " -kegg  : show KEGG anotation data (T or F)[%s]"%(OPT['kegg']) 
  print "<options only for pdf/eps (G='P')>"
  print " -opdf  : output pdf file [%s]"%(OPT['opdf'])
  print " -oeps  : output eps file [%s]"%(OPT['oeps'])
  print " -scp   : scale (pixel/angstrom) for PyX   [%f]"%(OPT['scp']) 
  print " -pfontsize  : fontsize for pyx [%s]"%(OPT['pfontsize'])
  print " -fQ    : file type for moleculeQ. 'S'df,'P'db 'K'cf,'2':mol2[%s]"%(OPT['fQ'])
  print " -fL    : file type for moleculeL. 'S'df,'P'db 'K'cf,'2':mol2[%s]"%(OPT['fL'])
  print " -ncol  : number of column for lay='M' [%d]"%(OPT['ncol']) 
  print " -nrow  : number of row    for lay='M' [%d]"%(OPT['nrow']) 
  sys.exit(1)

molfunc.read_option(sys.argv,OPT)
  
OPT['scp']  = float(OPT['scp'])
OPT['mgn']  = float(OPT['mgn'])
OPT['omgn'] = float(OPT['omgn'])
OPT['lw']   = float(OPT['lw'])

### [1] Read kcomu Search Result file (OPT['isc']) ##

if (OPT['isc'] != ''):
  libmol_list = [] 
  schdat = {}
  kcombu_func.read_kcombu_search_result(OPT['isc'],libmol_list,schdat)
else:
  print "#ERROR:input search result file (-isc) is obligatory."
  sys.exit(1)

print "#Nlib %d"%(len(libmol_list)) 
print "LIB_DIR  %s"%(schdat['LIBRARY_DIRECTORY']) 
print "LIB_FILE %s"%(schdat['LIBRARY_FILE']) 

### [2] Read query molecule as molQ ##
molQ = molecule.Molecule()
molfunc.read_molecule_in_various_formats(molQ,schdat['QUERY_FILE'],OPT['fQ'])
molQ.marked_atoms = [i for i in range(molQ.Natom)]
molfunc.set_rgb_to_atoms(molQ,OPT['col'],1.0)
if (OPT['pca']=='T'):
  pca.rotate_molecule_XY_flat_plane_by_PCA(molQ)

print "molQ Natom %d Nbond %d"%(molQ.Natom, molQ.Nbond) 

molQ.set_min_max_XY_molecule()
ave_WX =  molQ.maxX - molQ.minX 
ave_WY =  molQ.maxY - molQ.minY

### [3] Read library molecules ##
for lib in (libmol_list):
  lib.mol = molecule.Molecule()
  fname = "%s/%s"%(schdat['LIBRARY_DIRECTORY'],lib.molname)
  if (molfunc.read_molecule_in_various_formats(lib.mol,fname,OPT['fL'])==0):
    print "#ERROR:file type for '%s' cannot be determined."%(fname)
    sys.exit(1)
  print "#%s Natom %d len_num_in_fileL %d"%(lib.molname, lib.mol.Natom, len(lib.num_in_fileL))
  lib.mol.marked_atoms = ['' for i in range(lib.mol.Natom)]
  lib.mol.rotate_180_degree_around_Xaxis()
  lib.numQ = []
  lib.numL = []
  kcombu_func.make_atom_num_list_from_num_in_file(molQ,   lib.numQ,lib.num_in_fileQ)
  kcombu_func.make_atom_num_list_from_num_in_file(lib.mol,lib.numL,lib.num_in_fileL)
  for i in range(len(lib.numL)):
    lib.mol.marked_atoms[lib.numL[i]] = "%d"%(i+1)
  molfunc.set_rgb_to_atoms(lib.mol,OPT['col'],1.0)
  if (OPT['rms']=='T'):
    rmsd.rotate_molA_onto_molB_for_matched_atoms_with_minRMSD(lib.mol, molQ, lib.numL, lib.numQ)
  lib.mol.set_min_max_XY_molecule()
  ave_WX += lib.mol.maxX - lib.mol.minX
  ave_WY += lib.mol.maxY - lib.mol.minY

if (len(libmol_list)>0):
 ave_WX /= (len(libmol_list)+1)
 ave_WY /= (len(libmol_list)+1)


#### [4] Making PDF file ####

pagelist = [] 
PaperWidth  = 21.00 ## cm (for A4) ##
PaperHeight = 29.70 ## cm (for A4) ##

### 'O'ne-column layout ###

if (OPT['lay']=='O'):
  Nmol_in_page = 5 
  MolHeight = (PaperHeight - 2*OPT['mgn'])/Nmol_in_page
  MolWidth  =  PaperWidth  - 2*OPT['mgn']
  print "#Nmol_in_page %d MolHeight %f"%(Nmol_in_page, MolHeight)
  
  for rank in range(len(libmol_list)):
    x = libmol_list[rank]
    rank_page = rank%Nmol_in_page 
  
    [minX,maxX,minY,maxY] = x.mol.set_min_max_XY_molecule()
    scX = (MolWidth-2.0*OPT['omgn'])/(maxX-minX)
    scY = (MolHeight-2.0*OPT['omgn'])/(maxY-minY)
    if (scX<scY):
      sc = scX
    else:
      sc = scY
    if (rank_page==0):
      cvs = canvas.canvas()
      print "new canvas"
   
    gminX = OPT['mgn'] + OPT['omgn']
    gminY = PaperHeight - OPT['mgn'] - MolHeight -  MolHeight*rank_page + OPT['omgn']
   # sc = 0.5 
    print ">%d %s rank %d rank_page %d minX %f minY %f gminX %f gminY %f"%(rank,x.molname,x.rank, rank_page,minX,minY,gminX,gminY)
  #  print "#sc %f scX %f scY %f"%(sc,scX,scY) 
    pyx_mol.draw_bonds(OPT,x.mol,cvs,sc,gminX,gminY,OPT['lw'])
    pyx_mol.draw_atoms(OPT,x.mol,cvs,sc,gminX,gminY)
  
    cvs.stroke(path.line(gminX-OPT['omgn'],gminY-OPT['omgn'],PaperWidth-OPT['mgn'],gminY-OPT['omgn']),[style.linewidth(0.05),color.rgb(0.0,0.0,0.0)])
    cvs.text(gminX,gminY,x.molname,[color.rgb(0,0,0)])
  
    if (rank_page==(Nmol_in_page-1)) or (rank==(len(libmol_list)-1)):
      page = document.page(cvs,paperformat=document.paperformat.A4,centered=0)
      pagelist.append(page)
      cvs = ''
      print "#end_page rank %d Nmol_in_page-1 %d rank_page %d"%(rank,Nmol_in_page-1,rank_page)

### 'M'ulti-column layout ###

if (OPT['lay']=='M') or (OPT['lay']=='Q'):
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


  Ncol = int(OPT['ncol'])
  Nrow = int(OPT['nrow']) 
  pagelist = [] 
  PaperWidth  = 21.00 ## cm (for A4) ##
  PaperHeight = 29.70 ## cm (for A4) ##

  BoxWidth  = (PaperWidth  - 2*OPT['mgn'])/Ncol
  BoxHeight = (PaperHeight - 2*OPT['mgn'])/Nrow

  MolWidth   = BoxWidth  - 2*OPT['omgn']
  MolHeight  = BoxHeight - 2*OPT['omgn'] - OPT['thgt']
 
  ave_scX = MolWidth/ave_WX
  ave_scY = MolHeight/ave_WY
  if (ave_scX<ave_scY):
    ave_sc = ave_scX
  else:
    ave_sc = ave_scY

  print "Paper (%f %f) Box(%f %f) Mol(%f %f)"%(PaperWidth,PaperHeight,BoxWidth,BoxHeight,MolWidth,MolHeight)

  cvs = canvas.canvas()
  ## draw title ##
  text_list = [] 
  text_list.append("KCOMBU SEARCH RESULT")
  text_list.append("QUERY FILE: %s"%(schdat['QUERY_FILE']))
  text_list.append("LIBRARY FILE: %s"%(schdat['LIBRARY_FILE'])) 
  text_list.append("LIBRARY DIRECTORY: %s"%(schdat['LIBRARY_DIRECTORY'])) 
  text_list.append("THRE TANIMOTO FILTER:%s"%(schdat['THRE_TANIMOTO_FILTER'])) 
  text_list.append("THRE TANIMOTO OUT:   %s"%(schdat['THRE_TANIMOTO_OUT'])) 
  text_list.append("NMOL IN LIBRARY:%s"%(schdat['NMOL_IN_LIBRARY'])) 
  text_list.append("NMOL OVER THRE:%s"%(schdat['NMOL_OVER_THRE'])) 
  text_list.append("DATE START: %s"%(schdat['DATE_START']))
  Xbottom_left = OPT['mgn']
  Ybottom_left = PaperHeight - OPT['mgn'] - BoxHeight
  draw_text_box(cvs,text_list,Xbottom_left,Ybottom_left,BoxWidth, BoxHeight)
 
  ## draw the query molecule ##
  Xbottom_left = OPT['mgn']+BoxWidth*(Ncol-1)
  Ybottom_left = PaperHeight - OPT['mgn'] - BoxHeight
 
  draw_molecular_box(cvs,molQ,"[query]",Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,OPT['omgn'],OPT['thgt'],ave_sc,OPT['lw'])

  if (OPT['lay']=='M'):
    ## draw library molecules ##
    for rank in range(len(libmol_list)):
      lib = libmol_list[rank]
      page_pos = (rank + Ncol)%(Ncol*Nrow)
      ncol  = page_pos % Ncol 
      nrow  = page_pos / Ncol
      if (page_pos==0):
        cvs = canvas.canvas()
        print "#new canvas"

      text = "[%d]%s tanimoto=%.3f"%(rank+1,lib.molname,lib.tanimoto) 
      Xbottom_left = OPT['mgn'] + BoxWidth * ncol
      Ybottom_left = OPT['mgn'] + BoxHeight * (Nrow - 1 - nrow)
  
      draw_molecular_box(cvs,lib.mol,text,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,OPT['omgn'],OPT['thgt'],ave_sc,OPT['lw'])
  
      if (page_pos==(Ncol*Nrow-1)) or (rank==(len(libmol_list)-1)):
        text = "-%d-"%(len(pagelist)+1) 
        cvs.text(PaperWidth/2,OPT['mgn']/2,text,[color.rgb(0,0,0)])
        page = document.page(cvs,paperformat=document.paperformat.A4,centered=0)
        pagelist.append(page)
        print "#add page"
        cvs = ''

  if (OPT['lay']=='Q'):
    ## draw library molecules ##
    for rank in range(len(libmol_list)):
      lib = libmol_list[rank]
      page_pos = (rank + 1)%(Nrow)
      ncol  = 0
      nrow  = page_pos
      if (page_pos==0):
        cvs = canvas.canvas()
        print "#new canvas"

      if (OPT['kegg']=='T'):
        field = lib.molname.split('.')
        keggid = field[0]
        keggdat = {}
        extDBMkegg.extract_kegg_dbm_file_content(keggid,keggdat)
        print "keggid '%s' %s"%(keggid,keggdat['NAME'][0])
        text = "[%d](%.3f)%s %s"%(rank+1,lib.tanimoto,keggid,keggdat['NAME'][0]) 
      else:
        text = "[%d]%s tanimoto=%.3f"%(rank+1,lib.molname,lib.tanimoto) 
  
      Xbottom_left = OPT['mgn'] + BoxWidth * ncol
      Ybottom_left = OPT['mgn'] + BoxHeight * (Nrow - 1 - nrow)
  
      draw_molecular_box(cvs,lib.mol,text,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,OPT['omgn'],OPT['thgt'],ave_sc,OPT['lw'])
  
      molQ.marked_atoms = ['' for i in range(molQ.Natom)]
      for i in range(len(lib.numQ)):
        molQ.marked_atoms[lib.numQ[i]] = "%d"%(i+1)
      molfunc.set_rgb_to_atoms(molQ,OPT['col'],1.0)

      Xbottom_left = OPT['mgn'] + BoxWidth * 1 
      text = "[query] Nmatched atom %d/(%d,%d)"%(len(lib.numQ),lib.mol.Natom, molQ.Natom) 
      draw_molecular_box(cvs,molQ,text,Xbottom_left,Ybottom_left,BoxWidth,BoxHeight,OPT['omgn'],OPT['thgt'],ave_sc,OPT['lw'])

      if (page_pos==(Nrow-1)) or (rank==(len(libmol_list)-1)):
        text = "-%d-"%(len(pagelist)+1) 
        cvs.text(PaperWidth/2,OPT['mgn']/2,text,[color.rgb(0,0,0)])
        page = document.page(cvs,paperformat=document.paperformat.A4,centered=0)
        if (len(pagelist)==0):
          pagelist.append(page)
        #pagelist.append(page)
        print "#add page"
        cvs = ''

  pass






#### [5] Making document and output as the PDF file ####

print "#Npage %d"%(len(pagelist))
d = document.document(pages=pagelist)
if (OPT['opdf'] != ''):
  print "#-->'%s'"%(OPT['opdf'])
  d.writePDFfile(OPT['opdf'])
if (OPT['oeps'] != ''):
  print "#-->'%s'"%(OPT['oeps'])
  d.writeEPSfile(OPT['oeps'])
