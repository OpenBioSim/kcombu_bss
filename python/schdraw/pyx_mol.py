##
## <pyx_mol.py>
##
##  functions for molecular drawing using PyX
##
## << cordinates of PyX >> (chimera and ChemDraw employ this scheme) 
## Y-direction : 'U'p
##
## y
## /\ 
##  |
##  |
##  ------> x
##

import sys
import os 
import math
from pyx import *
import molecule 
import molfunc

LastModDate = "Nov 11, 2013"




def draw_bonds(OPT,mol,cvs,sc,gxo,gyo,lwidth):
  for b in (mol.bonds):
    sa = b.natom1
    ea = b.natom2
    slabel = mol.atoms[sa].label
    elabel = mol.atoms[ea].label
    drawbond = 0
    if ((mol.atoms[sa].element != 'H') and (mol.atoms[ea].element != 'H')):
      drawbond = 1
    elif ((mol.atoms[sa].element == 'H') and (mol.atoms[ea].element != 'H')):
      if (slabel != '') or (mol.atoms[sa].circlesize != ''):
        drawbond = 1
    elif ((mol.atoms[ea].element == 'H') and (mol.atoms[sa].element != 'H')):
      if (elabel != '') or (mol.atoms[ea].circlesize != ''):
        drawbond = 1


    if (drawbond==1):
    
      sx = mol.atoms[sa].x
      sy = mol.atoms[sa].y
      ex = mol.atoms[ea].x
      ey = mol.atoms[ea].y
  
      sr = mol.atoms[sa].rgb[0]
      sg = mol.atoms[sa].rgb[1]
      sb = mol.atoms[sa].rgb[2]
  
      er = mol.atoms[ea].rgb[0]
      eg = mol.atoms[ea].rgb[1]
      eb = mol.atoms[ea].rgb[2]
    
      gsx = float(sc*(sx - mol.minX) + gxo)
      gsy = float(sc*(sy - mol.minY) + gyo)
      gex = float(sc*(ex - mol.minX) + gxo)
      gey = float(sc*(ey - mol.minY) + gyo)
      [gsx,gsy,gex,gey] = molfunc.font_radius_modified_gxy(OPT,gsx,gsy,gex,gey,slabel,elabel,OPT['pfontsize'])
   #   print "gsx %f gsy %f gex %f gey %f"%(gsx,gsy,gex,gey)
  
      if (OPT['updw']=='T') and ((b.bondstereo=='1') or (b.bondstereo=='6')):
        width_dbond = 0.14
        [gsx1,gsy1,gex1,gey1,gsx2,gsy2,gex2,gey2] = molfunc.double_bond_gxy(sc,sx,sy,ex,ey,mol.minX,mol.minY,gxo,gyo,width_dbond)
        [gsx1,gsy1,gex1,gey1] = molfunc.font_radius_modified_gxy(OPT,gsx1,gsy1,gex1,gey1,slabel,elabel,OPT['pfontsize'])
        [gsx2,gsy2,gex2,gey2] = molfunc.font_radius_modified_gxy(OPT,gsx2,gsy2,gex2,gey2,slabel,elabel,OPT['pfontsize'])
        if (b.bondstereo == '1'):
          arrow = path.path(path.moveto(gsx,gsy),path.lineto(gex1,gey1),path.lineto(gex2,gey2),path.closepath())
          cvs.fill(arrow,[style.linewidth(lwidth),color.rgb(0.0,0.0,0.0)])
        if (b.bondstereo == '6'):
          draw_dotted_triangle(cvs,gsx,gsy,gex1,gey1,gex2,gey2)
          #cvs.stroke(path.line(gsx,gsy,gex,gey),[style.linewidth(lwidth),color.rgb(1.0,0.0,0.0)])
          pass
      else: 
        if (OPT['btype']!='I') and ((b.bondtype=='2') or (b.bondtype=='3')):
          width_dbond = 0.07
          [gsx1,gsy1,gex1,gey1,gsx2,gsy2,gex2,gey2] = molfunc.double_bond_gxy(sc,sx,sy,ex,ey,mol.minX,mol.minY,gxo,gyo,width_dbond)
          [gsx1,gsy1,gex1,gey1] = molfunc.font_radius_modified_gxy(OPT,gsx1,gsy1,gex1,gey1,slabel,elabel,OPT['pfontsize'])
          [gsx2,gsy2,gex2,gey2] = molfunc.font_radius_modified_gxy(OPT,gsx2,gsy2,gex2,gey2,slabel,elabel,OPT['pfontsize'])
          cvs.stroke(path.line(gsx1,gsy1,gex1,gey1),[style.linewidth(lwidth),color.rgb(0.0,0.0,0.0)])
          cvs.stroke(path.line(gsx2,gsy2,gex2,gey2),[style.linewidth(lwidth),color.rgb(0.0,0.0,0.0)])
  
        if (b.bondtype!='2') or (OPT['btype']=='I'):
          cvs.stroke(path.line(gsx,gsy,gex,gey),[style.linewidth(lwidth),color.rgb(0.0,0.0,0.0)])




def draw_atoms(OPT,mol,cvs,sc,gxo,gyo):
  #gr     = float(sc*float(OPT['rcir']))
  gr_large   = float(sc*float(OPT['rcirl']))
  gr_medium  = float(sc*float(OPT['rcirm']))
  gr_small   = float(sc*float(OPT['rcirs']))
#  print "#def draw_atoms(OPT,mol,cvs,sc,gxo,gyo) sc %s gr %f %f %f"%(sc,gr_large,gr_medium,gr_small)
 
  grtxt  = float(OPT['pfontsize'])/2.0
  index = [i for i in range(mol.Natom)] 
  sindex = sorted(index,lambda x,y:cmp(mol.atoms[x].z,mol.atoms[y].z))
  ## smaller 'z'-atom is  draw firster, larger 'z'-atom is draw last.

  for i in (sindex):
    a = mol.atoms[i]
    gx = float(sc*(a.x - mol.minX) + gxo)
    gy = float(sc*(a.y - mol.minY) + gyo)

    scr = 0.0
    gr = 0.0
    if (a.circlesize == 'M'):
#      scr = 1.0
#      scr = 0.8
      gr = gr_medium
    elif (a.circlesize == 'L'):
#      scr = 2.0
#      scr = 1.2
      gr = gr_large
    elif (a.circlesize == 'S'):
#      scr  = 0.5
#      scr  = 0.6
      gr = gr_small

#    if (scr > 0.0):
#      cvs.fill(path.circle(gx,gy,gr*scr),[color.rgb(a.rgb[0],a.rgb[1],a.rgb[2])])
#      cvs.stroke(path.circle(gx,gy,gr*scr),[color.rgb(0,0,0)])
    if (gr > 0.0):
      cvs.fill(path.circle(gx,gy,gr),[color.rgb(a.rgb[0],a.rgb[1],a.rgb[2])])
      cvs.stroke(path.circle(gx,gy,gr),[color.rgb(0,0,0)])

    if (OPT['tcol']== 'A'):
      (rt,gt,bt) = (a.rgb[0],a.rgb[1],a.rgb[2])
    else:
      (rt,gt,bt) = (0.0,0.0,0.0)

    if (a.label != ''):
      #cvs.text(gx+grtxt,gy+grtxt,a.label,[color.rgb(rt,gt,bt)])
      #cvs.text(gx+grtxt,gy+grtxt,a.label,[color.rgb(0.0,0.0,0.0)])
      cvs.text(gx-grtxt,gy-grtxt,a.label,[color.rgb(0.0,0.0,0.0)])
      #print "cvs.text('%s')"%(a.label)

    if (a.sidelabel != ''):
      (lx,ly) = molfunc.find_empty_2D_point(mol,a.num)
     # cvs.text(gx+grtxt+3.0*scr*lx,gy+grtxt+3.0*scr*ly,mol.marked_atoms[a.num],[color.rgb(rt,gt,bt)])
    #  cvs.text(gx-grtxt+3.0*scr*lx,gy-grtxt+3.0*scr*ly,a.sidelabel,[color.rgb(0,0,0)])
   #   cvs.text(gx-grtxt+1.0*scr*lx,gy-grtxt+1.0*scr*ly,a.sidelabel,[color.rgb(0,0,0)])
      cvs.text(gx-grtxt+gr*lx,gy-grtxt+gr*ly,a.sidelabel,[color.rgb(0,0,0)])
  pass


def draw_dotted_triangle(cvs,gsx,gsy,gex1,gey1,gex2,gey2):
#   ge1*-----*ge2
#        --- 
#         -
#         * gs 
  v1x = gex1 - gsx
  v1y = gey1 - gsy

  v2x = gex2 - gsx
  v2y = gey2 - gsy

  Ndiv = 12 

  for i in range(Ndiv):
    x1 = gsx + i*float(v1x)/Ndiv
    y1 = gsy + i*float(v1y)/Ndiv

    x2 = gsx + i*float(v2x)/Ndiv
    y2 = gsy + i*float(v2y)/Ndiv

    cvs.stroke(path.line(x1,y1,x2,y2),[color.rgb(0.0,0.0,0.0)])
  pass



