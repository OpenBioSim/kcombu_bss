##
## <svg_mol.py>
##
##  functions for molecular drawing using SVG (Scalable Vector Graphics)
##
## << cordinates of SVG >> (rasmol employ this scheme)
##  Ydirection: 'D'own
##   ------> x
##  | 
##  | 
##  \/ 
##  y
##

import sys
import os 
import math
import molecule 
import molfunc

LastModDate = "2019/11/14"



def draw_bonds(OPT,mol,lines,sc,gxo,gyo):

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
 
    #print "sa %s ea %s drawbond %d bondstereo '%s'"%(mol.atoms[sa].element,mol.atoms[ea].element,drawbond,b.bondstereo)
    if (drawbond==1):

      sx = mol.atoms[sa].x
      sy = mol.atoms[sa].y
      ex = mol.atoms[ea].x
      ey = mol.atoms[ea].y
  
      #gsx = int(sc*(sx - mol.minX) + gxo)
      #gsy = int(sc*(sy - mol.minY) + gyo)
      #gex = int(sc*(ex - mol.minX) + gxo)
      #gey = int(sc*(ey - mol.minY) + gyo)
      gsx = sc*(sx - mol.minX) + gxo
      gsy = sc*(sy - mol.minY) + gyo
      gex = sc*(ex - mol.minX) + gxo
      gey = sc*(ey - mol.minY) + gyo

      [gsx,gsy,gex,gey] = molfunc.font_radius_modified_gxy(OPT,gsx,gsy,gex,gey,slabel,elabel,OPT['fontsize'])


      if (OPT['updw']=='T') and ((b.bondstereo=='1') or (b.bondstereo=='6')):
        gwidth_dbond = 0.2 * sc
        [gsx1,gsy1,gex1,gey1,gsx2,gsy2,gex2,gey2] = molfunc.double_bond_gxy(sc,sx,sy,ex,ey,mol.minX,mol.minY,gxo,gyo,gwidth_dbond)
        [gsx1,gsy1,gex1,gey1] = molfunc.font_radius_modified_gxy(OPT,gsx1,gsy1,gex1,gey1,slabel,elabel,OPT['fontsize'])
        [gsx2,gsy2,gex2,gey2] = molfunc.font_radius_modified_gxy(OPT,gsx2,gsy2,gex2,gey2,slabel,elabel,OPT['fontsize'])
        ## 'Up'--> solid arrow, 'Down'--> dotted arrow. 
        if (b.bondstereo == '1'):
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black"/>'%(gsx,gsy,gex,gey))
          lines.append('<polyline points="%f %f , %f %f, %f %f" style="stroke black; fill:black"/>'%(gsx,gsy,gex1,gey1,gex2,gey2))
          #imdraw.polygon((gsx,gsy)+(gex1,gey1)+(gex2,gey2), outline=(mol.atoms[sa].rgb), fill=((mol.atoms[sa].rgb)))
        if (b.bondstereo == '6'):
          #lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black"/>'%(gsx,gsy,gex,gey))
          #lines.append('<polyline points="%f %f , %f %f, %f %f" style="stroke black; fill:black"/>'%(gsx,gsy,gex1,gey1,gex2,gey2))
          draw_dotted_triangle(lines,gsx,gsy,gex1,gey1,gex2,gey2,mol.atoms[sa].rgb)
          pass
      else: 
        if (b.bondtype=='2') or (b.bondtype=='3'):
          gwidth_dbond = 0.1 * sc
          [gsx1,gsy1,gex1,gey1,gsx2,gsy2,gex2,gey2] = molfunc.double_bond_gxy(sc,sx,sy,ex,ey,mol.minX,mol.minY,gxo,gyo,gwidth_dbond)
          [gsx1,gsy1,gex1,gey1] = molfunc.font_radius_modified_gxy(OPT,gsx1,gsy1,gex1,gey1,slabel,elabel,OPT['fontsize'])
          [gsx2,gsy2,gex2,gey2] = molfunc.font_radius_modified_gxy(OPT,gsx2,gsy2,gex2,gey2,slabel,elabel,OPT['fontsize'])
  
          #imdraw.line((gsx1,gsy1,int((gsx1+gex1)/2.0),int((gsy1+gey1)/2.0)),fill=mol.atoms[sa].rgb,width=1)
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" style="stroke: #%s;"/>'%(gsx1,gsy1,(gsx1+gex1)/2.0,(gsy1+gey1)/2.0,mol.atoms[sa].hexrgb))
          #imdraw.line((gex1,gey1,int((gsx1+gex1)/2.0),int((gsy1+gey1)/2.0)),fill=mol.atoms[ea].rgb,width=1)
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" style="stroke: #%s;"/>'%(gex1,gey1,(gsx1+gex1)/2.0,(gsy1+gey1)/2.0,mol.atoms[ea].hexrgb))
  
          #imdraw.line((gsx2,gsy2,int((gsx2+gex2)/2.0),int((gsy2+gey2)/2.0)),fill=mol.atoms[sa].rgb,width=1)
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" style="stroke: #%s;"/>'%(gsx2,gsy2,(gsx2+gex2)/2.0,(gsy2+gey2)/2.0,mol.atoms[sa].hexrgb))
          #imdraw.line((gex2,gey2,int((gsx2+gex2)/2.0),int((gsy2+gey2)/2.0)),fill=mol.atoms[ea].rgb,width=1)
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" style="stroke: #%s;"/>'%(gex2,gey2,(gsx2+gex2)/2.0,(gsy2+gey2)/2.0,mol.atoms[ea].hexrgb))
  
        if (b.bondtype!='2'):
          #imdraw.line((gsx,gsy,int((gsx+gex)/2.0),int((gsy+gey)/2.0)),fill=mol.atoms[sa].rgb,width=1)
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" style="stroke: #%s;"/>'%(gsx,gsy,(gsx+gex)/2.0,(gsy+gey)/2.0,mol.atoms[sa].hexrgb))
          #imdraw.line((gex,gey,int((gsx+gex)/2.0),int((gsy+gey)/2.0)),fill=mol.atoms[ea].rgb,width=1)
          lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" style="stroke: #%s;"/>'%(gex,gey,(gsx+gex)/2.0,(gsy+gey)/2.0,mol.atoms[ea].hexrgb))



def draw_atoms(OPT,mol,lines,sc,gxo,gyo):
  gr_large   = sc*float(OPT['rcirl'])
  gr_medium  = sc*float(OPT['rcirm'])
  gr_small   = sc*float(OPT['rcirs'])
  #txtshiftX = float(OPT['fontsize'])/2.0
  #txtshiftY = float(OPT['fontsize'])/2.0
  txtshiftX = float(OPT['fontsize'])/4.0
  txtshiftY = float(OPT['fontsize'])/4.0

  index = [i for i in range(mol.Natom)]
  sindex = sorted(index,lambda x,y:cmp(mol.atoms[y].z,mol.atoms[x].z))
  ## larger 'z'-atom is  draw firster, smaller 'z'-atom is draw last.

  for i in (sindex):
    a = mol.atoms[i]
    #gx = int(sc*(a.x - mol.minX) + gxo)
    #gy = int(sc*(a.y - mol.minY) + gyo)
    gx = sc*(a.x - mol.minX) + gxo
    gy = sc*(a.y - mol.minY) + gyo

    gr = 0
    if (a.circlesize == 'M'):
      gr = gr_medium
    elif (a.circlesize == 'L'):
      gr = gr_large
    elif (a.circlesize == 'S'):
      gr  = gr_small
    #gr_slbl    = int(gr * 3.0)
    #gr_slbl    = gr * 3.0
    gr_slbl    = gr * 2.0

    if (gr>0):
      #imdraw.ellipse((gx-gr,gy-gr,gx+gr,gy+gr),fill=a.rgb)
      lines.append('<circle cx="%f" cy="%f" r="%f" style="stroke: #%s; fill: #%s;"/>'%(gx,gy,gr,a.hexrgb,a.hexrgb))

    if (a.label != ''):
      #imdraw.text((gx-txtshiftX,gy-txtshiftY),a.label,font=myfont,fill=a.rgb)
      pass
    if (a.sidelabel != ''):
      (lx,ly) = molfunc.find_empty_2D_point(mol,a.num)
      #print "num_in_file %d lx %f ly %f"%(a.num_in_file,lx,ly)
      #imdraw.text((gx-txtshiftX+gr_slbl*lx,gy-txtshiftY+gr_slbl*ly),a.sidelabel,font=myfont,fill=a.rgb)
      #txtshiftX = 0
      #txtshiftY = 0
      #lx = 0
      #ly = 0
      ## text (x,y) is the bottom left postion of the font  
      lines.append('<text x="%f" y="%f" style="font-family: %s; font-size:%spx; stroke:none; fill:black">%s</text>'%(gx-txtshiftX+gr_slbl*lx, gy+txtshiftY+gr_slbl*ly, OPT['font'],OPT['fontsize'], a.sidelabel))
      #lines.append('<text x="%f" y="%f" style="font-family: %s; font-size:%spx; stroke:none; fill:black">%s</text>'%(gx+gr_slbl*lx, gy+gr_slbl*ly, OPT['font'],OPT['fontsize'], a.sidelabel))


def draw_dotted_triangle(lines,gsx,gsy,gex1,gey1,gex2,gey2,rgb):
#   ge1*-----*ge2
#        --- 
#         -
#         * gs 
  v1x = gex1 - gsx
  v1y = gey1 - gsy

  v2x = gex2 - gsx
  v2y = gey2 - gsy

  #Ndiv = 8
  Ndiv = 12

  for i in range(Ndiv):
    #x1 = int(gsx + i*float(v1x)/Ndiv) 
    #y1 = int(gsy + i*float(v1y)/Ndiv)
    #x2 = int(gsx + i*float(v2x)/Ndiv) 
    #y2 = int(gsy + i*float(v2y)/Ndiv)
    x1 = gsx + i*float(v1x)/Ndiv
    y1 = gsy + i*float(v1y)/Ndiv
    x2 = gsx + i*float(v2x)/Ndiv 
    y2 = gsy + i*float(v2y)/Ndiv

    lines.append('<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black"/>'%(x1,y1,x2,y2))
#    imdraw.line((x1,y1,x2,y2),fill=rgb,width=1)
  # imdraw.polygon((gsx,gsy)+(gex1,gey1)+(gex2,gey2), outline=(mol.atoms[sa].rgb), fill=((mol.atoms[sa].rgb)))
  pass


def save_lines(ofname,lines):
  if (ofname == '-'):
    of = sys.stdout
  else:
    print "#save_lines() --> '%s'"%(ofname)
    of = open(ofname,'w')
  for line in lines:
    of.write("%s\n"%(line))
  if (ofname != '-'):
    of.close()
  pass
