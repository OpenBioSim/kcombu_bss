#
# <molfunc.py>
#
# This is mainly for "moldraw.cgi"
#

import sys
import os
import math
from datetime import datetime
import webfunc
import ligandbox_func
import kcombuweb_func

LastModDate = "2019/11/27"


def read_atom_match_file_in_simcomp(ifname,ABtype,rank,dat):
#>> FILE FORMAT EXAMPLE <<
#Input1 (6): alanin
#Input2 (5): glycin
#Result: 4	0.571429
#Alignment: 4 node (100.00%) identity in 4 node overlap
#1:[O6a]	2:[O6a]	## SCCS 1
#2:[O6a]	1:[O6a]	## SCCS 1
#6:[C6a]	5:[C6a]	## SCCS 1	## L= 4:[C1c]	## R= 4:[C1b]
#3:[N1a]	3:[N1a]	## SCCS 2	## L= 4:[C1c]	## R= 4:[C1b]
  if os.access(ifname,os.R_OK):
    f = open(ifname)
  else:
    print "#ERROR:Can't open '%s'"%(filename)
  f = open(ifname)

  matchlist = []
  dat['anumA'] = []
  dat['anumB'] = []
  Nrank = 0
  status = ' '

  for line in f:
    line = line.rstrip('\n')
    #print line,status  
    X = line.split()

#Alignment: 4 node (100.00%) identity in 4 node overlap
    if (line.startswith('#Input1')):
#Input1 (6): alanin
#Input2 (5): glycin
      fields = line.split(' ')
      dat['molA'] = fields[2]
    elif (line.startswith('#Input2')):
      fields = line.split(' ')
      dat['molB'] = fields[2]
    elif (status=='A') and (len(line)>10):
#1:[O6a]	2:[O6a]	## SCCS 1
#1:[O6a]\t2:[O6a]\t## SCCS 1
#2:[O6a]	1:[O6a]	## SCCS 1
#2:[O6a]\t1:[O6a]\t## SCCS 1
#6:[C6a]	5:[C6a]	## SCCS 1	## L= 4:[C1c]	## R= 4:[C1b]
#6:[C6a]\t5:[C6a]\t## SCCS 1\t## L= 4:[C1c]\t## R= 4:[C1b]
#3:[N1a]	3:[N1a]	## SCCS 2	## L= 4:[C1c]	## R= 4:[C1b]
      fields = line.split('\t')
      (atomnum,atomtype) = fields[0].split(':')
      dat['anumA'].append(int(atomnum)-1)
      (atomnum,atomtype) = fields[1].split(':')
      dat['anumB'].append(int(atomnum)-1)
      if (ABtype == 'A'):
        matchlist.append(atomnum)
      if (ABtype == 'B'):
        matchlist.append(atomnum)
    elif (line.startswith('Alignment')):
      status = 'A'

  f.close()
  #print "#matchlist_length %d"%(len(matchlist))
  return(matchlist)



def make_rotation_matrix_by_xyz_angle(th_x,th_y,th_z):

 Cx = math.cos(th_x*math.pi/180.0)
 Sx = math.sin(th_x*math.pi/180.0)
 Cy = math.cos(th_y*math.pi/180.0)
 Sy = math.sin(th_y*math.pi/180.0)
 Cz = math.cos(th_z*math.pi/180.0)
 Sz = math.sin(th_z*math.pi/180.0)

 Rx = [[0.0 for i in range(3)] for j in range(3)]
 Ry = [[0.0 for i in range(3)] for j in range(3)]
 Rz = [[0.0 for i in range(3)] for j in range(3)]

 Rx[0][0] = 1.0
 Rx[1][1] =  Cx
 Rx[1][2] = -Sx
 Rx[2][1] =  Sx
 Rx[2][2] =  Cx;

 Ry[0][0] = Cy
 Ry[0][2] = -Sy
 Ry[1][1] =  1.0
 Ry[1][2] = 0.0
 Ry[2][0] = Sy
 Ry[2][1] =  0.0
 Ry[2][2] = Cy

 Rz[0][0] = Cz
 Rz[0][1] = -Sz
 Rz[1][0] = Sz
 Rz[1][1] =  Cz
 Rz[2][2] = 1.0

 Ryx  = multiply_matrix3D(Ry,Rx)
 Rzyx = multiply_matrix3D(Rz,Ryx)
 return(Rzyx)


def multiply_matrix3D(A,B):
  C= [[0.0 for i in range(3)] for j in range(3)]
  for i in range(3):
    for j in range(3):
      C[i][j] = 0.0;
      for k in range(3):
        C[i][j] += A[i][k] * B[k][j]
  return(C)



def value_to_RGB(x,type,maxval):
  ## variable x,R,G,B : 0.0 <= [x,R,G,B] <= 1.0

  if (type=='RWB') or (type=='BGR') or (type=='WB') or (type=='RB'):
    x = 1.0 - x

  if (type=='BW') or (type=='WB'):
    R = x
    G = x
    B = x

  if (type=='BR') or (type=='RB'):
    R = x
    G = 0.0
    B = 1.0 - x

  if (type=='BWR') or (type=='RWB'):
    if (x<=0.5):
      R = 2.0*x
      G = 2.0*x
      B = 1.0
    elif (x<=1.0):
      R = 1.0
      G = 1.0-2.0*(x-0.5)
      B = 1.0-2.0*(x-0.5)

  if (type=='RGB') or (type=='BGR'):
    if (x<=0.25):
      R = 1.0
      G = 4.0*x
      B = 0.0
    elif (x<=0.50):
      R = 1.0 -4.0*(x-0.25)
      G = 1.0
      B = 0.0
    elif (x<=0.75):
      R = 0.0
      G = 1.0
      B = 4.0*(x-0.5)
    elif (x<=1.0):
      R = 0.0
      G = 1.0-4.0*(x-0.75)
      B = 1.0

  return(int(maxval*R),int(maxval*G),int(maxval*B))


def element_rgb_color(ele,maxval,coltype):


  if (coltype=='rgb'):
   ##
   ## basically taken from RasMol cpk color scheme

    if (ele == 'C'):
      rgb = [127,127,127]
   #   rgb = [0,0,0]
    elif (ele == 'O'):
      rgb = [255,0,0]
    elif (ele == 'N'):
      rgb = [100,100,255]
    elif (ele == 'S'):
      rgb = [255,230,0]
    elif (ele == 'P'):
      rgb = [255,170,0]
    elif (ele == 'F'):
      rgb = [218,165,32]
    elif (ele == 'He')or (ele == 'HE'):
      rgb = [255,192,203]
    elif (ele == 'Si')or (ele == 'Si'):
      rgb = [218,165,32]
    elif (ele == 'Au')or (ele == 'AU'):
      rgb = [218,165,32]
    elif (ele == 'Cl')or (ele == 'CL'):
      rgb = [0,255,0]
    elif (ele == 'Br')or (ele == 'BR'):
      rgb = [128,40,40]
    elif (ele == 'I'):
      rgb = [160,32,240]
    elif (ele == 'Mg')or (ele == 'MG'):
      rgb = [34,139,34]
    elif (ele == 'Ca')or (ele == 'CA'):
      rgb = [105,105,105]
    elif (ele == 'Ti')or (ele == 'TI'):
      rgb = [105,105,105]
    elif (ele == 'Cr')or (ele == 'CR'):
      rgb = [105,105,105]
    elif (ele == 'Mn')or (ele == 'MG'):
      rgb = [105,105,105]
    elif (ele == 'Ag')or (ele == 'AG'):
      rgb = [105,105,105]
    elif (ele == 'Fe')or (ele == 'FE'):
      rgb = [255,170,0]
    elif (ele == 'Ba')or (ele == 'BA'):
      rgb = [255,170,0]
    elif (ele == 'Ni')or (ele == 'NI'):
      rgb = [128,40,40]
    elif (ele == 'Cu')or (ele == 'CU'):
      rgb = [128,40,40]
    elif (ele == 'Zn')or (ele == 'ZN'):
      rgb = [128,40,40]
    elif (ele == 'Na')or (ele == 'NA'):
      rgb = [0,0,255]

    elif (ele == 'H'):
      rgb = [200,200,200]
    else:
      rgb = [250,22,145]



  if (coltype=='gray'):
    if (ele == 'C'):
   #   rgb = (160,160,160)
      rgb = [255,255,255]
    elif (ele == 'O'):
      rgb = [0,0,0]
    elif (ele == 'N'):
      rgb = [200,200,200]
    else:
      rgb = [128,128,128]

  if (maxval!=255):
    rgb[0] = rgb[0]*maxval/255.0
    rgb[1] = rgb[1]*maxval/255.0
    rgb[2] = rgb[2]*maxval/255.0
  RGB = (rgb[0],rgb[1],rgb[2])
  #print "RGB %s %s %s"%(rgb[0],rgb[1],rgb[2])
  return(RGB)



def set_rgb_to_atoms(mol,coltype,maxval):
  if (coltype=='R'):
    max_mark = 0.0
    for a in (mol.atoms):
     if (mol.marked_atoms[a.num]!= '') and (float(mol.marked_atoms[a.num])>max_mark):
       max_mark = float(mol.marked_atoms[a.num])

  for a in (mol.atoms):
    a.rgb = (0,0,0)
    if (coltype=='C'):
      a.rgb =  element_rgb_color(a.element,maxval,'rgb')
    elif (coltype=='G'):
      a.rgb =  element_rgb_color(a.atomtype,maxval,'gray')
    elif (coltype=='M'):
      if (mol.marked_atoms[a.num] != ''):
        a.rgb = (maxval,0,0)
      else:
        a.rgb = (0,0,0)
    elif (coltype=='R'):
      if (mol.marked_atoms[a.num] != ''):
        a.rgb = value_to_RGB(float(mol.marked_atoms[a.num])/max_mark,"BGR",maxval)
      else:
        a.rgb = (0,0,0)
    elif (coltype=='B'):
        a.rgb = (0,0,0)
    else:
        a.rgb = (0,0,0)

    ### added for SVG (T.Kawabata, 2019/10/24) ###
    a.hexrgb = "%s%s%s"%( int255_to_hex2(a.rgb[0]), int255_to_hex2(a.rgb[1]), int255_to_hex2(a.rgb[2]))


def int255_to_hex2(x):
  X = "%X"%(x)
  if (len(X)==1):
    X = '0' + X
  return(X)

def set_label_to_atoms(mol,labeltype,updwtype):
  #print "#set_label_to_atoms(labeltype %s updwtype %s)"%(labeltype, updwtype)
  for a in (mol.atoms):
    a.label = ''
    if (labeltype == 'A'):
      a.label = a.atomtype
    elif (labeltype == 'N')or (labeltype == 'I'):
      if (a.element !='C') and (a.element !='H'):
        a.label = a.element
      if (a.element == 'H'):
        NnoCnei = 0
        for i in (a.neighbors):
          if (mol.atoms[i].element != 'C'):
            NnoCnei += 1
        if (NnoCnei>0):
          a.label = a.element

        for b in (mol.bonds):
          if (b.bondstereo=='1') or (b.bondstereo=='6'):
            if (b.natom1 == a.num) or (b.natom2==a.num):
               a.label = a.element
    elif (labeltype == 'M'):
      a.label = mol.marked_atoms[a.num]
    elif (labeltype == 'n'):
      a.label = "%d"%(a.num_in_file)

  if (labeltype == 'I'):
    for b in (mol.bonds):
      sa = b.natom1
      ea = b.natom2
      anumH = -1
      if (mol.atoms[sa].element == 'H'):
       anumH    = sa
       anumNonH = ea
      if (mol.atoms[ea].element == 'H'):
       anumH    = ea
       anumNonH = sa
      if (anumH>=0): 
        if (updwtype=='T') and ((b.bondstereo=='1') or (b.bondstereo=='6')):
          mol.atoms[anumH].label = mol.atoms[anumH].element
        if (mol.atoms[anumNonH].element != 'C'):
          mol.atoms[anumH].label = mol.atoms[anumH].element
  #for a in (mol.atoms):
  #  print "%s '%s'"%(a.element,a.label)
  pass


def set_sidelabel_to_atoms(mol,labeltype):
  for a in (mol.atoms):
    a.sidelabel = ''
    if (labeltype == 'A'):
      a.sidelabel = a.atomtype
    elif (labeltype == 'N'):
      if (a.element !='C'):
        a.sidelabel = a.element
    elif (labeltype == 'M'):
      a.sidelabel = mol.marked_atoms[a.num]
    elif (labeltype == 'n'):
      a.sidelabel = "%d"%(a.num_in_file)
  pass



def set_circlesize_to_atoms(mol,cirtype,updwtype,hydeHtype):
  # cirtype:   : draw circle on the atom. 'A'll, 'H':heavyatom, 'M'arked_atom",'I':implicit H and C
  # circlesize = 'L'arge,'M'edium,'S'mall, '':not drawing circle
  # If cirtype='A'll, 'H':heavyatom, 'M'arked_atom".
  #  ==> 'M'edium, or ''  
  # If cirtype='M'arked_atom".
  #  ==> 'L'arge, or 'S'mall  
  # If cirtype='I'mplicit H and C:
  #  ==> 'M'edium, or 'S'mall  

 
  #print "#set_circlesize_to_atoms(cirtype:%s updwtype:%s hydeHtype:%s)\n"%(cirtype,updwtype,hydeHtype)
  for a in (mol.atoms):
    a.circlesize = ''
    if (cirtype == 'A'):
      a.circlesize = 'M'  
    if (cirtype == 'H') and (a.element != 'H'):
      a.circlesize = 'M'  
    if (cirtype == 'I') and (a.element != 'H'):
      a.circlesize = 'M'  
    if (cirtype=='M'):
      if (mol.marked_atoms[a.num]!=''):
        a.circlesize = 'L'  
      else: 
        a.circlesize = 'S'  
    if (hydeHtype =='T') and (a.element == 'H'):
      a.circlesize = ''
    #print "%d '%s' '%s'"%(a.num,a.circlesize,mol.marked_atoms[a.num])

  if (cirtype == 'I'):
    for b in (mol.bonds):
      sa = b.natom1
      ea = b.natom2
      anumH = -1
      if (mol.atoms[sa].element == 'H'):
       anumH    = sa
       anumNonH = ea
      if (mol.atoms[ea].element == 'H'):
       anumH    = ea
       anumNonH = sa
      if (anumH>=0): 
        if (updwtype=='T') and ((b.bondstereo=='1') or (b.bondstereo=='6')):
          mol.atoms[anumH].circlesize = 'S'
        if (mol.atoms[anumNonH].element != 'C'):
          mol.atoms[anumH].circlesize = 'S'
#  for a in (mol.atoms):
#    print "#%d %s circlesize '%s'"%(a.num_in_file,a.element,a.circlesize) 
# sys.exit(1)
  pass







def set_mol_marked_atoms_by_atom_num(mol,anumlist):
  i = 0 
  for num in (anumlist):
    if (int(num) <= mol.Natom):
      mol.marked_atoms[int(num)-1] = "%d"%(i+1)
      #print "num %d"%(int(num)) 
    i += 1




def set_mol_marked_atoms_by_atom_num_in_file(mol,anumlist):
  i = 0 
  for num_in_file in (anumlist):
    if (mol.num_from_num_in_file.has_key(int(num_in_file))):
      anum = mol.num_from_num_in_file[int(num_in_file)]
      mol.marked_atoms[anum] = "%d"%(i+1)
    else:
      print "#WARNING:not found num_in_file '%d'"%(int(num_in_file))
    i += 1


def return_larger(x,y):
  if (x>y):
    return(x)
  else:
    return(y)

def return_smaller(x,y):
  if (x<y):
    return(x)
  else:
    return(y)


def read_molecule_in_various_formats(mol='',filename='',filetype='',dirtype='',subdir='',contents='',read_policy='local',config_file=''):
  #print "#def read_molecule_in_various_formats(filename='%s',filetype='%s' dirtype='%s'):"%(filename,filetype,dirtype)
  filename_full = ''
  WEB_CONF = {}
  if (read_policy=='ligandbox'): 
    webfunc.read_environment_file(WEB_CONF,filename=config_file)
    filename_full = ligandbox_func.library_2D_full_filename_from_ligand_id(filename,WEB_CONF)
 
  elif (read_policy=='kcombuweb'): 
    webfunc.read_environment_file(WEB_CONF,filename=config_file)
    filename_full = kcombuweb_func.full_filename_from_short_filename(filename,dirtype,WEB_CONF)
 
  elif (read_policy=='local'): 
    filename_full = filename 
  else: 
    webfunc.read_environment_file(WEB_CONF,filename=config_file)
    filename_full = WEB_CONF.get(dirtype,'') + '/' + filename
    if (subdir=='1'):
      filename_full = WEB_CONF.get(dirtype,'') + '/' + filename[0] + '/'+ filename
 
  ok = 0
  if (filename_full.endswith(".sdf")) or (filetype == 'S') :
    ok = mol.read_in_sdf(filename=filename_full,contents=contents)
  elif (filename.endswith(".pdb")) or (filetype == 'P') :
    ok = mol.read_in_pdb(filename=filename_full,contents=contents)
    mol.guess_bonds_from_atom_xyz()
  elif (filename.endswith(".kcf")) or (filetype == 'K') :
    ok = mol.read_in_kcf(filename=filename_full,contents=contents)
  elif (filename.endswith(".mol2")) or (filetype == '2') :
    ok = mol.read_in_mol2(filename=filename_full,contents=contents)
  else:
    return(0)
  return(ok)



def cal_scale_by_min_max_XY_molecule(self,gwidth,gheight):
  xmin = xmax = self.atoms[0].x
  ymin = ymax = self.atoms[0].y
  for a in range(1,self.Natom):
    if (xmin>self.atoms[a].x):
      xmin = self.atoms[a].x
    if (xmax<self.atoms[a].x):
      xmax = self.atoms[a].x
    if (ymin>self.atoms[a].y):
       ymin = self.atoms[a].y
    if (ymax<self.atoms[a].y):
      ymax = self.atoms[a].y

  scx = gwidth /(xmax-xmin)
  scy = gheight/(ymax-ymin)
  if (scx<scy):
    return(scx)
  else:
    return(scy)


def find_empty_2D_point(mol,o):
## find empty space around the focusing atom 'o'
## only connected neighbor atoms are considered.
##  return((empX,empY))   where empX^2 + empY^2 = 1.


  neiX = []
  neiY = []
#  print "o:%d num_in_file %d"%(o,mol.atoms[o].num_in_file)
  for i in (mol.atoms[o].neighbors):
    nx = mol.atoms[i].x - mol.atoms[o].x
    ny = mol.atoms[i].y - mol.atoms[o].y
    L  = math.sqrt(nx*nx + ny*ny)
    nx = nx/L
    ny = ny/L
#    print "nx %f ny %f L %f"%(nx,ny,L)
    neiX.append(nx)
    neiY.append(ny)

  Ndiv = 18
  max_minDD = -1.0
  empX = 0.0
  empY = 0.0
  empTh = 0.0
  for i in range(Ndiv):
    theta = 2.0 * math.pi * i /float(Ndiv)
    pX = math.cos(theta)
    pY = math.sin(theta)

    minDD = -1.0
    for j in range(len(neiX)):
      DD = (neiX[j]-pX)*(neiX[j]-pX) + (neiY[j]-pY)*(neiY[j]-pY)
      if (minDD<0.0) or (DD<minDD):
        minDD = DD
#    print "#theta %f pX %f pY %f minDD %f"%(180.0*theta/math.pi,pX,pY,minDD)
    if (max_minDD<0.0) or (minDD>max_minDD):
      max_minDD = minDD
      empX = pX
      empY = pY
      empTh = theta
#  print "#num_in_file %d emp (%f %f) theta %f max_minDD %f"%(mol.atoms[o].num_in_file,empX,empY,180*empTh/math.pi,max_minDD)
  return((empX,empY))
  pass


def double_bond_gxy(sc,sx,sy,ex,ey,xo,yo,gxo,gyo,gwidth_dbond):
# s1  ------  e1
#  s *------* e
# s2  ------  e2

  nx =  (ey-sy)
  ny = -(ex-sx)
  nl = math.sqrt(nx*nx + ny*ny)
  if (nl>0.0):
    nx /= nl
    ny /= nl

  gsx1 = sc*(sx - xo) + gxo + gwidth_dbond*nx
  gsy1 = sc*(sy - yo) + gyo + gwidth_dbond*ny
  gex1 = sc*(ex - xo) + gxo + gwidth_dbond*nx
  gey1 = sc*(ey - yo) + gyo + gwidth_dbond*ny
  gsx2 = sc*(sx - xo) + gxo - gwidth_dbond*nx
  gsy2 = sc*(sy - yo) + gyo - gwidth_dbond*ny
  gex2 = sc*(ex - xo) + gxo - gwidth_dbond*nx
  gey2 = sc*(ey - yo) + gyo - gwidth_dbond*ny
  return([gsx1,gsy1,gex1,gey1,gsx2,gsy2,gex2,gey2])



def font_radius_modified_gxy(OPT,gsx,gsy,gex,gey,slabel,elabel,fsize):
  vx = gex - gsx;
  vy = gey - gsy;
  vlen = math.sqrt(vx*vx + vy*vy);
  if (vlen > 0.0):
    vx /= vlen
    vy /= vlen

  if (slabel!=''):
    gsmx = gsx + vx * 1.2*fsize/2
    gsmy = gsy + vy * 1.2*fsize/2
  else:
    gsmx = gsx
    gsmy = gsy

  if (elabel!=''):
    gemx = gex - vx * 1.2*fsize/2
    gemy = gey - vy * 1.2*fsize/2
  else:
    gemx = gex
    gemy = gey

  #print "gsx %d gsy %d gsmx %d gsmy %d"%(gsx,gsy,gsmx,gsmy)
  return([gsmx,gsmy,gemx,gemy])


def rgbhex_to_rgb255(rgbhex):
  if (len(rgbhex)<6):
    print "#ERROR(rgbhex_to_rgb255):insufficient length rgbhex '%s'."%(rgbhex)
    sys.exit(1)

  hex = {'0':0,'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'A':10,'B':11,'C':12,'D':13,'E':14,'F':15}

  R1 = rgbhex[0]
  R0 = rgbhex[1]
  G1 = rgbhex[2]
  G0 = rgbhex[3]
  B1 = rgbhex[4]
  B0 = rgbhex[5]
  R = 16*hex[R1] + hex[R0]
  G = 16*hex[G1] + hex[G0]
  B = 16*hex[B1] + hex[B0]
  return((R,G,B)) 









#############
### MAIN ####
#############

def _main():
  if (len(sys.argv)<2):
    print "molfunc.py"
    sys.exit(1)

  if (len(sys.argv)>1):
    (R,G,B) = rgbhex_to_rgb255(sys.argv[1])
    print "%s --> R %d G %d B %d"%(sys.argv[1],R,G,B)
if __name__ == '__main__': _main()
