#!/usr/bin/env python

##
## <moldraw_svg.cgi>
##

import sys
import os 
import math
import cgi


LastModDate = "2020/03/18"

OPT = {}

##
## BEFORE USING "moldraw_svg.cgi", FOLLOWING THREE VARIABLES ([1],[2],[3]) SHOULD BE MODIFIED.
##
## [1] Assign the directory containing "molecule.py" "molfunc.py" "pca.py" "rmsd.py" "kcombu_func.py" ##
##  by sys.path.append([directory]).

#sys.path.append('/home/web/kcombu/cgi-bin/moldraw_svg')
#sys.path.append('/home/takawaba/work/kcombu/src/moldraw')
#sys.path.append('/var/www/html/ligandbox/cgi-bin/moldraw')
#sys.path.append('/var/www/html/kcombu/cgi-bin/moldraw')

### [2] The variable 'READ_POLICY' is very important          ###
### 'local':read local files. 'ligandbox':ligandbox server file,'kcombuweb': kcombu server file ###
### If READ_POLICY != 'local', the program can only read molecular files in the allowed directory,
### assigned by CONFIG_FILE.
#READ_POLICY = 'web'
#CONFIG_FILE = '/home/web/kcombu/conf/webkcombu.conf'

READ_POLICY = 'local'
CONFIG_FILE = ''


import molecule 
import webfunc 
import molfunc
import pca 
import rmsd 
import kcombu_func
import ligandbox_func
import svg_mol


#############
#### MAIN ###
#############

## for directory for import  ##
if (sys.argv[0].find('/')!=-1):
  if (os.path.islink(sys.argv[0])):
    fname = os.readlink(sys.argv[0])
  else:
    fname = sys.argv[0]

  fields = fname.split('/')
  filedir   = fname.replace(fields[len(fields)-1],'')
  sys.path.append(filedir)

 
## read options ##

OPT['ifA']  = '' 
OPT['ifB']  = '' 
OPT['isdf']  = '' 
OPT['ipdb']  = '' 
OPT['sc']    = '20.0'
OPT['amgn'] = '1.0'
OPT['mk'] = ''
OPT['mknf'] = ''
OPT['of'] = 'out.svg'
OPT['font'] = 'sans-serif'
OPT['cir'] = 'M'
OPT['lbl']  = 'F'
OPT['slbl'] = 'M'
OPT['col'] = 'C'
OPT['tcol'] = 'C'
OPT['bcol'] = 'C'
OPT['iam'] = ''
OPT['iams'] = ''
OPT['AB']  = 'A'
OPT['rk']  = '1' 
OPT['fA'] = ''
OPT['fB'] = ''
OPT['dA'] = ''
OPT['sA'] = ''
OPT['dB'] = ''
OPT['sB'] = ''
OPT['ttl'] = 'T'
OPT['rxyz'] = '0:0:0'
OPT['gxyz'] = '0:0:0'
OPT['abs'] = 'F'
OPT['pca'] = 'F'
OPT['rms'] = 'F'
OPT['M'] = 'S';
OPT['btype'] = 'N';
OPT['lw'] = '0.1'
OPT['rcirl'] = '0.3'
OPT['rcirm'] = '0.2'
OPT['rcirs'] = '0.15'
OPT['fontsize'] = '10'
OPT['rsc'] = 'T'
OPT['hydeH'] = 'F'
OPT['updw'] = 'F'
OPT['ydir'] = 'U'
OPT['style'] = 'C'
OPT['bgrgb'] = 'FFFFFF'

cgiform = cgi.FieldStorage()
OPT['ttl'] = 'T'
OPT['rxyz'] = '0:0:0'
OPT['gxyz'] = '0:0:0'
OPT['abs'] = 'F'
OPT['pca'] = 'F'
OPT['rms'] = 'F'
OPT['M'] = 'S';
OPT['btype'] = 'N';
OPT['lw'] = '0.1'
OPT['rcirl'] = '0.3'
OPT['rcirm'] = '0.2'
OPT['rcirs'] = '0.15'
OPT['rsc'] = 'T'
OPT['hydeH'] = 'F'
OPT['updw'] = 'F'
OPT['ydir'] = 'U'
OPT['style'] = 'C'
OPT['bgrgb'] = 'FFFFFF'

OPT['rA'] = '1:0:0:0'
OPT['tA'] = '0:0:0'

OPT['webconf'] = ''

cgiform = cgi.FieldStorage()

if (len(sys.argv)<2) and (len(cgiform.keys())==0):
  print "Content-type: text/html\n\n<HTML><BODY><PRE>";
  print "moldraw_svg.cgi <options>"
  print " make image/pdf/eps file for molecules."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " ex.) for command line: moldraw_svg.cgi ifA=[sdffile] fA=S of=[output.svg]"
  print "      for cgi-bin     : moldraw_svg.cgi?ifA=[sdffile]&fA=S&of=webimg"
#  print "      for double mols : moldraw_svg.cgi -iam [kcombu-atom_match_file] -of [output.png]"
  print " -M     : MODE. 'S'ingle,'D'ouble [%s]"%(OPT['M']) 
  print " -h     : show detailed help messages"
  print "[caution]:  'molecule.py',  'molfunc.py', 'pca.py',  'rmsd.py', 'kcombu_func.py', 'svg_mol.py' are required."
  print "</PRE></BODY></HTML>";
  sys.exit(1)

if (len(sys.argv)==2) and ((sys.argv[1]=='-h') or (sys.argv[1]=='-help') or (sys.argv[1]=='--h') or (sys.argv[1]=='--help')):
  print "<options>"
  print "<options for single molecule>"
  print " -ifA   : input file for moleculeA. [%s]"%(OPT['ifA'])
  print " -fA    : file type for moleculeA. 'S'df,'P'db 'K'cf,'2':mol2[%s]"%(OPT['fA'])
  print " -isdf  : input SDF file [%s]"%(OPT['isdf']) 
  print " -ipdb  : input PDB file [%s]"%(OPT['ipdb']) 
  print " -iam   : input atom match file [%s]"%(OPT['iam'])
  print " -mk    : marked atom number         (1,2,3,...) [%s]"%(OPT['mk'])
  print " -mknf  : marked atom number_in_file (1,2,3,...) [%s]"%(OPT['mknf'])
  print " -AB    : molecule A or B for '-iam' [%s]"%(OPT['AB'])
  print " -dA    : directory type for molA ('key'-value described in 'webconf') [%s]"%(OPT['dA'])
  print " -sA    : sub directory type for molA [%s]"%(OPT.get('sA',''))
  print " -pca   : xy-flat orientation by PCA ('T' or 'F') [%s]"%(OPT['pca'])
  print " -rsc   : rescale by ave_bond length into 1.5 A  ('T' or 'F') [%s]"%(OPT['rsc'])
  print " -rA    : rotation for molecule (axis_x):(axis_y):(axis_z):(angle)[%s]"%(OPT['rA'])
  print " -gxyz  : center for rotation[%s]"%(OPT['gxyz'])
  print " -webconf: predefined web type. 'kcombu' or 'homcos'[%s]"%(OPT['webconf'])
  print "<options for double molecule>"
  print " -iam   : input atom match file in kcombu format[%s]"%(OPT['iam'])
  print " -iams  : input atom match file in simcomp format[%s]"%(OPT['iams'])
  print " -ifB   : input file for moleculeB. [%s]"%(OPT['ifB'])
  print " -fB    : file type for moleculeB. 'S'df,'P'db 'K'cf,'2':mol2[%s]"%(OPT['fB'])
  print " -dA    : input directory type for molB [%s]"%(OPT['dB'])
  print " -abs   : keep absolute 3D orientation between two molecules (T or F) [%s]"%(OPT['abs'])
  print " -rms   : rmsd-rotation for molA onto molB ('T' or 'F') [%s]"%(OPT['rms'])
  print "<options for both single and double molecule>"
  print " -rk    : rank for '-iam' [%s]"%(OPT['rk'])
  print "<options for display in global>"
  print " -style : Display style. 'M'onochrome(implicit HandC label, no circle, UpDown bonds)"
  print "        :                'C'olor(circle with CPK color,UpDown bonds) 'F:no style-assigend.[%s]"%(OPT['style'])
  print "  ** this option controls -cir,-lbl,-updw,-col options"
  print "<options for colors>"
  print " -col   : coloring: 'B':black_all 'C'pk_color, 'G'ray,'M'arked color(red),'R'ainbow_marked_color [%s]"%(OPT['col'])
  print " -bcol  : bond coloring: 'b':black_all 'A'tom_color [%s]"%(OPT['bcol'])
  print " -tcol  : text coloring: 'b':black_all 'A'tom_color [%s]"%(OPT['tcol'])
  print "<options for circle>"
  print " -cir   : draw circle on the atom. 'A'll, 'H':heavyatom, 'M'arked_atom"
  print "        :  'I':implicit H and heavyatoms,'F':not draw circle.[%s]"%(OPT['cir'])
  print " -rcirl : scale for radius of circle for 'large'  atoms(pixel/angstrom) [%s]"%(OPT['rcirl'])
  print " -rcirm : scale for radius of circle for 'medium' atoms(pixel/angstrom) [%s]"%(OPT['rcirm'])
  print " -rcirs : scale for radius of circle for 'small'  atoms(pixel/angstrom) [%s]"%(OPT['rcirs'])
  print "<options for labels>"
  print " -lbl       : label text in the center of atom. 'A'll_atom_element, 'N':non-C and non-H"
  print "            :  'I':implicit H and C 'M'atch_num,'n'um_in_file 'F'alse.[%s]"%(OPT['lbl'])
  print " -slbl      : side label text at the side of the atom. 'A'll_atom_element, 'N':non-C and non-H"
  print "            :  'I':implicit H and C 'M'atch_num,'n'um_in_file 'F'alse.[%s]"%(OPT['slbl'])
  print " -font      : fontname [%s]"%(OPT['font'])
  print " -fontsize  : fontsize for image (point) [%s]"%(OPT['fontsize'])
  print " -ttl       : show title (T or F)[%s]"%(OPT['ttl']) 
  print "<options for bonds>"
  print " -btype : bond type(single,double,...) display. 'N'ormal, 'I'gnore [%s]"%(OPT['btype'])
  print " -updw  : draw bond 'up'/'down' representation (T or F) [%s]"%(OPT['updw'])
  print " -lw    : line width [%s]"%(OPT['lw'])
  print "<other options>"
  print " -ydir  : Y-axis direction 'D'own, 'U'p [%s]"%(OPT['ydir'])
  print " -amgn  : angstrom-margin around molecule (angstrom) [%s]"%(OPT['amgn']) 
  print " -rxyz  : rotation angle(degree) for x,y,z axis[%s]"%(OPT['rxyz'])
  print " -gxyz  : center for rotation[%s]"%(OPT['gxyz'])
  print " -bgrgb : back ground color [%s]"%(OPT['bgrgb'])
  print " -hydeH : hyde Hydrogen (T or F) [%s]"%(OPT['hydeH']) 
  print " -of   : output image file name (stdout, webimg) [%s]"%(OPT['of'])
  print " -sc   : scale (pixel/angstrom) for Image [%s]"%(OPT['sc']) 
  print "</PRE></BODY></HTML>"
  sys.exit(1)

CONFIG_FILE = ''
READ_POLICY = 'cgi'

(Noption, READ_POLICY) = webfunc.get_cgi_values(OPT,sys.argv)

#print "#Noption %d READ_POLICY %s"%(Noption, READ_POLICY)
if (OPT.get('webconf','') == 'kcombu'):
  CONFIG_FILE = '/home/web/kcombu/conf/webkcombu.conf'
  sys.path.append('/var/www/html/kcombu/cgi-bin/moldraw_svg')
  READ_POLICY = 'cgi'
elif (OPT.get('webconf','') == 'ligandbox'):
  CONFIG_FILE = '/var/www/html/ligandbox/conf/webligandbox.conf'
  sys.path.append('/var/www/html/ligandbox/cgi-bin/moldraw_svg')
  READ_POLICY = 'cgi'
elif (OPT.get('webconf','') == 'homcos'):
  CONFIG_FILE = '/var/www/html/homcos/conf/webhomcos.conf'
  sys.path.append('/var/www/html/homcos/cgi-bin/moldraw_svg')
  READ_POLICY = 'cgi'


#for key in OPT.keys():
#  print "#key '%s' option '%s'"%(key,OPT[key])


OPT['sc']  = float(OPT['sc'])
OPT['amgn'] = float(OPT['amgn'])


if (OPT['isdf']!='') or (OPT['ipdb']!='') or ((OPT['ifA']!='') and (OPT['ifB']=='')):
  OPT['M'] = 'S'
elif (OPT['iam']!='') or (OPT['iams']!=''):
  OPT['M'] = 'D'

if (OPT['rxyz'] != '0:0:0'):
  thxyz = OPT['rxyz'].split(':')
  Rmat = molfunc.make_rotation_matrix_by_xyz_angle(float(thxyz[0]),float(thxyz[1]),float(thxyz[2]))
if (OPT['gxyz'] != '0:0:0'):
  gcen = OPT['gxyz'].split(':')
  gcen[0] = float(gcen[0])
  gcen[1] = float(gcen[1])
  gcen[2] = float(gcen[2])

#if (OPT['M']=='S'):
#  if (OPT['style']=='M'):
#    OPT['cir'] = 'F'
#    OPT['lbl']   = 'I'
#    OPT['updw'] = 'T'
#    OPT['col'] = 'B'
#  elif (OPT['style']=='C'):
#    OPT['cir']  = 'I'
#    OPT['lbl']  = 'F'
#    OPT['col']  = 'C'
#    OPT['updw'] = 'T'

SVGLINES = []



###################################
## ['S'] : Single Molecule Mode ### 
###################################
if (OPT['M']=='S'):
  
  ### [1] Read Molecule ##
  molA = molecule.Molecule()
 
  if (OPT['isdf'] != ''):
    OPT['ifA'] = OPT['isdf']
    OPT['fA'] = 'S'
  if (OPT['ipdb'] != ''):
    OPT['ifA'] = OPT['ipdb']
    OPT['fA'] = 'P'

  #print "#ifA '%s'"%(OPT['ifA']) 
  #print "#fA '%s'"%(OPT['fA']) 
  #print "#dA '%s'"%(OPT['dA']) 
  #print "#sA '%s'"%(OPT['sA']) 
  #print "#read_policy '%s'"%(READ_POLICY) 
  ok = molfunc.read_molecule_in_various_formats(mol=molA,filename=OPT['ifA'],filetype=OPT['fA'],dirtype=OPT['dA'],subdir=OPT['sA'],read_policy=READ_POLICY,config_file=CONFIG_FILE)
  #molA.write_in_sdf("molA.sdf")
  if (ok==0) or (molA.Natom == 0):
    print "Content-type: text/html\n\n<HTML><BODY><PRE>";
    print "#ERROR:Improper options: -isdf %s and -ipdb %s -ifA %s"%(OPT.get('isdf',''),OPT.get('ipdb',''), OPT.get('ifA')) 
    print "</PRE></BODY></HTML>"
    sys.exit(1)

  if (OPT['rA'] != '1:0:0:0') or (OPT['tA'] != '0:0:0'):
    molA.rotate_and_tranlate(OPT['rA'],OPT['tA'])

  if (OPT['ydir']=='U'):
    molA.rotate_180_degree_around_Xaxis() 


  if (OPT['rxyz'] != '0:0:0'):
    if (OPT['gxyz'] == '0:0:0'):
      gcen = molA.cal_gcenter()
    molA.rotate_gcenter(Rmat,gcen)
    #molA.write_in_sdf("rot.sdf")

  if (OPT['rsc']=='T'):
    molA.rescale_by_average_bond_length()

  if (OPT['pca']=='T'):
    if (OPT['ifB']!=''):
      pca.rotate_two_molecules_XY_flat_plane_by_PCA(molA,molB)
    else:
      pca.rotate_molecule_XY_flat_plane_by_PCA(molA)

  #print "#file '%s' m.Natom %d m.Nheavyatom %d"%(ifname,mol.Natom, mol.Nheavyatom)
 
 
  ### [2] Prepare marked_atoms[] array  ##
  molA.marked_atoms = ['' for i in range(molA.Natom)]
 
  if (OPT['mk'] != ''):
    mklist = OPT['mk'].split(',')
    molfunc.set_mol_marked_atoms_by_atom_num(molA,mklist)

  if (OPT['mknf'] != ''):
    mklist = OPT['mknf'].split(',')
    molfunc.set_mol_marked_atoms_by_atom_num_in_file(molA,mklist)

  if (OPT['iam'] != ''):
    amlist = [] 
    amdat = {}
    mklist = kcombu_func.read_atom_match_file(OPT['iam'],OPT['AB'],int(OPT['rk']),amdat)
    if (len(amlist)>0):
      for i in range(len(amlist[int(OPT['rk'])]['A'])):
        if (OPT['AB']=='A'):
          n = amlist[int(OPT['rk'])-1]['A'][i]-1 
        else:
          n = amlist[int(OPT['rk'])-1]['B'][i]-1
        molA.marked_atoms[int(n)] =  "%d"%(i+1)
     
  molfunc.set_label_to_atoms(molA,OPT['lbl'],OPT['updw'])
  molfunc.set_sidelabel_to_atoms(molA,OPT['slbl'])
  molfunc.set_circlesize_to_atoms(molA,OPT['cir'],OPT['updw'],OPT['hydeH'])
  
  #### [3I] Construct Image ###
  molfunc.set_rgb_to_atoms(molA,OPT['col'],255)
  [minX,maxX,minY,maxY] = molA.set_min_max_XY_molecule()
  
  XSIZE = int(math.ceil(OPT['sc']*(maxX-minX + 2*OPT['amgn'])))
  YSIZE = int(math.ceil(OPT['sc']*(maxY-minY + 2*OPT['amgn'])))
  GXO = int(math.ceil(OPT['sc']*OPT['amgn']))
  GYO = int(math.ceil(OPT['sc']*OPT['amgn']))
  #print "#minX %f maxX %f minY %f maxY %f XSIZE %d YSIZE %d"%(minX,maxX,minY,maxY,XSIZE,YSIZE)
  [R,G,B] =  molfunc.rgbhex_to_rgb255(OPT['bgrgb'])
  #im = Image.new("RGB",(XSIZE,YSIZE),(R,G,B))
  #imdraw = ImageDraw.Draw(im)

  SVGLINES.append('<svg width="%d" height="%d"'%(XSIZE,YSIZE))
  SVGLINES.append('xmlns="http://www.w3.org/2000/svg"')
  SVGLINES.append('xmlns:xlink="http://www.w3.org/1999/xlink"')
  SVGLINES.append('>')
  SVGLINES.append('<!-- COMMAND:%s -->'%(OPT.get('COMMAND','')))
  SVGLINES.append('<!-- DATE   :%s -->'%(OPT.get('START_DATE','')))
  
  #### [4] Draw bonds and atoms ####
  #image_mol.draw_bonds(OPT,molA,imdraw,       OPT['sc'],GXO,GYO)
  svg_mol.draw_bonds(OPT,molA,SVGLINES,OPT['sc'],GXO,GYO)
  svg_mol.draw_atoms(OPT,molA,SVGLINES,OPT['sc'],GXO,GYO)
  #image_mol.draw_atoms(OPT,molA,imdraw,myfont,OPT['sc'],GXO,GYO)



 
###################################
## ['D'] : Double Molecule Mode ### 
###################################
if (OPT['M']=='D'):

  ## [1] Read Molecules ##
  molA = molecule.Molecule()
  molB = molecule.Molecule()

  amdat = {}
  if (OPT['iam']!=''):
    amlist = []
    amdat = {}
    kcombu_func.read_atom_match_file(OPT['iam'],amlist,amdat)
    if (amdat.has_key('FiletypeA')):
      OPT['fA'] = amdat['FiletypeA']
    if (amdat.has_key('FiletypeB')):
      OPT['fB'] = amdat['FiletypeB']
    molfunc.read_molecule_in_various_formats(mol=molA,filename=amdat['MoleculeA'],filetype=OPT['fA'])
    molfunc.read_molecule_in_various_formats(mol=molB,filename=amdat['MoleculeB'],filetype=OPT['fB'])


  elif (OPT['iams']!=''):
    mklist = molfunc.read_atom_match_file_in_simcomp(OPT['iams'],OPT['AB'],int(OPT['rk']),amdat)
    amdat['molA'] = OPT['ifA']
    amdat['molB'] = OPT['ifB']
  else:
    molfunc.read_molecule_in_various_formats(mol=molA,filename=OPT['ifA'],filetype=OPT['fA'])
    molfunc.read_molecule_in_various_formats(mol=molB,filename=OPT['ifB'],filetype=OPT['fB'])

  if (OPT['rsc']=='T'):
    molA.rescale_by_average_bond_length()
    molB.rescale_by_average_bond_length()

  molA.rotate_180_degree_around_Xaxis() 
  molB.rotate_180_degree_around_Xaxis() 
 
  if (OPT['pca']=='T') and (OPT['rms']=='T'):
    pca.rotate_molecule_XY_flat_plane_by_PCA(molB)
    rmsd.rotate_molA_onto_molB_for_matched_atoms_with_minRMSD(molA, molB, amdat['anumA'], amdat['anumB']) 
  elif (OPT['rms']=='T'):
    #rmsd.rotate_molA_onto_molB_for_matched_atoms_with_minRMSD(molA, molB, amdat['anumA'], amdat['anumB']) 
    rmsd.rotate_molA_onto_molB_for_matched_atoms_with_minRMSD(molA, molB, amlist[0].numA, amlist[0].numB) 
  elif (OPT['pca']=='T'):
    pca.rotate_two_molecules_XY_flat_plane_by_PCA(molA,molB)

  if (OPT['rxyz'] != '0:0:0'):
    if (OPT['gxyz'] == '0:0:0'):
      gcen = molA.cal_gcenter()
    molA.rotate_gcenter(Rmat,gcen)
    molB.rotate_gcenter(Rmat,gcen)

 
  #print "#file '%s' m.Natom %d m.Nheavyatom %d"%(ifname,mol.Natom, mol.Nheavyatom)
  ### [2] make marked_atoms[] 
  molA.marked_atoms = ['' for i in range(molA.Natom)]
  molB.marked_atoms = ['' for i in range(molB.Natom)]

  if (len(amlist)>0):
    for i in range(len(amlist[int(OPT['rk'])-1]['A'])):
      a = amlist[int(OPT['rk'])-1]['A'][i] - 1
      b = amlist[int(OPT['rk'])-1]['B'][i] - 1
      if (int(a)<len(molA.marked_atoms)):
        molA.marked_atoms[int(a)] =  "%d"%(i+1)
      if (int(b)<len(molB.marked_atoms)):
        molB.marked_atoms[int(b)] =  "%d"%(i+1)
  
  molfunc.set_label_to_atoms(molA,OPT['lbl'],OPT['updw'])
  molfunc.set_sidelabel_to_atoms(molA,OPT['slbl'])
  molfunc.set_label_to_atoms(molB,OPT['lbl'],OPT['updw'])
  molfunc.set_sidelabel_to_atoms(molB,OPT['slbl'])

  molfunc.set_circlesize_to_atoms(molA,OPT['cir'],OPT['updw'],OPT['hydeH'])
  molfunc.set_circlesize_to_atoms(molB,OPT['cir'],OPT['updw'],OPT['hydeH'])
 
  #[minX_A,maxX_A,minY_A,maxY_A] = molA.min_max_XY_molecule()
  #[minX_B,maxX_B,minY_B,maxY_B] = molB.min_max_XY_molecule()
  molA.set_min_max_XY_molecule()
  molB.set_min_max_XY_molecule()

  minX = molfunc.return_smaller(molA.minX,molB.minX)
  maxX = molfunc.return_larger( molA.maxX,molB.maxX)
  minY = molfunc.return_smaller(molA.minY,molB.minY)
  maxY = molfunc.return_larger( molA.maxY,molB.maxY)

  wX_A = molA.maxX - molA.minX  
  wX_B = molB.maxX - molB.minX 
  wY_A = molA.maxY - molA.minY  
  wY_B = molB.maxY - molB.minY 

  max_wX = molfunc.return_larger(wX_A,wX_B)
  max_wY = molfunc.return_larger(wY_A,wY_B)

  #### [3I] Construct Image ####

  molfunc.set_rgb_to_atoms(molA,OPT['col'],255)
  molfunc.set_rgb_to_atoms(molB,OPT['col'],255)

  if (OPT['abs']=='T'):
    XSIZE = int(math.ceil(OPT['sc']*(max_wX*2 + 4*OPT['amgn'])))
    YSIZE = int(math.ceil(OPT['sc']*(max_wY   + 2*OPT['amgn'])))
  else:
    XSIZE = int(math.ceil(OPT['sc']*(wX_A + wX_B + 4*OPT['amgn'])))
    YSIZE = int(math.ceil(OPT['sc']*(max_wY + 2*OPT['amgn'])))

  #print "#XSIZE %d YSIZE %d"%(XSIZE, YSIZE)
  #print "#minX %f maxX %f minY %f maxY %f XSIZE %d YSIZE %d"%(minX,maxX,minY,maxY,XSIZE,YSIZE)
  [R,G,B] =  molfunc.rgbhex_to_rgb255(OPT['bgrgb'])
  im = Image.new("RGB",(XSIZE,YSIZE),(R,G,B))
  imdraw = ImageDraw.Draw(im)
 
  gminX_A = OPT['sc']*OPT['amgn']
  gminY_A = OPT['sc']*OPT['amgn']
  gminX_B = OPT['sc']*(3*OPT['amgn'] + wX_A)
  gminY_B = OPT['sc']*OPT['amgn']

  if (OPT['abs']=='T'):
    image_mol.draw_bonds(OPT,molA,imdraw,       OPT['sc'],gminX_A,gminX_A)
    image_mol.draw_atoms(OPT,molA,imdraw,myfont,OPT['sc'],gminX_A,gminX_A)
 
    image_mol.draw_bonds(OPT,molB,imdraw,       OPT['sc'],gminX_A,gminX_A)
    image_mol.draw_atoms(OPT,molB,imdraw,myfont,OPT['sc'],gminX_A,gminX_A)
  else:
    image_mol.draw_bonds(OPT,molA,imdraw,       OPT['sc'],gminX_A,gminY_A)
    image_mol.draw_atoms(OPT,molA,imdraw,myfont,OPT['sc'],gminX_A,gminY_A)
 
    image_mol.draw_bonds(OPT,molB,imdraw,       OPT['sc'],gminX_B,gminY_B)
    image_mol.draw_atoms(OPT,molB,imdraw,myfont,OPT['sc'],gminX_B,gminY_B)

  if (OPT['ttl']=='T'):
    titlestr = OPT['iam']
    #titlestr = amdat['molA'] + '_vs_' + amdat['molB']
    titlestr = OPT['iam'] + OPT['iams'] + ':' + amdat.get('molA','') + '_vs_' + amdat.get('molB','')
    #imdraw.text((OPT['ifontsize']/2,OPT['ifontsize']/4),titlestr,font=myfont,fill=(0,0,0))



########################### 
####  Output the image ####
########################### 

SVGLINES.append('</svg>\n')

if (OPT['of'] == 'stdout'):
  svg_mol.save_lines('-',SVGLINES)
elif (OPT['of'] == 'webimg'):
  print "Content-type: image/svg+xml\n";
  svg_mol.save_lines('-',SVGLINES)
elif (OPT['of']!= '') and ((OPT['of'].endswith('.svg')) or (OPT['of'].endswith('.SVG'))):
  svg_mol.save_lines(OPT['of'],SVGLINES)
else:
  print "Content-type: text/html\n\n<HTML><BODY><PRE>";
  print "#ERROR:Improper '-of' option" 
  print "</PRE></BODY></HTML>"
  sys.exit(1)
