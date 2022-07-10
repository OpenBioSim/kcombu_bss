#
# <molecule.py>
#
# Class for dealing various format of molecule.
# mainly for 'moldraw.cgi'
#


import sys
import os
import re 
import math

LastModDate = "2019/12/01"

class Molecule:

  def __init__(self):
    self.filename    = ''
    self.Natom       = 0 ## number of atoms (int)
    self.Nheavyatom  = 0 ## number of heavy atoms (int)
    self.Ncarbon     = 0 ## number of carbon atoms (int)
    self.formula     = 0 ## molecular formula
    self.Nbond    = 0    ## number of bonds (int)
    self.atoms     = []   ## list of "class Atom"
    self.bonds     = []   ## list of "class Bond"
    self.contable = []   ## conection table matrix (2D matrix of "str")
    self.Nele     = {} 
    self.num_from_num_in_file = {}
    self.num_from_num_ic      = {}
    self.minX = 0.0 
    self.maxX = 0.0 
    self.minY = 0.0 
    self.maxY = 0.0 
    self.annotation = {}
    pass

  def read_in_sdf(self,filename='',contents=''):
#>>> sample of SDF file (benzene.sdf) <<
#241
#  -OEChem-08070801562D
#
# 12 12  0     0  0  0  0  0  0999 V2000
#    2.8660    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    3.7321    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.0000   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    3.7321   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.8660   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.8660    1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    1.4631    0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    4.2690    0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    1.4631   -0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    4.2690   -0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    2.8660   -1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#  1  2  2  0  0  0  0
#  1  3  1  0  0  0  0
#  1  7  1  0  0  0  0
#  2  4  1  0  0  0  0
#  2  8  1  0  0  0  0
#  3  5  2  0  0  0  0
#  3  9  1  0  0  0  0
#  4  6  2  0  0  0  0
#  4 10  1  0  0  0  0
#  5  6  1  0  0  0  0
#  5 11  1  0  0  0  0
#  6 12  1  0  0  0  0
#M  END
#
#>>> sample of SDF file (KEGG_DRUG D00021) <<
#
#
#
# 12 12  0  0  1  0  0  0  0  0999 V2000
#   23.1700   -8.6100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   23.1700  -10.0100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   24.3824  -10.7100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   25.5949  -10.0100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   25.5949   -8.6100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   24.3824   -7.9100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   21.9576   -7.9100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   20.7621   -8.6004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   19.5747   -7.9149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   18.3835   -8.6029    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#   19.5745   -6.5102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#   20.7620  -10.0097    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
#  1  2  2  0     0  0
#  2  3  1  0     0  0
#  3  4  2  0     0  0
#  4  5  1  0     0  0
#  5  6  2  0     0  0
#  1  6  1  0     0  0
#  1  7  1  0     0  0
#  7  8  1  0     0  0
#  8  9  1  0     0  0
#  9 10  1  0     0  0
#  9 11  2  0     0  0
#  8 12  1  6     0  0
#M  END



#>>> sample of SDF file (pubchem CID_5950;alanine) <<
#
#5950
#  -OEChem-10111103342D
#
# 13 12  0     1  0  0  0  0  0999 V2000
#    5.1350   -0.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#    4.2690    1.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#    2.5369    0.2500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
#    3.4030   -0.2500    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
#    3.4030   -1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    4.2690    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    3.4030    0.3700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    2.7830   -1.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    3.4030   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    4.0230   -1.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    2.0000   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    2.5369    0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    5.6720    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#  1  6  1  0  0  0  0
#  1 13  1  0  0  0  0
#  2  6  2  0  0  0  0
#  4  3  1  6  0  0  0
#  3 11  1  0  0  0  0
#  3 12  1  0  0  0  0
#  4  5  1  0  0  0  0
#  4  6  1  0  0  0  0
#  4  7  1  0  0  0  0
#  5  8  1  0  0  0  0
#  5  9  1  0  0  0  0
#  5 10  1  0  0  0  0
#M  END
#> <PUBCHEM_COMPOUND_CID>
#5950
#
#> <PUBCHEM_COMPOUND_CANONICALIZED>
#1
#
#> <PUBCHEM_CACTVS_COMPLEXITY>
#61.8
#
#> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>
#3
#
#> <PUBCHEM_CACTVS_HBOND_DONOR>
#2
#
#> <PUBCHEM_CACTVS_ROTATABLE_BOND>
#1
#
#> <PUBCHEM_CACTVS_SUBSKEYS>
#AAADcYBCMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHgAQCAAACCjBgAQCCABAAgAIAACQCAAAAAAAAAAAAIGAAAACAAAAAAAAQAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==
#
#> <PUBCHEM_IUPAC_OPENEYE_NAME>
#(2S)-2-aminopropanoic acid
#> <PUBCHEM_COORDINATE_TYPE>
#1
#5
#255

#> <PUBCHEM_BONDANNOTATIONS>
#4  3  6




#
#>> mol file from JME editor <<
##c1ccccc1
#JME 2010.01.2 Fri Aug 05 10:37:36 JST 2011
# 
#  6  6  0  0  0  0  0  0  0  0999 V2000
#    0.0000    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    0.0000    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    1.2124    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    1.2124    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.4249    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#  1  2  2  0  0  0  0
#  1  3  1  0  0  0  0
#  2  4  1  0  0  0  0
#  3  5  2  0  0  0  0
#  4  6  2  0  0  0  0
#  5  6  1  0  0  0  0
#M  END
#
##>> SDF file from namiki ###
#
#  Symyx   10231017462D 1   1.00000     0.00000     0
#
# 19 19  0     0  0            999 V2000
#   -1.0042    0.0000    0.0000 C   0  0  0  0  0  0           0  0  0
#   -0.5042    0.8708    0.0000 C   0  0  0  0  0  0           0  0  0
#   -0.5042   -0.8667    0.0000 C   0  0  0  0  0  0           0  0  0
#   -2.0083    0.0000    0.0000 O   0  0  0  0  0  0           0  0  0
#:
#M  END
#> <NAMIKI_ID> (NS-00004048)
#NS-00004048
#
#> <SUPPLIERNAME_1> (Zelinsky)
#Zelinsky
#
#> <SUPPLIERID_1> (UZI/9254967)
#UZI/9254967
#
#> <SUPPLIERNAME_2> (Vitas-M)
#Vitas-M
#
#> <SUPPLIERID_2> (STK398194)
#STK398194
#
#> <SUPPLIERNAME_3> (Princeton)
#Princeton
#
#> <SUPPLIERID_3> (OSSL_150646)
#OSSL_150646
#
#
#> <LIGANDBOX_ID>
#00003242-01


    if (filename!='') and (os.access(filename,os.R_OK)):
      f = open(filename)
      contents = f.read()
      f.close()

    if (contents==''):
      print "#ERROR:Can't open sdf file '%s'"%(filename)
      return(0)


    lines = contents.split('\n')
    self.Nline = len(lines)
    self.filename = filename
    self.Nele = {} 
    Nline = 0
    status = ' ' 
    natom = 0
    self.Nheavyatom = 0
    nbond = 0
    key_annot   = ''
    value_annot = ''

    for line in (lines): 
      line = line.rstrip('\n')
      line = line.rstrip('\r')
      #print "<PRE>%s\n</PRE>"%(line)
      #if (line.startswith("M  END")):
      if ((status==' ') and (line.startswith("M  "))):
        status =  'E' 
      ### read head Natom,Nbond lines ###
      ###012345678901234567890123456789012345678901234567890 
      ###           1         2         3         4          
      ### 13 12  0  0  0  0  0  0  0  0  0 
      ### 10  9  0     0  0  0  0  0  0999 V2000 
      ###107110  0     0  0  0  0  0  0999 V2000 

      if (status==' ') and (Nline==3) and (len(line)>=6):
        if (re.match(r'^[0-9\s]+$',line[0:3])):
          self.Natom = int(line[0:3])
        if (re.match(r'^[0-9\s]+$',line[3:6])):
          self.Nbond = int(line[3:6])
        for i in range(self.Natom):
          a = Atom()
          self.atoms.append(a)
        for i in range(self.Nbond):
          b = Bond()
          self.bonds.append(b)
        self.contable = [['0' for i in range(self.Natom)] for j in range(self.Natom)]
      ### read atom lines ###
      if (status== ' ') and (Nline>3) and (line.startswith('#')==0) and (len(line)>30) and (self.Natom>0):
#    2.8660   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    2.8660    1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#0123456789012345678901234567890123456789012345678901234567890
#          1         2         3         4
        self.atoms[natom].num = natom 
        self.atoms[natom].num_in_file = natom + 1
        self.num_from_num_in_file[self.atoms[natom].num_in_file] = self.atoms[natom].num 
        self.atoms[natom].x = float(line[0:10]) 
        self.atoms[natom].y = float(line[10:20]) 
        self.atoms[natom].z = float(line[20:30])
        ele  = line[31:33]
        if (ele[-1]==' '):
          ele = ele[0:1]
        self.atoms[natom].element  = ele
        self.atoms[natom].atomname  = ele
        self.Nele[ele] = self.Nele.get(ele,0) + 1
        self.atoms[natom].atomtype = self.atoms[natom].element
        #print x,y,z,atomtype
        if (self.atoms[natom].element != 'H'):
          self.Nheavyatom  += 1 
        natom += 1
      ### read bond lines ###
      if (status == ' ') and (Nline>3) and (line.startswith('#')==0) and (len(line)<30) and (self.Nbond>0):
#  1  2  2  0  0  0  0
#  1  3  1  0  0  0  0
#  9 10  1  0     0  0
#  9 11  2  0     0  0
#  8 12  1  6     0  0
#0123456789012345678901234567890123456789012345678901234567890
#          1         2         3         4
        natom1   =    int(line[0:3]) - 1 
        natom2      = int(line[3:6]) - 1 
        bondtype    = line[8:9] 
        bondstereo = line[11:12] 
        self.bonds[nbond].num      = nbond 
        self.bonds[nbond].natom1   = natom1
        self.bonds[nbond].natom2   = natom2
        self.bonds[nbond].bondtype = bondtype
        self.bonds[nbond].bondstereo = bondstereo
        self.contable[natom1][natom2] = bondtype
        self.contable[natom2][natom1] = bondtype
        nbond += 1
        pass   

      ### read annotation lines ###
      if (status=='E') and (line.startswith('>')==1):
        #print line
        field = line.split(' ')
        if (len(field)>1):
          key_annot  = field[1]
          key_annot  = key_annot.lstrip('<')
          key_annot  = key_annot.rstrip('>')
      if (status=='E') and (len(key_annot)>1) and (line.startswith('>')==0) and (line.startswith('#')==0) and (len(line)>1):
        #value_annot = line
        self.annotation[key_annot] = line
        pass

      Nline += 1
    self.formula =  make_molecular_formula(self.Nele)
    self.count_neighbor_atoms()
    if (self.Nele.has_key('C')):
      self.Ncarbon = self.Nele['C']
    else:
      self.Ncarbon = 0
    return(1)


  def read_in_pdb(self,filename='',contents=''):

    if (filename!='') and (os.access(filename,os.R_OK)):
      f = open(filename)
      contents = f.read()
      f.close()

    if (contents==''):
      print "#ERROR:Can't open pdb file '%s'"%(filename)
      return(0)
    lines = contents.split('\n')

    self.filename = filename
    self.Nele = {} 
    Nline = 0
    natom = 0
    self.Nheavyatom = 0
    nbond = 0
    for line in (lines): 
      line = line.rstrip('\n')
      if (line.startswith("ATOM")) or (line.startswith("HETATM")):
    #          1         2         3         4         5         6         7
    #01234567890123456789012345678901234567890123456789012345678901234567890123456789
    #ATOM    676  CB  GLU A  85      10.440  29.552  12.788  6.00 16.96           C
    #ATOM    680  OE2 GLU A  85      10.230  30.451  16.374  8.00 41.03           O
    #ATOM    682  CA  LEU A  86       7.618  29.487   9.238  6.00 12.23           C
    #HETATM 1236  O4  SO4   154      33.810  28.815  -4.624  8.00 14.90           O
    #HETATM 1237 FE   HEM   155      15.271  27.962   0.622 24.00  7.86          FE
    #:   
    #HETATM 2036 ZN    ZN A 261      37.506  16.588 -14.749  0.97  7.81          ZN
    #HETATM 2039 HG    HG A 266      31.074  10.149  -2.419  0.34 13.26          HG
    #HETATM 2040 HG    HG A 267      30.947  10.315  -0.870  0.39 14.45          HG
    #HETATM 2041 HG    HG A 269      36.267  14.504  -1.501  0.10 13.42          HG
    #:   
    #HETATM 4028  H91 EQP A 500      30.115  19.716  64.216  1.00  0.00
    #HETATM 4029  H92 EQP A 500      29.902  21.176  63.270  1.00  0.00
    #HETATM 4030 H111 EQP A 500      28.343  14.492  59.955  1.00  0.00
    #HETATM 4031 H112 EQP A 500      26.751  13.647  60.120  1.00  0.00
    #HETATM 4032 H113 EQP A 500      28.237  12.939  60.857  1.00  0.00
    #HETATM 4033  HN5 EQP A 500      26.029  15.764  61.013  1.00  0.00
    #HETATM 4034 HOP3 EQP A 500      22.616  20.686  66.538  1.00  0.00
    #ATOM      6  H   ALA     2      -0.131  -2.162  -0.491  1.00  0.00           H
    #ATOM      7  HA  ALA     2      -0.269   0.603  -1.418  1.00  0.00           H
    #ATOM      8 1HB  ALA     2      -1.605   1.006   0.691  1.00  0.00           H
    #ATOM      9 2HB  ALA     2      -0.285   0.342   1.681  1.00  0.00           H
    #ATOM     10 3HB  ALA     2      -0.053   1.861   0.784  1.00  0.00           H
    #ATOM     21 1HG1 VAL     3       5.994  -1.541  -1.334  1.00  0.00           H
    #ATOM     22 2HG1 VAL     3       5.904  -1.816  -3.093  1.00  0.00           H
    #ATOM     23 3HG1 VAL     3       6.311  -0.221  -2.485  1.00  0.00           H
    #ATOM     24 1HG2 VAL     3       2.300  -1.925  -2.264  1.00  0.00           H
    #ATOM     25 2HG2 VAL     3       3.658  -2.724  -3.064  1.00  0.00           H
    #ATOM     26 3HG2 VAL     3       3.611  -2.626  -1.284  1.00  0.00           H

        anum  = line[6:11]
        atom  = line[12:16]
        resi  = line[17:20]
        chain = line[21:22]
        rnum  = line[22:27]
        a = Atom()
        self.atoms.append(a)
        self.atoms[natom].num = natom 
        self.atoms[natom].num_in_file = int(anum.replace(' ','')) 
        self.num_from_num_in_file[self.atoms[natom].num_in_file] = self.atoms[natom].num 
        self.atoms[natom].x = float(line[31:38]) 
        self.atoms[natom].y = float(line[38:46]) 
        self.atoms[natom].z = float(line[46:54])
        ele = line[12:14].replace(' ','')
        self.atoms[natom].atomname = line[12:16].replace(' ','') 
        if (ele[0]=="H") and (line[12:15] != "HG"):
          ele = "H" 
        if ((line[12]=="1")or(line[12]=="2")or(line[12]=="3")) and (line[13]=='H'):
          ele = "H" 
        self.atoms[natom].element  = ele
        #print "atomname '%s' ele '%s'"%(self.atoms[natom].atomname,self.atoms[natom].element)
        self.Nele[ele] = self.Nele.get(ele,0) + 1
        self.atoms[natom].atomtype = self.atoms[natom].element
        if (self.atoms[natom].element != 'H'):
          self.Nheavyatom  += 1 
        natom += 1 
    self.Natom = natom
    self.formula =  make_molecular_formula(self.Nele)
    return(1)

  def write_in_pdb(self,ofname):
    if (ofname=='stdout'):
      of = sys.stdout
    else:
      print "#write_in_pdb()-->'%s'"%(ofname)
      of = open(ofname,'w')
    for i in range(self.Natom):
      prenum1 = prenum2 = prenum3 = '0'
      if (self.atoms[i].prenum1>=0):
        prenum1 = self.atoms[self.atoms[i].prenum1].num_in_file
      if (self.atoms[i].prenum2>=0):
        prenum2 = self.atoms[self.atoms[i].prenum2].num_in_file
      if (self.atoms[i].prenum3>=0):
        prenum3 = self.atoms[self.atoms[i].prenum3].num_in_file
      of.write("ATOM  %5d %3s  %3s %5s    %8.3f%8.3f%8.3f%6.3f%6.1f%6.1f%5s%5s%5s\n"%(self.atoms[i].num_in_file,self.atoms[i].element,"MOL","1",self.atoms[i].x,self.atoms[i].y,self.atoms[i].z,self.atoms[i].bond_len,self.atoms[i].bond_ang,self.atoms[i].torsion_ang,prenum1,prenum2,prenum3))
    if (ofname!='stdout'):
      of.close()

  def read_in_kcf(self,filename='',contents=''):
    if (filename!='') and (os.access(filename,os.R_OK)):
      f = open(filename)
      contents = f.read()
      f.close()

    if (contents==''):
      print "#ERROR:Can't open kcf file '%s'"%(filename)
      return(0)

    lines = contents.split('\n')

    self.filename = filename
    self.Nele = {} 
    Nline = 0
    status = ' ' 
    natom = 0
    self.Nheavyatom = 0
    nbond = 0
    for line in (lines): 
      line = line.rstrip('\n')
      field = line.split()
      #if (line.startswith("M  END")):
      if (line.startswith("///")):
        status = 'E'
      ## read head Natom,Nbond lines ###
      ##012345678901234567890123456789012345678901234567890 **/
      ##           1         2         3         4          **/
      ##ENTRY       iressa                      Compound
      ##ATOM        31
      ##            1   X   Cl   12.3739   -3.5272
      ##            2   X   F    14.1174   -2.5473
      ##            3   O2x O     2.0000   -0.4728
      ##            4   O2a O     7.1962    0.5272
      ##            29  C8y C    12.3854   -2.5273
      ##            :
      ##            30  C8x C    13.2687   -1.0374
      ##            31  C8y C    13.2572   -2.0373
      ##BOND        34
      ##            1     1  29 1
      ##            2     2  31 1
      ##            3     3  14 1
      ##            4     3  15 1

      ### read ATOM lines ###
      elif (status == 'A') and  (line.startswith("      ")):
        self.atoms[natom].num = natom 
        self.atoms[natom].num_in_file = natom + 1
        self.num_from_num_in_file[self.atoms[natom].num_in_file] = self.atoms[natom].num 
        self.atoms[natom].x = float(field[3]) 
        self.atoms[natom].y = float(field[4]) 
        self.atoms[natom].z = 0.0
        self.atoms[natom].element  = field[2] 
        self.atoms[natom].atomname = field[1]
        self.atoms[natom].atomtype = field[1]
        if (self.atoms[natom].element != 'H'):
          self.Nheavyatom  += 1 
        natom += 1

      elif (line.startswith("ATOM")):
        self.Natom = int(field[1]) 
        for i in range(self.Natom):
          a = Atom()
          self.atoms.append(a)
        self.contable = [['0' for i in range(self.Natom)] for j in range(self.Natom)]
        status = 'A'
      elif (line.startswith("BOND")):
        self.Nbond = int(field[1]) 
        for i in range(self.Nbond):
          b = Bond()
          self.bonds.append(b)
        status = 'B'
      ### read BOND lines ###
      elif (status == 'B') and  (line.startswith("      ")):
        natom1   = int(field[1]) - 1 
        natom2   = int(field[2]) - 1 
        bondtype = field[3]
        self.bonds[nbond].natom1   = natom1
        self.bonds[nbond].natom2   = natom2
        self.bonds[nbond].bondtype = bondtype
        self.contable[natom1][natom2] = bondtype
        self.contable[natom2][natom1] = bondtype
        nbond += 1
        pass   
      Nline += 1
    self.formula =  make_molecular_formula(self.Nele)
    if (self.Nele.has_key('C')):
      self.Ncarbon = self.Nele['C']
    else:
      self.Ncarbon = 0
    return(1)

  def read_in_mol2(self,filename='',contents=''):
  #>>EXAMPLE OF MOL2 FILE FORMAT <<
  #@<TRIPOS>MOLECULE
  #5950
  # 13 12 0 0 0
  #SMALL
  #GASTEIGER
  #Energy = 0
  #
  #@<TRIPOS>ATOM
  #      1 O           5.1350   -0.2500    0.0000 O.3     1  ALA1       -0.4795
  #      2 OXT         4.2690    1.2500    0.0000 O.2     1  ALA1       -0.2493
  #      3 N           2.5369    0.2500    0.0000 N.3     1  ALA1       -0.3186
  #      4 CA          3.4030   -0.2500    0.0000 C.3     1  ALA1        0.1004
  #      5 CB          3.4030   -1.2500    0.0000 C.3     1  ALA1       -0.0395
  #      6 C           4.2690    0.2500    0.0000 C.2     1  ALA1        0.3214
  #      7 HA          3.4030    0.3700    0.0000 H       1  ALA1        0.0571
  #      8 HB1         2.7830   -1.2500    0.0000 H       1  ALA1        0.0250
  #      9 HB2         3.4030   -1.8700    0.0000 H       1  ALA1        0.0250
  #     10 HB3         4.0230   -1.2500    0.0000 H       1  ALA1        0.0250
  #     11 H1          2.0000   -0.0600    0.0000 H       1  ALA1        0.1189
  #     12 H2          2.5369    0.8700    0.0000 H       1  ALA1        0.1189
  #     13 H           5.6720    0.0600    0.0000 H       1  ALA1        0.2951
  #@<TRIPOS>BOND
  #     1     1     6    1
  #     2     1    13    1
  #     3     2     6    2
  #     4     4     3    1
  #     5     3    11    1
  #     6     3    12    1
  #     7     4     5    1
  #     8     4     6    1
  #     9     4     7    1
  #    10     5     8    1
  #    11     5     9    1
  #    12     5    10    1
    if (filename!='') and (os.access(filename,os.R_OK)):
      f = open(filename)
      contents = f.read()
      f.close()

    if (contents==''):
      print "#ERROR:Can't open mol2 file '%s'"%(filename)
      return(0)

    lines = contents.split('\n')

    self.filename = filename
    self.Nele = {} 
    Nline = 0
    status = ' ' 
    natom = 0
    self.Nheavyatom = 0
    nbond = 0
    for line in (lines): 
      line = line.rstrip('\n')
      field = line.split()
      if (status=='M') and (Nline==2):
        self.Natom = int(field[0]) 
        for i in range(self.Natom):
          a = Atom()
          self.atoms.append(a)
        self.contable = [['0' for i in range(self.Natom)] for j in range(self.Natom)]
        self.Nbond = int(field[1]) 
        #print "#Natom %d Nbond %d"%(self.Natom, self.Nbond) 
        for i in range(self.Nbond):
          b = Bond()
          b.num = len(self.bonds) + 1
          self.bonds.append(b)
 
      if (line.startswith("@<TRIPOS>MOLECULE")):
        status = 'M'
        Nline = 0
      elif (line.startswith("@<TRIPOS>ATOM")):
        status = 'A'
      elif (line.startswith("@<TRIPOS>BOND")):
        status = 'B'
      #if (line.startswith("M  END")):

      ### read ATOM lines ###
      elif (status == 'A') and (line.startswith("@")==0) and (len(field)>=5):
  #      4 CA          3.4030   -0.2500    0.0000 C.3     1  ALA1        0.1004
        self.atoms[natom].num = natom 
        self.atoms[natom].num_in_file = int(field[0])
        self.num_from_num_in_file[self.atoms[natom].num_in_file] = self.atoms[natom].num 
        self.atoms[natom].atomtype = field[1]
        self.atoms[natom].x = float(field[2]) 
        self.atoms[natom].y = float(field[3]) 
        self.atoms[natom].z = float(field[4])
        self.atoms[natom].atomname = field[5]
        self.atoms[natom].ligname  = field[7]
        self.atoms[natom].charge = float(field[8])
        elefield = []
        elefield = field[5].split('.') 
        self.atoms[natom].element  = elefield[0] 
        #print "'%s' --> element '%s'"%(field[5],elefield[0])
        if (self.atoms[natom].element != 'H'):
          self.Nheavyatom  += 1 
        natom += 1

      ### read BOND lines ###
      elif (status == 'B') and  (line.startswith("@")==0) and (len(field)>=4):
  #     1     1     6    1
  #     2     1    13    1
        natom1   = int(field[1]) - 1 
        natom2   = int(field[2]) - 1 
        bondtype = field[3]
        if (bondtype=='ar'):
          bondtype = '1'

        self.bonds[nbond].natom1   = natom1
        self.bonds[nbond].natom2   = natom2
        self.bonds[nbond].bondtype = bondtype
        self.contable[natom1][natom2] = bondtype
        self.contable[natom2][natom1] = bondtype
        nbond += 1
        pass   
      Nline += 1
    self.formula =  make_molecular_formula(self.Nele)
    if (self.Nele.has_key('C')):
      self.Ncarbon = self.Nele['C']
    else:
      self.Ncarbon = 0
    return(1)





  def write_in_sdf(self,ofname):
    if (ofname=='stdout'):
      of = sys.stdout
    else:
      print "#write_in_sdf()-->'%s'"%(ofname)
      of = open(ofname,"w")
    #of.write("0\n")
    of.write("%s\n"%(self.filename))
    #of.write("  -OEChem-09210803422D\n\n")
    #of.write("  -ISIS-            3D\n\n")
    of.write("  input_file:'%s'\n\n"%(self.filename))
    #of.write(" 13 12  0     1  0  0  0  0  0999 V2000\n")
    of.write("%3d%3d  0     1  0  0  0  0  0999 V2000\n"%(len(self.atoms),len(self.bonds)))
    # (1) write atoms #
    for a in (self.atoms):
      of.write("%10.4f%10.4f%10.4f %-2s  0  0  0  0  0  0  0  0  0  0  0  0\n"%(a.x,a.y,a.z,a.element))
    # (2) write bonds #
    for b in (self.bonds):
      of.write("%3d%3d%3s  0  0  0  0\n"%(b.natom1+1,b.natom2+1,b.bondtype))
    of.write("M  END\n")
    if (ofname!='stdout'):
      of.close()


  def write_in_mol2(self,ofname):
    if (ofname=='stdout'):
      of = sys.stdout
    else:
      print "#write_in_mol2()-->'%s'"%(ofname)
      of = open(ofname,"w")
    of.write("@<TRIPOS>MOLECULE\n")   
    of.write("%s\n"%(self.filename))
    of.write("%3d%3d 0 0 0\n"%(self.Natom, self.Nbond))
    of.write("\n")
    of.write("\n")
    of.write("@<TRIPOS>ATOM\n")   
    for a in (self.atoms):
#      1 C          -6.9684    0.8200    0.0000 C.3     1  LIG1        0.1046
#      2 C          -8.3999    1.2681    0.0000 C.3     1  LIG1        0.1595
#      3 O          -8.7276    2.7319    0.0000 O.3     1  LIG1       -0.3508
      a.ligname = 'LIG1'
      of.write("%7d %-2s      %10.4f%10.4f%10.4f %3s %5d  %4s    %10.4f\n"%(a.num+1,a.element,a.x,a.y,a.z,a.atomname,1,a.ligname,a.charge)) 
    # (2) write bonds #
    of.write("@<TRIPOS>BOND\n")   
#     9     8     9    1
#    10     9    10    1
#    11    10    11    1
#    12    11    12   ar
    for b in (self.bonds):
      of.write("%6d%6d%6d   %-2s\n"%(b.num+1,b.natom1+1,b.natom2+1,b.bondtype))
    if (ofname!='stdout'):
      of.close()

  def guess_bonds_from_atom_xyz(self):
# <Typical Bond Length>
#
# aromatic C-C : 1.40 A
# sp3      C-C : 1.50--1.55 A
#          S-C : 1.94 A
# aromtic  N-C : 1.34--1.36 A
#        C-O-C : 1.45 A
#          C=O : 1.24 A
#        C-S-C : 1.80 A
#          S-C : 1.80 A
#        O-P-O : 1.77 A
# others:
#
#        I-C   : 2.31 A
#        F-C   : 1.19 A
#        F-C   : 0.97 A
# hydrogen:
#        H-C   : 1.11 A
#        H-O   : 0.96 A
#        H-N   : 1.08 A
#
# Therefore, we set normally, Dbond_min = 1.1 A, Dbond_max = 1.9 A.
# If one of the atoms are not 'C','N','O','S','P',
#   Dbond_min = 0.9 A, Dbond_max = 2.4 A.


    #print "#guess_bonds_from_atom_xyz(self):"
    self.contable = [['0' for i in range(self.Natom)] for j in range(self.Natom)]
    dmap = [[0.0 for i in range(self.Natom)] for j in range(self.Natom)]
    Nbond_max = {}
    Nbond_max['C'] = 4;
    Nbond_max['O'] = 2;
    Nbond_max['N'] = 4;
    Nbond_max['S'] = 4;
    Nbond_max['P'] = 5;
    Nbond_max['F'] = 1;
    Nbond_max['Cl'] = 1;
    Nbond_max['CL'] = 1;
    Nbond_max['Br'] = 1;
    Nbond_max['BR'] = 1;
    Nbond_max['I'] = 1;
    standard = {}
    standard['C'] =  standard['O'] =  standard['N'] =  standard['S'] =  standard['P'] =  1

    ## (1) calculate distance map ##
    for a in range (self.Natom):
      self.atoms[a].Nnei_heavy = 0
      self.atoms[a].Nneighbor  = 0
      for b in range (a,self.Natom):
        dx = self.atoms[a].x - self.atoms[b].x
        dy = self.atoms[a].y - self.atoms[b].y
        dz = self.atoms[a].z - self.atoms[b].z
        dmap[a][b] = dmap[b][a] = math.sqrt(dx*dx + dy*dy + dz*dz)

    ## (2) calculate contable[][] ##
    end = 0
    while (end==0):
      minD = -1.0
      for a in range (self.Natom):
        ele_a = self.atoms[a].element
        if (ele_a != 'H'):
          for b in range (self.Natom):
            ele_b = self.atoms[a].element
            if (ele_b != 'H'):
               if (self.contable[a][b] == '0'): 
                 if (self.atoms[a].Nnei_heavy < Nbond_max.get(ele_a,4)): 
                   if (self.atoms[b].Nnei_heavy < Nbond_max.get(ele_b,4)):
                     D = dmap[a][b]
                     if (standard.has_key(ele_a)==1) and (standard.has_key(ele_b)==1):
                       Dbond_min = 1.1
                       Dbond_max = 1.9
                     else: 
                       Dbond_min = 0.9
                       Dbond_max = 2.4
                     if (D>Dbond_min) and (D<Dbond_max) and ((minD<0.0) or (D<minD)):
                       minD = D
                       min_a = a
                       min_b = b
      if (minD>=0.0):
        #print "a %d b %d"%(min_a,min_b)
        self.contable[min_a][min_b] = self.contable[min_b][min_a] = '1'
        self.atoms[min_a].Nnei_heavy += 1 
        self.atoms[min_b].Nnei_heavy += 1 
        self.atoms[min_a].Nneighbor += 1 
        self.atoms[min_b].Nneighbor += 1 
      else:
        end = 1

    ## (3) calculate contable[][] for hydrogen ##
    Dbond_min = 0.90
    Dbond_max = 1.20
    end = 0
    while (end==0):
      minD = -1.0
      for a in range (self.Natom):
        ele_a = self.atoms[a].element
        if (ele_a == 'H') and (self.atoms[a].Nneighbor < 1):
          #print a,ele_a,self.atoms[a].Nneighbor
          for b in range (self.Natom):
            ele_b = self.atoms[b].element
            if (ele_b != 'H'):
              if (self.contable[a][b] == '0'): 
                #print a,ele_a,b,ele_b,dmap[a][b],self.atoms[a].Nneighbor,self.atoms[b].Nneighbor
                D = dmap[a][b]
                if (D>Dbond_min) and (D<Dbond_max) and ((minD<0.0) or (D<minD)):
                  minD = D
                  min_a = a
                  min_b = b
      if (minD>=0.0):
        self.contable[min_a][min_b] = self.contable[min_b][min_a] = '1'
        self.atoms[min_a].Nneighbor += 1 
        self.atoms[min_b].Nneighbor += 1 
        #print "make hydrogen bonding",min_a,min_b
      else:
        end = 1

    ## (4) make bond[]  ##
    self.Nbond = 0
    for a in range (self.Natom):
      for b in range (self.Natom):
        if (self.contable[a][b] != '0'):
          B = Bond()
          self.bonds.append(B)
          self.bonds[self.Nbond].natom1   = a
          self.bonds[self.Nbond].natom2   = b
          self.bonds[self.Nbond].bondtype = self.contable[a][b]
          self.Nbond += 1


  def count_neighbor_atoms(self):
    for atm in (self.atoms):
      self.atoms[atm.num].neighbors = []
      atm.Nneighbor = 0 
      for i in range(self.Natom):
        if (self.contable[atm.num][i]!='0'):
          atm.Nneighbor += 1
          self.atoms[atm.num].neighbors.append(i)
          if (self.atoms[i].element=='C'):
            atm.NneiC += 1
          elif (self.atoms[i].element=='O'):
            atm.NneiO += 1
          elif (self.atoms[i].element=='N'):
            atm.NneiN += 1
          elif (self.atoms[i].element=='H'):
            atm.NneiH += 1
          elif (self.atoms[i].element=='X'):
            atm.NneiX += 1
        if (self.contable[atm.num][i]=='1'):
          atm.Nbond_single += 1
        if (self.contable[atm.num][i]=='2'):
          atm.Nbond_double += 1


  def set_atomtype(self,AtomTypeMethod):
    for atm in (self.atoms):
        
      if (AtomTypeMethod=='B'):
        atm.atomtype  = atm.element
        atm.atomtype += 'C' * atm.NneiC 
        atm.atomtype += 'O' * atm.NneiO 
        atm.atomtype += 'N' * atm.NneiN 
        atm.atomtype += 'H' * atm.NneiH 
        atm.atomtype += 'X' * atm.NneiX 
      
      if (AtomTypeMethod=='BE'):
        atm.atomtype  = atm.element
        atm.atomtype += 'C' * atm.NneiC 
        atm.atomtype += 'O' * atm.NneiO 
        atm.atomtype += 'N' * atm.NneiN 
        atm.atomtype += 'H' * atm.NneiH 
        atm.atomtype += 'X' * atm.NneiX 
        atm.atomtype += "%d"%(atm.EC1)
      
      if (AtomTypeMethod=='E'):
        atm.atomtype  = atm.element
        atm.atomtype += "%d"%(atm.EC1)


      if (atm.element=='C') and (AtomTypeMethod=='2'):
        if (Nbond_double==1): 
          atm.atomtype = 'c'
        else:
          atm.atomtype = 'C'

      if (atm.element=='C') and (AtomTypeMethod=='4'):
        if (Nbond_double==1) and (atm.NneiC==2):
          atm.atomtype = 'c'
        if (Nbond_double==1) and (atm.NneiO>0):
          atm.atomtype = 'cO'
        if (Nbond_double==1) and (atm.NneiN>0):
          atm.atomtype = 'cN'

      if (atm.element=='C') and (AtomTypeMethod=='7'):
        if (Nbond_double==1):
          if (atm.NneiH==0): 
            atm.atomtype = 'c0'
          if (atm.NneiH==1): 
            atm.atomtype = 'c1'
          if (atm.NneiH==2): 
            atm.atomtype = 'c2'
        if (Nbond_double==0):
          if (atm.NneiH==0): 
            atm.atomtype = 'C0'
          if (atm.NneiH==1): 
            atm.atomtype = 'C1'
          if (atm.NneiH==2): 
            atm.atomtype = 'C2'
          if (atm.NneiH==3): 
            atm.atomtype = 'C3'

        print "#atom %d ele %s BOND s %d d %d NEI C %d O %d N %d H %d --> %s"%(atm.num,atm.element,Nbond_single,Nbond_double,NneiC,NneiO,NneiN,NneiH,atm.atomtype)  


  def set_DABRLtype(self):
    for atm in (self.atoms):
      if (atm.element=='C'):
        if (atm.Nneighbor == 4):
          atm.DABRLtype = 'L'  ##  a'L'iphatic  
        elif (atm.Nneighbor == 3):
          if (atm.NneiO>0):
            atm.DABRLtype = 'L'  ##  a'L'iphatic  
          else:
            atm.DABRLtype = 'R'  ## aromatic   
        else:
          atm.DABRLtype = 'L'  ## a'L'iphatic   
      elif (atm.element=='O'):
        if (atm.NneiH>0):
          atm.DABRLtype = 'B'  ##  'B'oth acceptor and doner  
        else: 
          atm.DABRLtype = 'A'  ##  'A'cceptor 
      elif (atm.element=='N'):
        if (atm.Nneighbor==2):
          atm.DABRLtype = 'A'  ##  'A'cceptor  
        elif (atm.NneiH>0):
          atm.DABRLtype = 'D'  ##  'D'oner 
        elif (atm.NneiH==0):
          atm.DABRLtype = 'R'  ##   sp2 (aromatic) 
      else:
        atm.DABRLtype = 'X'  ##    (others) 
      

  def set_min_max_XY_molecule(self):
    if (len(self.atoms)==0):
      return([0.0, 0.0, 0.0, 0.0])

    self.minX = self.maxX = self.atoms[0].x
    self.minY = self.maxY = self.atoms[0].y
    for a in range(1,self.Natom):
      if (self.minX > self.atoms[a].x):
        self.minX = self.atoms[a].x
      if (self.maxX < self.atoms[a].x):
        self.maxX = self.atoms[a].x
      if (self.minY > self.atoms[a].y):
        self.minY = self.atoms[a].y
      if (self.maxY < self.atoms[a].y):
        self.maxY = self.atoms[a].y
    return([self.minX,self.maxX,self.minY,self.maxY])


  def cal_gcenter(self):
    gcen = [0.0,0.0,0.0]
    for a in range(self.Natom):
      gcen[0] += self.atoms[a].x
      gcen[1] += self.atoms[a].y
      gcen[2] += self.atoms[a].z
    for k in range(3):
      gcen[k] /= self.Natom
    return(gcen)
 

  def rotate_gcenter(self,R,gcen):
    pos  = [0.0,0.0,0.0]
    npos = [0.0,0.0,0.0]

    for a in range(self.Natom):
      pos[0] = self.atoms[a].x
      pos[1] = self.atoms[a].y
      pos[2] = self.atoms[a].z
      for i in range(3):
        npos[i] = 0
        for k in range(3):
          npos[i] +=  R[i][k] * (pos[k]-gcen[k])

      self.atoms[a].x = npos[0] + gcen[0]
      self.atoms[a].y = npos[1] + gcen[1]
      self.atoms[a].z = npos[2] + gcen[2]


  def rotate_180_degree_around_Xaxis(self):
    for a in range(self.Natom):
      self.atoms[a].y = -self.atoms[a].y
      self.atoms[a].z = -self.atoms[a].z


  def rotate_and_tranlate(self,raxisang_str,trans_vec_str):
    raxisang = raxisang_str.split(':')
    for i in range (len(raxisang)):
      raxisang[i] = float(raxisang[i])
    #print "raxisang '%s' %f %f %f %f"%(raxisang_str,raxisang[0],raxisang[1],raxisang[2],raxisang[3])
    trans_vec = trans_vec_str.split(':')
    for i in range (len(trans_vec)):
      trans_vec[i] = float(trans_vec[i])

    R =  make_rot_matrix_from_rot_axis_angle(raxisang)
    pos   = [0.0,0.0,0.0]
    npos  = [0.0,0.0,0.0]
    #print "#R0 %f %f %f"%(R[0][0],R[0][1],R[0][2])
    #print "#R1 %f %f %f"%(R[1][0],R[1][1],R[1][2])
    #print "#R2 %f %f %f"%(R[2][0],R[2][1],R[2][2])
    for a in range(self.Natom):
      pos[0] = self.atoms[a].x
      pos[1] = self.atoms[a].y
      pos[2] = self.atoms[a].z

      for i in range(3):
        npos[i] = 0.0
        for k in range(3):
          npos[i] +=  R[i][k] * pos[k]
        npos[i] += trans_vec[i]

      self.atoms[a].x = npos[0]
      self.atoms[a].y = npos[1]
      self.atoms[a].z = npos[2]
      #print "#pos %f %f %f npos %f %f %f"%(pos[0],pos[1],pos[2],npos[0],npos[1],npos[2])

  def rescale_by_average_bond_length(self):
    #print "#def rescale_by_average_bond_length(self):"
    Npair = 0
    ave_dis = 0.0
    standard_dis = 1.5
    for a in range(self.Natom):
      if (self.atoms[a].element != 'H'):
        for b in range(self.Natom):
          if (self.atoms[b].element != 'H'):
            if (self.contable[a][b] != '0'):
              dx = self.atoms[a].x - self.atoms[b].x 
              dy = self.atoms[a].y - self.atoms[b].y 
              dz = self.atoms[a].z - self.atoms[b].z 
              dis = math.sqrt(dx*dx + dy*dy + dz*dz)
              ave_dis += dis
              Npair += 1
    if (Npair>0) and (ave_dis>0.0):
      ave_dis = ave_dis/Npair 
      scale = standard_dis/ave_dis
      #print "#ave_dis %f scale %f"%(ave_dis,scale)
      for a in range(self.Natom):
        self.atoms[a].x *= scale 
        self.atoms[a].y *= scale 
        self.atoms[a].z *= scale 
 
               

  def set_internal_coodinate_topology(self):
    print "#set_internal_coodinate_topology()"
    for a in range(self.Natom):
      self.atoms[a].mark = 0
    ## (1) set the initial three atoms ##
    hit = 0
    a = 0
    for a in range(self.Natom):
      for b in (self.atoms[a].neighbors):
        for c in (self.atoms[b].neighbors):
          if (c != a):
            hit = 1
            break
        if (hit==1):
          break 
      if (hit==1):
        break 
    if (hit==0):
      print "#ERROR:Can't find the initial three atoms" 
      sys.exit(1)

    self.atoms[a].mark    = 1
    self.atoms[a].num_ic  = 0
    self.num_from_num_ic[0]  = a

    self.atoms[b].mark    = 1
    self.atoms[b].num_ic  = 1
    self.num_from_num_ic[1]  = b 
    self.atoms[b].prenum1 = a

    self.atoms[c].mark    = 1
    self.atoms[c].num_ic  = 2
    self.num_from_num_ic[2]  = c 
    self.atoms[c].prenum1 = b
    self.atoms[c].prenum2 = a
    print "init_three %d %d %d"%(self.atoms[a].num_in_file,self.atoms[b].num_in_file,self.atoms[c].num_in_file)
    ## (2) set prenum1, prenum2, prenum3  for other atoms ##
    Nnum_ic = 3 
    end = 0   
    while (end!=1): 
      Naddmark = 0
      Nmarkzero = 0
      for a in range(self.Natom):
        if (self.atoms[a].mark==0):
          Nmarkzero += 1
          hit = 0 
          for b in (self.atoms[a].neighbors):
            if (self.atoms[b].mark==1):
              for c in (self.atoms[b].neighbors):
                if (c!=a) and (self.atoms[c].mark==1):
                  for d in (self.atoms[c].neighbors):
                    if (d!=a) and (d!=b) and (self.atoms[d].mark==1):
                      hit = 1
                      Naddmark += 1
                      self.atoms[a].mark   = 1
                      self.atoms[a].num_ic = Nnum_ic 
                      self.num_from_num_ic[Nnum_ic]  = a
                      self.atoms[a].prenum1 = b
                      self.atoms[a].prenum2 = c
                      self.atoms[a].prenum3 = d
                      Nnum_ic += 1
                    if (hit==1):
                      break
                if (hit==1):
                  break
            if (hit==1):
              break
      #print "ROUND Naddmark %d Nmarkzero %d"%(Naddmark,Nmarkzero) 
      if (Naddmark==0):
        end = 1

  def set_internal_coodinate_angle(self):
    d    = [0.0, 0.0, 0.0]
    p    = [0.0, 0.0, 0.0]
    p1   = [0.0, 0.0, 0.0]
    p2   = [0.0, 0.0, 0.0]
    p3   = [0.0, 0.0, 0.0]
    p1p  = [0.0, 0.0, 0.0]
    p1p2 = [0.0, 0.0, 0.0]
 
    for a in range(self.Natom):
      set_vec3(p, self.atoms[a].x, self.atoms[a].y, self.atoms[a].z)
      prenum1 = self.atoms[a].prenum1
      prenum2 = self.atoms[a].prenum2
      prenum3 = self.atoms[a].prenum3
      if (prenum1>=0):
        set_vec3(p1, self.atoms[prenum1].x, self.atoms[prenum1].y, self.atoms[prenum1].z)
        sub_vec3(d,p,p1)
        self.atoms[a].bond_len = len_vec3(d)
      if (prenum1>=0) and (prenum2>=0):
        set_vec3(p2, self.atoms[prenum2].x, self.atoms[prenum2].y, self.atoms[prenum2].z)
        self.atoms[a].bond_ang = bond_angle(p,p1,p2)
      if (prenum1>=0) and (prenum2>=0) and (prenum3>=0):
        set_vec3(p3, self.atoms[prenum3].x, self.atoms[prenum3].y, self.atoms[prenum3].z)
        self.atoms[a].torsion_ang = torsion_angle(p,p1,p2,p3)

      #sys.stdout.write("len %8.3f bang %8.3f tang %8.3f"%(self.atoms[a].bond_len, self.atoms[a].bond_ang, self.atoms[a].torsion_ang))
      #sys.stdout.write(" %2s "%(self.atoms[a].num_in_file))
      #if (prenum1>=0):
      #  sys.stdout.write("%2d "%(self.atoms[prenum1].num_in_file))
      #if (prenum2>=0):
      #  sys.stdout.write("%2d "%(self.atoms[prenum2].num_in_file))
      #if (prenum3>=0):
      #  sys.stdout.write("%2d "%(self.atoms[prenum3].num_in_file))
      #sys.stdout.write("\n")


  def set_xyz_coodinates_from_internal(self):
    print "#set_xyz_coodinates_from_internal()"
    ex   = [0.0, 0.0, 0.0]
    ey   = [0.0, 0.0, 0.0]
    ez   = [0.0, 0.0, 0.0]
    p1  = [0.0, 0.0, 0.0]
    p2  = [0.0, 0.0, 0.0]
    p3  = [0.0, 0.0, 0.0]
    p23 = [0.0, 0.0, 0.0]
   
    for i in range(self.Natom):
      num  = self.num_from_num_ic[i]
      if (i<3):
        self.atoms[num].mark = 1
      else:
        self.atoms[num].mark = 0
        self.atoms[num].x = 0.0
        self.atoms[num].y = 0.0
        self.atoms[num].z = 0.0
 
    for i in range(self.Natom):
      num  = self.num_from_num_ic[i]
      atm  = self.atoms[num]
      pre1 = self.atoms[num].prenum1
      pre2 = self.atoms[num].prenum2
      pre3 = self.atoms[num].prenum3
  
      if (pre1>=0) and (pre2>=0) and (pre3>=0):
        if (self.atoms[pre1].mark==0) or (self.atoms[pre2].mark==0) or (self.atoms[pre3].mark==0):
          print "woops!!"

        #print "%d %d %d %d"%(num,self.atoms[pre1].mark,self.atoms[pre2].mark,self.atoms[pre3].mark)
        set_vec3(p1,self.atoms[pre1].x, self.atoms[pre1].y, self.atoms[pre1].z)
        set_vec3(p2,self.atoms[pre2].x, self.atoms[pre2].y, self.atoms[pre2].z)
        set_vec3(p3,self.atoms[pre3].x, self.atoms[pre3].y, self.atoms[pre3].z)
        sub_vec3(ez,p2,p1)
        normalize(ez)
        sub_vec3(p23,p3,p2)
        cross_prod(ey,p23,ez)
        normalize(ey)
        cross_prod(ex,ez,ey)
        cos_t = math.cos(atm.bond_ang*math.pi/180.0) 
        sin_t = math.sin(atm.bond_ang*math.pi/180.0) 
        cos_p = math.cos(atm.torsion_ang*math.pi/180.0) 
        sin_p = math.sin(atm.torsion_ang*math.pi/180.0) 
        x = atm.bond_len * sin_t * cos_p
        y = atm.bond_len * sin_t * sin_p
        z = atm.bond_len * cos_t
        atm.x = self.atoms[pre1].x + x*ex[0] + y*ey[0] + z*ez[0] 
        atm.y = self.atoms[pre1].y + x*ex[1] + y*ey[1] + z*ez[1] 
        atm.z = self.atoms[pre1].z + x*ex[2] + y*ey[2] + z*ez[2] 
        atm.mark = 1
        #print "ex %f %f %f ey %f %f %f ez %f %f %f"%(ex[0],ex[1],ex[2],ey[0],ey[1],ey[2],ez[0],ez[1],ez[2])

  def write_in_mopin(self,ofname):
    ## mopin : MOPAC internal coordinate
    print "#write_in_mopin()-->'%s'"%(ofname)
    of = open(ofname,"w")
    of.write("PUT KEYWORDS HERE\n")
    of.write("%s\n"%(self.filename))
    of.write("\n")
    for i in range(self.Natom):
      if (self.num_from_num_ic.has_key(i)==0):
        print "#key error for '%s'"%(i)
      num  = self.num_from_num_ic[i]
      atm  = self.atoms[num]
      pre1 = self.atoms[num].prenum1
      pre2 = self.atoms[num].prenum2
      pre3 = self.atoms[num].prenum3
      pre1_ic = pre2_ic = pre3_ic = 0 
      if (pre1>=0):
        pre1_ic = self.atoms[pre1].num_ic + 1
      if (pre2>=0):
        pre2_ic = self.atoms[pre2].num_ic + 1
      if (pre3>=0):
        pre3_ic = self.atoms[pre3].num_ic + 1
      ##N    0.000000  1    0.000000  1    0.000000  1     0   0   0
      ##N    0.000000  1    0.000000  1    0.000000  1     0   0   0
      ##C    1.461655  1    0.000000  1    0.000000  1     1   0   0
      ##C    1.522128  1  113.626770  1    0.000000  1     2   1   0
      ##O    1.232197  1  121.101292  1  337.712286  1     3   2   1
      ##01 0123456789012  0123456789012  0123456789012  012301230123 
      of.write("%-2s%11.6f  1 %11.6f  1 %11.6f  1  %4d%4d%4d\n"%(atm.element,atm.bond_len,atm.bond_ang,atm.torsion_ang,pre1_ic,pre2_ic,pre3_ic))
    of.close()

  def cal_scale_for_gwidth_gheight(self,gwidth,gheight):
    scx = gwidth /(self.maxX - self.minX)
    scy = gheight/(self.maxY - self.minY)
    if (scx<scy):
      sc = scx
    else:
      sc = scy
    return(sc)

  def __str__(self):
    return "#sdf.Molecule (Natom %d Nbond %d '%s')"%(self.Natom,self.Nbond,self.formula)
    pass




#########################
### Vector functions ####
#########################

def set_vec3(p,x,y,z):
  p[0] = x
  p[1] = y
  p[2] = z

def sub_vec3(c,a,b):
  c[0] = a[0] - b[0]
  c[1] = a[1] - b[1]
  c[2] = a[2] - b[2]

def len_vec3(p):
  dd = p[0]*p[0] + p[1]*p[1] + p[2]*p[2]
  return(math.sqrt(dd))

def dot_prod(a,b):
  return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

def normalize(p):
  d = math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
  p[0] /= d
  p[1] /= d
  p[2] /= d

def cross_prod(c,a,b):
  c[0] = a[1]*b[2] - a[2]*b[1]
  c[1] = a[2]*b[0] - a[0]*b[2]
  c[2] = a[0]*b[1] - a[1]*b[0]

def bond_angle(po,pa,pb):
  ## po--pa--pb ##
  ao = [0.0,0.0,0.0]
  ab = [0.0,0.0,0.0]
  sub_vec3(ao,po,pa)
  sub_vec3(ab,pb,pa)
  cos_bang = dot_prod(ao,ab)/len_vec3(ao)/len_vec3(ab)
  return(math.acos(cos_bang)*180/math.pi)


def torsion_angle(o,a,b,c):
  ## o--a--b--c ##
  ## c--b--a--o ##
  ao = [0.0,0.0,0.0]
  ab = [0.0,0.0,0.0]
  ba = [0.0,0.0,0.0]
  bc = [0.0,0.0,0.0]
  bo = [0.0,0.0,0.0]
  noab = [0.0,0.0,0.0]
  nabc = [0.0,0.0,0.0]
  sub_vec3(ao,o,a)
  sub_vec3(ab,b,a)
  sub_vec3(ba,a,b)
  sub_vec3(bc,c,b)
  sub_vec3(bo,o,b)
  cross_prod(noab,ao,ab)
  cross_prod(nabc,ba,bc)
  cos_torsion = dot_prod(noab,nabc)/len_vec3(noab)/len_vec3(nabc)
  dprod = dot_prod(nabc,bo)
  if (cos_torsion>1.0):
    cos_torsion = 1.0 
  if (cos_torsion<-1.0):
    cos_torsion = -1.0 
  angle = math.acos(cos_torsion)*180/math.pi
  if (dprod>=0.0):
    return(angle)
  if (dprod<0.0):
    return(-angle)
  



class Atom:
  def __init__(self):
    self.x = 0.0
    self.y = 0.0
    self.z = 0.0
    self.num          = '' # Atom number (0,1,2,...) 
    self.num_in_file  = '' # Atom number specified in file (1,2,3,...)
    self.element  = ''
    self.atomname = ''
    self.atomtype = ''
    self.ligname  = ''
    self.stereo_parity = '0' # '0':not stereo,'1':clockwise, '2':counter-clockwise '3':either or unmarked
    self.charge   = 0.0
    self.DABRLtype = '' # 'D'oner, 'A'cceptor, 'B'oth, sp'2', sp'3'      
    self.Nneighbor = 0 
    self.neighbors = [] # atom numbers for nerighboring atoms 
    self.NneiC = 0  # number of connected 'C' atoms
    self.NneiO = 0  # number of connected 'O' atoms
    self.NneiN = 0  # number of connected 'N' atoms
    self.NneiH = 0  # number of connected 'H' atoms
    self.NneiX = 0  # number of connected  other atoms

    self.Nbond_single = 0
    self.Nbond_double = 0
    self.num_ic  = -1 # atom_num for internal coodinates
    self.prenum1 = -1 # atom_num for one   before this atom (for internal coodinates)
    self.prenum2 = -1 # atom_num for two   before this atom (for internal coodinates)
    self.prenum3 = -1 # atom_num for three before this atom (for internal coodinates)
    self.bond_len = 0.0    # bond    length                      (for internal coodinates)
    self.bond_ang = 0.0    # bond    angle  (degree)             (for internal coodinates)
    self.torsion_ang = 0.0 # torsion angle  (degree)             (for internal coodinates)

    self.mark = 0     # mark for various purposes
    pass
  def __str__(self):
    return("sdf.Atom(xyz %d %d %d element %s)"%(self.x,self.y,self.z,self.element))



class Bond:
  def __init__(self):
    self.num         = 0  # bonds number (0,1,2,...)
    self.natom1      = 0  # 0,1,...,mol.Natom-1
    self.natom2      = 0  # 0,1,...,mol.Natom-1
    self.bondtype    = '' 
    self.bondstereo = '0'  # '0':not stereo,'1':Up, '6':Down, '4':either #
    pass
  def __str__(self):
    return("sdf.Bond(natom %d %d bondtype %s)"%(self.natom1,self.natom2,self.bondtype))


def make_molecular_formula(Nele):
  ele_list = ['C','H','HE','LI','BE','B','N','O','F','NE','NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE'] 

  molform = '' 
  for e in (ele_list):
    if (Nele.has_key(e)):
      if (molform !=''):
        molform += ' ' 
      molform += e 
      if (Nele[e]>1):
        molform = molform + "%d"%(Nele[e]) 
  return(molform)
  pass

def cal_Nele_from_molecular_formulra(Nele,exist,molform):
  # ex) 'C O4',  'C23 H35 N10 O18 P3',  'C19 H39 N O5 S','C34 H32 N4 O4 FE'
  Natom = 0
  num = molform.split()
  for x in (num):
    e = ''
    n = ''

    for i in range(len(x)):
      if (x[i].isalpha()):
        e = e + x[i]
      else:
        n = n + x[i]

    if (n==''):
      n = '1'
    if (e != 'H'):
      exist[e] = 1
      Nele[e] = int(n)
      Natom += int(n) 
  return(Natom)


def compare_two_molecular_formulas(molformA,molformB):
  # ex) 'C O4',  'C23 H35 N10 O18 P3',  'C19 H39 N O5 S','C34 H32 N4 O4 FE'

  exist = {}
  NeleA = {}
  NatomA =  cal_Nele_from_molecular_formulra(NeleA,exist,molformA)
  NeleB = {}
  NatomB =  cal_Nele_from_molecular_formulra(NeleB,exist,molformB)

  Natom_common = 0
  for e in (exist.keys()):
    if (NeleA.has_key(e)==0) or (NeleB.has_key(e)==0): 
      nmin = 0  
    elif (NeleA[e] < NeleB[e]): 
      nmin = NeleA[e]
    else:
      nmin = NeleB[e]
    Natom_common += nmin

  print "Natom_common %d NatomA %d NatomB %d"%(Natom_common,NatomA,NatomB)
  tanimoto = float(Natom_common)/float(NatomA+NatomB-Natom_common)
  return([Natom_common,tanimoto]) 
  pass


def make_rot_matrix_from_rot_axis_angle(w):
  ## w[4]:  rot axis [0:1:2] and rot_angle [3]  (input) 

  Ux = [[0.0,-w[2],w[1]], [w[2],0.0,-w[0]], [-w[1],w[0],0.0]]

  #Ux[0][0] =   0.0  Ux[0][1] = -w[2] Ux[0][2] =  w[1] 
  #Ux[1][0] =  w[2]  Ux[1][1] =   0.0 Ux[1][2] = -w[0];
  #Ux[2][0] = -w[1]; Ux[2][1] =  w[0];  Ux[2][2] =  0.0;
  Uxx = [[w[0]*w[0],w[0]*w[1],w[0]*w[2]],[w[1]*w[0],w[1]*w[1],w[1]*w[2]],[w[2]*w[0],w[2]*w[1],w[2]*w[2]]]

  #Uxx[0][0] =  w[0]*w[0]; Uxx[0][1] = w[0]*w[1]; Uxx[0][2] = w[0]*w[2];
  #Uxx[1][0] =  w[1]*w[0]; Uxx[1][1] = w[1]*w[1]; Uxx[1][2] = w[1]*w[2];
  #Uxx[2][0] =  w[2]*w[0]; Uxx[2][1] = w[2]*w[1]; Uxx[2][2] = w[2]*w[2];

  cos_t = math.cos(w[3])
  sin_t = math.sin(w[3])
  one_cos_t = 1.0 - cos_t
  R = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  for i in range(3): 
    for j in range(3): 
      R[i][j] = sin_t * Ux[i][j] + one_cos_t*Uxx[i][j]
      if (i==j):
        R[i][j] += cos_t;
  #print "<PRE>"
  #print "w %f %f %f %f"%(w[0],w[1],w[2],w[3])
  #print "R0 %f %f %f"%(R[0][0],R[0][1],R[0][2])
  #print "R1 %f %f %f"%(R[1][0],R[1][1],R[1][2])
  #print "R2 %f %f %f"%(R[2][0],R[2][1],R[2][2])
  #print "</PRE>"
  return(R)



def _main():
  if (len(sys.argv)<2):
    print "molecule.py [sdf|mol|pdb|kcf|mol2 file]"
    sys.exit(1)
  m = Molecule()
  if (sys.argv[1].endswith(".sdf")): 
    m.read_in_sdf(sys.argv[1])
  elif (sys.argv[1].endswith(".mol")): 
    m.read_in_sdf(sys.argv[1])
  elif (sys.argv[1].endswith(".pdb")): 
    m.read_in_pdb(sys.argv[1])
    m.guess_bonds_from_atom_xyz()
  elif (sys.argv[1].endswith(".kcf")): 
    m.read_in_kcf(sys.argv[1])
  elif (sys.argv[1].endswith(".mol2")): 
    m.read_in_mol2(sys.argv[1])
  else:
    print "#no proper input file"
    sys.exit(1)
 
  print "Natom %d"%(m.Natom) 
  m.count_neighbor_atoms()
  m.set_DABRLtype()
  print "formula '%s'"%(m.formula)

  for key in m.annotation.keys():
    print key,m.annotation[key]
  



  sys.exit(1)
  Nele = {}
  exist = {} 
  cal_Nele_from_molecular_formulra(Nele,exist,m.formula)
  for e in (Nele.keys()):
    print "%s %d"%(e,Nele[e]) 
  #[Ncatm,tanimoto] = compare_two_molecular_formulas(m.formula,m.formula)
  #print "Ncatm %d tanimoto %f"%(Ncatm,tanimoto)


  m.set_internal_coodinate_topology()
  m.set_internal_coodinate_angle()

  #m.write_in_sdf("out0.sdf")
  #m.write_in_sdf('stdout')
  #m.write_in_pdb('stdout')
  m.write_in_mol2('stdout')
  m.write_in_mol2('out.mol2')
  sys.exit(1)
  m.write_in_pdb("out0.pdb")
  m.write_in_mopin("out.mopin")
  m.set_xyz_coodinates_from_internal()
  m.write_in_sdf("out1.sdf")
  m.write_in_pdb("out1.pdb")

  #m.atoms[12].torsion_ang += 10.0
  #m.atoms[20].torsion_ang += 10.0
  #m.atoms[12].torsion_ang = -60 
  #m.atoms[20].torsion_ang = -60 
  #m.set_xyz_coodinates_from_internal()
  #m.set_internal_coodinate_angle()
  #m.write_in_sdf("out2.sdf")
  #m.write_in_pdb("out2.pdb")
  #m.write_in_mopin("out2.mopin")

  sys.exit(1)
  
  if (len(sys.argv)>2):
    m.set_atomtype(sys.argv[2])
  else:
    m.set_atomtype('B')
  print m
  
  if (m.Nele.has_key('C')) and (m.Nele.has_key('H')) and (m.Nele.has_key('O')):
    print "%f %f %s %s KREVELEN"%(float(m.Nele['O'])/m.Nele['C'],float(m.Nele['H'])/m.Nele['C'],m.filename,m.formula)
  sys.exit(1)
  for a in range(m.Natom):
    print "ATOM %2d x elem '%s' DABRL %s atom_type '%s'"%(a+1,m.atom[a].element, m.atom[a].DABRLtype,m.atom[a].atomtype)
  sys.exit(1) 
  for b in range(m.Nbond):
    print "BOND %d %d bondtype '%s'"%(m.bond[b].natom1,m.bond[b].natom2,m.bond[b].bondtype)
  pass

if __name__ == '__main__': _main()
