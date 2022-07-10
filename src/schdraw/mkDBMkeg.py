#!/usr/bin/env python
import sys
import os
import anydbm

LastModDate = "July 9, 2011"

#>> FILE EXAMPLE OF KEGG "drug" <<
#ENTRY       D08872                      Drug
#NAME        Besifloxacin hydrochloride (USAN);
#            Besivance (TN)
#FORMULA     C19H21ClFN3O3. HCl
#EXACT_MASS  429.1022
#MOL_WEIGHT  430.3007
#ACTIVITY    Anti-infective
#REMARK      ATC code: S01AX23
#DBLINKS     CAS: 405165-61-9
#            PubChem: 96025555
#ATOM        28
#            1   C8y C    23.5900  -15.6100
#            2   C8y C    23.5900  -17.0100
#            3   C8y C    24.7800  -17.7100
#:
#            4   C8y C    26.0400  -17.0100
#            26  C1x C    21.1967  -20.6045
#            27  N1a N    18.8385  -17.7124
#            28  X   Cl   33.6000  -17.5000
#BOND        30
#            1     1   2 2
#            2     2   3 1
#            3     3   4 2
#            4     4   5 1
#            5     5   6 2
#
#:
#            27   20  21 1
#            28   25  26 1
#            29    2  22 1
#            30   20  27 1 #Up
#///
#ENTRY       D08873                      Drug
#NAME        Betrixaban (USAN)
#FORMULA     C23H22ClN5O3
#EXACT_MASS  451.1411
#MOL_WEIGHT  451.9055
#ACTIVITY    Antithrombotic;
#            Prevention of deep vein thrombosis and pulmonary embolism after surgery
#DBLINKS     CAS: 330942-05-7
#:


if (len(sys.argv)<2):
  print "mkDBMkegg.py [input_KEGG_datafile] [output_dbmfile]"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print " for making DBM file for accessing a following format file"
  print ">[key]"
  print "[content]"
  print "//"  
  sys.exit(1)

ifile    = sys.argv[1]
odbmfile = sys.argv[2] 

f = open(ifile)
fdbm_dic = anydbm.open(odbmfile,'n')

line = 'woops'
while (line != ''):
  fpos_head_of_line = f.tell()
  line = f.readline()
  fpos_end_of_line = f.tell()
  #print "%d '%s'"%(f.tell(),line)
  if (line.startswith('ENTRY')):
    line  = line.rstrip('\n')
    field = line.split()
    key = field[1]
    start = fpos_head_of_line
  if (line.startswith('//')):
    end = fpos_end_of_line
    fdbm_dic[key] = "%d %d"%(start,end)

f.close()
fdbm_dic.close()
