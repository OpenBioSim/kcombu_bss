#!/usr/bin/env python
import os
import sys

LastModDate = "2020/04/11"

if (len(sys.argv)<2):
  print "ext_chain.py [pdbfile] [chain_id]"
  sys.exit(1)

ifname = sys.argv[1]
chain_id = sys.argv[2]
if (os.access(ifname,os.R_OK)==0):
  print "#ERROR:Can't open '%s'"%(ifname)
  sys.exit(1)

f = open(ifname)

#          1         2         3
#0123456789012345678901234567890123456789
#ATOM      1  N   ASP A   3     -11.380 -17.161  44.850  1.00 36.55           N  
#ATOM      2  CA  ASP A   3     -10.861 -17.534  43.496  1.00 35.06           C  
#ATOM      3  C   ASP A   3     -11.465 -16.737  42.339  1.00 33.73           C  
#ATOM      4  O   ASP A   3     -10.752 -16.368  41.464  1.00 32.38           O  
#:
#HETATM 1847  N   SEP P   5     -17.798   2.808  16.569  1.00 20.71           N  
#HETATM 1848  CA  SEP P   5     -18.799   2.054  17.299  1.00 19.24           C  
#HETATM 1849  CB  SEP P   5     -18.729   0.588  16.824  1.00 17.23           C  
#:
#HETATM 1871  O   HOH A2001     -10.109 -20.756  20.635  1.00 57.40           O
#HETATM 1872  O   HOH A2002     -12.936  -8.790  42.591  1.00 36.47           O
#HETATM 1873  O   HOH A2003      -9.945 -18.455  20.908  1.00 46.93           O
#HETATM 1874  O   HOH A2004     -10.139  -9.246  42.759  1.00 30.88           O


for line in f:
  if (line.startswith('ATOM') or line.startswith('HETATM')):
    if (line[21]==chain_id) and (line[17:20]!="HOH"):
      sys.stdout.write("%s"%(line))

f.close()
