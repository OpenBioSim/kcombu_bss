#!/usr/bin/env python
import sys
import os
import kcombu

LastModDate = "June 13, 2011"

############
### MAIN ###
############

OPT = {}
OPT['L'] = ''
OPT['idl'] = 'SupLIG'
OPT['A'] = 'F'
OPT['idR'] = 'AtmMch';
OPT['idT'] = 'AtmMch';
OPT['idS'] = '';
OPT['cmp'] = 'aa'
OPT['dset'] = 'lig_HIV:lig_THR:lig_CDK2:lig_CAH2:lig_NEU'
OPT['A'] = 'F'
OPT['cmpcom'] = '/home/takawaba/work/kcombu/src/python/cmp_atommatch.py'
OPT['stacom'] = '/home/takawaba/work/kcombu/src/python/stat_range.py'
OPT['cmp'] = 'ta' 
OPT['L'] = 'ligand.list'
OPT['idT'] = 'AM_BuC'
OPT['idR'] = 'AM_ExC'
OPT['oc']  = 'BuC_ExC'
OPT['M'] = 'C'
OPT['x'] = '8'
OPT['y'] = '-1'
OPT['ymin'] = '0.0'
OPT['ymax'] = '1.0'

if (len(sys.argv)<2):
  print "allcmp_am.py <options>"
  print " for doing 'cmp_atommatch.py' or 'stat_range.py' for various dataset."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>" 
  print " -M      : Mode. 'C':make atom_match comparision file 'S':show statisitical agreement [%s]"%(OPT['M']) 
  print " -dset   : dataset list [%s]"%(OPT['dset'])  
  print " -cmpcom : command for 'cmp_atommatch.py' [%s]"%(OPT['cmpcom'])  
  print " -stacom : command for 'stat_range.py' [%s]"%(OPT['stacom'])  
  print " -L      :list of library molecules[%s]"%(OPT['L'])
  print " -idl    :input directory for library molecules[%s]"%(OPT['idl'])
  print " -idT    :input directory (Test)      for atom_maching file [%s]"%(OPT['idT'])
  print " -idR    :input directory (Reference) for atom_maching file [%s]"%(OPT['idR'])
  print " -oc     :output atom_match comparison result file [%s]"%(OPT['oc'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])  
  print "<options for 'cmp_atommatch.py'>" 
  print " -cmp    :comparison(test_vs_ref) 'tt':top-vs-top,'aa':average-vs-average" 
  print "         :  'at':average-vs-top,'ta':top-vs-average 'am':ave-vs-max [%s]"%(OPT['cmp']) 
  print "<options for 'stat_range.py'>" 
  print " -x      :column for 'x'(target) (1,2,...) [%s]"%(OPT['x']) 
  print " -y   : column for 'y'(range restriction) (1,2,...) -1:accept everything[%s]"%(OPT['y'])
  print " -ymin: minimum value for y (>ymin)[%s]"%(OPT['ymin'])
  print " -ymax: maximum value for y (<=ymax)[%s]"%(OPT['ymax'])
  sys.exit(1)

### [1] read option ###

kcombu.read_option(sys.argv,OPT)

dsetlist = OPT['dset'].split(':')

allfile = ''
for dset in (dsetlist):
  #print ">%s"%(dset)
  if (OPT['M']=='C'):
    command = "%s -L %s/%s -idT %s/%s -idR %s/%s -cmp %s > %s/COMP_AM/%s"%(OPT['cmpcom'],dset,OPT['L'],dset,OPT['idT'],dset,OPT['idR'],OPT['cmp'],dset,OPT['oc'])
    print "#%s"%(command)

  if (OPT['M']=='S'):
    command = "%s %s/COMP_AM/%s -x %s -y %s -ymin %s -ymax %s"%(OPT['stacom'],dset,OPT['oc'],OPT['x'],OPT['y'],OPT['ymin'],OPT['ymax'])
    if (OPT['A']!='T'):
      print "#%s"%(command)

  if (OPT['A']=='T'):
    os.system(command)

  allfile = allfile + "%s/COMP_AM/%s"%(dset,OPT['oc']) 
  allfile = allfile + ':' 


if (OPT['M']=='S'):
  command = "%s %s -x %s -y %s -ymin %s -ymax %s"%(OPT['stacom'],allfile,OPT['x'],OPT['y'],OPT['ymin'],OPT['ymax'])
  if (OPT['A']!='T'):
    print "#%s"%(command)
  elif (OPT['A']=='T'):
    os.system(command)
