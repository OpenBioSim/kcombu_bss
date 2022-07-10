#!/usr/bin/env python
##
## <split_ccd.py>
##
## for splitting a chamical component dictionary cif files into many one-compound cif files



import sys
import re
import math
import os
from datetime import datetime

LastModDate = '2020/04/01'


def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        opt_dic[argv[i][1:]] = argv[i+1]


def make_output_filename(comp_id,opt_odir,opt_subdir):
  ofname = ''
  if (opt_subdir  == 'T'): 
     odir  = opt_odir +  '/' + comp_id[0] 
     ofname = opt_odir +  '/' + comp_id[0] + '/' + comp_id  + '.cif'
  else:
     odir = opt_odir
     ofname = opt_odir + '/' + comp_id  + '.cif'
  return(ofname)

def write_content(ofname,contents,opt_action):
  if (opt_action =='T'):
    fields = ofname.split('/')
    ofilename = fields[-1]
    ofiledir = ofname.rstrip(ofilename)
    print "#write_to --> '%s'"%(ofname)
    if (os.path.exists(ofiledir)==0):
      print "#mkdir '%s'"%(ofiledir)
      os.mkdir(ofiledir)      
    fo = open(ofname,'w')
    fo.write("%s"%(contents))
    fo.close()
  else:
    print "#Contents will be written to --> '%s'"%(ofname)

def read_content(ifname):
  if (os.access(ifname, os.R_OK)==0):
    return('')
  contents = ''
  f = open(ifname)
  for line in f:
    contents += line
  f.close() 
  return(contents)


def write_one_content(comp_id,contents,opt_odir,opt_subdir,opt_action):
  if (opt_subdir  == 'T'): 
     odir  = opt_odir +  '/' + comp_id[0] 
     ofname = opt_odir +  '/' + comp_id[0] + '/' + comp_id  + '.cif'
  else:
     odir = opt_odir
     ofname = opt_odir + '/' + comp_id  + '.cif'
  if (opt_action =='T'):
    print "#'%s' write_to --> '%s'"%(comp_id,ofname)
    if (os.path.exists(odir)==0):
      print "#mkdir '%s'"%(odir)
      os.mkdir(odir)      
    fo = open(ofname,'w')
    fo.write("%s"%(contents))
    fo.close()
  else:
    print "#'%s' will be written to --> '%s'"%(comp_id,ofname)


###############
#### MAIN #####
###############

OPT = {}
OPT['icif'] = ''
OPT['subdir'] = 'T'
OPT['odir'] = 'tmpout'
OPT['A'] = 'F'
OPT['olog'] = 'split_ccd.log'
OPT['onew'] = ''

if (len(sys.argv)<3):
  print "split_ccd.py <options>"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print " for splitting a chamical component dictionary (CCD) cif files into many one-compound cif files"
  print "<options>"
  print "  -icif : input CCD 'components.cif' file [%s]"%(OPT.get('icif',''))
  print "  -odir : output directory [%s]"%(OPT.get('odir','')) 
  print "  -subdir : use first char subdir ex) [odir] A/ATP.sdf  ('T' or 'F') [%s]"%(OPT.get('subdir','')) 
  print "  -onew : output list file for new comp_ids [%s]"%(OPT.get('onew',''))
  print "  -olog : output log file [%s]"%(OPT.get('olog',''))
  print "  -A    : Action. 'T' or 'F' [%s]"%(OPT.get('A',''));  
  sys.exit(1)

read_option(sys.argv,OPT)

iciffile = OPT['icif']
if not os.access(iciffile,os.R_OK):
  print "#ERROR:Can't open ciffile '%s'"%(iciffile)
  sys.exit(1)

if (iciffile.endswith('.gz')):
  command = "zcat %s"%(iciffile)
  fi = os.popen(command)
  print "#os.popen('%s')"%(command)
else:
  fi = open(iciffile)


if (OPT.get('odir','') != '') and (os.path.exists(OPT['odir'])==0) and (OPT.get('A','')=='T'):
  os.mkdir(OPT['odir'])



contents = ''
comp_id  = ''
Nout_mol = 0
new_comp_ids = []

for line in fi:
  #sys.stdout.write("%s"%(line))
  dline = line.rstrip('\n')
  dline = dline.rstrip('\r')

  if (dline.startswith("data_")):
    if (contents != '') and (comp_id != ''):
      Nout_mol += 1
      ofname = make_output_filename(comp_id,OPT.get('odir',''),OPT.get('subdir',''))
      contents_pre = read_content(ofname)
      if (contents_pre != contents):
        new_comp_ids.append(comp_id)     
      write_content(ofname,contents,OPT.get('A',''))
      #write_one_content(comp_id,contents,OPT.get('odir',''),OPT.get('subdir',''),OPT.get('A',''))
    (head,comp_id) = dline.split('_')
    contents = ''
 
  contents = contents + line 

if (contents != '') and (comp_id != ''):
  Nout_mol += 1
  ofname = make_output_filename(comp_id,OPT.get('odir',''),OPT.get('subdir',''))
  write_content(ofname,contents,OPT.get('A',''))
  #write_one_content(comp_id,contents,OPT.get('odir',''),OPT.get('subdir',''),OPT.get('A',''))

fi.close()

print "#Nout_mol %d"%(Nout_mol)

### write new_comp_ids[] --> OPT['onew'] ###
if (OPT.get('onew','') != ''):
  of = open(OPT['onew'],'w')
  print "#write_new_comp_ids(Nnew_comp_ids:%d) --> '%s'"%(len(new_comp_ids),OPT['onew'])
  of.write("#NEW COMP_IDs LIST\n")
  of.write("#COMMAND %s\n"%(OPT.get('COMMAND','')))
  of.write("#DATE    %s\n"%(OPT.get('START_DATE','')))
  of.write("#Nnew_comp_ids %d\n"%(len(new_comp_ids)))
  for comp_id in (new_comp_ids):
    of.write("%s\n"%(comp_id))
  of.close() 


if (OPT.get('olog','') != '') and (OPT.get('A','') == 'T'):
  print "#write_to_logfile --> '%s'"%(OPT['olog'])
  olog = open(OPT['olog'],'w')
  olog.write("#COMMAND %s\n"%(OPT.get('COMMAND','')))
  olog.write("#DATE    %s\n"%(OPT.get('START_DATE','')))
  olog.write("#Nnew_comp_ids %d\n"%(len(new_comp_ids)))
  olog.write("#Nout_mol %d\n"%(Nout_mol))
  for key in (OPT.keys()):
    olog.write("#OPTION '%s' '%s'\n"%(key,OPT[key]))
  olog.close()


