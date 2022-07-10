#
# <kcombuweb_func.py>
#

import sys
import os 
import re
import math


LastModDate = "June 9, 2014"

def full_filename_from_short_filename(shortfile,dirtype,CONF):
## >> For temporary uploaded Query Molecule <<
## if lig_id = '12345.Qmol' --> /var/www/html/ligandbox_out/TMPOUT/12345.Qmol ##

  ## If shortfile contains '..' or '/' or ' or  ", the script stops.
  if (shortfile.find('..')>=0) or (shortfile.find('/')>=0) or (shortfile.find('\'')>=0) or (shortfile.find('\"')>=0):
      print "Content-type: text/html\n";
      print "<HTML><BODY><PRE>"
      print "#ERROR:Improper short_filename('%s')"%(shortfile)
      print "</PRE></BODY></HTML>"
      sys.exit(1)

  if (dirtype=='2'):
    full_filename = CONF.get('PDB_3D_LIGAND_DIR','') + '/' + shortfile 
  elif (dirtype=='3'):
    full_filename = CONF.get('PDB_3D_LIGAND_DIR','') + '/' + shortfile 
  elif (dirtype=='T'):
    full_filename = CONF.get('TMPOUT_DIR','') + '/' + shortfile 
  else:
    print "Content-type: text/html\n";
    print "<HTML><BODY><PRE>"
    print "#ERROR:Improper short_filename('%s')"%(shortfile)
    print "</PRE></BODY></HTML>"
    sys.exit(1)

  return(full_filename)


def read_environment_file(ienvfile,envdic):
# print "#read_environment_file(ienvfile '%s')"%(ienvfile)
# # FORMAT EXAMPLE ##
# #COMMENT
# BASE_DIR  /usr/people/takawaba/work/STRMAT/Matras10
# SCORE_DIR /usr/people/takawaba/work/STRMAT/Matras10/ROM-00Oct6
# BSSP_DIR  /lab/DB/BSSP
# PDB_DIR   /lab/PDB
# LBOXDIR    /home/WEBDB/201110
# DOWNLOAD_DIR $LBOXDIR/DOWNLOAD
# PROPERTY_DIR $LBOXDIR/PROPERTIES

  if (os.access(ienvfile,os.R_OK)==False):
    return 0
  f = open(ienvfile)
  for line in f:
   line = line.rstrip('\n')
   if (line.startswith('#')==0) and (len(line)>5):
     field = line.split()
     key = field[0]
     value = field[1]
     for k in (envdic.keys()):
       pattern = '$'+k
       if (value.find(pattern)>=0):
         #print "key '%s' includes '%s'."%(k,pattern)
         value = value.replace(pattern,envdic[k])
     envdic[key] = value

  f.close()

