#
# <kcombuweb_func.py>
#

import sys
import os 
import re
import math


LastModDate = "2019/11/14"

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



