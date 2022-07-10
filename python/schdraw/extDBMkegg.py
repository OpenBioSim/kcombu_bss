#
# <extDBMkegg.py>
#

import sys
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
#012345678901234567890
#          1         2



def translate_kegg_content(dat,content):
  header = ''
  for line in (content.split('\n')):
    if (len(line)>12) and (line.startswith('/')==0) and (line.startswith('#')==0):
      line = line.rstrip('\n')
      headstr = line[0:12].replace(' ','') 
      tailstr = line[12:]
      if (headstr!='') and (len(headstr)>1):
        header = headstr 
      #print "header '%s' headstr '%s'"%(header,headstr)
      if (dat.has_key(header)==0):
        dat[header] = []
      dat[header].append(tailstr)
  pass


def extract_dbm_file_content(entry,idatfile,idbmfile):
  ## (1) read dbmfile ##
  fse_dic = anydbm.open(idbmfile,"r")
  [file_start,file_end] = fse_dic[entry].split(' ')
  fse_dic.close()
  file_length =  int(file_end) - int(file_start)
  #print "#start %s end %s len %d"%(file_start,file_end,file_length)
  ## (2) read datfile ##
  f = open(idatfile)
  f.seek(int(file_start))
  content = f.read(file_length)
  f.close()
  return(content)


def extract_kegg_dbm_file_content(entry,dat,idatfile='/DB/kegg/medicus/drug/drug',idbmfile='/DB/kegg/medicus/drug/drug.db'):
  content = extract_dbm_file_content(entry,idatfile,idbmfile)
  translate_kegg_content(dat,content)
  

def _main():
  if (len(sys.argv)<2):
    print 'extDBMkegg.py [entry] [database_file] [dbm_file]'
    sys.exit(1)

  entry = sys.argv[1]
#  print "#entry %s"%(entry)

  dat = {}
  if (len(sys.argv)>3):
    extract_kegg_dbm_file_content(entry,dat,idatfile=sys.argv[2],idbmfile=sys.argv[3])
  else:
    extract_kegg_dbm_file_content(entry,dat)

#  for head in (dat.keys()):
#    print ">%s"%(head)
#    for x in (dat[head]):
#      print x


if __name__ == '__main__': _main()



