#
# <webfunc.py>
#
# sub rourintines for web server

import os
import sys 
import re
from datetime import datetime
import cgi

LastModDate = '2020/03/18'


def read_environment_file(envdic,filename='../conf/webhomcos.conf',version=''):
##  version should be such as '201211' or '201505'.
  if (version != ''):
    if (re.match(r'^[0-9]+$',version))   and (len(version)<7):
      filename = '../conf/webhomcos_' + version + '.conf'
    else:
      error_exit("version='%s' is not proper."%(version))


## FORMAT EXAMPLE ##
# BASE_DIR  /usr/people/takawaba/work/STRMAT/Matras10
# SCORE_DIR /usr/people/takawaba/work/STRMAT/Matras10/ROM-00Oct6
# BSSP_DIR  /lab/DB/BSSP
# PDB_DIR   /lab/PDB
# DATE_PDB          20150506
# DATE_UNIPROT      2015_05
# RDB_DBNAME        pdb_mmcif_sheep
# JMOL_WEB          http://ipproo.protein.osaka-u.ac.jp/jmol
# BLASTDBDIR        /work/takawaba/mmCIF_SHEEP_MOL/mmCIF_SEQ/$DATE_PDB
# BLASTDBFILE       unitmol_$DATE_PDB.fasta
# UNIPROT_DATFILE   /work/takawaba/DB/UNIPROT/$DATE_UNIPROT/uniprot_sprot.dat
# PRECALC_AC_LIST          /work/takawaba/HOMCOS_PRECALC/$DATE_UNIPROT_HUMAN/AC.list
# PRECALC_AASEQ_DIR        /work/takawaba/HOMCOS_PRECALC/$DATE_UNIPROT_HUMAN/AASEQ
 
  if (os.access(filename,os.F_OK)==False):
    return 0
  f = open(filename)
  for line in f:
   line = line.rstrip('\n')
   if (line.startswith('#')==0) and (len(line)>5):
     field = line.split()
     key   = field[0] 
     value = field[1] 
     for X in (envdic.keys()):
       value = value.replace('$'+X, envdic[X]) 
     envdic[key] = value
       
  f.close()




def script_directory(script_file):
## return the absolute path directory of the given script file.
## >> example <<
## script_file         = 'prot_sch_conbars_tsv.cgi'
## abspath_script_file = '/var/www/html/homcos/cgi-bin/prot_sch_conbars_tsv.cgi'
## script_dir          = '/var/www/html/homcos/cgi-bin'
  abspath_script_file = os.path.abspath(script_file)
  dirs = abspath_script_file.split('/')
  script_dir = ''
  for i in range(0,len(dirs)-1):
    script_dir += '/' + dirs[i]
  return(script_dir)




def get_cgi_values(optdic,argv):
  #[key1]=[value1]&[key2]=[value2]&[key3]=[value3]
  # or
  #[key1]=[value1] [key2]=[value2] [key3]=[value3]
  #
  # **CAUTION**. If key is duplicated, cgi.FieldStorage() stops by the error.
  #   example of errors:  precalc_blt=T&precalc_blt=T 
  #
  # return((Noption,option_type))
  #   option_type is : 'cgi' or 'local'
  # 

  Noption = 0 
  form = cgi.FieldStorage()

  #print "<PRE>#Nkey %d</PRE>"%len(form.keys())
  option_type = 'cgi'

  for key in (form.keys()):
    if (key != ''):
      option_type = 'cgi'
      if (key.startswith('QUPLOAD_PDB')) or (key == 'QLIG_UPLOADFILE') or (key == 'jsme_mol_string') or (key == 'QLIG_JSME_MOL_STRING'):
        optdic[key] = form[key].value
      else:
        optdic[key] = form[key].value.replace(' ','')
      Noption += 1
    #print "<PRE>key '%s' value '%s'</PRE>"%(key,optdic[key]) 
  
 
  ## added by T.Kawabata (2016/12/04) to prevent "directory traversal" and "cross site scripting".##
  ## modified for exceptions for "upload_*" variables (2017/02/08) and QLIG_UPLOADFILE (2017/06/26)##
  ## modified for exceptions for "jsme_mol_string" and 'QLIG_JSME_MOL_STRING' (2017/09/23)##
  ## modified for exceptions for "jsme_mol_string" and 'QLIG_JSME_MOL_STRING' and 'QLIG_SMILES' (2017/09/27)##
  ## -- from here -- ##
  for key in (optdic.keys()):
    val = optdic[key]
    #print "#key %s val %s"%(key,val)
    if (key == 'jsme_mol_string') or (key == 'QLIG_JSME_MOL_STRING') or (key == 'QLIG_SMILES'):
      if (val.find('<')>=0) or (val.find('>')>=0) or (val.find('..')>=0) or (val.find('\'')>=0) or (val.find('\"')>=0) or (val.find('"')>=0) or (val.find("'")>=0) or (val.find('*')>=0) or (val.find('%')>=0) or (val.find(';')>=0) or (re.match(r'passwd',val)) or (len(val)>3000):
        error_exit("IMPROPER INPUT OPTION !! (%s=%s)"%(key,optdic[key]))
    elif (key == 'uploadqsites'):
      pass
      #if (val.find('<')>=0) or (val.find('>')>=0) or (val.find('..')>=0) or (val.find('\'')>=0) or (val.find('\"')>=0) or (val.find('"')>=0) or (val.find("'")>=0) or (val.find('*')>=0) or (val.find('%')>=0) or (val.find(';')>=0) or (re.match(r'passwd',val)) or (len(val)>3000):
      #  error_exit("IMPROPER INPUT OPTION !! (%s=%s)"%(key,optdic[key]))

    elif (key.startswith('QUPLOAD_PDB')==0) and (key != 'QLIG_UPLOADFILE'):
      if (val.find('<')>=0) or (val.find('>')>=0) or (val.find('..')>=0) or (val.find('\'')>=0) or (val.find('\"')>=0) or (val.find('/')>=0) or (val.find('"')>=0) or (val.find("'")>=0) or (val.find('+')>=0) or  (val.find('*')>=0) or (val.find('/')>=0) or (val.find('%')>=0) or (val.find(';')>=0) or (val.find('\\')>=0) or (re.match(r'passwd',val)) or (len(val)>3000):
        error_exit("IMPROPER INPUT OPTION !! (%s=%s)"%(key,optdic[key]))



  ## -- to here -- ##



  now = datetime.now()
  optdic['START_DATE'] = now.strftime("%Y_%m_%d_%H_%M_%S")

  if (Noption<2):
    #optdic['COMMAND'] = argv[0]
    for i in range(1,len(argv)):
      #optdic['COMMAND'] += ' '
      #optdic['COMMAND'] += argv[i]
      if (argv[i].find('=')>=0): 
        option_type = 'local'
        [head,val] = argv[i].split('=')
        val = val.replace(' ','')
        optdic[head] = val
        Noption += 1

  return((Noption,option_type))




def read_argv_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

def read_cgi_bin_option(cgiform,opt_dic):
  N = 0
  for key in (cgiform.keys()):
    opt_dic[key] = cgiform[key].value
    N += 1



def write_logfile(ologfile,comment_list=[],comment_str=''):
 #print "#write_logfile(ologfile,comment_list=[],comment_str=''):"
 if (os.path.isfile(ologfile)) and (os.access(ologfile,os.W_OK)==0):
   #error_exit("#LOGFILE '%s' CAN NOT BE WRITTEN."%(ologfile))
   return(0)
 of = open(ologfile,"a")
 now     = datetime.now()
 datestr = now.strftime("%Y/%m/%d %H:%M:%S")
 PID  = os.getpid()
 of.write(">%s %s %d %s\n"%(datestr,os.getenv("REMOTE_ADDR"),PID,comment_str))
 for content in (comment_list):
   lines = content.split('\n')
   for line in (lines):
     of.write("#%s\n"%(line))
 of.close()
 return(1)


def message(message):
  print "Content-type: text/html\n";
  print "<HTML><BODY><PRE>"
  print "#%s"%(message)
  print "</PRE></BODY></HTML>"

def error_exit(message):
  print "Content-type: text/html\n";
  print "<HTML><BODY><PRE>"
  print "#ERROR:%s"%(cgi.escape(message))
  print "</PRE></BODY></HTML>"
  sys.exit(1)




def _main():

  if (len(sys.argv)<2):
    print "python webfunc.py"
    print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
    sys.exit(1)

  #write_logfile("temp.log",comment_str=sys.argv[1])
  CONF = {} 
  read_environment_file(sys.argv[1],CONF)
  for key in (CONF.keys()):
    print "'%s' '%s'"%(key,CONF[key])


if __name__ == '__main__': _main()
