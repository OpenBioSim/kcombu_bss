#
# <kcombu.py>
#
# subrourintines for kcombu 

import sys
import os
from datetime import datetime

LastModDate = 'Jan 10, 2012'

def read_kcombu_search_result(ifname,hitmol_list):
#>> FILE_FORMAT_EXAMPLE (OLD) <<
# [SIMILAR_MOLECULE_LIST]
# #[rank(1)] [tanimoto(2)] [Nmcs(3)] [Natom(4)] [select_dis(5)]]
# #[molform(6)] [molname(7)] [aligned_atoms(8)]
# #[rmsd(9)] [rmsd_frm_init(10)] [property(11)]
# #                                           [element_for_query]:NCCCOCCOCCCOPOOO
# #                                              [ring_for_query]:aa a a  aa      
# #                                       [atom_number_for_query]:         1111111
# #                                                              :1234567890123456
# #[rk][molname] [tani][Nmcs][Natm] [sel_dis] [molform]          :----------------
# 1    PLP      1.000000  16  16    0.000 C8_N_O6_P               NCCCOCCOCCCOPOOO
# 2    PXP      1.000000  16  16    0.000 C8_N_O6_P               NCCCOCCOCCCOPOOO
# 3    PMP      0.882353  15  16    2.000 C8_N2_O5_P              NCCCOCC-CCCOPOOO
# 4    X04      0.833333  15  17    5.000 C8_N2_O6_P              NCCCOCCOC-COPOOO
# 5    MPL      0.833333  15  17    5.000 C9_N_O6_P               -CCCOCCOCCCOPOOO
# 6    NOP      0.833333  15  17    5.000 C8_N_O7_P               -CCCOCCOCCCOPOOO
#:
# 55   ILP      0.600000  15  24    6.000 C14_N2_O7_P             NCCCOCC-CCCOPOOO
# 56   KE4      0.600000  15  24    6.000 C13_N3_O7_P             NCCCOCC-CCCOPOOO
# [ATOM_MATCHING]
# 1    |0 0|1 1|2 2|3 3|4 4|5 5|6 6|7 7|8 8|9 9|10 10|11 11|12 12|13 13|14 14|15 15
# 2    |0 0|1 1|2 2|3 3|4 4|5 5|6 6|7 7|8 8|9 9|10 10|11 11|12 12|13 13|14 14|15 15
# 3    |0 0|1 1|2 2|3 3|4 4|5 5|6 6|8 8|9 9|10 10|11 11|12 12|13 13|14 14|15 15
# 4    |0 9|1 10|2 11|3 12|4 13|5 14|6 15|7 16|8 6|10 5|11 4|12 0|13 1|14 2|15 3
# 5    |1 10|2 11|3 12|4 13|5 14|6 15|7 16|8 6|9 9|10 5|11 4|12 0|13 1|14 2|15 3
#:

#>> FILE_FORMAT_EXAMPLE <<
# [SIMILAR_MOLECULE_LIST]
# #[rank(1)] [molecular_num(2)] [file_num(3)] [file_offset(4)] [molecular_name(5)]
# #[similarity MCS(6)] [similarity ONEATOM(7)] [similarity ATOMPAIR(8)]
# #[Natompair(9)] [Nheavyatom(10)] [select_dis(11)] [molecular formula(12)]
# #                                                                             [element_for_query]:OCCNCCCCOCCCCCCNCCCCCC
# #                                                                                [ring_for_query]:         aaaaaaaaaaaaa
# #                                                                         [atom_number_for_query]:         1111111111222
# #                                                                                                :1234567890123456789012
# #[rk][num][name]                       [sMCS][sONE][sPAIR][Npair][Nhvyatm][seldis][molform]      :----------------------
# 1    1162122 19 24003076  02590432-01  0.870 0.955 0.940  20  21   35.000 C18_H24_N_O2            OCCNCCCCOCCCCC--CCCCCC
# 2    1162123 19 24012778  02590433-01  0.833 0.913 0.904  20  22   47.000 C19_H26_N_O2            OCCNCCCCOCCCCC--CCCCCC
# 3    1162121 19 23993396  02590431-01  0.833 0.913 0.909  20  22   42.000 C19_H26_N_O2            OCCNCCCCOCCCCC--CCCCCC
# 4    2795962 46 77720418  00018120-01  0.818 0.818 0.831  18  18   29.000 C14_H21_N2_O2           OCCNCCCCOCCCCCCNCC----
# 5    876651  14 224437020 00997701-01  0.792 0.870 0.857  19  21   39.000 C17_H18_N2_O2           OC-NCCCCOCCCCC--CCCCCC
# 6    855782  14 68348542  00976832-01  0.760 0.913 0.851  19  22   46.000 C18_H20_N2_O2           OC-NCCCCOCCCCC--CCCCCC
# 7    1162120 19 23978495  02590430-01  0.760 0.913 0.894  19  22   43.000 C19_H26_N_O2            OCCN-CCCOCCCCC--CCCCCC
# 8    1632085 26 419193924 00694344-01  0.750 0.750 0.800  18  20   33.000 C16_H25_N2_O2           OCCNCCCCOCCCCCCNCC----
# 9    1079776 17 421920618 03273831-01  0.750 0.909 0.864  18  20   47.000 C17_H19_N_O2            O--NCCCCOCCCCC--CCCCCC
# 
  if os.access(ifname,os.R_OK):
    f = open(ifname)
  else:
    print "#WARNING:Can't open '%s'"%(ifname)
    return(0)
  state = ''
  for line in f:  
    line = line.rstrip('\n')
    #print "%s:%s"%(line,state)
    if (line.startswith('[SIMILAR_MOLECULE_LIST]')==1):
      state = 'list'
    elif (line.startswith('[ATOM_MATCHING]')==1):
      state = 'match'
    elif (line.startswith('#')==0) and (len(line)>10):
      if (state=='list'):
        field = line.split()
        dic = {}
        dic['rank']         = int(field[0])
        dic['num_mol']      = int(field[1])
        dic['num_file']     = int(field[2])
        dic['file_offset']  = int(field[3])
        dic['molname']      = field[4]
        dic['sim_mcs']      = float(field[5])
        dic['sim_one']      = float(field[6])
        dic['sim_pair']     = float(field[7])
        dic['Npair']        = int(field[8])
        dic['Nhvyatm']      = int(field[9])
        dic['seldis']       = float(field[10])
        dic['molform']      = field[11]
        hitmol_list.append(dic)
      if (state=='match'):
        field = line.split('|')
        rank = int(field[0])
        hitmol_list[rank-1]['anumA'] = []
        hitmol_list[rank-1]['anumB'] = []
        for i in range(1,len(field)):
          [a,b] = field[i].split()
          hitmol_list[rank-1]['anumA'].append(int(a))
          hitmol_list[rank-1]['anumB'].append(int(b))


  return(1)

def read_atom_match_file(ifname,matchlist,dat):
# matchlist[i]['A'][j] : atomnumber for molecule A for j-th pair in i-th matchlist.
# matchlist[i]['B'][j] : atomnumber for molecule B for j-th pair in i-th matchlist.

#  print "#read_atom_match_file('%s'):"%(ifname)

# #>> FILE FORMAT EXAMPLE <<
# #COMP_TIME  0.004088 seconds
# #MoleculeA SupLIG/SIA_1mweA
# #MoleculeB SupLIG/G39_2qwhA
# #NatomA 21
# #NatomB 20
# #NheavyatomA 21
# #NheavyatomB 20
# #TotalNatompair 112
# #[numA] [num_in_fileA] [atomnameA] --- [numB] [num_in_fileB] [atomnameB]
# >1
# #Npair_atom  13
# #tanimoto    0.464286
# #select_dis  31.000000
# #Ncomponent  1
# #Maxdiff_topodis  0
# 1     3216   C1  --- 1     3194   C1
# 2     3217   C2  --- 4     3197   C2
# 3     3218   C3  --- 13    3206   C7
# 4     3219   C4  --- 12    3205   C6
# 5     3220   C5  --- 7     3200   C5
# 6     3221   C6  --- 6     3199   C4
# 10    3225   C10 --- 9     3202   C10
# 11    3226   C11 --- 11    3204   C11
# 12    3227   N5  --- 8     3201   N5
# 13    3228   O1A --- 2     3195   O1A
# 14    3229   O1B --- 3     3196   O1B
# 16    3231   O4  --- 14    3207   O7
# 21    3236   O10 --- 10    3203   O10
# //
# >2
# >2
# #Npair_atom  13
# #tanimoto    0.464286
# #select_dis  31.000000
# #Ncomponent  1
# #Maxdiff_topodis  0
# 1     3216   C1  --- 1     3194   C1
# 2     3217   C2  --- 4     3197   C2
# 3     3218   C3  --- 13    3206   C7
# 4     3219   C4  --- 12    3205   C6
# 5     3220   C5  --- 7     3200   C5
# 6     3221   C6  --- 6     3199   C4
# 10    3225   C10 --- 9     3202   C10
# 11    3226   C11 --- 11    3204   C11
# 12    3227   N5  --- 8     3201   N5
# 13    3228   O1A --- 3     3196   O1B
# :

  if not os.access(ifname,os.R_OK):
    if  os.access(ifname+'.gz',os.R_OK):
      f = os.popen('gunzip -c ' + ifname + '.gz')
    else:
      print "#WARNING:Can't open filename '%s' and '%s'" %(ifname,ifname+'.gz')
      return(0)
  else:
    if (ifname.endswith('.gz')):
      f = os.popen('gunzip -c ' + ifname)
    else:
      f = open(ifname)
  readmatch = 0

  for line in f:
    line = line.rstrip('\n')
    X = line.split()

    if (line.startswith('/')):
      matchlist.append(match)
      readmatch = 0
    elif (line.startswith('>')):
      match = {}
      match['A'] = []
      match['B'] = []
      match['Ncomponent'] = 0
      readmatch = 1
    elif (line.startswith('#NheavyatomA')):
      dat['NheavyatomA'] = int(X[1])
    elif (line.startswith('#NheavyatomB')):
      dat['NheavyatomB'] = int(X[1])
    elif (line.startswith('#TotalNatompair')):
      dat['TotalNatompair'] = float(X[1])
    elif (line.startswith('#COMP_TIME')):
      dat['COMP_TIME'] = float(X[1])
    elif (line.startswith('#tanimoto')):
      match['tanimoto'] = float(X[1])
    elif (line.startswith('#select_dis')):
      match['selct_dis'] = float(X[1])
    elif (line.startswith('#Npair_atom')):
      match['Npair_atom'] = float(X[1])
    elif (line.startswith('#Ncomponent')):
      match['Ncomponent'] = int(X[1])
    elif (line.startswith('#Maxdiff_topodis')):
      match['Maxdiff_topodis'] = int(X[1])
    elif (line.startswith('#')==0) and (len(line)>10) and (readmatch==1):
      match['A'].append(int(X[0]))
      match['B'].append(int(X[4]))
  f.close()
  if (dat.has_key('NheavyatomA')==0) or (dat.has_key('NheavyatomB')==0):
    return(0)

  return(1)

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



def _main():
  if (len(sys.argv)<2):
    print "python kcombu.py [kcombu_search_result]"
    sys.exit(1)
  ifname = sys.argv[1]
  hitmol_list = []
  read_kcombu_search_result(ifname,hitmol_list)
  print "#Nlist %d"%(len(hitmol_list))
  for m in (hitmol_list):
    sys.stdout.write("%d %s %.3f"%(m['rank'],m['molname'],m['sim_mcs']))
    for i in range(m['Npair']):
      sys.stdout.write(":%d %d"%(m['anumA'][i],m['anumB'][i]))
    sys.stdout.write("\n")   
    
if __name__ == '__main__': _main()
