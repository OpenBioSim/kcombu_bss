#!/usr/bin/env python
import sys
import os
from datetime import datetime
import time
import socket

LastModDate = "Apr 22, 2014"

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


def read_list_file(ifname,list):
  print "#read_list_file('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     list.append(fields[0])
  f.close()

def read_ligand_protein_list_file(ifname,list):
  print "#read_ligand_protein_list_file('%s')"%(ifname)
#>> FILE EXAMPLE 
# #[ligand_file] [protein_file]
# G20_2qwf_1_G_1.pdb 2qwf_1_A_1.pdb
# GNA_2qwe_1_G_1.pdb 2qwe_1_A_1.pdb
# DAN_2qwc_1_G_1.pdb 2qwc_1_A_1.pdb
# 49A_1f8e_1_H_1.pdb 1f8e_1_A_1.pdb
# 9AM_1f8d_1_H_1.pdb 1f8d_1_A_1.pdb
# G39_2qwh_1_E_1.pdb 2qwh_1_A_1.pdb
# SIA_2c4a_1_B_1.pdb 2c4a_1_A_1.pdb

  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR(read_ligand_protein_list_file): Can't open '%s'."%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0) and (len(line)>1):
      fields = line.split()
      ligand  = fields[0]
      protein = fields[1]
      dat = {}
      dat['lig'] = ligand
      dat['pro'] = protein
      list.append(dat)
  f.close()


def write_chimera_command_file_for_ligand(ofname,isdffile,omol2file,formal_charge,method):
#open SupPRO/1bji_1_A_1.pdb
#addh
#addcharge
#write format mol2 0 out.mol2
  print "#write_chimera_command_file_ligand(formal_charge:%d) -->'%s'"%(formal_charge,ofname)
  if (ofname=='-'):
    of = sys.stdout
  else:
    of = open(ofname,'w')
  of.write("open %s\n"%(isdffile))
#  of.write("addh\n")
  of.write("addcharge nonstd :#0 %d method %s\n"%(formal_charge,method))
  of.write("write format mol2 0 %s\n"%(omol2file))
  of.write("stop\n")
  if (ofname!='-'):
    of.close()


def write_chimera_command_file_for_receptor(ofname,ipdbfile,omol2file,onohfile):
#open SupPRO/1bji_1_A_1.pdb
#addh
#addcharge
#write format mol2 0 out.mol2
  print "#write_chimera_command_file_receptor() -->'%s'"%(ofname)
  if (ofname=='-'):
    of = sys.stdout
  else:
    of = open(ofname,'w')
  of.write("open %s\n"%(ipdbfile))
  of.write("delete @/element=H\n")
  of.write("addh\n")
  of.write("addcharge\n")
  of.write("write format mol2 0 %s\n"%(omol2file))
  of.write("delete @/element=H\n")
  of.write("write format pdb 0 %s\n"%(onohfile))
  of.write("stop\n")
  if (ofname!='-'):
    of.close()



def write_sphgen_input_file(ofname,imsfile,osphfile):
#rec.ms
#R
#X
#0.0
#4.0
#1.4
#rec.sph
  print "#write_sphgen_input_file() -->'%s'"%(ofname)
  if (ofname=='-'):
    of = sys.stdout
  else:
    of = open(ofname,'w')
  of.write("%s\n"%(imsfile))
  of.write("R\n")
  of.write("X\n")
  of.write("0.0\n")
  of.write("4.0\n")
  of.write("1.4\n")
  of.write("%s\n"%(osphfile))
  if (ofname!='-'):
    of.close()


def write_box_input_file(ofname,isphfile,oboxfile):
#-- from here (box.in)--
#Y
#5
#selected_spheres.sph
#1
#rec_box.pdb
#--- to here ---
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_box_input_file() -->'%s'"%(ofname)
    of = open(ofname,'w')
  of.write("Y\n")
  of.write("5\n")
  of.write("%s\n"%(isphfile))
  of.write("1\n")
  of.write("%s\n"%(oboxfile))
  if (ofname!='-'):
    of.close()

def write_grid_input_file(ofname,irecmol2file,iboxpdbfile,grid_prefix):
#- from here (grid.in) --
#compute_grids                  yes
#grid_spacing                   0.3
#output_molecule                no
#contact_score                  no
#energy_score                   yes
#energy_cutoff_distance         9999
#atom_model                     a
#attractive_exponent            6
#repulsive_exponent             12
#distance_dielectric            yes
#dielectric_factor              4
#bump_filter                    yes
#bump_overlap                   0.75
#receptor_file                  rec_charged.mol2
#box_file                       rec_box.pdb
#vdw_definition_file            /home/takawaba/tool/dock6.6/parameters/vdw_AMBER_parm99.defn
#score_grid_prefix              grid
#--to here --
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_grid_input_file() -->'%s'"%(ofname)
    of = open(ofname,'w')
  of.write("compute_grids                  yes\n")
  of.write("grid_spacing                   0.3\n")
  of.write("output_molecule                no\n")
  of.write("contact_score                  no\n")
  of.write("energy_score                   yes\n")
  of.write("energy_cutoff_distance         9999\n")
  of.write("atom_model                     a\n")
  of.write("attractive_exponent            6\n")
  of.write("repulsive_exponent             12\n")
  of.write("distance_dielectric            yes\n")
  of.write("dielectric_factor              4\n")
  of.write("bump_filter                    yes\n")
  of.write("bump_overlap                   0.75\n")
  of.write("receptor_file                  %s\n"%(irecmol2file))
  of.write("box_file                       %s\n"%(iboxpdbfile))
  of.write("vdw_definition_file            %s/parameters/vdw_AMBER_parm99.defn\n"%(OPT['dock_dir']))
  of.write("score_grid_prefix              %s\n"%(grid_prefix))
  if (ofname!='-'):
    of.close()

def write_rigid_docking_input_file(ofname,iligmol2file,isphfile,igrid_prefix,oligmol2_prefix):
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_rigid_docking_input_file() -->'%s'"%(ofname)
    of = open(ofname,'w')

  of.write("ligand_atom_file                                             %s\n"%(iligmol2file))
  of.write("limit_max_ligands                                            no\n")
  of.write("skip_molecule                                                no\n")
  of.write("read_mol_solvation                                           no\n")
  of.write("calculate_rmsd                                               no\n")
  of.write("use_database_filter                                          no\n")
  of.write("orient_ligand                                                yes\n")
  of.write("automated_matching                                           yes\n")
  of.write("receptor_site_file                                           %s\n"%(isphfile))
  of.write("max_orientations                                             1000\n")
  of.write("critical_points                                              no\n")
  of.write("chemical_matching                                            no\n")
  of.write("use_ligand_spheres                                           no\n")
  of.write("use_internal_energy                                          yes\n")
  of.write("internal_energy_rep_exp                                      12\n")
  of.write("flexible_ligand                                              no\n")
  of.write("bump_filter                                                  no\n")
  of.write("score_molecules                                              yes\n")
  of.write("contact_score_primary                                        no\n")
  of.write("contact_score_secondary                                      no\n")
  of.write("grid_score_primary                                           yes\n")
  of.write("grid_score_secondary                                         no\n")
  of.write("grid_score_rep_rad_scale                                     1\n")
  of.write("grid_score_vdw_scale                                         1\n")
  of.write("grid_score_es_scale                                          1\n")
  of.write("grid_score_grid_prefix                                       %s\n"%(igrid_prefix))
  of.write("multigrid_score_secondary                                    no\n")
  of.write("dock3.5_score_secondary                                      no\n")
  of.write("continuous_score_secondary                                   no\n")
  of.write("descriptor_score_secondary                                   no\n")
  of.write("gbsa_zou_score_secondary                                     no\n")
  of.write("gbsa_hawkins_score_secondary                                 no\n")
  of.write("SASA_descriptor_score_secondary                              no\n")
  of.write("amber_score_secondary                                        no\n")
  of.write("minimize_ligand                                              yes\n")
  of.write("simplex_max_iterations                                       1000\n")
  of.write("simplex_tors_premin_iterations                               0\n")
  of.write("simplex_max_cycles                                           1\n")
  of.write("simplex_score_converge                                       0.1\n")
  of.write("simplex_cycle_converge                                       1.0\n")
  of.write("simplex_trans_step                                           1.0\n")
  of.write("simplex_rot_step                                             0.1\n")
  of.write("simplex_tors_step                                            10.0\n")
  of.write("simplex_random_seed                                          0\n")
  of.write("simplex_restraint_min                                        no\n")
  of.write("atom_model                                                   all\n")
  of.write("vdw_defn_file                                                %s/parameters/vdw_AMBER_parm99.defn\n"%(OPT['dock_dir']))
  of.write("flex_defn_file                                               %s/parameters/flex.defn\n"%(OPT['dock_dir']))
  of.write("flex_drive_file                                              %s/parameters/flex_drive.tbl\n"%(OPT['dock_dir']))
  of.write("ligand_outfile_prefix                                        %s\n"%(oligmol2_prefix))
  of.write("write_orientations                                           no\n")
  of.write("num_scored_conformers                                        1\n")
  of.write("rank_ligands                                                 no\n")

  if (ofname!='-'):
    of.close()





def write_flexible_docking_input_file(ofname,iligmol2file,isphfile,igrid_prefix,oligmol2_prefix,max_orientations='500',pruning_max_orients='100',score_cutoff='25.0'):
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_flexible_docking_input_file() -->'%s'"%(ofname)
    of = open(ofname,'w')
  of.write("ligand_atom_file                                             %s\n"%(iligmol2file))
  of.write("limit_max_ligands                                            no\n")
  of.write("skip_molecule                                                no\n")
  of.write("read_mol_solvation                                           no\n")
  of.write("calculate_rmsd                                               no\n")
  of.write("use_database_filter                                          no\n")
  of.write("orient_ligand                                                yes\n")
  of.write("automated_matching                                           yes\n")
  of.write("receptor_site_file                                           %s\n"%(isphfile))
  #of.write("max_orientations                                             500\n")
  of.write("max_orientations                                             %s\n"%(max_orientations))
  of.write("critical_points                                              no\n")
  of.write("chemical_matching                                            no\n")
  of.write("use_ligand_spheres                                           no\n")
  of.write("use_internal_energy                                          yes\n")
  of.write("internal_energy_rep_exp                                      12\n")
  of.write("flexible_ligand                                              yes\n")
  of.write("user_specified_anchor                                        no\n")
  of.write("limit_max_anchors                                            no\n")
  of.write("min_anchor_size                                              40\n")
  of.write("pruning_use_clustering                                       yes\n")
#  of.write("pruning_max_orients                                          100\n")
  of.write("pruning_max_orients                                          %s\n"%(pruning_max_orients))
  of.write("pruning_clustering_cutoff                                    100\n")
## >> The default value for 'pruning_conformer_score_cutoff' is 25.0.                           ##
##    However, the default sometimes fail to "complete growth'. A larger value (~100) is better ##
#  of.write("pruning_conformer_score_cutoff                               25.0\n")
  of.write("pruning_conformer_score_cutoff                               %s\n"%(score_cutoff))
  of.write("use_clash_overlap                                            no\n")
  of.write("write_growth_tree                                            no\n")
  of.write("bump_filter                                                  no\n")
  of.write("score_molecules                                              yes\n")
  of.write("contact_score_primary                                        no\n")
  of.write("contact_score_secondary                                      no\n")
  of.write("grid_score_primary                                           yes\n")
  of.write("grid_score_secondary                                         no\n")
  of.write("grid_score_rep_rad_scale                                     1\n")
  of.write("grid_score_vdw_scale                                         1\n")
  of.write("grid_score_es_scale                                          1\n")
  of.write("grid_score_grid_prefix                                       %s\n"%(igrid_prefix))
  of.write("multigrid_score_secondary                                    no\n")
  of.write("dock3.5_score_secondary                                      no\n")
  of.write("continuous_score_secondary                                   no\n")
  of.write("descriptor_score_secondary                                   no\n")
  of.write("gbsa_zou_score_secondary                                     no\n")
  of.write("gbsa_hawkins_score_secondary                                 no\n")
  of.write("SASA_descriptor_score_secondary                              no\n")
  of.write("amber_score_secondary                                        no\n")
  of.write("minimize_ligand                                              yes\n")
  of.write("minimize_anchor                                              yes\n")
  of.write("minimize_flexible_growth                                     yes\n")
  of.write("use_advanced_simplex_parameters                              no\n")
  of.write("simplex_max_cycles                                           1\n")
  of.write("simplex_score_converge                                       0.1\n")
  of.write("simplex_cycle_converge                                       1.0\n")
  of.write("simplex_trans_step                                           1.0\n")
  of.write("simplex_rot_step                                             0.1\n")
  of.write("simplex_tors_step                                            10.0\n")
  of.write("simplex_anchor_max_iterations                                500\n")
  of.write("simplex_grow_max_iterations                                  500\n")
  of.write("simplex_grow_tors_premin_iterations                          0\n")
  of.write("simplex_random_seed                                          0\n")
  of.write("simplex_restraint_min                                        no\n")
  of.write("atom_model                                                   all\n")
  of.write("vdw_defn_file                                                %s/parameters/vdw_AMBER_parm99.defn\n"%(OPT['dock_dir']))
  of.write("flex_defn_file                                               %s/parameters/flex.defn\n"%(OPT['dock_dir']))
  of.write("flex_drive_file                                              %s/parameters/flex_drive.tbl\n"%(OPT['dock_dir']))
  of.write("ligand_outfile_prefix                                        %s\n"%(oligmol2_prefix))
  of.write("write_orientations                                           no\n")
  of.write("num_scored_conformers                                        1\n")
  of.write("rank_ligands                                                 no\n")

  if (ofname!='-'):
    of.close()





def read_formal_charge_from_CIF_file(iciffile):
# data_CT6
# #
# _chem_comp.id                                    CT6
# _chem_comp.name                                  "(5Z)-5-(3-BROMOCYCLOHEXA-2,5-DIEN-1-YLIDENE)-N-(PYRIDIN-4-YLMETHYL)-1,5-DIHYDROPYRAZOLO[1,5-A]PYRIMIDIN-7-AMINE"
# _chem_comp.type                                  NON-POLYMER
# _chem_comp.pdbx_type                             HETAIN
# _chem_comp.formula                               "C18 H15 Br N5"
# _chem_comp.mon_nstd_parent_comp_id               ?
# _chem_comp.pdbx_synonyms                         ?
# _chem_comp.pdbx_formal_charge                    1
# _chem_comp.pdbx_initial_date                     2005-11-08
# _chem_comp.pdbx_modified_date                    2011-06-04
# _chem_comp.pdbx_ambiguous_flag                   N
  if (os.access(iciffile,os.R_OK)==0):
    return(1000)
  f = open(iciffile)
  for line in f:
    if (line.startswith('_chem_comp.pdbx_formal_charge')):
      line = line.rstrip('\n')
      (key,value) = line.split()
      return(int(value))
  f.close()
  return(1000)


def read_mol2_file_to_check_net_charge(imol2file):

#@<TRIPOS>ATOM
#          1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#      1 C1         -3.3340    0.1440   -1.1920 C.2       1 UNK    0.6361
#      2 O1         -4.4350   -0.1550   -0.7750 O.2       1 UNK   -0.5410
#      3 O2         -3.1920    0.6070   -2.4490 O.3       1 UNK   -0.5811
#      4 C2         -2.1530    0.0010   -0.3220 C.2       1 UNK    0.0391
#      5 C3         -2.3650   -0.4580    0.8960 C.2       1 UNK   -0.1502
#:
#     36 H16         2.9460    1.3680   -2.6450 H         1 UNK    0.0437
#     37 H17         3.6540   -0.2210   -2.2710 H         1 UNK    0.0437
#     38 H18         3.7690    0.3020   -4.5000 H         1 UNK    0.4110
#@<TRIPOS>BOND
#     1    1    2 2
#     2    1    3 1
  netcharge = 0.0

  if (os.access(imol2file,os.R_OK)==0):
    print "#WARNING: Can't open '%s'."%(imol2file)
    return("file_not_found")
  f = open(imol2file)

  atomline = 0
  for line in f:
    if (line.startswith('@<TRIPOS>BOND')):
      atomline = 0 
  
    if (atomline==1) and (len(line)>40):
      line = line.rstrip('\n')
      field = line.split()
      charge = float(field[8])
      netcharge += charge       

    if (line.startswith('@<TRIPOS>ATOM')):
      atomline = 1 

  f.close()
  return("%f"%(netcharge))


def add_comp_time_to_mol2(ifname,comp_time):
  print "#add_comp_time_to_mol2('%s',comp_time %s)"%(ifname,comp_time)
  lines = []
  if (os.access(ifname,os.R_OK)==0):
    print "#WARNING(add_comp_time_to_mol2):Can't open '%s'."%(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    lines.append(line)
  f.close()
  of = open(ifname,'w')
# ########## Name:                G20
# ##########    Grid Score:          -57.195641
# ##########      Grid_vdw:          -39.115467
# ##########       Grid_es:          -18.080172
# ##########    Int_energy:           14.405195

  now = datetime.now()
  date_string = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(len(lines)):
    line = lines[i]
    of.write("%s"%(line))
    if (line.startswith('#')) and (i<(len(lines)-1)) and (lines[i+1].startswith('#')==0):
      of.write("##########    COMP_TIME:           %f seconds (added by T.Kawabata)\n"%(comp_time)) 
      of.write("##########    DATE     :           %s (added by T.Kawabata)\n"%(date_string)) 
      of.write("##########    COMMAND  :           %s (added by T.Kawabata)\n"%(OPT['COMMAND']))
      of.write("##########    HOST     :           %s (added by T.Kawabata)\n"%(socket.gethostname())) 
  of.close()
  return(1)



############
### MAIN ###
############

OPT = {}
OPT['L']   = 'ligand.list'
OPT['idref'] = 'SupLIG'

OPT['M'] = 'L'
OPT['A'] = 'F'
OPT['E'] = 'A'
OPT['S'] = 'F'
OPT['div'] = '0/1'
OPT['dsdf']  = 'SDF_3D'
OPT['dcif']  = 'CIF_3D'
OPT['dmol2'] = 'MOL2_CHIMERA_AM1'
OPT['rms'] = 'T'
OPT['chsp3'] = 'F'
OPT['stprng'] = 'T'
OPT['gra'] = 'T'
OPT['mcs'] = 'T'
OPT['ligch'] = 'am1'
OPT['ochg'] = ''
OPT['drpro'] = 'SupPRO'
OPT['drmol2'] = 'SupPRO_MOL2'
OPT['drnoh']  = 'SupPRO_NOH'
OPT['drms']  = 'SupPRO_MS'
OPT['dsph']  = 'SupPRO_SPH'
OPT['drgrid'] = 'SupPRO_GRID'
OPT['dms'] = 'dms'

OPT['sphgen'] = 'sphgen'
OPT['sphsele'] = 'sphere_selector'
OPT['drefmol2'] = 'MOL2_SupLIG'
OPT['Dsphsele'] = '10.0'
OPT['showbox'] = 'showbox'
OPT['grid'] = 'grid'
OPT['dock'] = 'dock6'
OPT['ddock'] = 'DOCK'
OPT['rgdflx'] = 'F'
OPT['subdir'] = 'T'
OPT['avex'] = 'F'
OPT['self'] = 'F'
OPT['score_cutoff'] = '100.0'
OPT['dock_dir'] = '/home/takawaba/tool/dock6.6'
OPT['max_orientations']='500'
OPT['pruning_max_orients']='100'
OPT['offchg'] = '0'
OPT['c'] = '0'
OPT['C'] = '1'
OPT['n'] = '0'
OPT['N'] = '1'



PID = os.getpid()
HOSTNAME = socket.gethostname()


if (len(sys.argv)<2):
  print "prep_dock.py <options>"
  print " for preparing UCSF DOCK using UCSF Chimera."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -M      : MODE. 'L'igand preparation, 'R'eceptor preparation [%s]"%(OPT['M'])
  print "         :       'P'robe generation 'S'election of probe 'G'rid generation 'C'heck generated files."
  print "         :       'D'ocking." 
  print " -L      : list of library molecules[%s]"%(OPT['L'])
  print " -div    : Job division (bunshi/bunbo) [%s]"%(OPT['div'])
  print " -c      : core number                       [%s]"%(OPT['c'])
  print " -C      : number of total cores for one node[%s]"%(OPT['C'])
  print " -n      : node number           [%s]"%(OPT['n'])
  print " -N      : number of total nodes [%s]"%(OPT['N'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print "<options for MODE='L'igand preparation>"
  print " -dsdf  : input dir. for 3D SDF molecules [%s]"%(OPT['dsdf'])
  print " -dcif  : input dir. for 3D SDF molecules for formal charge[%s]"%(OPT['dcif'])
  print " -dmol2 : output directory for ligand 3D MOL2 molecules[%s]"%(OPT['dmol2'])
  print " -ligch  : add charge method for ligand 'am1' or 'gas' [%s]"%(OPT['ligch']) 
  print " -ochg   : output file to  check netcharge in the dmol2 [%s]"%(OPT['ochg'])
  print " -offchg : offset for formal charge [%s]"%(OPT['offchg']) 
  print "<options for MODE='R'ecptor preparation>"
  print " -drpro : input  dir. for receptor protein file in PDB [%s]"%(OPT['drpro'])
  print " -drmol2: output dir. for receptor protein in mol2 [%s]"%(OPT['drmol2'])
  print " -drnoh : output dir. for receptor protein in PDB  [%s]"%(OPT['drnoh'])
  print " -drms  : output dir. for receptor molecular surface  [%s]"%(OPT['drms'])
  print " -dms    : program for surface calculation [%s]"%(OPT['dms'])
  print "<options for MODE='P'robe genration>"
  print " -drms     : input dir.  for receptor molecular surface  [%s]"%(OPT['drms'])
  print " -dsph     : output dir. for sphere probes    [%s]"%(OPT['dsph'])
  print " -sphgen   : program for sphere generation [%s]"%(OPT['sphgen'])
  print "<options for MODE='S'election of the probes>"
  print " -dsph     : input dir. for sphere probes    [%s]"%(OPT['dsph'])
  print " -sphsele  : program for sphere selector [%s]"%(OPT['sphsele'])
  print " -drefmol2 : reference ligand mol2 file [%s]"%(OPT['drefmol2'])
  print " -Dsphsele : threshold distance from reference [%s]"%(OPT['Dsphsele'])
  print "<options for MODE='G'rid generation>"
  print " -dsph     : input dir. for sphere probes    [%s]"%(OPT['dsph'])
  print " -drmol2   : output dir. for receptor protein in mol2 [%s]"%(OPT['drmol2'])
  print " -drgrid   : output dir. for grid file [%s]"%(OPT['drgrid'])
  print " -showbox  : program for showbox [%s]"%(OPT['showbox'])
  print " -grid     : program for grid [%s]"%(OPT['grid'])
  print "<options for MODE='D'ocking>"
  print " -dmol2   : input dir. for ligand 3D MOL2 molecules[%s]"%(OPT['dmol2'])
  print " -drmol2  : input dir. for receptor protein in mol2 [%s]"%(OPT['drmol2'])
  print " -dsph    : input dir. for sphere probes    [%s]"%(OPT['dsph'])
  print " -drgrid  : input dir. for grid file [%s]"%(OPT['drgrid'])
  print " -ddock   : output dir. for docking    [%s]"%(OPT['ddock'])
  print " -rgdflx  : 'R'igid or 'F'lexible [%s]"%(OPT['rgdflx'])
  print " -dock    : program for dock [%s]"%(OPT['dock'])
  print " -subdir  : use target ligand name as output subdirectory. (T or F)[%s]"%(OPT['subdir'])
  print " -avex    : avoid existing model file (T or F)[%s]"%(OPT['avex'])
  print " -self    : do only self modeling (T or F) [%s]"%(OPT['self'])
  print " -score_cutoff     : value for 'pruning_conformer_score_cutoff' [%s]"%(OPT['score_cutoff'])
  print " -max_orientations : value for 'max_orientations' [%s]"%(OPT['max_orientations'])
  print " -pruning_max_orients : value for 'pruning_max_orients' [%s]"%(OPT['pruning_max_orients'])
  print " -dock_dir: directory for the dock program [%s]"%(OPT['dock_dir'])
  sys.exit(1)

### [1] read option ###
read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)
if (OPT['C']!='1') or (OPT['N']!='1'):
  bunshi = int(OPT['n'])*int(OPT['C']) + int(OPT['c'])
  bunbo  = int(OPT['N'])*int(OPT['C'])


### [2] read liglist (-L) ###
ligpro_list = []
read_ligand_protein_list_file(OPT['L'],ligpro_list)

Nlig  = len(ligpro_list)



####################################
### MODE 'L': Ligand preparation ###
####################################
if (OPT['M']=='L'):
  Nlig  = len(ligpro_list)
  Nsta = Nlig*bunshi/bunbo
  Nend = (Nlig*(bunshi+1))/bunbo

  ocomfile = '%d.com'%(PID)
  FORMAL_CHARGE = {}
  NET_CHARGE    = {}
  
  #for ligpro in (ligpro_list):
  for i in range(Nsta,Nend):
    ligpro = ligpro_list[i]
  # G39_2qwh_1_E_1.pdb 2qwh_1_A_1.pdb
  # SIA_2c4a_1_B_1.pdb 2c4a_1_A_1.pdb
    print ligpro
    ligT = ligpro['lig']
    (ligname3T,pdbT,assembly_idT,asym_idT,operT) = ligT.split('_')
    isdffile = OPT['dsdf']  + '/' + ligname3T + '.sdf'
    iciffile = OPT['dcif']  + '/' + ligname3T + '.cif'
    if (os.access(isdffile,os.R_OK) and os.access(iciffile,os.R_OK)):
      omol2file = OPT['dmol2'] + '/' + ligname3T + '.mol2'
      print isdffile
      print iciffile
      FORMAL_CHARGE[ligT] = read_formal_charge_from_CIF_file(iciffile)
      print "#%s formal_charge %d"%(ligname3T,FORMAL_CHARGE[ligT])
      if (OPT['A']=='T'):
        write_chimera_command_file_for_ligand(ocomfile,isdffile,omol2file,int(FORMAL_CHARGE[ligT])+int(OPT['offchg']),OPT['ligch'])
      else:
        write_chimera_command_file_for_ligand('-',isdffile,omol2file,int(FORMAL_CHARGE[ligT])+int(OPT['offchg']),OPT['ligch'])
      command = "chimera --nogui %s"%(ocomfile)
      print "#command %s"%(command)
      if (OPT['A']=='T'):
        os.system(command)
  
      if (OPT['ochg']!=''):
        NET_CHARGE[ligT] = read_mol2_file_to_check_net_charge(omol2file)
  
  
  
  if (OPT['ochg']!=''):
    of = open(OPT['ochg'],'w')
    print "#write_netcharges --> '%s'"%(OPT['ochg']) 
    of.write("#COMMAND %s\n"%(OPT['COMMAND']))
    of.write("#DATE    %s\n"%(OPT['START_DATE']))
    of.write("#[ligname3T] [formal_charge] [net_charge]\n")
    for ligpro in (ligpro_list):
      ligT = ligpro['lig']
      of.write("%s  %s  %s\n"%(ligT,FORMAL_CHARGE[ligT],NET_CHARGE[ligT]))
    of.close()

######################################
### MODE 'R': Receptor preparation ###
######################################
if (OPT['M']=='R'):
  ocomfile = '%d.com'%(PID)
  for ligpro in (ligpro_list):
  # G39_2qwh_1_E_1.pdb 2qwh_1_A_1.pdb
  # SIA_2c4a_1_B_1.pdb 2c4a_1_A_1.pdb
    print ligpro
    pro_pdbfile = ligpro['pro']
#    (pdb,assembly_id,asym_id,oper) = pro.split('_')
    (prohead,tail) = pro_pdbfile.split('.')
    ipdbfile  = OPT['drpro'] + '/' + prohead + '.pdb'
    omol2file = OPT['drmol2'] + '/' + prohead + '.mol2'
    onohfile  = OPT['drnoh'] + '/' + prohead + '.pdb'
    omsfile   = OPT['drms'] + '/' + prohead + '.ms'
    if (OPT['A']=='T'):
      write_chimera_command_file_for_receptor(ocomfile,ipdbfile,omol2file,onohfile)
    else:
      write_chimera_command_file_for_receptor('-',ipdbfile,omol2file,onohfile)
    #command     = "chimera --nogui %s"%(ocomfile)
    command_dms = "dms %s -a -n -w 1.4 -v -o %s"%(onohfile,omsfile)
    #print "#command     %s"%(command)
    print "#command_dms %s"%(command_dms)
    if (OPT['A']=='T'):
      #os.system(command) <-- only for 20140417
      os.system(command_dms)

##################################
### MODE 'P': Probe generation ###
##################################
if (OPT['M']=='P'):
  insphfile = 'INSPH'
  Nlig  = len(ligpro_list)
  Nsta = Nlig*bunshi/bunbo
  Nend = (Nlig*(bunshi+1))/bunbo
  for ligpro in (ligpro_list):
  #for i in range(Nsta,Nend):
    #ligpro = ligpro_list[i]
    print ligpro
    pro_pdbfile = ligpro['pro']
    lig_pdbfile  = ligpro['lig']
    (prohead,tail) = pro_pdbfile.split('.')
    (lighead,tail) = lig_pdbfile.split('.')
    imsfile       = OPT['drms'] + '/' + prohead + '.ms'
    osphfile_all  = OPT['dsph'] + '/' + prohead + '.sph.all'
    osphfile      = OPT['dsph'] + '/' + prohead + '.sph'
    command_sphgen  = "%s"%(OPT['sphgen'])
    if (OPT['A']=='T'):
      write_sphgen_input_file(insphfile,imsfile,osphfile_all)
      if (os.path.isfile('OUTSPH')):
        os.system('rm OUTSPH')
      if (os.path.isfile('temp1.ms')):
        os.system('rm temp1.ms')
      if (os.path.isfile('temp3.atc')):
        os.system('rm temp3.atc')
      os.system(command_sphgen)

    else:
      write_sphgen_input_file('-',imsfile,osphfile_all)
      print "#command_sphgen  %s"%(command_sphgen)
  pass


#####################################
### MODE 'S': Selection of Probes ###
#####################################
if (OPT['M']=='S'):
  insphfile = 'INSPH'
  for ligpro in (ligpro_list):
    print ligpro
    pro_pdbfile = ligpro['pro']
    lig_pdbfile  = ligpro['lig']
    (prohead,tail) = pro_pdbfile.split('.')
    (lighead,tail) = lig_pdbfile.split('.')
    osphfile_all  = OPT['dsph'] + '/' + prohead + '.sph.all'
    osphfile      = OPT['dsph'] + '/' + prohead + '.sph'
    refmol2file   = OPT['drefmol2'] + '/' + lighead + '.mol2'
    command_sphsele = "%s %s %s %s"%(OPT['sphsele'],osphfile_all,refmol2file,OPT['Dsphsele'])
    if (OPT['A']=='T'):
      os.system(command_sphsele)
      if (os.path.isfile('selected_spheres.sph')):
        os.system("mv selected_spheres.sph %s"%(osphfile))

    else:
      print "#command_sphsele %s"%(command_sphsele)
      print "#command_mv      mv selected_spheres.sph %s"%(osphfile)
  pass












##################################
### MODE 'G': Grid generation ###
##################################
if (OPT['M']=='G'):
  Nlig  = len(ligpro_list)
  Nsta = Nlig*bunshi/bunbo
  Nend = (Nlig*(bunshi+1))/bunbo
  print "#Nlig %d bunshi %d bunbo %d Nsta %d Nend %d"%(Nlig,bunshi,bunbo,Nsta,Nend)


  oboxinfile  = '%d.boxin'%(PID)
  ogridinfile = '%d.gridin'%(PID)
  for i in range (Nsta,Nend):
    ligpro = ligpro_list[i]
    pro_pdbfile = ligpro['pro']
    (prohead,tail) = pro_pdbfile.split('.')
    isphfile   = OPT['dsph'] + '/' + prohead + '.sph'
    oboxfile   = OPT['dsph'] + '/' + prohead + '.box.pdb'
    imol2file = OPT['drmol2'] + '/' + prohead + '.mol2'
    ogrid_prefix = OPT['drgrid'] + '/' + prohead 
    command_showbox  = "%s < %s"%(OPT['showbox'],oboxinfile)
    command_grid  = "%s -i %s"%(OPT['grid'],ogridinfile)
    if (OPT['A']=='T'):
      write_box_input_file(oboxinfile,isphfile,oboxfile)
      os.system(command_showbox)
      write_grid_input_file(ogridinfile,imol2file,oboxfile,ogrid_prefix)
      os.system(command_grid)
    else: 
      write_box_input_file('-',isphfile,oboxfile)
      print "#command_showbox %s"%(command_showbox)
      write_grid_input_file('-',imol2file,oboxfile,ogrid_prefix)
      print "#command_grid %s"%(command_grid)

######################################
### MODE 'C': Check generated file ###
######################################
if (OPT['M']=='C'):
  receptor_error_dic = {}
  for ligpro in (ligpro_list):
    pro_pdbfile = ligpro['pro']
    lig_pdbfile  = ligpro['lig']
    (prohead,tail) = pro_pdbfile.split('.')
    (lighead,tail) = lig_pdbfile.split('.')
    omol2file = OPT['drmol2'] + '/' + prohead + '.mol2'
    onohfile  = OPT['drnoh'] + '/' + prohead + '.pdb'
    omsfile   = OPT['drms'] + '/' + prohead + '.ms'
    osphfile  = OPT['dsph'] + '/' + prohead + '.sph'
    ogrid_file = OPT['drgrid'] + '/' + prohead  + '.bmp'
    ok_string = '' 
    if (os.path.isfile(omol2file)==0):
      ok_string += ' no_mol2' 
    if (os.path.isfile(onohfile)==0):
      ok_string += ' no_noh' 
    if (os.path.isfile(omsfile)==0):
      ok_string += ' no_ms' 
    if (os.path.isfile(osphfile)==0):
      ok_string += ' no_sph' 
    if (os.path.isfile(ogrid_file)==0):
      ok_string += ' no_grid' 

    (ligname3,pdb,assembly_id,asym_id,oper) = lig_pdbfile.split('_')
    oligmol2file = OPT['dmol2'] + '/' + ligname3 + '.mol2'
    if (os.path.isfile(oligmol2file)==0):
      ok_string += ' no_lig_mol2' 
      print "no:",oligmol2file
 
    if (ok_string != ''):
      receptor_error_dic[lig_pdbfile + ':' + pro_pdbfile]=ok_string


  for key in (receptor_error_dic.keys()):
    (lig_pdbfile,pro_pdbfile) =  key.split(':')
    print "%s %s %s"%(lig_pdbfile,pro_pdbfile,receptor_error_dic[key])






#########################
### MODE 'D': Docking ###
#########################
if (OPT['M']=='D'):
  oinfile = '%s_%d.dockin'%(HOSTNAME,PID)
  Nligpro = len(ligpro_list)

#  Npair = Nligpro * Nligpro
  ### count Npair ###
  Npair = 0
  for t in range(Nligpro):
    ligproT = ligpro_list[t]
    ligT = ligproT['lig']
    (ligname3T,pdbT,assembly_idT,asym_idT,operT) = ligT.split('_')
    for r in range(Nligpro):
      if ((OPT['self']!='T') or (t == r)):
        ligproR = ligpro_list[r]
        ligR = ligproR['lig']
        (ligheadR,tail) = ligR.split('.')
        proR = ligproR['pro']
        (proheadR,tail) = proR.split('.')
        if (OPT['subdir']=='T'):
          odir_dock = "%s/%s"%(OPT['ddock'],ligname3T)
          oligmol2_prefix   = odir_dock + '/' + ligname3T + '_' + ligheadR
        else:
          oligmol2_prefix   = OPT['ddock'] + '/' + ligname3T + '_' + ligheadR

        if (OPT['avex']!='T') or (os.path.isfile(oligmol2_prefix+'.mol2')==0) or (os.path.getsize(oligmol2_prefix+'.mol2')<10):
          Npair += 1



  Nsta = Npair*bunshi/bunbo
  Nend = (Npair*(bunshi+1))/bunbo
  print "#Nlig %d Npair %d bunshi %d bunbo %d Nsta %d Nend %d"%(Nlig,Npair,bunshi,bunbo,Nsta,Nend)




  ### Do docking ###
  npair = 0

  for t in range(Nligpro):
    ligproT = ligpro_list[t]
    ligT = ligproT['lig']
    (ligname3T,pdbT,assembly_idT,asym_idT,operT) = ligT.split('_')
    iligmol2file = OPT['dmol2'] + '/' + ligname3T + '.mol2'
    for r in range(Nligpro):
      if ((OPT['self']!='T') or (t == r)):
        ligproR = ligpro_list[r]
        ligR = ligproR['lig']
        (ligheadR,tail) = ligR.split('.')
        proR = ligproR['pro']
        (proheadR,tail) = proR.split('.')
        irecmol2file   = OPT['drmol2'] + '/' + proheadR + '.mol2'
        isphfile       = OPT['dsph'] + '/'   + proheadR + '.sph'
        igrid_prefix   = OPT['drgrid'] + '/' + proheadR


        if (OPT['subdir']=='T'):
          odir_dock = "%s/%s"%(OPT['ddock'],ligname3T)
          oligmol2_prefix   = odir_dock + '/' + ligname3T + '_' + ligheadR
          if (os.path.isdir(odir_dock)==0):
            print "subdir '%s' does not exist."%(odir_dock)
            mkdir_command = "mkdir %s"%(odir_dock)
            if (OPT['A']=='T'):
              os.system(mkdir_command)
        else:
          #outpdbfile = "%s/%s_%s"%(OPT['odtar'],ligname3T,ligR)
          oligmol2_prefix   = OPT['ddock'] + '/' + ligname3T + '_' + ligheadR


        if (OPT['avex']!='T') or (os.path.isfile(oligmol2_prefix+'.mol2')==0) or (os.path.getsize(oligmol2_prefix+'.mol2')<10):
          if (Nsta<=npair) and (npair<Nend):
            print "#>%s"%(oligmol2_prefix)
            command_dock = "%s -i %s"%(OPT['dock'],oinfile)
            if (OPT['A']=='T'):
              if (OPT['rgdflx']=='F'):
                write_flexible_docking_input_file(oinfile,iligmol2file,isphfile,igrid_prefix,oligmol2_prefix,score_cutoff=OPT['score_cutoff'])
              else:
                write_rigid_docking_input_file(oinfile,iligmol2file,isphfile,igrid_prefix,oligmol2_prefix)
              sta_sec = time.time()
              os.system(command_dock)
              end_sec = time.time()
              mvmol2_command = "mv %s_scored.mol2 %s.mol2"%(oligmol2_prefix,oligmol2_prefix)
              os.system(mvmol2_command)
              add_comp_time_to_mol2(oligmol2_prefix + '.mol2',end_sec-sta_sec)
              cpdock_command = "cp %s_%d.dockin %s.dockin"%(HOSTNAME,PID,oligmol2_prefix)
              os.system(cpdock_command)
            else:
              if (OPT['rgdflx']=='F'):
                write_flexible_docking_input_file('-',iligmol2file,isphfile,igrid_prefix,oligmol2_prefix,score_cutoff=OPT['score_cutoff'])
              else:
                write_rigid_docking_input_file('-',iligmol2file,isphfile,igrid_prefix,oligmol2_prefix)
              print "#command_dock %s"%(command_dock)
          npair += 1
