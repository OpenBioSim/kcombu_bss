/*

 <MOLECULE_from_CCD_CIF.c>

=============================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.
This software is released under
the GNU Lesser General Public License (LGPL), see LICENSE.txt.
=============================================================


 This file is for getting a linked list of "struct ATOM" 
 from the "struct PDB_DATA" (generated from "struct BLOCK").

 This file should be modified considering the definition 
 of "struct ATOM".

 This file is modified for the program "kcombu".

 LastModDate: 2019/11/12.

*/
/* 
 * Example of CCD(Chemical Components Dictionary) CIF file 
-------------------------------------------------------------------
data_ALA
# 
_chem_comp.id                                    ALA 
_chem_comp.name                                  ALANINE 
_chem_comp.type                                  "L-PEPTIDE LINKING" 
_chem_comp.pdbx_type                             ATOMP 
_chem_comp.formula                               "C3 H7 N O2" 
_chem_comp.mon_nstd_parent_comp_id               ? 
_chem_comp.pdbx_synonyms                         ? 
_chem_comp.pdbx_formal_charge                    0 
_chem_comp.pdbx_initial_date                     1999-07-08 
_chem_comp.pdbx_modified_date                    2011-06-04 
_chem_comp.pdbx_ambiguous_flag                   N 
_chem_comp.pdbx_release_status                   REL 
_chem_comp.pdbx_replaced_by                      ? 
_chem_comp.pdbx_replaces                         ? 
_chem_comp.formula_weight                        89.093 
_chem_comp.one_letter_code                       A 
_chem_comp.three_letter_code                     ALA 
_chem_comp.pdbx_model_coordinates_details        ? 
_chem_comp.pdbx_model_coordinates_missing_flag   N 
_chem_comp.pdbx_ideal_coordinates_details        ? 
_chem_comp.pdbx_ideal_coordinates_missing_flag   N 
_chem_comp.pdbx_model_coordinates_db_code        ? 
_chem_comp.pdbx_subcomponent_list                ? 
_chem_comp.pdbx_processing_site                  RCSB 
# 
loop_
_chem_comp_atom.comp_id 
_chem_comp_atom.atom_id 
_chem_comp_atom.alt_atom_id 
_chem_comp_atom.type_symbol 
_chem_comp_atom.charge 
_chem_comp_atom.pdbx_align 
_chem_comp_atom.pdbx_aromatic_flag 
_chem_comp_atom.pdbx_leaving_atom_flag 
_chem_comp_atom.pdbx_stereo_config 
_chem_comp_atom.model_Cartn_x 
_chem_comp_atom.model_Cartn_y 
_chem_comp_atom.model_Cartn_z 
_chem_comp_atom.pdbx_model_Cartn_x_ideal 
_chem_comp_atom.pdbx_model_Cartn_y_ideal 
_chem_comp_atom.pdbx_model_Cartn_z_ideal 
_chem_comp_atom.pdbx_component_atom_id 
_chem_comp_atom.pdbx_component_comp_id 
_chem_comp_atom.pdbx_ordinal 
ALA N   N   N 0 1 N N N 2.281  26.213 12.804 -0.966 0.493  1.500  N   ALA 1  
ALA CA  CA  C 0 1 N N S 1.169  26.942 13.411 0.257  0.418  0.692  CA  ALA 2  
ALA C   C   C 0 1 N N N 1.539  28.344 13.874 -0.094 0.017  -0.716 C   ALA 3  
ALA O   O   O 0 1 N N N 2.709  28.647 14.114 -1.056 -0.682 -0.923 O   ALA 4  
ALA CB  CB  C 0 1 N N N 0.601  26.143 14.574 1.204  -0.620 1.296  CB  ALA 5  
ALA OXT OXT O 0 1 N Y N 0.523  29.194 13.997 0.661  0.439  -1.742 OXT ALA 6  
ALA H   H   H 0 1 N N N 2.033  25.273 12.493 -1.383 -0.425 1.482  H   ALA 7  
ALA H2  HN2 H 0 1 N Y N 3.080  26.184 13.436 -0.676 0.661  2.452  H2  ALA 8  
ALA HA  HA  H 0 1 N N N 0.399  27.067 12.613 0.746  1.392  0.682  HA  ALA 9  
ALA HB1 1HB H 0 1 N N N -0.247 26.699 15.037 1.459  -0.330 2.316  HB1 ALA 10 
ALA HB2 2HB H 0 1 N N N 0.308  25.110 14.270 0.715  -1.594 1.307  HB2 ALA 11 
ALA HB3 3HB H 0 1 N N N 1.384  25.876 15.321 2.113  -0.676 0.697  HB3 ALA 12 
ALA HXT HXT H 0 1 N Y N 0.753  30.069 14.286 0.435  0.182  -2.647 HXT ALA 13 
# 
loop_
_chem_comp_bond.comp_id 
_chem_comp_bond.atom_id_1 
_chem_comp_bond.atom_id_2 
_chem_comp_bond.value_order 
_chem_comp_bond.pdbx_aromatic_flag 
_chem_comp_bond.pdbx_stereo_config 
_chem_comp_bond.pdbx_ordinal 
ALA N   CA  SING N N 1  
ALA N   H   SING N N 2  
ALA N   H2  SING N N 3  
ALA CA  C   SING N N 4  
ALA CA  CB  SING N N 5  
ALA CA  HA  SING N N 6  
ALA C   O   DOUB N N 7  
ALA C   OXT SING N N 8  
ALA CB  HB1 SING N N 9  
ALA CB  HB2 SING N N 10 
ALA CB  HB3 SING N N 11 
ALA OXT HXT SING N N 12 
# 
loop_
_pdbx_chem_comp_descriptor.comp_id 
_pdbx_chem_comp_descriptor.type 
_pdbx_chem_comp_descriptor.program 
_pdbx_chem_comp_descriptor.program_version 
_pdbx_chem_comp_descriptor.descriptor 
ALA SMILES           ACDLabs              10.04 "O=C(O)C(N)C"                                                 
ALA SMILES_CANONICAL CACTVS               3.341 "C[C@H](N)C(O)=O"                                             
ALA SMILES           CACTVS               3.341 "C[CH](N)C(O)=O"                                              
ALA SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "C[C@@H](C(=O)O)N"                                            
ALA SMILES           "OpenEye OEToolkits" 1.5.0 "CC(C(=O)O)N"                                                 
ALA InChI            InChI                1.03  "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1" 
ALA InChIKey         InChI                1.03  QNAYBMKLOCPYGJ-REOHCLBHSA-N                                   
# 
loop_
_pdbx_chem_comp_identifier.comp_id 
_pdbx_chem_comp_identifier.type 
_pdbx_chem_comp_identifier.program 
_pdbx_chem_comp_identifier.program_version 
_pdbx_chem_comp_identifier.identifier 
ALA "SYSTEMATIC NAME" ACDLabs              10.04 L-alanine                    
ALA "SYSTEMATIC NAME" "OpenEye OEToolkits" 1.5.0 "(2S)-2-aminopropanoic acid" 
# 
loop_
_pdbx_chem_comp_audit.comp_id 
_pdbx_chem_comp_audit.action_type 
_pdbx_chem_comp_audit.date 
_pdbx_chem_comp_audit.processing_site 
ALA "Create component"  1999-07-08 RCSB 
ALA "Modify descriptor" 2011-06-04 RCSB 
# 
-------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "molecule.h" 
#include "globalvar.h" 
#include "molprop.h" 
#include "kcombu_io_CIF.h" 

/** FUNCTIONS (GLOBAL) **/
int Translate_Lines_into_CCD_CIF_Molecule();

/** FUNCTIONS (LOCAL) **/
static int find_atom_num_from_atomname();


int Translate_Lines_into_CCD_CIF_Molecule(HeadLineNode,mol)
/*
 for reading chemical components dictionary file (CCD; components.cif ) 
*/ 
  struct LINENODE *HeadLineNode;
  struct MOLECULE *mol;
{
  struct BLOCK BLK;
  struct WORD HeadWord;
  struct TABLE *tb;
  struct ROWLINE *rn;
  struct ATOM *atm;
  struct BOND *bnd;
  int    i;
  char   buff[MAX_LENGTH_LINE];

  printf("#Translate_Lines_into_CCD_CIF_Molecule(HeadLineNode,mol)\n"); fflush(stdout);
  /*** (1) Translate LINEs into WORD and BLOCK  ***/
  Translate_Lines_into_WORDs(HeadLineNode,&HeadWord);
  translate_CIF_WORDs_into_BLOCK_data(&HeadWord,&BLK);
  make_rowline_from_temp_head_word_for_one_row_table(&BLK);
  /* printf("#Ntable %d\n",BLK.Ntable); */

  mol->Natom = 0;  
  mol->Nbond = 0;  
  mol->resi[0] = '\0';
  mol->chain   = ' ';
  mol->Natom = mol->Nheavyatom = 0;
  mol->Nbond = 0;
  mol->conmap.malloced = mol->dismap.malloced = 0;
  mol->Nmatch = 0;
  mol->chiral_flag = '0';
  mol->Npermu = 0;
  mol->permuhead.next = NULL;


  sprintf(mol->name,"%s",mol->core_molname);
  mol->comment[0] = '\0';
 /*
  tb = &(BLK.head_table);
  while (tb->next != NULL){
    tb = tb->next;
    printf("TABLE table_name '%s' Nrow %d Ncolumn %d\n",tb->table_name,tb->Nrow, tb->Ncolumn);
  }
  */


  /*** (2)  get Natom and Nbond and malloc MOLECULE ***/
  tb = TABLE_from_key(BLK.table_hash_tab,"chem_comp_atom");
  if (tb != NULL){ mol->Natom = tb->Nrow; }

  tb = TABLE_from_key(BLK.table_hash_tab,"chem_comp_bond");
  if (tb != NULL){ mol->Nbond = tb->Nrow; }

  if (mol->Natom == 0){
    printf("#ERROR[Translate_Lines_into_CCD_CIF_Molecule()]: No atom is found.\n");
    exit(1);
  }
 /* printf("#Natom %d Nbond %d\n",mol->Natom, mol->Nbond); */
  Malloc_MOLECULE(mol,mol->Natom,mol->Nbond);

  /*** (3) Setup ATOMs from 'chem_comp_atom' table ***/
  tb = TABLE_from_key(BLK.table_hash_tab,"chem_comp_atom");
  if (tb != NULL){ 
    rn = &(tb->head_rowline);
    i = 0;
    while (rn->next != NULL){
      rn = rn->next;
      atm = &(mol->atoms[i]);
      atm->num = i;
      sprintf(atm->atomname,"%s",data_from_rowline(tb,rn,"atom_id"));
      sprintf(atm->element,"%s",data_from_rowline(tb,rn,"type_symbol"));
      atm->Pos[0] = atof(data_from_rowline(tb,rn,"pdbx_model_Cartn_x_ideal"));
      atm->Pos[1] = atof(data_from_rowline(tb,rn,"pdbx_model_Cartn_y_ideal"));
      atm->Pos[2] = atof(data_from_rowline(tb,rn,"pdbx_model_Cartn_z_ideal"));
      atm->charge = atof(data_from_rowline(tb,rn,"charge"));
      atm->num_in_file = atoi(data_from_rowline(tb,rn,"pdbx_ordinal"));
      /*
      printf("#ATOM %d %s %s %f %f %f %f\n",atm->num, atm->atomname,atm->element,atm->Pos[0], atm->Pos[1],atm->Pos[2],atm->charge);
      */ 
      i += 1; 
    }
  }
  /*** (4) Setup BONDs from 'chem_comp_bond' table ***/
  tb = TABLE_from_key(BLK.table_hash_tab,"chem_comp_bond");
  if (tb != NULL){ 
  rn = &(tb->head_rowline);
    i = 0;
    while (rn->next != NULL){
      rn = rn->next;
      bnd = &(mol->bonds[i]);
      bnd->num = i;
      sprintf(buff,"%s",data_from_rowline(tb,rn,"value_order"));
      if (strcmp(buff,"SING")==0){ bnd->bondtype = '1'; }
      if (strcmp(buff,"DOUB")==0){ bnd->bondtype = '2'; }
      if (strcmp(buff,"TRIP")==0){ bnd->bondtype = '3'; }
      sprintf(buff,"%s",data_from_rowline(tb,rn,"atom_id_1"));
      bnd->anum1 =  find_atom_num_from_atomname(mol->atoms, mol->Natom,buff);
      sprintf(buff,"%s",data_from_rowline(tb,rn,"atom_id_2"));
      bnd->anum2 =  find_atom_num_from_atomname(mol->atoms, mol->Natom,buff);
     /*
      printf("#BOND %d '%c' %d %d %s %s\n",bnd->num, bnd->bondtype,bnd->anum1, bnd->anum2,
          mol->atoms[bnd->anum1].atomname, mol->atoms[bnd->anum2].atomname);
      */
      bnd->bondstereo =  '0';
      mol->conmap.map[bnd->anum1][bnd->anum2] = bnd->bondtype; 
      mol->conmap.map[bnd->anum2][bnd->anum1] = bnd->bondtype; 
      i += 1; 
    }
  }
  /*** (5) Ending ***/
  free_WORDs(&HeadWord);
  free_BLOCK(&BLK);

  return(1);
} /* end of Translate_Lines_into_CCD_CIF_Molecule() */


int find_atom_num_from_atomname(atoms, Natom, query_atomname)
  struct ATOM *atoms;
  int    Natom;
  char   *query_atomname;
{
  int a;
  for (a=0;a<Natom;++a){
    if (strcmp(atoms[a].atomname,query_atomname)==0){
      return(a);
    }
  }
  return(-1);
} /* end of find_atom_num_from_atomname() */




