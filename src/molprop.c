/*

 <molprop.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  for assigning various properties for the molecule,
  such as extended connectivity, ring, connected graphs.

*/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "molprop.h"
#include "molring.h"
#include "moltopodis.h"
#include "molpermu.h"
#include "ioLINE.h"
#include "ioSDF.h"
#include "ioPDB.h"
#include "ioKCF.h"
#include "ioMOL2.h"
#include "ioSMILES.h"
#include "kcombu_io_CIF.h"
#include "MOLECULE_from_CCD_CIF.h"
#include "ringblock.h"
#include "stereo_check.h"

/** FUNCTIONS (GLOBAL) **/
char MoleculeFileType();
int  Read_MOLECULE();
int  Set_MOLECULE();
void Initialize_MOLECULE();
void Initialize_MOLECULE_string();
void Malloc_MOLECULE();
void Free_MOLECULE();
void Copy_MOLECULE();
void Copy_MOLECULE_XYZ();
void Copy_ATOM();
void Copy_BOND();
void Set_Center_of_Gravity();
void Set_one_char_ele();
void Set_Nheavyatom_Nheavybond();
void Set_Nneighbor_CHNOSP();
void Set_Nneighbor_atomtype();
void Cal_Distance_Map();
void Cal_Distance_Map_for_Heavy_atoms();
void Set_Molform();
void Set_atomtype_of_Atoms();
void Count_atomtype_in_Natomtype();
float Max_Tanimoto_from_Nelement();
float Max_Tanimoto_from_Natomtype();
float Tanimoto_bwn_oneatom_descriptors();
int  Substructure_condition_oneatom_descriptors();
int  Isomorphic_condition_oneatom_descriptors();
int  Set_Extended_Connectivity();
void Set_DAMRL_from_2D_structure_with_hydrogen();
void Show_conmap();
void Set_max_atomtype_by_atomtype_class();
int Number_of_atomtype();
float Score_bwn_atomtypes();
char *atomtype_string_from_number();
char molecular_file_type_from_tail_of_file();
void Get_Part_Of_Line();
void Split_to_Words();
void Split_to_Words_With_Double_Splsym();
int First_Appear_Pos_in_String();
void Remove_Symbol_from_String();
void Remove_HeadTail_Space_from_String();
void Change_String_ToUpper();
void Change_Space_To_Underbar();
int Match_Tail_String();
int Number_of_Element();
int Check_Element_String_Or_Not();
void Write_atomtype_in_rasmol_script();
char* string_getenv();

/** FUNCTIONS (LOCAL) **/
static void Set_Aromatic_Ring_by_Sp2_atoms();
static void Set_Atom_Aromatic_by_Aromatic_Ring();
static void Make_Core_Filename();



static char element_list[MAX_NELEMENT][3] =
{"H",
 "HE","LI","BE","B","C","N","O","F","NE",
 "NA","MG","AL","SI","P","S","CL","AR","K",
 "CA","SC","TI","V","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR",
 "KR","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I","XE",
 "CS","BA",
 "LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU",
 "HF","TA","W", "RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN",
 "?"
};


static char element_list_lower[MAX_NELEMENT][3] =
{"H",
 "He","Li","Be","B","C","N","O","F","Ne",
 "Na","Mg","Al","Si","P","S","Cl","Ar","K",
 "Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br",
 "Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
 "Cs","Ba",
 "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
 "Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
 "?"
};



char MoleculeFileType(str)
  char *str;
{
  int L;

  L = strlen(str);
  if (L==1){
         if (str[0]=='P') return('P'); /* PDB */
    else if (str[0]=='2') return('2'); /* MOL2 */
    else if (str[0]=='S') return('S'); /* SDF */
    else if (str[0]=='K') return('K'); /* KCF */
    else if (str[0]=='C') return('C'); /* CIF */
    else if (str[0]=='M') return('M'); /* SMILES */
    else if (str[0]=='m') return('m'); /* SMILES STDIN */
    else if (str[0]=='c') return('c'); /* CIF */
    else return('?');
  }
  else{
         if (strcmp(str,"PDB")==0){  return('P');} /* PDB */
    else if (strcmp(str,"pdb")==0){  return('P');} /* PDB */
    else if (strcmp(str,"MOL2")==0){ return('2');} /* MOL2 */
    else if (strcmp(str,"mol2")==0){ return('2');} /* MOL2 */
    else if (strcmp(str,"KCF")==0){  return('K');} /* KCF */
    else if (strcmp(str,"kcf")==0){  return('K');} /* KCF */
    else if (strcmp(str,"SDF")==0){  return('S');} /* SDF */
    else if (strcmp(str,"sdf")==0){  return('S');} /* SDF */
    else if (strcmp(str,"SMI")==0){  return('M');} /* SMILES */
    else if (strcmp(str,"smi")==0){  return('M');} /* SMILES */
    else if (strcmp(str,"SMILES")==0){  return('M');} /* SMILES */
    else if (strcmp(str,"smiles")==0){  return('M');} /* SMILES */
    else if (strcmp(str,"CIF")==0){  return('C');} /* CIF */
    else if (strcmp(str,"cif")==0){  return('C');} /* CIF */
    else if (strcmp(str,"mmCIF")==0){  return('C');} /* CIF */
    else if (strcmp(str,"mmcif")==0){  return('C');} /* CIF */
    else if (strcmp(str,"MMCIF")==0){  return('C');} /* CIF */
    else return('?');
  }

} /* end of MoleculeFileType() */


int Read_MOLECULE(ifname,mol,FileType,AtmHetType,BondType)
  char *ifname;
  struct MOLECULE *mol;
  char   FileType;    /* '-' not assigned. 'S':sdf, 'P':pdb  'C':cif */
  char   AtmHetType;  /* 'A'tom, 'H'etatm, 'B'oth (only for PDB) */
  char   BondType;    /* 'B': consider bonds (cal topo-distances and ring)  */
{
  int ok;
  struct LINENODE HeadLineNode;

  printf("#Read_MOLECULE('%s',FileType:'%c' AtmHetType:'%c' BondType:'%c' chain:'%c')\n",
      ifname,FileType,AtmHetType,BondType,mol->chain); fflush(stdout);

  Initialize_MOLECULE(mol);
  Make_Core_Filename(mol->core_molname,ifname);
 
 
  if (FileType=='-'){
    mol->filetype =  molecular_file_type_from_tail_of_file(ifname);
  }
  else { mol->filetype = FileType;}
  ok  = 0;

  if (mol->filetype=='m'){ 
    HeadLineNode.next = NULL;
    HeadLineNode.prev = NULL;
    Add_string_to_LINENODEs(ifname,&HeadLineNode);
  }
  else{
    Read_LINENODEs(ifname,&HeadLineNode);
  }


  printf("#mol->filetype '%c'\n",mol->filetype);
  /* Show_LINENODEs(&HeadLineNode); */

       if (mol->filetype=='P'){ ok = Translate_Lines_into_PDB_Molecule(&HeadLineNode,mol,AtmHetType,BondType);}
  else if (mol->filetype=='S'){ ok = Translate_Lines_into_SDF_Molecule(&HeadLineNode,mol);}
  else if (mol->filetype=='K'){ ok = Translate_Lines_into_KCF_Molecule(&HeadLineNode,mol); } 
  else if (mol->filetype=='2'){ ok = Translate_Lines_into_MOL2_Molecule(&HeadLineNode,mol);}
  else if (mol->filetype=='M'){ ok = Translate_Lines_into_SMILES_Molecule(HeadLineNode.next->line,mol);}
  else if (mol->filetype=='m'){ ok = Translate_Lines_into_SMILES_Molecule(HeadLineNode.next->line,mol);}
  else if (mol->filetype=='C'){ ok = Translate_Lines_into_CCD_CIF_Molecule(&HeadLineNode,mol);}
 
  Free_LINENODEs(&HeadLineNode);
 /*
  printf("#ok %d\n",ok); fflush(stdout);  
  if (ok==0) exit(1);
 */
  if ((mol->maxNatom > 0) && (mol->Natom > mol->maxNatom)){
    printf("#ERROR:Natom %d is larger than maxNatom %d.\n",mol->Natom, mol->maxNatom);
    return(0); 
  }

  if (ok==0){return(0);}

  Set_MOLECULE(mol,BondType);

  return(ok);

} /* end of Read_MOLECULE() */




int Set_MOLECULE(mol,BondType)
  struct MOLECULE *mol;
  char   BondType;   /* 'B': consider bonds (cal topo-distances and ring)  */
{
 /*
  printf("#Set_MOLECULE(BondType '%c')\n",BondType); fflush(stdout);
 */
  Set_one_char_ele(mol);
  Set_Nheavyatom_Nheavybond(mol);

 if (BondType =='B'){
    if (PAR.Wdistan>0.0){
      Malloc_FLOAT2DMAP(&(mol->dismap),mol->Natom);
      /* Cal_Distance_Map_for_Heavy_atoms(mol); */
      Cal_Distance_Map(mol);
    }
    if (mol->Nbond==0){
      Guess_Bonds_from_Atom_XYZ(mol);
    }
 }

  Set_Center_of_Gravity(mol);

  if (mol->filetype != 'P'){Set_Unique_atomname_to_Molecule(mol);}

  Set_radius(mol);
  Set_Nneighbor_CHNOSP(mol);
  Set_Extended_Connectivity(mol);
  Set_Molform(mol);
  
  if (BondType == 'B'){
    Set_Topological_Distance(mol,'-');
    Set_Ring_Block_and_SSSR_Ring(mol); 
 /*
    Show_Ring_Structure(mol,&(mol->HeadRing),"mol.HeadRing");
    Show_Ring_Structure(mol,&(mol->HeadBlock),"mol.HeadBlock");
  */
  }

  Set_atomtype_of_Atoms(mol,PAR.atomtype_class);
  Count_atomtype_in_Natomtype(mol); 
  Set_Nneighbor_atomtype(mol);
  
       if (PAR.bondtype_class == 'A'){ Set_Aromatic_Bond_Atoms(mol);}
  else if (PAR.bondtype_class == 'R'){ Set_Rotatable_Bond(mol);}

  if ((PAR.AtomStereoFromBondStereo!='F') && (mol->filetype=='S')){
    Set_Atomic_Stereo_Parity_from_2D_BondStereo(mol,PAR.AtomStereoFromBondStereo);
    /* printf("#Set_Stereo_Parity_of_Molecule_from_2D_BondStereo('%s');\n",mol->filename); */
  }

  if (mol->SetStereo3D == 'T'){ 
    Set_stereo_parity_of_3D_molecule(mol); 
    printf("#Set_stereo_parity_of_3D_molecule('%s');\n",mol->filename);
  }

  if (PAR.Wrank > 0.0){
    Set_Aromatic_Bond_Atoms(mol);
    Set_Connected_Structure(mol);
    Set_Unique_Extended_Connectivity(mol,-1,0);
  }

  printf("#Natom %d Nbond %d Nheavyatom %d Nheavybond %d Nring %d\n",mol->Natom, mol->Nbond, mol->Nheavyatom, mol->Nheavybond, mol->Nring);

  return(1);

} /* end of Set_MOLECULE() */








void Initialize_MOLECULE(mol)
  struct MOLECULE *mol;
{
  mol->atoms = NULL; 
  mol->bonds = NULL; 
  mol->Natom = 0;
  mol->Nbond = 0;
  mol->Nheavyatom = 0; 
  mol->Nheavybond = 0;
  mol->num_frm_heavy = NULL;
  mol->heavy_frm_num = NULL;
  mol->molform[0] = '\0';
  mol->chiral_flag = '0';
  mol->Nring = 0; 
  mol->Nringblock = 0; 
  mol->Nconstr = 0;
  mol->Nmatch = 0;
  mol->Npermu  = 0;
  mol->conmap.N = 0;
  mol->conmap.malloced = 0;
 
  mol->dismap.N = 0;
  mol->dismap.malloced = 0;
  
  mol->topodismap.N = 0;
  mol->topodismap.malloced = 0;
  
  mol->pathmap.N = 0;
  mol->pathmap.malloced = 0;
  
  mol->HeadRing.next  = NULL; mol->HeadRing.prev = NULL;
  mol->HeadRing.num   = -1; 
  mol->HeadRing.Natom = 0; 
  mol->HeadRing.num_atoms = NULL;
  mol->HeadBlock.next = NULL; mol->HeadBlock.prev = NULL;
  mol->HeadBlock.num   = -1; 
  mol->HeadBlock.Natom = 0; 
  mol->HeadBlock.num_atoms = NULL;
  
  mol->permuhead.next = NULL; mol->permuhead.prev = NULL;
  mol->permuhead.N    = -1;
  mol->permuhead.index_heavy = NULL;

  mol->HeadProperty.next = NULL;
  mol->HeadProperty.prev = NULL;
  
  mol->atompair_descriptor =NULL;
  mol->ring_descriptor =NULL;

} /* end of Initialize_MOLECULE() */


void Initialize_MOLECULE_string(mol)
  struct MOLECULE *mol;
{
  mol->name[0] = '\0';           /* name of molecule */
  mol->core_molname[0] = '\0';
  mol->filename[0] = '\0';       /* filename */
  mol->osdffile[0] = '\0';       /* output file in sdf */
  mol->opdbfile[0] = '\0';       /* output file in pdb */
  mol->omol2file[0] = '\0';      /* output file in mol2 */
  mol->comment[0] = '\0';        /* comment of molecule (2nd line of SDF file)*/


} /* end of Initialize_MOLECULE_string() */

void Malloc_MOLECULE(mol,Natom,Nbond)
 struct MOLECULE *mol;
 int Natom, Nbond;
{
  int a,b;
  /* printf("#Malloc_MOLECULE(mol,Natom %d Nbond %d)\n",Natom,Nbond); fflush(stdout); */
  mol->atoms = NULL; 
  mol->bonds = NULL; 
  mol->HeadProperty.next = NULL;
  mol->atompair_descriptor = NULL;
  mol->ring_descriptor = NULL; 

 
  /* (1) Malloc atoms */
  if (Natom>0){
    mol->Natom = Natom;
    mol->atoms = (struct ATOM*)malloc(sizeof(struct ATOM)*Natom);
    mol->num_frm_heavy = (int*)malloc(sizeof(int)*Natom);
    mol->heavy_frm_num = (int*)malloc(sizeof(int)*Natom);

    Malloc_CHAR2DMAP(&(mol->conmap),Natom); 
    Malloc_FLOAT2DMAP(&(mol->dismap),Natom);
    Malloc_INT2DMAP(&(mol->topodismap),Natom);
    Malloc_INT2DMAP(&(mol->pathmap),Natom);

    for (a=0;a<Natom;++a){
      mol->atoms[a].num = a;
      mol->atoms[a].num_in_file = 0;
      mol->atoms[a].num_heavy = 0;
      mol->atoms[a].Pos[0] = 0.0;
      mol->atoms[a].Pos[1] = 0.0;
      mol->atoms[a].Pos[2] = 0.0;
      mol->atoms[a].element[0] = '\0';
      mol->atoms[a].one_char_ele = ' ';
      mol->atoms[a].atomtype[0] = '\0';
      mol->atoms[a].atomname[0] = '\0';
      mol->atoms[a].resi[0] = '\0';
      mol->atoms[a].rnum[0] = '\0';
      mol->atoms[a].chain = ' ';
      mol->atoms[a].atmhet = ' ';
      mol->atoms[a].ring = ' ';
      mol->atoms[a].ringblock = ' ';
      mol->atoms[a].aromatic = ' ';
      mol->atoms[a].constr_num = 0;
      mol->atoms[a].DAMRL = ' ';
      mol->atoms[a].orbital = ' ';
      mol->atoms[a].tFactor = 0.0;
      mol->atoms[a].mass_diff  = 0;
      mol->atoms[a].char_charge  = '0';
      mol->atoms[a].stereo_parity  = '0';
      mol->atoms[a].mol2atomtype[0]  = '\0';
      mol->atoms[a].radius = 0.0;
      mol->atoms[a].charge  = 0.0;
      mol->atoms[a].Nneighbor   = 0;
      mol->atoms[a].Nnei_heavy = 0;
      mol->atoms[a].NneiC = 0;
      mol->atoms[a].NneiN = 0;
      mol->atoms[a].NneiO = 0;
      mol->atoms[a].NneiS = 0;
      mol->atoms[a].NneiP = 0;
      mol->atoms[a].NneiH = 0;
      mol->atoms[a].EC0 = 0;
      mol->atoms[a].EC1 = 0;
      mol->atoms[a].EC2 = 0;
      mol->atoms[a].EC3 = 0;
      mol->atoms[a].EC4 = 0;
      mol->atoms[a].rank = 0;
      mol->atoms[a].Nmatch = 0;
      mol->atoms[a].subnum = 0;
      mol->atoms[a].subfloat = 0.0;
      mol->atoms[a].mark = 0;
 
      mol->atoms[a].neighbors = NULL;
    }  
  }

  /* (2) Malloc bonds */
  if (Nbond>0){
    mol->Nbond = Nbond;
    mol->bonds = (struct BOND*)malloc(sizeof(struct BOND)*Nbond);
    for (b=0;b<Nbond;++b){
      mol->bonds[b].num   = b;
      mol->bonds[b].anum1 = -1;
      mol->bonds[b].anum2 = -1;
      mol->bonds[b].bondtype   = ' ';
      mol->bonds[b].bondstereo = ' ';
   }
  }
  
  
} /* end of Malloc_MOLECULE() */


void Free_MOLECULE(mol)
  struct MOLECULE *mol;
{
 int a;
 
 /*
 printf("#Free_MOLECULE(mol)Natom %d Nbond %d\n",mol->Natom, mol->Nbond); fflush(stdout);
 */

 if (mol->Nbond > 0){
   /* printf("#mol->Nbond %d\n",mol->Nbond); fflush(stdout); */
   if (mol->bonds !=NULL) {
     free(mol->bonds); 
     mol->bonds = NULL;
   }
   mol->Nbond = 0; 
 }

 if (mol->Natom>0){ 
   if (mol->dismap.malloced==1) {  Free_FLOAT2DMAP(&(mol->dismap));}
   if (mol->conmap.malloced==1) { Free_CHAR2DMAP(&(mol->conmap)); }
   if (mol->topodismap.malloced==1) { Free_INT2DMAP(&(mol->topodismap)); }
   if (mol->pathmap.malloced==1)    { Free_INT2DMAP(&(mol->pathmap)); }

   for (a=0;a<mol->Natom;++a){
     if (mol->atoms[a].Nneighbor>0) { free(mol->atoms[a].neighbors); mol->atoms[a].neighbors = NULL;}
   }
   if (mol->num_frm_heavy != NULL) {free(mol->num_frm_heavy); mol->num_frm_heavy= NULL;}
   if (mol->heavy_frm_num != NULL) {free(mol->heavy_frm_num); mol->heavy_frm_num= NULL;}
   if (mol->atoms != NULL) {free(mol->atoms);}
 }

/*
 if (mol->Nbond > 0){
   printf("#mol->Nbond %d\n",mol->Nbond); fflush(stdout);
   if (mol->bonds !=NULL) {
     free(mol->bonds); 
     mol->bonds = NULL;
   }
 }
 */

 if (mol->HeadRing.next != NULL)     { Free_RING(mol);}
 if (mol->HeadBlock.next != NULL)    { Free_BLOCK(mol);}
 if (mol->permuhead.next != NULL)    { Free_PERMUTATION(&(mol->permuhead));}
 if (mol->HeadProperty.next != NULL) { Free_KEYVALUEs(&(mol->HeadProperty));}
 mol->filename[0] = mol->core_molname[0] = mol->name[0] = '\0';
 mol->Natom = mol->Nbond = mol->Nheavyatom = 0;
 if (mol->atompair_descriptor !=NULL){ free(mol->atompair_descriptor); mol->atompair_descriptor = NULL;}
 if (mol->ring_descriptor     !=NULL){ free(mol->ring_descriptor);     mol->ring_descriptor = NULL;}

} /* end of Free_MOLECULE() */




void Copy_MOLECULE(molN,molO) /* molN <- molO */
 struct MOLECULE *molN, *molO;
{
  int i,j;
  /* 
    printf("molN Natom %d molO Natom %d\n",molN->Natom, molO->Natom); 
   */


  sprintf(molN->filename,"%s",molO->filename);
  sprintf(molN->name,"%s",molO->name);
  sprintf(molN->resi,"%s",molO->resi);
  molN->chain = molO->chain;
  molN->Natom      = molO->Natom;
  molN->Nbond      = molO->Nbond;
  molN->Nheavyatom = molO->Nheavyatom;
  molN->Nring = molO->Nring;
  molN->Nconstr = molO->Nconstr;


  for (i=0;i<3;++i){
    molN->Gcen[i] = molO->Gcen[i];
    molN->eigval[i] = molO->eigval[i];
    for (j=0;j<3;++j){
      molN->eigvec[i][j] = molO->eigvec[i][j];
    }
  }


  for (i=0;i<molN->Natom;++i){ Copy_ATOM(&(molN->atoms[i]),&(molO->atoms[i])); }

  for (i=0;i<molN->Nheavyatom;++i){  molN->num_frm_heavy[i] = molO->num_frm_heavy[i];}


  for (i=0;i<molN->Nbond;++i){
    Copy_BOND(&(molN->bonds[i]),&(molO->bonds[i]));
    /*
    molN->bonds[i].num = molO->bonds[i].num;
    molN->bonds[i].anum1 = molO->bonds[i].anum1;
    molN->bonds[i].anum2 = molO->bonds[i].anum2;
    molN->bonds[i].bondtype = molO->bonds[i].bondtype;
    */
  }


  if ((molN->conmap.malloced==1)&&(molO->conmap.malloced==1)){
    for (i=0;i<molN->Natom;++i){
      for (j=0;j<molN->Natom;++j){
        molN->conmap.map[i][j] = molO->conmap.map[i][j];
      }
    }
  }

  if ((molN->dismap.malloced==1)&&(molO->dismap.malloced==1)){
    for (i=0;i<molN->Natom;++i){
      for (j=0;j<molN->Natom;++j){  
         molN->dismap.map[i][j] = molO->dismap.map[i][j];
      }
    }
  }


  if ((molN->topodismap.malloced==1)&&(molO->topodismap.malloced==1)){
    for (i=0;i<molN->Natom;++i){
      for (j=0;j<molN->Natom;++j){  
        molN->topodismap.map[i][j] = molO->topodismap.map[i][j];
      }
    }
  }


  if ((molN->pathmap.malloced==1)&&(molO->pathmap.malloced==1)){
    for (i=0;i<molN->Natom;++i){
      for (j=0;j<molN->Natom;++j){  
        molN->pathmap.map[i][j] = molO->pathmap.map[i][j];
      }
    }
  }

  Malloc_Copy_RINGs(&(molN->HeadRing),&(molO->HeadRing));
  Malloc_Copy_RINGs(&(molN->HeadBlock),&(molO->HeadBlock));

  /* permuhead is NOT copied !! */
  molN->permuhead.next = NULL; 
  molN->permuhead.prev = NULL; 
  molN->permuhead.N    = 0; 

} /* end of Copy_MOLECULE() */


void Copy_MOLECULE_XYZ(molN,molO) /* molN <- molO */
 struct MOLECULE *molN, *molO;
{
  int i;

  if (molN->Natom != molO->Natom){
    printf("#ERROR(Copy_MOLECULE_XYZ):Natom (%d, %d) are not same.\n",molN->Natom,molO->Natom);  
    exit(1); 
  }

  for (i=0;i<molN->Natom;++i){
    molN->atoms[i].Pos[0] = molO->atoms[i].Pos[0];
    molN->atoms[i].Pos[1] = molO->atoms[i].Pos[1];
    molN->atoms[i].Pos[2] = molO->atoms[i].Pos[2];
  }


} /* end of Copy_MOLECULE_XYZ() */



void Copy_ATOM(atomN,atomO)
  struct ATOM *atomN; /* New      ATOM */
  struct ATOM *atomO; /* Original ATOM */
{
  int j;
/*
  printf("#Copy_ATOM(atomN,atomO %d '%s' element '%s' atomtype '%s' mol2atomtype '%s' resi '%s' Nneighbor %d)\n",atomO->num,atomO->atomname,atomO->element,atomO->atomtype,atomO->mol2atomtype,atomO->resi,atomO->Nneighbor); 
  for (j=0;j<atomO->Nneighbor;++j){
    printf("#neighbors[%d] %d\n",j,atomO->neighbors[j]);
  }
  fflush(stdout);
*/
  atomN->num = atomO->num;
  atomN->num_in_file = atomO->num_in_file;
  atomN->num_heavy = atomO->num_heavy;
  atomN->one_char_ele = atomO->one_char_ele;
  atomN->tFactor = atomO->tFactor;
  atomN->radius = atomO->radius;
  sprintf(atomN->element,"%s",atomO->element);
  sprintf(atomN->atomname,"%s",atomO->atomname);
  sprintf(atomN->atomtype,"%s",atomO->atomtype);
  sprintf(atomN->mol2atomtype,"%s",atomO->mol2atomtype);
  sprintf(atomN->resi,"%s",atomO->resi);
  sprintf(atomN->rnum,"%s",atomO->rnum);
  atomN->chain  = atomO->chain;
  atomN->atmhet = atomO->atmhet;
  atomN->ring   = atomO->ring;
  atomN->ringblock     = atomO->ringblock;
  atomN->aromatic      = atomO->aromatic;
  atomN->constr_num    = atomO->constr_num;
  atomN->DAMRL         = atomO->DAMRL;
  atomN->orbital       = atomO->orbital;
  atomN->tFactor       = atomO->tFactor;
  atomN->mass_diff     = atomO->mass_diff;
  atomN->char_charge   = atomO->char_charge;
  atomN->stereo_parity = atomO->stereo_parity;
  atomN->charge        = atomO->charge;
  sprintf(atomN->mol2atomtype,"%s",atomO->mol2atomtype);
  atomN->radius        = atomO->radius;


  atomN->mark   = atomO->mark;  
  atomN->subnum = atomO->subnum;  
 
  for (j=0;j<3;++j){ atomN->Pos[j] = atomO->Pos[j];}
  if (atomN->Nneighbor>0){ free(atomN->neighbors); }
  atomN->neighbors = NULL;
  atomN->Nneighbor = 0;

 /*
  printf("#atomN->Nneighbor %d atomO->Nneighbor %d\n", atomN->Nneighbor,atomO->Nneighbor); fflush(stdout); 
 */
  atomN->Nneighbor = atomO->Nneighbor;

  if (atomN->Nneighbor>0) {
    /* printf("#malloc atomN->Nneighbor %d\n", atomN->Nneighbor); fflush(stdout);  */
    atomN->neighbors = (int *)malloc(sizeof(int)*atomN->Nneighbor); 
    for (j=0;j<atomN->Nneighbor;++j){ 
      atomN->neighbors[j] = atomO->neighbors[j]; 
    }
  }
  atomN->Nnei_heavy = atomO->Nnei_heavy;
  atomN->NneiC = atomO->NneiC;
  atomN->NneiN = atomO->NneiN;
  atomN->NneiO = atomO->NneiO;
  atomN->NneiS = atomO->NneiS;
  atomN->NneiP = atomO->NneiP;
  atomN->NneiH = atomO->NneiH;
  atomN->EC0   = atomO->EC0;
  atomN->EC1   = atomO->EC1;
  atomN->EC2   = atomO->EC2;

  for (j=0;j<PAR.max_atomtype;++j){ 
    atomN->Nnei_atomtype[j] = atomO->Nnei_atomtype[j];
  }

  atomN->ring  = atomO->ring;
  atomN->DAMRL  = atomO->DAMRL;
  atomN->constr_num = atomO->constr_num;

} /* end of Copy_ATOM() */


void Copy_BOND(bondN,bondO)
  struct BOND *bondN; /* New      BOND */
  struct BOND *bondO; /* Original BOND */
{
   bondN->num        = bondO->num;
   bondN->anum1      = bondO->anum1;
   bondN->anum2      = bondO->anum2;
   bondN->bondtype   = bondO->bondtype;
   bondN->bondstereo = bondO->bondstereo;
} /* end of Copy_BOND() */



void Set_Center_of_Gravity(mol)
  struct MOLECULE *mol;
{
 int a,k;

 for (k=0;k<3;++k){mol->Gcen[k] = 0.0;}

 for (a=0;a<mol->Natom;++a){
   for (k=0;k<3;++k)  mol->Gcen[k] += mol->atoms[a].Pos[k];
 }

 for (k=0;k<3;++k){mol->Gcen[k] /= mol->Natom;}

} /* end of Set_Center_of_Gravity() */



void Set_one_char_ele(mol)
 struct MOLECULE *mol;
{
 int a,L;

/*
HETATM 2676  C4  X37 A1480      10.497   8.882   3.134  1.00  0.00
HETATM 2677  N5  X37 A1480      10.448   7.376   3.115  1.00  0.00
HETATM 2678  HN5 X37 A1480      10.881   7.037   3.831  1.00  0.00
HETATM 2679 HN5A X37 A1480       9.574   7.078   3.136  1.00  0.00
HETATM 2944 CL1  GVP A1351      13.670  13.571  -0.385  1.00  0.00
HETATM 2945  C2  GVP A1351      12.956  12.366   0.639  1.00  0.00
HETATM 2946  C3  GVP A1351      12.776  11.079   0.166  1.00  0.00
HETATM 2947  C4  GVP A1351      12.209  10.116   0.972  1.00  0.00
HETATM 1452 CU    CU A 100      -8.925   3.189   2.681  1.00  0.00          CU  
*/

 mol->Nheavyatom = 0;
 for (a=0;a<mol->Natom;++a){
   L = strlen(mol->atoms[a].element);
   /** Hydrogen **/
   if ((mol->atoms[a].element[0]=='H') && (L==1)){
     mol->atoms[a].one_char_ele = 'H';
     mol->atoms[a].atomtype[0] = '\0';
     mol->atoms[a].num_heavy = -1;

   }
   /** Hydrogen **/
   else if ((mol->atoms[a].element[0]=='H') && 
       (strncmp(mol->atoms[a].element,"HG",2)!=0)&&(strncmp(mol->atoms[a].element,"Hg",2)!=0)&&
       (strncmp(mol->atoms[a].element,"HE",2)!=0)&&(strncmp(mol->atoms[a].element,"He",2)!=0)&&
       (strncmp(mol->atoms[a].element,"HF",2)!=0)&&(strncmp(mol->atoms[a].element,"Hf",2)!=0)&&
       (strncmp(mol->atoms[a].element,"Ho",2)!=0)&&(strncmp(mol->atoms[a].element,"HO",2)!=0)
         ) {
     mol->atoms[a].one_char_ele = 'H';
   }
   else{
   /** Heavy Atom **/
     mol->atoms[a].num_heavy = mol->Nheavyatom;
     mol->num_frm_heavy[mol->atoms[a].num_heavy] = mol->atoms[a].num;
     mol->heavy_frm_num[mol->atoms[a].num] = mol->atoms[a].num_heavy;
     mol->Nheavyatom += 1;
     if (L==0){ mol->atoms[a].one_char_ele = '?';}
     if (L==1){ mol->atoms[a].one_char_ele = mol->atoms[a].element[0];}
     else if (L>1){
            if (strncmp(mol->atoms[a].element,"CL",2)==0){ mol->atoms[a].one_char_ele = 'c';}
       else if (strncmp(mol->atoms[a].element,"Cl",2)==0){ mol->atoms[a].one_char_ele = 'c';}
       else if (strncmp(mol->atoms[a].element,"BR",2)==0){ mol->atoms[a].one_char_ele = 'b';}
       else if (strncmp(mol->atoms[a].element,"Br",2)==0){ mol->atoms[a].one_char_ele = 'b';}
       else if (strncmp(mol->atoms[a].element,"MG",2)==0){ mol->atoms[a].one_char_ele = 'm';}
       else if (strncmp(mol->atoms[a].element,"Mg",2)==0){ mol->atoms[a].one_char_ele = 'm';}
       else if (strncmp(mol->atoms[a].element,"NA",2)==0){ mol->atoms[a].one_char_ele = 'n';}
       else if (strncmp(mol->atoms[a].element,"Na",2)==0){ mol->atoms[a].one_char_ele = 'n';}
       else if (strncmp(mol->atoms[a].element,"ZN",2)==0){ mol->atoms[a].one_char_ele = 'z';}
       else if (strncmp(mol->atoms[a].element,"Zn",2)==0){ mol->atoms[a].one_char_ele = 'z';}
       else if (strncmp(mol->atoms[a].element,"FE",2)==0){ mol->atoms[a].one_char_ele = 'f';}
       else if (strncmp(mol->atoms[a].element,"Fe",2)==0){ mol->atoms[a].one_char_ele = 'f';}
       else if (strncmp(mol->atoms[a].element,"CU",2)==0){ mol->atoms[a].one_char_ele = 'u';}
       else if (strncmp(mol->atoms[a].element,"Cu",2)==0){ mol->atoms[a].one_char_ele = 'u';}
       else if (strncmp(mol->atoms[a].element,"CA",2)==0){ mol->atoms[a].one_char_ele = 'a';}
       else if (strncmp(mol->atoms[a].element,"Ca",2)==0){ mol->atoms[a].one_char_ele = 'a';}
       else if (strncmp(mol->atoms[a].element,"HG",2)==0){ mol->atoms[a].one_char_ele = 'h';}
       else if (strncmp(mol->atoms[a].element,"Hg",2)==0){ mol->atoms[a].one_char_ele = 'h';}
       else if (strncmp(mol->atoms[a].element,"RU",2)==0){ mol->atoms[a].one_char_ele = 'R';}
       else if (strncmp(mol->atoms[a].element,"Ru",2)==0){ mol->atoms[a].one_char_ele = 'R';}
       else if (strncmp(mol->atoms[a].element,"RE",2)==0){ mol->atoms[a].one_char_ele = 'r';}
       else if (strncmp(mol->atoms[a].element,"Re",2)==0){ mol->atoms[a].one_char_ele = 'r';}
       else {mol->atoms[a].one_char_ele = 'X';}
     }
   }
/*
   printf("#a %d %d atomname '%s' element '%s' one_char_ele '%c'\n",
    a,mol->atoms[a].num_in_file,mol->atoms[a].atomname, mol->atoms[a].element, mol->atoms[a].one_char_ele);
 */
 }

} /* end of Set_one_char_ele() */


void Set_Nheavyatom_Nheavybond(mol)
  struct MOLECULE *mol;
{
 int a,b,L;

/*
 printf("#Set_Nheavyatom_Nheavybond(Natom:%d conmap N %d mal %d cal %d)\n",
   mol->Natom, mol->conmap.N,mol->conmap.malloced, mol->conmap.calculated); fflush(stdout);
*/
/*
 for (a=0;a<mol->Natom;++a){
   for (b=a+1;b<mol->Natom;++b){
      printf("#BOND %d %d '%c'\n",a,b,mol->conmap.map[a][b]); fflush(stdout);
   }
 }
 */
 mol->Nheavyatom = 0;
 for (a=0;a<mol->Natom;++a){
   L = strlen(mol->atoms[a].element);
   if ((mol->atoms[a].element[0]!='H') || (L>1)){
     mol->atoms[a].num_heavy = mol->Nheavyatom;
     mol->num_frm_heavy[mol->atoms[a].num_heavy] = mol->atoms[a].num;
     mol->heavy_frm_num[mol->atoms[a].num] = mol->atoms[a].num_heavy;
     mol->Nheavyatom += 1;
   }
 }

 mol->Nheavybond = 0;
 for (a=0;a<mol->Natom;++a){
   for (b=0;b<mol->Natom;++b){
     /* printf("#BOND %d %d '%c'\n",a,b,mol->conmap.map[a][b]); */
     if ((a!=b)&&(mol->conmap.map[a][b] != '0')){
       if ((a<b)&&(mol->atoms[a].one_char_ele !='H')&&(mol->atoms[b].one_char_ele !='H')){ 
          mol->Nheavybond += 1;
        }
     }
   }
 }

/*
 printf("#Set_Nheavyatom_Nheavybond(Nheavyatom %d Nheavybond %d)\n",mol->Nheavyatom, mol->Nheavybond);
*/
} /* end of Set_Nheavyatom_Nheavybon() */





void Set_radius(mol)
 struct MOLECULE *mol;
{
  int a,L;
  struct ATOM *an;
/*
>> Van der Waals radii taken from Bondi's compilation (1964).[1] <<
Hydrogen(H) 1.20 
Carbon(C) 1.70 
Nitrogen(N) 1.55 
Oxygen(O) 1.52 
Fluorine(F) 1.47 
Phosphorus(P) 1.80 
Sulfur(S) 1.80 
Chlorine(Cl) 1.75 
Copper(Cu) 1.4 
Bondi, A. (1964). "Van der Waals Volumes and Radii". J. Phys. Chem. 68 (3): 441?51. doi:10.1021/j100785a001
*/
  for (a=0;a<mol->Natom;++a){
    an = &(mol->atoms[a]);
    L = strlen(an->element);
       if ((L==1)&&(strncmp(an->element,"H",1)==0)){  mol->atoms[a].radius = 1.20;}
  else if ((L==1)&&(strncmp(an->element,"C",1)==0)){  mol->atoms[a].radius = 1.70;}
  else if ((L==1)&&(strncmp(an->element,"N",1)==0)){  mol->atoms[a].radius = 1.55;}
  else if ((L==1)&&(strncmp(an->element,"O",1)==0)){  mol->atoms[a].radius = 1.52;}
  else if ((L==1)&&(strncmp(an->element,"F",1)==0)){  mol->atoms[a].radius = 1.47;}
  else if ((L==1)&&(strncmp(an->element,"S",1)==0)){  mol->atoms[a].radius = 1.80;}
  else if ((L==1)&&(strncmp(an->element,"P",1)==0)){  mol->atoms[a].radius = 1.80;}
  else if ((L==1)&&(strncmp(an->element,"I",1)==0)){  mol->atoms[a].radius = 1.98;}
  else if ((L==1)&&(strncmp(an->element,"K",1)==0)){  mol->atoms[a].radius = 2.75;}
  else if ((L==2)&&(strncmp(an->element,"Cl",2)==0)){ mol->atoms[a].radius = 1.75;}
  else if ((L==2)&&(strncmp(an->element,"CL",2)==0)){ mol->atoms[a].radius = 1.75;}
  else if ((L==2)&&(strncmp(an->element,"Cu",2)==0)){ mol->atoms[a].radius = 1.40;}
  else if ((L==2)&&(strncmp(an->element,"CU",2)==0)){ mol->atoms[a].radius = 1.40;}
  else if ((L==2)&&(strncmp(an->element,"Br",2)==0)){ mol->atoms[a].radius = 1.85;}
  else if ((L==2)&&(strncmp(an->element,"BR",2)==0)){ mol->atoms[a].radius = 1.85;}
  else if ((L==2)&&(strncmp(an->element,"Zn",2)==0)){ mol->atoms[a].radius = 1.39;}
  else if ((L==2)&&(strncmp(an->element,"ZN",2)==0)){ mol->atoms[a].radius = 1.39;}
  else if ((L==2)&&(strncmp(an->element,"Mg",2)==0)){ mol->atoms[a].radius = 1.73;}
  else if ((L==2)&&(strncmp(an->element,"MG",2)==0)){ mol->atoms[a].radius = 1.73;}
  else if ((L==2)&&(strncmp(an->element,"Na",2)==0)){ mol->atoms[a].radius = 2.27;}
  else if ((L==2)&&(strncmp(an->element,"NA",2)==0)){ mol->atoms[a].radius = 2.27;}
  else{
     mol->atoms[a].radius= 1.70; /* radius for Carbon */
  } 
/*
  printf("#%3d '%s' '%s'%c %f\n",an->num,an->atomname,an->element,an->one_char_ele,mol->atoms[a].radius);
*/
 }


} /* end of Set_radius() */





void Set_Nneighbor_CHNOSP(mol)
 struct MOLECULE *mol;
{
 int a,b,j;
 int neighbors[100];

 mol->Nheavybond = 0;

 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].Nneighbor  = 0;
   mol->atoms[a].Nnei_heavy = 0;
   mol->atoms[a].NneiC =  mol->atoms[a].NneiH =  mol->atoms[a].NneiN = 0;
   mol->atoms[a].NneiO =  mol->atoms[a].NneiS =  mol->atoms[a].NneiP = 0;

   for (b=0;b<mol->Natom;++b){
     if ((a!=b)&&(mol->conmap.map[a][b] != '0')){
       /* mol->atoms[a].neighbors[mol->atoms[a].Nneighbor] = b; */
       neighbors[mol->atoms[a].Nneighbor] = b; 
       mol->atoms[a].Nneighbor += 1;

       if (mol->atoms[b].one_char_ele !='H') mol->atoms[a].Nnei_heavy += 1;
            if (mol->atoms[b].one_char_ele=='C') mol->atoms[a].NneiC += 1;
       else if (mol->atoms[b].one_char_ele=='H') mol->atoms[a].NneiH += 1;
       else if (mol->atoms[b].one_char_ele=='N') mol->atoms[a].NneiN += 1;
       else if (mol->atoms[b].one_char_ele=='O') mol->atoms[a].NneiO += 1;
       else if (mol->atoms[b].one_char_ele=='S') mol->atoms[a].NneiS += 1;
       else if (mol->atoms[b].one_char_ele=='P') mol->atoms[a].NneiP += 1;
       if ((a<b)&&(mol->atoms[a].one_char_ele !='H')&&(mol->atoms[b].one_char_ele !='H')){ mol->Nheavybond += 1;}
     }
   }
 
   mol->atoms[a].neighbors = (int *)malloc(sizeof(int)*mol->atoms[a].Nneighbor);
   for (j=0;j<mol->atoms[a].Nneighbor;++j) mol->atoms[a].neighbors[j] = neighbors[j];
   mol->atoms[a].orbital = 'x';
   if (mol->atoms[a].Nneighbor==3) mol->atoms[a].orbital = '2';
   if (mol->atoms[a].Nneighbor==4) mol->atoms[a].orbital = '3';
 }

} /* end of Set_Nneighbor_CHNOSP() */




void Set_Nneighbor_atomtype(mol)
 struct MOLECULE *mol;
{
 int a,b,j,m;

 for (a=0;a<mol->Natom;++a){
   for (j=0;j<PAR.max_atomtype;++j) mol->atoms[a].Nnei_atomtype[j] = 0;
    mol->atoms[a].Nnei_heavy = 0;

   for (b=0;b<mol->Natom;++b){
     if (mol->conmap.map[a][b] != '0'){
       if (mol->atoms[b].one_char_ele !='H'){ 
          mol->atoms[a].Nnei_heavy += 1;
          m = Number_of_atomtype(mol->atoms[b].atomtype);
          mol->atoms[a].Nnei_atomtype[m] += 1;
        }
     }
   }
 }

} /* end of Set_Nneighbor_atomtype() */






void Cal_Distance_Map(mol)
 struct MOLECULE *mol;
{
 int a,b;
 float dx,dy,dz,dis;

 for (a=0;a<mol->Natom;++a){
     for (b=a;b<mol->Natom;++b){
         dx = mol->atoms[a].Pos[0] - mol->atoms[b].Pos[0];
         dy = mol->atoms[a].Pos[1] - mol->atoms[b].Pos[1];
         dz = mol->atoms[a].Pos[2] - mol->atoms[b].Pos[2];
         dis = sqrt(dx*dx + dy*dy + dz*dz);
         /* printf("a %d b %d N %d\n",a,b,mol->Natom); fflush(stdout); */
         mol->dismap.map[a][b] = mol->dismap.map[b][a] = dis;
     }
 }
 mol->dismap.calculated = 1;

} /* end of Cal_Distance_Map() */









void Cal_Distance_Map_for_Heavy_atoms(mol)
 struct MOLECULE *mol;
{
 int a,b;
 float dx,dy,dz,dis;

/*
 printf("#Cal_Distance_Map_for_Heavy_atoms(mol %s Natom %d %d)\n",mol->filename,mol->Natom,mol->dismap.N);
*/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].one_char_ele != 'H'){
     for (b=a;b<mol->Natom;++b){
       if (mol->atoms[b].one_char_ele != 'H'){
         dx = mol->atoms[a].Pos[0] - mol->atoms[b].Pos[0];
         dy = mol->atoms[a].Pos[1] - mol->atoms[b].Pos[1];
         dz = mol->atoms[a].Pos[2] - mol->atoms[b].Pos[2];
         dis = sqrt(dx*dx + dy*dy + dz*dz);
         /* printf("a %d b %d N %d\n",a,b,mol->Natom); fflush(stdout); */
         mol->dismap.map[a][b] = mol->dismap.map[b][a] = dis;
       }
     }
    }
 }
 mol->dismap.calculated = 1;

} /* end of Cal_Distance_Map_for_Heavy_atoms() */



void Set_Molform(mol)
 struct MOLECULE *mol;
{
 int i,a,e,m;
 int mform_len;
 char buff[16];

 /** [1] Count Number of element Nelement[] */
 for (a=0;a<MAX_NELEMENT;++a) mol->Nelement[a] = 0;
 for (i=0;i<mol->Natom;++i){
   e = Number_of_Element(mol->atoms[i].element);
   mol->Nelement[e] += 1;
 }

 /** [2] malloc molform[] **/
 mform_len = 2;
 for (a=0;a<MAX_NELEMENT;++a){
  if (mol->Nelement[a]>0)  mform_len += 3;
  if (mol->Nelement[a]>=10) mform_len += 2;
 }

 /** [3] make molform[] **/
 m = 0;
 for (a=0;a<MAX_NELEMENT;++a){
  if (mol->Nelement[a]>0){
    /* printf("Nelement[%d] %d\n",a,mol->Nelement[a]); */
    if (m>0) {mol->molform[m] = '_'; ++m;}
    if (mol->Nelement[a]==1) sprintf(buff,"%s",element_list[a]);
    else  sprintf(buff,"%s%d",element_list[a],mol->Nelement[a]);
    for (i=0;i<strlen(buff);++i) { 
      mol->molform[m] = buff[i]; 
      ++m;
      if (m==128){i= a = 10000000;} 
      }
  }
 }
 mol->molform[m] = '\0';

} /* end of Set_Molform() */




void Set_atomtype_of_Atoms(mol,atomtype_class)
 struct MOLECULE *mol;
 char atomtype_class;
{
 int a;
 char aromastr[2],element[8];
 /*
  ** Basically, Hydrogen atoms are not considered for atomtype classification.
     The atomtype for all the hydrogens should be 'H'. 

  */

 /** [1] Processes only for 'DAMRL' types **/
 if (atomtype_class=='D'){ 
    if ((PAR.DAMRLassign=='3') || ((PAR.DAMRLassign=='X') && (mol->filetype=='P'))){
      Set_Aromatic_Ring_by_3D_Planarity(mol);
      Set_Atom_Aromatic_by_Aromatic_Ring(mol);
      Guess_Orbital_from_Atom_XYZ(mol);
      Set_DAMRL_from_3D_structure(mol);
      Set_DAMRL_from_3D_RING(mol);
    }
    else if ((PAR.DAMRLassign=='H') || ((PAR.DAMRLassign=='X') && (mol->filetype!='P'))){
      Set_Aromatic_Ring_by_Sp2_atoms(mol);
      Set_Atom_Aromatic_by_Aromatic_Ring(mol);
      Set_DAMRL_from_2D_structure_with_hydrogen(mol);
    }
  }


 /*** [2] Assign atomtype ***/

 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].one_char_ele == 'H') sprintf(mol->atoms[a].atomtype,"H");
   else{ 
    
   sprintf(element,"%s",mol->atoms[a].element);
   Change_String_ToUpper(element);

   if (mol->atoms[a].aromatic!='-') sprintf(aromastr,"@"); else aromastr[0] = '\0';

          if (atomtype_class == 'E'){
       sprintf(mol->atoms[a].atomtype,"%s",element);
     }
     else if (atomtype_class == 'B'){
       sprintf(mol->atoms[a].atomtype,"%s%d",element,mol->atoms[a].Nnei_heavy);
     }
     else if (atomtype_class == 'R'){
       sprintf(mol->atoms[a].atomtype,"%s%s",element,aromastr);
     }
     else if (atomtype_class == 'K'){
       sprintf(mol->atoms[a].atomtype,"%s%s",element,aromastr);
       /*
            if (mol->atoms[a].one_char_ele=='C'){ sprintf(mol->atoms[a].atomtype,"C%s",aromastr); }
       else if (mol->atoms[a].one_char_ele=='N'){ sprintf(mol->atoms[a].atomtype,"N%s",aromastr); }
       else if (mol->atoms[a].one_char_ele=='O'){ sprintf(mol->atoms[a].atomtype,"O%s",aromastr); }
       else if (mol->atoms[a].one_char_ele=='P'){ sprintf(mol->atoms[a].atomtype,"P%s",aromastr); }
       else if (mol->atoms[a].one_char_ele=='S'){ sprintf(mol->atoms[a].atomtype,"S%s",aromastr); }
       else { sprintf(mol->atoms[a].atomtype,"X%s",aromastr); }
       */ 
       if ((mol->atoms[a].one_char_ele=='O') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"O1");
       if ((mol->atoms[a].one_char_ele=='N') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"N1");
     }
     else if (atomtype_class == 'H'){
       if ((strcmp(element,"F")==0)||(strcmp(element,"CL")==0)||(strcmp(element,"Cl")==0)
         ||(strcmp(element,"BR")==0)||(strcmp(element,"Br")==0)||(strcmp(element,"I")==0)){
         sprintf(mol->atoms[a].atomtype,"C%s",aromastr);
       }
       else{ sprintf(mol->atoms[a].atomtype,"%s%s",element,aromastr);}

       if ((mol->atoms[a].one_char_ele=='O') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"O1");
       if ((mol->atoms[a].one_char_ele=='N') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"N1");
     }
     else if (atomtype_class == 'F'){
       if ((strcmp(element,"F")==0)||(strcmp(element,"CL")==0)||(strcmp(element,"Cl")==0)
         ||(strcmp(element,"BR")==0)||(strcmp(element,"Br")==0)||(strcmp(element,"I")==0)){
         sprintf(mol->atoms[a].atomtype,"F%s",aromastr);
       }
       else{ sprintf(mol->atoms[a].atomtype,"%s%s",element,aromastr);}

       if ((mol->atoms[a].one_char_ele=='O') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"O1");
       if ((mol->atoms[a].one_char_ele=='N') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"N1");
     }

     else if (atomtype_class == 'x'){
       sprintf(mol->atoms[a].atomtype,"%s%s",element,aromastr);
       if ((mol->atoms[a].one_char_ele=='O') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"O1");
       if ((mol->atoms[a].one_char_ele=='N') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"N1");
       if ((mol->atoms[a].one_char_ele=='C') && (mol->atoms[a].Nnei_heavy==1)) sprintf(mol->atoms[a].atomtype,"C1");
     }
     else if (atomtype_class == 'T'){
       sprintf(mol->atoms[a].atomtype,"%s%d%s",element,mol->atoms[a].Nnei_heavy,aromastr);
     }
     else if (atomtype_class == 'D'){ 
       sprintf(mol->atoms[a].atomtype,"%c",mol->atoms[a].DAMRL);
     }
     else if (atomtype_class == 'X'){
       sprintf(mol->atoms[a].atomtype,"X");
     }
     else if (atomtype_class == 'k'){
       sprintf(mol->atoms[a].atomtype,"%s",mol->atoms[a].atomname);
     }
     else{
       sprintf(mol->atoms[a].atomtype,"X");
     }
 /*
     printf("element %s ring %c atomtype '%s'\n",element,mol->atoms[a].ring,mol->atoms[a].atomtype);
 */
   }
 }

} /* end of Set_atomtype_of_Atoms() */







float Max_Tanimoto_from_Nelement(molA, molB)
  struct MOLECULE *molA, *molB;
{
 int e,Nmcs_max;
 float tanimoto_max;

 Nmcs_max = 0;
 for (e=0;e<MAX_NELEMENT;++e){
   if (e!=1){
     /* printf("e %d NeleA %d NeleB %d\n",e,NelementA[e],NelementB[e]); */
     if (molA->Nelement[e]<molB->Nelement[e]) Nmcs_max += molA->Nelement[e];
                             else   Nmcs_max += molB->Nelement[e];
   }
 }

  if ((molA->Nheavyatom==0) && (molB->Nheavyatom==0)) return(1.0);
 tanimoto_max = (float)Nmcs_max/(float)(molA->Nheavyatom + molB->Nheavyatom - Nmcs_max);
 return(tanimoto_max);
} /* end of Max_Tanimoto_from_Nelement() */






float Max_Tanimoto_from_Natomtype(molA, molB)
  struct MOLECULE *molA, *molB;
{
  int m,Nmcs_max;
  float tanimoto_max;
  Nmcs_max = 0;
  /* printf("#Max_Tanimoto_from_Natomtype(max_atomtype:%d)\n",PAR.max_atomtype); */
  for (m=0;m<PAR.max_atomtype;++m){
    if (molA->oneatom_descriptor[m]<molB->oneatom_descriptor[m]) Nmcs_max += molA->oneatom_descriptor[m];
                                          else Nmcs_max += molB->oneatom_descriptor[m];
  }
  if ((molA->Nheavyatom==0) && (molB->Nheavyatom==0)) return(1.0);
  tanimoto_max = (float)Nmcs_max/(float)(molA->Nheavyatom + molB->Nheavyatom - Nmcs_max);
  return(tanimoto_max);
} /* end of Max_Tanimoto_from_Natomtype() */




float Tanimoto_bwn_oneatom_descriptors(odA, odB)
  unsigned char *odA,*odB; /* [MAX_ATOMTPE]: one atom descriptor */ 
{
  int a,Ncommon,NA,NB;
  float tanimoto; 
  NA = NB = Ncommon = 0;

  for (a=0;a<PAR.max_atomtype;++a){
    NA += odA[a];
    NB += odB[a];
    if (odA[a]<odB[a]) Ncommon += odA[a];
                 else  Ncommon += odB[a];
  }

  if ((NA==0) && (NB==0)) return(0.0);  
  tanimoto = (float)Ncommon/(float)(NA+NB-Ncommon);
  /* 
   printf(">max_atomtype %d Ncommon %d tanimoto %f\n",PAR.max_atomtype,Ncommon,tanimoto);
  for (a=0;a<PAR.max_atomtype;++a){
    printf("%d %d\n",odA[a],odB[a]); 
  }
  */
  return(tanimoto);
}





int Substructure_condition_oneatom_descriptors(odA, odB)
  unsigned char *odA, *odB;  /* one atom descriptors */
  /**  A is included by B **/
{
  int m;
  for (m=0;m<PAR.max_atomtype;++m){
    if (odA[m] > odB[m]) return(0);
  }
  return(1);
} /* end of Substructure_condition_oneatom_descriptors() */




int Isomorphic_condition_oneatom_descriptors(odA,odB)
  unsigned char *odA, *odB;  /* one atom descriptors */
{
  int m;
  for (m=0;m<PAR.max_atomtype;++m){
    if (odA[m] != odB[m]) return(0);
  }
  return(1);
} /* end of Isomorphic_condition_oneatom_descriptors() */




int Set_Extended_Connectivity(mol)
 struct MOLECULE *mol;
{
 int a,b;

 /*
 printf("#int Set_Extended_Connectivity('%s')\n",mol->filename);
 */

 /** [0] Calculate EC0 **/
 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].EC0 =  mol->atoms[a].EC1 =  mol->atoms[a].EC2 =  mol->atoms[a].EC3 = mol->atoms[a].EC4 = 0;
   if (mol->atoms[a].one_char_ele!='H'){
     for (b=0;b<mol->Natom;++b){
       if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
         mol->atoms[a].EC0 += 1;
       }
     }
   }
 }

 if (PAR.levelEC_Dextcon=='0') return(1);
 
 /** [1] Calculate EC1 **/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].element[0]!='H'){
     for (b=0;b<mol->Natom;++b){
       if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
         mol->atoms[a].EC1 += mol->atoms[b].EC0;
       }
     }
   }
 }

 if (PAR.levelEC_Dextcon=='1') return(1);
 
 /** [2] Calculate EC2 **/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].element[0]!='H'){
     for (b=0;b<mol->Natom;++b){
       if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
         mol->atoms[a].EC2 += mol->atoms[b].EC1;
       }
     }
   }
 }

 if (PAR.levelEC_Dextcon=='2') return(1);
 
 /** [3] Calculate EC3 **/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].element[0]!='H'){
     for (b=0;b<mol->Natom;++b){
       if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
         mol->atoms[a].EC3 += mol->atoms[b].EC2;
       }
     }
   }
 }

 if (PAR.levelEC_Dextcon=='3') return(1);
 
 /** [4] Calculate EC4 **/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].element[0]!='H'){
     for (b=0;b<mol->Natom;++b){
       if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
         mol->atoms[a].EC4 += mol->atoms[b].EC3;
       }
     }
   }
 }

 return(1);

} /* end of Set_Extended_Connectivity() */




void Set_Aromatic_Ring_by_Sp2_atoms(mol)
  struct MOLECULE *mol;
{
  struct RING *rn;
  int i,n,Nsp2;

  /*
  Aromatic ring :
    (1) ring->Natom == 5 or 6, and all of the atoms are carbons or nitrogen.
    (2) Number of sp2 atoms is >= 3  
    
    ** sp2 atoms : carbon with Nneighbor=3 and nitrogen with Nneighbor=2.
  */

  rn = &(mol->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;
    Nsp2 = 0;

    if ((rn->Natom == 5) || (rn->Natom == 6)){ 

      for (i=0;i<rn->Natom;++i){mol->atoms[rn->num_atoms[i]].mark = 1;}

      for (i=0;i<rn->Natom;++i){
        n = rn->num_atoms[i];
        if (((mol->atoms[n].one_char_ele=='C') && (mol->atoms[n].Nneighbor==3)) ||
            ((mol->atoms[n].one_char_ele=='N') && (mol->atoms[n].Nneighbor==2)) ) Nsp2 += 1;
        } 
     }

    if (Nsp2>=3) rn->aromatic = 'A'; else rn->aromatic = '-'; 
  }

} /* end of Set_Aromatic_Ring_by_Sp2_atoms() */



void Set_Atom_Aromatic_by_Aromatic_Ring(mol)
  struct MOLECULE *mol;
{
  int a,i,n;
  struct RING *rn;

  for (a=0;a<mol->Natom;++a){
    mol->atoms[a].aromatic = '-';
  }

  rn = &(mol->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;
    if (rn->aromatic=='A'){
      for (i=0;i<rn->Natom;++i){
        n = rn->num_atoms[i];
        mol->atoms[n].aromatic = 'A'; 
      }
    }
  }

} /* end of Set_Atom_Aromatic_by_Aromatic_Ring() */



void Set_DAMRL_from_2D_structure_with_hydrogen(mol)
 struct MOLECULE *mol;
{
 /*
   >> CAUTION << 
   The function Set_Nneighbor_CHNOSP(mol) should be executed before executing this funcion.
   Nneighbor, and NneiH and NneiO are neccesary.
 */ 
 int a,i,c;

 /** [1] Main loop for assigning DAMRL **/

 for (a=0;a<mol->Natom;++a){
   /** ['C'] **/
   if (mol->atoms[a].one_char_ele == 'C'){
     if (mol->atoms[a].aromatic=='A') mol->atoms[a].DAMRL = 'R';
                               else   mol->atoms[a].DAMRL = 'L'; 
   }

   /** ['O'] **/
   else if (mol->atoms[a].one_char_ele == 'O'){
     if (mol->atoms[a].NneiH >0) mol->atoms[a].DAMRL = 'M';
                            else mol->atoms[a].DAMRL = 'A';
     if (mol->atoms[a].NneiP>0) mol->atoms[a].DAMRL = 'A';
     if (mol->atoms[a].NneiS>0) mol->atoms[a].DAMRL = 'A';
     /* for 'C(=O)OH' */
     if (mol->atoms[a].NneiC==1){
       c = -1;
       for (i=0;i<mol->atoms[a].Nneighbor;++i){
         if (mol->atoms[mol->atoms[a].neighbors[i]].one_char_ele=='C') c = mol->atoms[a].neighbors[i];
       }
       if ((c>=0) && (mol->atoms[c].Nneighbor==3) && (mol->atoms[c].NneiO==2)) mol->atoms[a].DAMRL = 'A';
     } 
   }

   /** ['N'] **/
   else if (mol->atoms[a].one_char_ele == 'N'){

     if (mol->atoms[a].NneiH>0) mol->atoms[a].DAMRL = 'D';
     else{
       if (mol->atoms[a].Nneighbor==2) mol->atoms[a].DAMRL = 'A';
       else {
         if (mol->atoms[a].aromatic == 'A')   mol->atoms[a].DAMRL = 'R';
         else if (mol->atoms[a].Nneighbor==1) mol->atoms[a].DAMRL = 'L'; /* C-triple bond-N */
         else mol->atoms[a].DAMRL = 'L';
       }
     }

   }
   /** ['P','S','F','c'(Cl),'b'(Br),'I']  **/
   else if ((mol->atoms[a].one_char_ele == 'P') || (mol->atoms[a].one_char_ele == 'S')||
            (mol->atoms[a].one_char_ele == 'F') || (mol->atoms[a].one_char_ele == 'c') ||
            (mol->atoms[a].one_char_ele == 'b') || (mol->atoms[a].one_char_ele == 'I') )
      mol->atoms[a].DAMRL = 'L';

   /** ['H']  **/
   else if (mol->atoms[a].one_char_ele == 'H') mol->atoms[a].DAMRL = 'H';
   /** [others]  **/
   else mol->atoms[a].DAMRL = 'X';

  /*
    printf("%d %c Nnei %d NneiH %d aromatic '%c' ring '%c' dabrl '%c'\n",
     mol->atoms[a].num_in_file,mol->atoms[a].one_char_ele,mol->atoms[a].Nneighbor,mol->atoms[a].NneiH,
     mol->atoms[a].aromatic,mol->atoms[a].ring,mol->atoms[a].DAMRL);
  */
 }


} /* end of Set_DAMRL_from_2D_structure_with_hydrogen() */




int Number_of_Element(ele)
  char *ele;
{
 int i,Nele;
 char hit;
 Nele = 54;
 i = 0;
 hit = 0;
 while ((i<Nele) && (hit==0)){
   if (strcmp(ele,element_list[i])==0)       {hit = 1; return(i);}
   if (strcmp(ele,element_list_lower[i])==0) {hit = 1; return(i);}
  i += 1;
 }
 return(Nele); 
} /* end of Number_of_Element() */


int Check_Element_String_Or_Not(ele)
  char *ele;
{
 int i,Nele;
 char hit;
 Nele = 54;
 i = 0;
 hit = 0;
 while ((i<Nele) && (hit==0)){
  if (strcmp(ele,element_list[i])==0)       {hit = 1; return(1);}
  if (strcmp(ele,element_list_lower[i])==0) {hit = 1; return(1);}
  i += 1;
 }
 return(0); 
} /* end of Check_Element_String_Or_Not() */




void Count_atomtype_in_Natomtype(mol)
  struct MOLECULE *mol;
{
  int a,m;
  for (m=0;m<MAX_ATOMTYPE;++m) mol->oneatom_descriptor[m] = 0;
  for (a=0;a<mol->Natom;++a){
    if ((mol->atoms[a].one_char_ele != 'H')&&(mol->atoms[a].atomtype[0] != '\0')){
        m = Number_of_atomtype(mol->atoms[a].atomtype);
        if (mol->oneatom_descriptor[m]<255) mol->oneatom_descriptor[m] += 1;
    }
  }
} /* end of Count_atomtype_in_Natomtype() */



void Set_max_atomtype_by_atomtype_class(){
       if (PAR.atomtype_class == 'E') PAR.max_atomtype =  6;
  else if (PAR.atomtype_class == 'R') PAR.max_atomtype =  10;
  else if (PAR.atomtype_class == 'K') { PAR.max_atomtype =  12; PAR.ring_atomtype = 4;}
  else if (PAR.atomtype_class == 'H') { PAR.max_atomtype =  12; PAR.ring_atomtype = 4;}
  else if (PAR.atomtype_class == 'F') { PAR.max_atomtype =  13; PAR.ring_atomtype = 4;}
  else if (PAR.atomtype_class == 'x') PAR.max_atomtype =  13;
  else if (PAR.atomtype_class == 'D') PAR.max_atomtype =  6;
  else if (PAR.atomtype_class == 'k') PAR.max_atomtype =  68;
  else if (PAR.atomtype_class == 'X') PAR.max_atomtype =  1;
  else  PAR.max_atomtype =  1;
  /* printf("%c %d\n",PAR.atomtype_class,PAR.max_atomtype); */
} /* end of Set_max_atomtype_by_atomtype_class() */


int Number_of_atomtype(type)
  char *type;
{
  /* 
  L = strlen(type);
  if (L==0) {printf("#ERROR:Number_of_atomtype. Length of type is zero\n"); exit(1);}
  */

     if (PAR.atomtype_class == 'E'){
       if (strcmp(type,"C")==0) return(0);
  else if (strcmp(type,"O")==0) return(1);
  else if (strcmp(type,"N")==0) return(2);
  else if (strcmp(type,"S")==0) return(3);
  else if (strcmp(type,"P")==0) return(4);
  else return(5);
 }
 else if (PAR.atomtype_class == 'R'){
        if (strcmp(type,"C")==0)  return(0);
   else if (strcmp(type,"O")==0)  return(1);
   else if (strcmp(type,"N")==0)  return(2);
   else if (strcmp(type,"S")==0)  return(3);
   else if (strcmp(type,"P")==0)  return(4);
   else if (strcmp(type,"C@")==0) return(5);
   else if (strcmp(type,"O@")==0) return(6);
   else if (strcmp(type,"N@")==0) return(7);
   else if (strcmp(type,"S@")==0) return(8);
   else return(9);
 }
 else if (PAR.atomtype_class == 'K'){
        if (strcmp(type,"C@")==0) return(0);
   else if (strcmp(type,"O@")==0) return(1);
   else if (strcmp(type,"N@")==0) return(2);
   else if (strcmp(type,"S@")==0) return(3);
   else if (strcmp(type,"C")==0)  return(4);
   else if (strcmp(type,"O")==0)  return(5);
   else if (strcmp(type,"N")==0)  return(6);
   else if (strcmp(type,"S")==0)  return(7);
   else if (strcmp(type,"P")==0)  return(8);
   else if (strcmp(type,"O1")==0) return(9);
   else if (strcmp(type,"N1")==0) return(10);
   else return(11);
 }
 else if (PAR.atomtype_class == 'H'){
        if (strcmp(type,"C@")==0) return(0);
   else if (strcmp(type,"O@")==0) return(1);
   else if (strcmp(type,"N@")==0) return(2);
   else if (strcmp(type,"S@")==0) return(3);
   else if (strcmp(type,"C")==0)  return(4);
   else if (strcmp(type,"O")==0)  return(5);
   else if (strcmp(type,"N")==0)  return(6);
   else if (strcmp(type,"S")==0)  return(7);
   else if (strcmp(type,"P")==0)  return(8);
   else if (strcmp(type,"O1")==0) return(9);
   else if (strcmp(type,"N1")==0) return(10);
   else return(11);
 }

 else if (PAR.atomtype_class == 'F'){
        if (strcmp(type,"C@")==0) return(0);
   else if (strcmp(type,"O@")==0) return(1);
   else if (strcmp(type,"N@")==0) return(2);
   else if (strcmp(type,"S@")==0) return(3);
   else if (strcmp(type,"C")==0)  return(4);
   else if (strcmp(type,"O")==0)  return(5);
   else if (strcmp(type,"N")==0)  return(6);
   else if (strcmp(type,"S")==0)  return(7);
   else if (strcmp(type,"P")==0)  return(8);
   else if (strcmp(type,"O1")==0) return(9);
   else if (strcmp(type,"N1")==0) return(10);
   else if (strcmp(type,"F")==0)  return(11);
   else return(12);
 }

 else if (PAR.atomtype_class == 'x'){
        if (strcmp(type,"C")==0)  return(0);
   else if (strcmp(type,"O")==0)  return(1);
   else if (strcmp(type,"N")==0)  return(2);
   else if (strcmp(type,"S")==0)  return(3);
   else if (strcmp(type,"P")==0)  return(4);
   else if (strcmp(type,"C@")==0) return(5);
   else if (strcmp(type,"O@")==0) return(6);
   else if (strcmp(type,"N@")==0) return(7);
   else if (strcmp(type,"S@")==0) return(8);
   else if (strcmp(type,"O1")==0) return(9);
   else if (strcmp(type,"N1")==0) return(10);
   else if (strcmp(type,"C1")==0) return(11);
   else return(12);
 }
 else if (PAR.atomtype_class == 'D'){ 
        if (strcmp(type,"D")==0)  return(0);
   else if (strcmp(type,"A")==0)  return(1);
   else if (strcmp(type,"M")==0)  return(2);
   else if (strcmp(type,"R")==0)  return(3);
   else if (strcmp(type,"L")==0)  return(4);
   else return(5);
 }

 else if (PAR.atomtype_class == 'k'){ 
        if (strcmp(type,"C1a")==0)  return(0);
   else if (strcmp(type,"C1b")==0)  return(1);
   else if (strcmp(type,"C1c")==0)  return(2);
   else if (strcmp(type,"C1d")==0)  return(3);
   else if (strcmp(type,"C1x")==0)  return(4);
   else if (strcmp(type,"C1y")==0)  return(5);
   else if (strcmp(type,"C1z")==0)  return(6);
   else if (strcmp(type,"C2a")==0)  return(7);
   else if (strcmp(type,"C2b")==0)  return(8);
   else if (strcmp(type,"C2c")==0)  return(9);
   else if (strcmp(type,"C2x")==0)  return(10);
   else if (strcmp(type,"C2y")==0)  return(11);
   else if (strcmp(type,"C3a")==0)  return(12);
   else if (strcmp(type,"C3b")==0)  return(13);
   else if (strcmp(type,"C4a")==0)  return(14);
   else if (strcmp(type,"C5a")==0)  return(15);
   else if (strcmp(type,"C5x")==0)  return(16);
   else if (strcmp(type,"C6a")==0)  return(17);
   else if (strcmp(type,"C7a")==0)  return(18);
   else if (strcmp(type,"C7x")==0)  return(19);
   else if (strcmp(type,"C8x")==0)  return(20);
   else if (strcmp(type,"C8y")==0)  return(21);
   else if (strcmp(type,"C0")==0)   return(22);
   else if (strcmp(type,"N1a")==0)  return(23);
   else if (strcmp(type,"N1b")==0)  return(24);
   else if (strcmp(type,"N1c")==0)  return(25);
   else if (strcmp(type,"N1d")==0)  return(26);
   else if (strcmp(type,"N1x")==0)  return(27);
   else if (strcmp(type,"N1y")==0)  return(28);
   else if (strcmp(type,"N2a")==0)  return(29);
   else if (strcmp(type,"N2b")==0)  return(30);
   else if (strcmp(type,"N2x")==0)  return(31);
   else if (strcmp(type,"N2y")==0)  return(32);
   else if (strcmp(type,"N3a")==0)  return(33);
   else if (strcmp(type,"N4x")==0)  return(34);
   else if (strcmp(type,"N4y")==0)  return(35);
   else if (strcmp(type,"N5x")==0)  return(36);
   else if (strcmp(type,"N5y")==0)  return(37);
   else if (strcmp(type,"N0")==0)   return(38);
   else if (strcmp(type,"O1a")==0)  return(39);
   else if (strcmp(type,"O1b")==0)  return(40);
   else if (strcmp(type,"O1c")==0)  return(41);
   else if (strcmp(type,"O1d")==0)  return(42);
   else if (strcmp(type,"O2a")==0)  return(43);
   else if (strcmp(type,"O2b")==0)  return(44);
   else if (strcmp(type,"O2c")==0)  return(45);
   else if (strcmp(type,"O2x")==0)  return(46);
   else if (strcmp(type,"O3a")==0)  return(47);
   else if (strcmp(type,"O3b")==0)  return(48);
   else if (strcmp(type,"O3c")==0)  return(49);
   else if (strcmp(type,"O4a")==0)  return(50);
   else if (strcmp(type,"O5a")==0)  return(51);
   else if (strcmp(type,"O5x")==0)  return(52);
   else if (strcmp(type,"O6a")==0)  return(53);
   else if (strcmp(type,"O7a")==0)  return(54);
   else if (strcmp(type,"O7x")==0)  return(55);
   else if (strcmp(type,"O0")==0)   return(56);
   else if (strcmp(type,"S1a")==0)  return(57);
   else if (strcmp(type,"S2a")==0)  return(58);
   else if (strcmp(type,"S2x")==0)  return(59);
   else if (strcmp(type,"S3a")==0)  return(60);
   else if (strcmp(type,"S3x")==0)  return(61);
   else if (strcmp(type,"S4a")==0)  return(62);
   else if (strcmp(type,"S0")==0)   return(63);
   else if (strcmp(type,"P1a")==0)  return(64);
   else if (strcmp(type,"P1b")==0)  return(65);
   else if (strcmp(type,"X")==0)    return(66);
   else if (strcmp(type,"Z")==0)    return(67);
   else return(67);
 }

 else if (PAR.atomtype_class == 'X'){
   return(0);
 }
 else { printf("WARNING:atomtype_class '%c' is not proper.\n", PAR.atomtype_class); }
 return(0);
} /* end of Number_of_atomtype() */





float Score_bwn_atomtypes(typeA,typeB)
  char *typeA,*typeB;
{

  if (PAR.atomtype_class == 'D'){
     if (strcmp(typeA,typeB)==0) return(1.0);
   else{
          if ((strcmp(typeA,"D")==0)&&(strcmp(typeB,"M")==0))  return(1.0);
     else if ((strcmp(typeA,"M")==0)&&(strcmp(typeB,"D")==0))  return(1.0);
     else if ((strcmp(typeA,"A")==0)&&(strcmp(typeB,"M")==0))  return(1.0);
     else if ((strcmp(typeA,"M")==0)&&(strcmp(typeB,"A")==0))  return(1.0);
     else return(0.0);
   } 
  }
  else{
     if (strcmp(typeA,typeB)==0) return(1.0);
     else                        return(0.0);
  }

} /* end of Score_bwn_atomtypes() */









char *atomtype_string_from_number(n)
  int n;
{
 if (PAR.atomtype_class == 'K'){
        if (n==0) return("C@");
   else if (n==1) return("O@");
   else if (n==2) return("N@");
   else if (n==3) return("S@");
   else if (n==4) return("C");
   else if (n==5) return("O");
   else if (n==6) return("N");
   else if (n==7) return("S");
   else if (n==8) return("P");
   else if (n==9) return("O1");
   else if (n==10) return("N1");
   else if (n==11) return("X");
 }
 else if (PAR.atomtype_class == 'H'){
        if (n==0) return("C@");
   else if (n==1) return("O@");
   else if (n==2) return("N@");
   else if (n==3) return("S@");
   else if (n==4) return("C");
   else if (n==5) return("O");
   else if (n==6) return("N");
   else if (n==7) return("S");
   else if (n==8) return("P");
   else if (n==9) return("O1");
   else if (n==10) return("N1");
   else if (n==11) return("X");
 } 
 else if (PAR.atomtype_class == 'E'){
        if (n==0) return("C");
   else if (n==1) return("O");
   else if (n==2) return("N");
   else if (n==3) return("S");
   else if (n==4) return("P");
   else return("X");
 }

 return("-");
}


char molecular_file_type_from_tail_of_file(ifname)
  char *ifname;
{
  int Nword,Wsta[10],Wend[10];
  char tail[MAX_FILENAME];
 /*
  D00141.mol
  TRP.cif.gz 
 */
  Split_to_Words(ifname,'.',&Nword,Wsta,Wend,10);

  sprintf(tail,"%s",ifname);
  if (Nword>=2){ 
    Get_Part_Of_Line(tail,ifname,Wsta[Nword-1],Wend[Nword-1]);
  }
    
  if ((strcmp(tail,"Z")==0) || (strcmp(tail,"gz")==0)){
    Get_Part_Of_Line(tail,ifname,Wsta[Nword-2],Wend[Nword-2]);
  }

       if (strcmp(tail,"sdf")==0){  return('S');}
  else if (strcmp(tail,"pdb")==0){  return('P');}
  else if (strcmp(tail,"ent")==0){ return('P');}
  else if (strcmp(tail,"kcf")==0){ return('K');}
  else if (strcmp(tail,"mol2")==0){ return('2');}
  else if (strcmp(tail,"smi")==0){  return('M');}
  else if (strcmp(tail,"cif")==0){  return('C');}
  else if (strcmp(tail,"SDF")==0){  return('S');}
  else if (strcmp(tail,"PDB")==0){  return('P');}
  else if (strcmp(tail,"KCF")==0){  return('K');}
  else if (strcmp(tail,"MOL2")==0){ return('2');}
  else if (strcmp(tail,"SMI")==0){ return('M');}
  else if (strcmp(tail,"CIF")==0){  return('C');}
  else{ return('S');}  /* the default is 'S'DF !! */
  
} /* end of molecular_file_type_from_tail_of_file() */



void Show_conmap(mol,comment)
  struct MOLECULE *mol;
  char *comment;
{
 int i,j;

 printf(">%s\n",comment);
 printf("     ");
 for (j=0;j<mol->Natom;++j) {
      if (mol->atoms[j].one_char_ele != 'H'){ printf("%c",mol->atoms[j].one_char_ele);}
 }
 printf("\n");
 
 for (i=0;i<mol->Natom;++i){
   if (mol->atoms[i].one_char_ele != 'H'){
    printf("%c%3d:",mol->atoms[i].one_char_ele,i+1);
    for (j=0;j<mol->Natom;++j){
      if (mol->atoms[j].one_char_ele != 'H'){
        printf("%c",mol->conmap.map[i][j]);
      }
    }
    printf("\n"); 
   }
 }


} /* end of Show_conmap() */


void Get_Part_Of_Line(part,line,s,e)
  char *part;
  char *line;
  int  s,e;
{
 int i,E,L;
 L = strlen(line)-1;
 /* printf("#Get_Part_Of_Line(line:'%s',s:%d e:%d L:%d)\n",line,s,e,L); fflush(stdout); */
 if (line[L] == '\n') L -= 1;
 if (s<0) s = 0;
 if (e>L) E = L; else E = e;
 for (i=s;i<=E;++i) part[i-s] = line[i];
 part[E-s+1] = '\0';
} /* end of Get_Part_of_Line() */



void Split_to_Words(str,splsym,Nword,Wsta,Wend,Nwordmax)
 char *str;          /* Input String */
 char splsym;        /* Symbol for split */
 int *Nword;         /* Number of words  */
 int Wsta[];         /* Start point of str for a word */
 int Wend[];         /* End point of str for a word */
 int Nwordmax;       /* Maxium number of Nword  */
{
 /* [Example]
      str = "//abc/d/ef//ghi//"
      splsym = '/'
        -->
      Nword 4
      (Wsta,Wend) = {(0,2), (4,4), (6,7), (10,12)}
  */
 int i,L;
 L = strlen(str);
 *Nword = 0; i = 0;
 while ((i<L)&&(*Nword < Nwordmax)){
  if (str[i]!=splsym){
     Wsta[*Nword] = i;
     while ((str[i]!=splsym)&&(i<=(L-1))) { ++i; }
     Wend[*Nword] = i-1;
     ++(*Nword);
   }
  ++i;
 }
} /* end of Split_to_Words() */



void Split_to_Words_With_Double_Splsym(str,splsym1,splsym2,Nword,Wsta,Wend,Nwordmax)
 char *str;          /* Input String */
 char splsym1;       /* Symbol1 for split */
 char splsym2;       /* Symbol2 for split */
 int *Nword;         /* Number of words  */
 int Wsta[];         /* Start point of str for a word */
 int Wend[];         /* End point of str for a word */
 int Nwordmax;       /* Maxium number of Nword  */
{
 /* [Example]
      str = " / /abc/d/ef/ ghi//"
      splsym1 = '/' splsym2 = '/'
        -->
      Nword 4
      (Wsta,Wend) = {(0,2), (4,4), (6,7), (10,12)}
  */
 int i,L;
 L = strlen(str);
 *Nword = 0; i = 0;
 while ((i<L)&&(*Nword < Nwordmax)){
  if ((str[i]!=splsym1)&&(str[i]!=splsym2)){
     Wsta[*Nword] = i;
     while ((str[i]!=splsym1) && (str[i]!=splsym2)&&(i<=(L-1))) { ++i; }
     Wend[*Nword] = i-1;
     ++(*Nword);
   }
  ++i;
 }
} /* end of Split_to_Words_With_Double_Splsym() */



int First_Appear_Pos_in_String(line,sym)
 char *line,sym;
{
 int L,i,hit;
 L = strlen(line);
 hit = i = 0;
 while ((i<L)&&(hit==0))
  if (line[i]==sym) hit = 1; else ++i;
 if (hit==0) return(-1); else return(i);

} /* end of First_Appear_Pos_in_String() */


void Remove_Symbol_from_String(rsp,orig,target_symbol)
 char *rsp,*orig;
 char target_symbol;
{
 int i,j;
 j =0;
 for (i=0;i<strlen(orig);++i){
   if (orig[i]!=target_symbol) { rsp[j] = orig[i];  ++j;}
  }
 rsp[j] = '\0';
} /* end of Remove_Symbol_from_String() */


void Remove_HeadTail_Space_from_String(rsp,orig)
 char *rsp,*orig;
{
 int i,j,sta_i,end_i,Lorig;
 sta_i = end_i = -1;
 Lorig = strlen(orig);

 i = 0;
 while (i<Lorig){
   if (orig[i]!=' '){sta_i = i; i = Lorig+1;}
   else {i += 1;}
 }
 i = Lorig -1;
 while (i>=0){
   if (orig[i]!=' '){end_i = i; i = 0;}
   else { i -= 1;}
 }

 j = 0;
 for (i=sta_i;i<=end_i;++i){ rsp[j] = orig[i];  ++j;}
 rsp[j] = '\0';

} /* end of Remove_HeadTail_Space_from_String() */



void Change_String_ToUpper(str)
 char *str;
{
 int i,L;
 L = strlen(str);
 for (i=0;i<L;++i){
   str[i] = toupper(str[i]);
 }

} /* end of Change_String_ToUpper() */



void Change_Space_To_Underbar(str)
 char *str;
{
 int i,L;
 L = strlen(str);
 for (i=0;i<L;++i){
   if (str[i] == ' ') str[i] = '_';
 }

} /* end of Change_Space_to_Underbar() */


int Match_Tail_String(str,tailpat)
 char *str,*tailpat;
{
 int i,Lstr,Lpat;
 Lstr = strlen(str);
 Lpat = strlen(tailpat);
 for (i=0;i<Lpat;++i){
   if (tailpat[i] != str[Lstr-Lpat+i]) return(0);
 }
 return(1);

} /* end of Match_Tail_String() */


void Make_Core_Filename(core,filename)
  char *core;    /* core of filename ( to be calculated) */
  char *filename; /* filename (input) */
{
  int Nword,Wsta[100],Wend[100];

  printf("#void Make_Core_Filename(core,filename '%s')\n",filename);
  Split_to_Words(filename,'/',&Nword,Wsta,Wend,10);
  printf("#Nword %d \n",Nword);
  Get_Part_Of_Line(core,filename,Wsta[Nword-1],Wend[Nword-1]);
  printf("#core '%s'\n",core);

} /* end of Make_Core_Filename() */



void Write_atomtype_in_rasmol_script(ofname,mol)
  char *ofname;
  struct MOLECULE *mol;
{
  FILE *fp;
  int i;

  fp = fopen(ofname,"w");
  printf("#Write_atomtype_in_rasmol_script()-->'%s'\n",ofname);
  if (fp==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
  }

  fprintf(fp,"background white\n");
  fprintf(fp,"select all\n");
  fprintf(fp,"color white\n");

  for (i=0;i<mol->Natom;++i){
    if (PAR.atomtype_class=='D'){
      fprintf(fp,"select atomno=%d\n",mol->atoms[i].num_in_file);
      if ((mol->atoms[i].DAMRL=='D')||(mol->atoms[i].DAMRL=='A')||(mol->atoms[i].DAMRL=='M')||(mol->atoms[i].DAMRL=='R')||(mol->atoms[i].DAMRL=='L')){
        if (mol->atoms[i].DAMRL=='D'){fprintf(fp,"color blue\n");}
        if (mol->atoms[i].DAMRL=='A'){fprintf(fp,"color red\n");}
        if (mol->atoms[i].DAMRL=='M'){fprintf(fp,"color purple\n");}
        if (mol->atoms[i].DAMRL=='R'){fprintf(fp,"color darkgreen\n");}
        if (mol->atoms[i].DAMRL=='L'){fprintf(fp,"color gray\n");}
        fprintf(fp,"label \"%c\"\n",mol->atoms[i].DAMRL);
         fprintf(fp,"echo \"atomno %d element %s atomtype %s\"\n",mol->atoms[i].num_in_file,mol->atoms[i].element,mol->atoms[i].atomtype);
      }
    }
  }
  fprintf(fp,"select all\n");
  fclose(fp);
} /* end of Write_atomtype_in_rasmol_script() */


char* string_getenv(item)
  char *item;
{
  if (getenv(item)==NULL){
     return("NOT_REGISTERED");
  }
  else{
     return(getenv(item));
  }
}

