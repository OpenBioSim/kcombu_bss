/*

 <molring.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  for assigning ring structures for the molecule,

*/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "molprop.h"

/** FUNCTIONS (GLOBAL) **/
void Set_Ring_Structure_with_maxNatom_in_ring();
void Set_Property_using_RING();
void Show_Ring_Structure();
void Set_DAML_ring_ATOMs();
void Check_Planarity_of_RING_from_Atom_XYZ();
int  Delete_Disconnected_Parts_Of_Molecule();
void Write_DAMLatom_and_RING_in_PDB();
void Read_DAMLring_Molecule_in_PDB();
int Number_of_RING();
void Free_RING();
void Free_BLOCK();
void Malloc_Copy_RINGs();
void Set_Rotatable_Bond();


/** FUNCTIONS (LOCAL) **/
static int Mark_Neighbors_for_Cycle_with_maxNatom_in_ring();
static int Mark_Neighbors_for_Connected();
static void Sub_Vec(); 
static void Cross_Vec();
static void Normalize_Vec();
static float Angle_bwn_Vec();


void Set_Ring_Structure_with_maxNatom_in_ring(mol,HeadRing,maxNatom_in_ring)
  struct MOLECULE *mol;
  struct RING *HeadRing;
  int maxNatom_in_ring;
{
 int a,b;
 struct ATOM *an;
 struct RING *rn;

 /** [1] initialize (ring = '-', mark = 0) **/
 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].ring     = '-';
   mol->atoms[a].ringblock = ' ';
   mol->atoms[a].mark = 0;
 }
 HeadRing->num  = -1;
 HeadRing->prev = HeadRing->next = NULL;

 /** [2] Mark_Neighbors_for_Cycle_with_maxNatom_in_ring for each atom **/

 for (a=0;a<mol->Natom;++a){
   an = &(mol->atoms[a]);
/*
   printf(">>>>>a %d %s mark %d\n",a,mol->atoms[a].atomname,mol->atoms[a].mark); fflush(stdout); 
 */ 
  if ((an->ring == '-')&&(an->one_char_ele != 'H')&&(an->Nnei_heavy>1)){
     for (b=0;b<mol->Natom;++b) mol->atoms[b].mark = 0;
     Mark_Neighbors_for_Cycle_with_maxNatom_in_ring(mol,HeadRing,a,a,1,maxNatom_in_ring);
    }
 }
 rn = HeadRing;
 while (rn->next != NULL){
   rn = rn->next;
   mol->Nring += 1;
 }

}  /* end of Set_Ring_Structure_with_maxNatom_in_ring() */





void Set_Property_using_RING(mol)
 struct MOLECULE *mol;
{
 int i,j,a,b;
 struct RING *rn;
 char abcstr[100];
 sprintf(abcstr,"abcdefghijklmnopqrstuvwxyz1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZ");
 
 /** [1] Set 'r'ing bondtype **/
   rn = &(mol->HeadRing);
   while (rn->next != NULL){
     rn = rn->next;
     for (i=0;i<rn->Natom;++i){
       a = rn->num_atoms[i];
       for (j=i;j<rn->Natom;++j){
         b = rn->num_atoms[j];
         if (mol->conmap.map[a][b]!='0') mol->conmap.map[a][b] = mol->conmap.map[b][a] = 'r'; /* 'r'ing bond */

       } 
    }  
   }

 /*** [2] Set "block" for each atom ***/
 rn = &(mol->HeadRing);
 while (rn->next != NULL){
   rn = rn->next;
   for (a=0;a<rn->Natom;++a){
     if (mol->atoms[rn->num_atoms[a]].ringblock==' '){
       mol->atoms[rn->num_atoms[a]].ringblock = abcstr[rn->num%62];
     }
     else mol->atoms[rn->num_atoms[a]].ringblock = 'X';
   } 
 }
} /* end of Set_Property_using_RING() */





void Show_Ring_Structure(mol, HeadRing,comment)
  struct MOLECULE *mol;
  struct RING *HeadRing;
  char *comment;
{
  struct RING *rn;
  struct ATOM *an;
  int a;

  if (comment[0]!='\0') printf(">%s\n",comment); 
 
  rn = HeadRing;
  printf("#HeadRing %d Natom %d\n",rn->num,rn->Natom);
  while (rn->next != NULL){
    rn = rn->next;
    printf("#ring %d Natom %d [",rn->num,rn->Natom);
    for (a=0;a<rn->Natom;++a){
      an = &(mol->atoms[rn->num_atoms[a]]);
      printf(" %c%d",an->one_char_ele,an->num_in_file);
   }
    printf("]\n");
  }
  printf("\n");

} /* end of Show_Ring_Structure() */






void Set_DAML_ring_ATOMs(mol)
  struct MOLECULE *mol;
/*
 Making Pseudo Atom centers  only containig 
 "DON", "ACC", "MIX", "CAL", "RGC", "RGU", "RGL".
*/
{
 struct ATOM *oatms;
 struct RING *rn;
 int i,k,Natom_new,n;
 char DAMRL;

 /* printf("#void Set_DAML_ring_ATOMs(mol %s)\n",mol->name); */
 oatms = (struct ATOM*)malloc(sizeof(struct ATOM)*mol->Natom);
 for (i=0;i<mol->Natom;++i){
   Copy_ATOM(&(oatms[i]),&(mol->atoms[i]));
 }

 /** (1) Count  'D','A','B','L' 'RGC','RGL','RGU, and remalloc mol->atoms[] **/ 
 n = 0;
 for (i=0;i<mol->Natom;++i){
   /* printf("%s %s DAMRL %c oatms %c\n",mol->atoms[i].atomname,mol->atoms[i].resi,mol->atoms[i].DAMRL, oatms[i].DAMRL); */
   DAMRL = oatms[i].DAMRL;
   if ((DAMRL=='D')||(DAMRL=='A')||(DAMRL=='B')||(DAMRL=='L')){ n += 1;}
 }
 printf("#n %d\n",n); fflush(stdout);
 rn = &(mol->HeadRing);
 while (rn->next != NULL){
   rn = rn->next;
   /* printf("ring [%d] planarity %c\n",rn->num,rn->planarity); */
   if (rn->planarity=='p'){ n += 3;}
 }

 printf("#n %d\n",n); fflush(stdout);

 free(mol->atoms);
 mol->atoms = (struct ATOM*)malloc(sizeof(struct ATOM)*n);
 Free_FLOAT2DMAP(&(mol->dismap));  
 Malloc_FLOAT2DMAP(&(mol->dismap),n); 

 /** (2) Keep only 'D','A','B','L' atoms **/ 
 Natom_new = 0; 
 for (i=0;i<mol->Natom;++i){
   DAMRL = oatms[i].DAMRL;
   if ((DAMRL=='D')||(DAMRL=='A')||(DAMRL=='B')||(DAMRL=='L')){
     Copy_ATOM(&(mol->atoms[Natom_new]),&(oatms[i]));
     mol->atoms[Natom_new].one_char_ele = 'C';
     if (DAMRL=='D') sprintf(mol->atoms[Natom_new].atomname,"DON");
     if (DAMRL=='A') sprintf(mol->atoms[Natom_new].atomname,"ACC");
     if (DAMRL=='M') sprintf(mol->atoms[Natom_new].atomname,"MIX");
     if (DAMRL=='L') sprintf(mol->atoms[Natom_new].atomname,"CAL");
/*
     printf("%d %s %s %f %f %f\n",Natom_new,
mol->atoms[Natom_new].atomname,mol->atoms[Natom_new].resi,mol->atoms[Natom_new].Pos[0],mol->atoms[Natom_new].Pos[1],mol->atoms[Natom_new].Pos[2]);
 */ 
    Natom_new += 1;
   }
 }

 /** (3) Add "RGC","RGL","RGU" atoms **/
 rn = &(mol->HeadRing);
 while (rn->next != NULL){
   rn = rn->next;
   if (rn->planarity=='p'){
     sprintf(mol->atoms[Natom_new].atomname,"RGC");
     mol->atoms[Natom_new].one_char_ele = 'C';
     mol->atoms[Natom_new].DAMRL   = 'R';
     for (k=0;k<3;++k) mol->atoms[Natom_new].Pos[k] = rn->cpos[k];
     Natom_new += 1;

     sprintf(mol->atoms[Natom_new].atomname,"RGU");
     mol->atoms[Natom_new].one_char_ele = 'C';
     mol->atoms[Natom_new].DAMRL   = '-';
     for (k=0;k<3;++k)
       mol->atoms[Natom_new].Pos[k] = rn->cpos[k] + rn->rg * rn->nvec[k];
     Natom_new += 1;

     sprintf(mol->atoms[Natom_new].atomname,"RGL");
     mol->atoms[Natom_new].one_char_ele = 'C';
     mol->atoms[Natom_new].DAMRL   = '-';
     for (k=0;k<3;++k)
       mol->atoms[Natom_new].Pos[k] = rn->cpos[k] - rn->rg * rn->nvec[k];
     Natom_new += 1;
   } 
 }
 printf("#Natom_orig %d --> Natom_new %d\n",mol->Natom,Natom_new);
 mol->Natom = Natom_new;

 printf("Natom %d dismap_N %d\n",mol->Natom,mol->dismap.N); fflush(stdout);
 Cal_Distance_Map_for_Heavy_atoms(mol);
 /*
 for (i=0;i<mol->Natom;++i){
   printf("reset %5d %s %c %s %f %f %f\n",
     i,mol->atoms[i].atomname,mol->atoms[i].DAMRL,mol->atoms[i].resi,mol->atoms[i].Pos[0],mol->atoms[i].Pos[1],mol->atoms[i].Pos[2]);
 } 
 */


 free(oatms);

} /* end of Set_DAML_ring_ATOMs() */






int Mark_Neighbors_for_Cycle_with_maxNatom_in_ring(mol,HeadRing,s,a,Nmember,maxNatom_in_ring)
  struct MOLECULE *mol;
  struct RING *HeadRing;
  int s;       /* index of the start atom */
  int a;       /* index of the focused atom */
  int Nmember; /* Number of member atoms for the ring (start atom:1) */
  int maxNatom_in_ring;
{
 /*
 >> CAUTION << 
   The function Set_Nneighbor_CHNOSP(mol) should be executed before executing this funcion.
      Nneighbor, neighbors[] are neccesary.
 */

 int i,b,c,min_Nmember;
 struct RING *rn;
 
 min_Nmember = 3;

 /*
 printf("Mark_Neighbors_for_Cycle_with_maxNatom_in_ring  s %d ring %c a %d  Nmember %d\n",
   mol->atoms[s].num_in_file, mol->atoms[s].ring, mol->atoms[a].num_in_file,Nmember); fflush(stdout); 
 */

 /** [1] Stop if Nmember > maxNatom_in_ring , or the start atom is already marked as 'ring'. **/

 if (Nmember > maxNatom_in_ring)    return(1);
 if (mol->atoms[s].ring=='r')  return(1);

 mol->atoms[a].mark = 1;

 /** [2] Scan neighbors(b) around the specified atom (a) **/
 for (i=0;i<mol->atoms[a].Nneighbor;++i){
   b  = mol->atoms[a].neighbors[i];
   if ((a!=b) && (mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
     if (mol->atoms[b].mark==0){ 
        Mark_Neighbors_for_Cycle_with_maxNatom_in_ring(mol,HeadRing,s,b,Nmember+1,maxNatom_in_ring); 
     }
     else if ((Nmember>=min_Nmember) && (b==s)){
       /** (C) Procedures for after finding cycle  **/
       /* 
       printf("#Find cycle!! s[%d] a[%d] b[%d] Nmember %d\n",
        mol->atoms[s].num_in_file, mol->atoms[a].num_in_file,mol->atoms[b].num_in_file,Nmember); 
       */
       /**  (C-1) Malloc new RING structure **/
       rn = HeadRing;
       while (rn->next !=NULL) rn = rn->next;
       rn->next = (struct RING*)malloc(sizeof(struct RING));
       rn->next->prev = rn;
       rn->next->next = NULL;
       rn->next->num = rn->num + 1;
       rn = rn->next;
       rn->Natom = 0;      
       rn->num_atoms = (int *)malloc(sizeof(int)*maxNatom_in_ring);
        /** (C-2) store atom numbers of ring members  **/
       for (c=0;c<mol->Natom;++c){
         if (mol->atoms[c].mark > 0){
           if (rn->Natom >= maxNatom_in_ring) 
             {printf("#ERROR:rn->Natom is over maxNatom_in_ring %d!!\n",maxNatom_in_ring); exit(1);}
           rn->num_atoms[rn->Natom] = c;
           rn->Natom += 1;
           mol->atoms[c].ring = 'r';
         } 
       }
       /** (C-3) escape with return(1) **/
       return(1); 
    }
  }

  } /* i */

  /** [3] Reset mark for the focused atom a **/
  mol->atoms[a].mark = 0;
 
  return(0);

} /* end of Mark_Neighbors_for_Cycle_with_maxNatom_in_ring() */



void Check_Planarity_of_RING_from_Atom_XYZ(mol)
  struct MOLECULE *mol;
{
 struct RING *rn;
 int o,a,b,c,i,k;
 float OA[3],OB[3],OC[3],Nvec[3],angle,angle_tole,d;

 /*
 printf("#Check_Planarity_of_RING_from_Atom_XYZ()\n"); 
 */
 angle_tole = 8.0;  /** degree **/

 rn = &(mol->HeadRing);
 while (rn->next != NULL){
   rn = rn->next;

   /** (1) Calculate center position (cpos) ***/
   for (k=0;k<3;++k) rn->cpos[k] = rn->nvec[k] = 0.0;
   for (i=0;i<rn->Natom;++i){
     for (k=0;k<3;++k) rn->cpos[k] += mol->atoms[rn->num_atoms[i]].Pos[k];
   } 
   for (k=0;k<3;++k) rn->cpos[k] /= rn->Natom;

  rn->rg = 0.0;  
  for (i=0;i<rn->Natom;++i){
    for (k=0;k<3;++k){
      d = mol->atoms[rn->num_atoms[i]].Pos[k] - rn->cpos[k];
      rn->rg += d*d;
    }
  } 
  rn->rg /= rn->Natom;
  if (rn->rg > 0.0) rn->rg = sqrt(rn->rg); 

  /** (2) Check Planarity ***/
        if (rn->Natom < 4){ rn->planarity = 'N';}
   else if (rn->Natom >= 4){
     o = rn->num_atoms[0];
     a = rn->num_atoms[1];
     b = rn->num_atoms[2];
     Sub_Vec(OA,mol->atoms[a].Pos, mol->atoms[o].Pos);
     Sub_Vec(OB,mol->atoms[b].Pos, mol->atoms[o].Pos);
     Cross_Vec(Nvec,OA,OB);
     
     rn->planarity = 'p';
     for (i=3;i<rn->Natom;++i){
       c = rn->num_atoms[i];
       Sub_Vec(OC,mol->atoms[c].Pos, mol->atoms[o].Pos);
       angle =  Angle_bwn_Vec(Nvec,OC);
       if (fabs(angle-90.0) > angle_tole) rn->planarity = 'N';
/*
       printf("RING %d angle %f diff90 %f plnarity %c\n",rn->num,angle,fabs(angle-90.0),rn->planarity); 
*/
     }
 
   /** (3) Set normal vector(nvec) for "plane" ring  ***/
     if (rn->planarity=='p'){
       for (i=0;i<rn->Natom;++i){
         c = rn->num_atoms[i];
         mol->atoms[c].ring = 'p';
       }
      Normalize_Vec(Nvec);
     for (k=0;k<3;++k) rn->nvec[k] = Nvec[k]; 
     }
 
   }
 }

} /* end of Check_Planarity_of_RING_from_Atom_XYZ() */








int Delete_Disconnected_Parts_Of_Molecule(molA)
  struct MOLECULE *molA;
 /* Leave the largest connected chemical, and delete other connected chemicals */
{
 struct MOLECULE molB;
 int a,b,A,B;
 int *a_from_b,Natom_delete;

 /* printf("#Delete_Disconnected_Parts_Of_Molecule(molA Natom %d)\n",molA->Natom); */
 Malloc_MOLECULE(&molB,molA->Natom,molA->Nbond);
 Copy_MOLECULE(&molB,molA);
 a_from_b = (int *)malloc(sizeof(int)*(molB.Natom+1)); 
 molA->Natom = 0; 
 /** [1] Copy Atoms **/
 Natom_delete = 0; 
 a = 0;
 for (b=0;b<molB.Natom;++b){
   a_from_b[b] = -1;
   if (molB.atoms[b].constr_num==1){
     Copy_ATOM(&(molA->atoms[a]),&(molB.atoms[b]));
     molA->atoms[a].num = a; 
     a_from_b[b] = a;
     a += 1;
     molA->Natom += 1;
   }
   else Natom_delete += 1;
 }


 /** [2] Copy Bonds **/
 molA->Nbond = 0;
 A = 0; 
 for (B=0;B<molB.Nbond;++B){
   if ((a_from_b[molB.bonds[B].anum1]>=0) && (a_from_b[molB.bonds[B].anum2]>=0)){
     molA->bonds[A].anum1    = a_from_b[molB.bonds[B].anum1]; 
     molA->bonds[A].anum2    = a_from_b[molB.bonds[B].anum2]; 
     molA->bonds[A].bondtype = molB.bonds[B].bondtype;
     molA->bonds[A].num = A; 
     A += 1;
     molA->Nbond += 1;     
   } 

 }
 free(a_from_b);

 Free_MOLECULE(&molB);
 return(Natom_delete);
} /* end of Delete_Disconnected_Parts_Of_Molecule() */




void Write_DAMLatom_and_RING_in_PDB(ofname,mol,mode,ChainID,model_num,comment)
 char  *ofname;
 struct MOLECULE *mol;
 char mode;       /* 'w'rite,'a'ppend */
 char ChainID;    /* ChainID for output. '-':use original ChainID */
 int  model_num;  /* number of MODEL */
 char *comment;
{
 int i,k;
 FILE *fpo;
 char atom_str[5],chain;
 float tFactor,pos[3];
 struct RING *rn;
 int    natom;

 printf("#Write_DAMLatom_and_RING_in_PDB(mode '%c')-->'%s'\n",mode,ofname);
 if (mode=='a') fpo = fopen(ofname,"a");
 else           fpo = fopen(ofname,"w");

 fprintf(fpo,"REMARK name_of_molecule '%s'\n",mol->name);
 fprintf(fpo,"REMARK output_filename '%s'\n",mol->filename);
 fprintf(fpo,"REMARK Natom %4d Nbond %4d\n",mol->Natom,mol->Nbond);
 fprintf(fpo,"REMARK COMMAND \"%s\"\n",PAR.COMMAND);
 fprintf(fpo,"REMARK DATE    %s\n",PAR.START_DATE);
 if (comment[0]!='\0') fprintf(fpo,"REMARK COMMENT '%s'\n",comment);
 fprintf(fpo,"MODEL    %5d\n",model_num);
 fprintf(fpo,"REMARK   atomname 'DON':donor,         'ACC':acceptor,          'MIX':mixed donor/acceptor\n");
 fprintf(fpo,"REMARK   atomname 'CAL':carbon_aliphatic\n");
 fprintf(fpo,"REMARK   atomname 'RGC':center_of_ring,'RGU':upper_side_of_ring 'RGL':lower_side_of_ring\n");
 fprintf(fpo,"REMARK                                                  Rg          nvec_x  nvec_y  nvec_z anum atom ele\n");

 natom = 0;
 /*** output DON, ACC, MIX, CAL ***/ 
 for (i=0;i<mol->Natom;++i){
    if (strlen(mol->atoms[i].atomname)!=4){
      sprintf(mol->atoms[i].atomname," %2s ",mol->atoms[i].element);
    }
    tFactor = (float)mol->atoms[i].Nmatch;
  
    if ((mol->atoms[i].DAMRL=='D') || (mol->atoms[i].DAMRL=='A') || (mol->atoms[i].DAMRL=='B') || (mol->atoms[i].DAMRL=='L')){
      if (ChainID=='-') chain = mol->atoms[i].chain; else chain = ChainID;

      if (mol->atoms[i].DAMRL=='D') sprintf(atom_str,"DON");
      if (mol->atoms[i].DAMRL=='A') sprintf(atom_str,"ACC");
      if (mol->atoms[i].DAMRL=='L') sprintf(atom_str,"CAL");
      if (mol->atoms[i].DAMRL=='M') sprintf(atom_str,"MIX");

      natom += 1;
      fprintf(fpo,"HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.3f%6.2f%8.3f%8.3f%8.3f %5d %4s %s %c\n",
          natom,atom_str,mol->atoms[i].resi,chain,mol->atoms[i].rnum,
          mol->atoms[i].Pos[0], mol->atoms[i].Pos[1],mol->atoms[i].Pos[2],0.0, 0.0, 0.0, 0.0, 0.0, 
          mol->atoms[i].num_in_file,mol->atoms[i].atomname, mol->atoms[i].element,mol->atoms[i].DAMRL);

    }
 }
 /*** output RING (RGC,RGU,RGL)***/ 
 rn = &(mol->HeadRing);
 while (rn->next != NULL){
   rn = rn->next;
   if (rn->planarity=='p'){
     natom += 1;
     fprintf(fpo,"HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.3f%6.2f%8.3f%8.3f%8.3f\n",
        natom,"RGC",mol->atoms[0].resi,mol->atoms[0].chain,mol->atoms[0].rnum, 
        rn->cpos[0],rn->cpos[1],rn->cpos[2],rn->rg,0.0,rn->nvec[0],rn->nvec[1],rn->nvec[2]);

     for (k=0;k<3;++k) pos[k] = rn->cpos[k] + rn->rg * rn->nvec[k];
     natom += 1;
     fprintf(fpo,"HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.3f%6.2f\n",
        natom,"RGU",mol->atoms[0].resi,mol->atoms[0].chain,mol->atoms[0].rnum,pos[0],pos[1],pos[2],rn->rg,0.0);

     for (k=0;k<3;++k) pos[k] = rn->cpos[k] - rn->rg * rn->nvec[k];
     natom += 1;
     fprintf(fpo,"HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.3f%6.2f\n",
        natom,"RGL",mol->atoms[0].resi,mol->atoms[0].chain,mol->atoms[0].rnum,pos[0],pos[1],pos[2],rn->rg,0.0);

   }
 }

 fprintf(fpo,"TER\n");
 fprintf(fpo,"ENDMDL\n");
 fclose(fpo);

} /* end of Write_DAMLatom_and_RING_in_PDB() */






int Mark_Neighbors_for_Connected(mol,a,Nconstr)
  struct MOLECULE *mol;
  int a;        /* index of the focused atom */
  int Nconstr;  /* number of connected graph */
{
 int b;

 mol->atoms[a].constr_num = Nconstr;

 for (b=0;b<mol->Natom;++b){
  if ((a!=b) && (mol->conmap.map[a][b] != '0')){
    if (mol->atoms[b].constr_num==0) Mark_Neighbors_for_Connected(mol,b,Nconstr); 
    }

 } /* b */

  return(1);

} /* end of Mark_Neighbors_for_Connected() */



void Sub_Vec(C,A,B)    /* C = A - B */
 float C[3],A[3],B[3];
{  C[0] = A[0]-B[0];
   C[1] = A[1]-B[1];
   C[2] = A[2]-B[2]; } 


void Cross_Vec(C,A,B)   /* C = A x B */
 float C[3],A[3],B[3];
{ C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0]; }


void Normalize_Vec(A)
 float A[3];
{
 float norm;
 norm = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
 norm = sqrt(norm);
 if (norm>0.0){
   A[0] /= norm; A[1] /= norm; A[2] /= norm;
 }
}

float Angle_bwn_Vec(A,B)
 float A[3],B[3];
{ float dot_prod,normA,normB,cosAB;
  normA    = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  normB    = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
  dot_prod = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  cosAB = dot_prod/sqrt(normA)/sqrt(normB);
  if (cosAB < 0.0) cosAB *= -1.0;
  return(acos(cosAB)*180.0/M_PI); 
}



void Read_DAMLring_Molecule_in_PDB(ifname,mol)
 char *ifname;
 struct MOLECULE *mol;
{
 FILE *fp;
 char line[256];
 char end,B[128],atmhet;
 int Nline,Lline,natom,nbond;

 printf("#Read_DAMLring_Molecule_in_PDB('%s')\n",ifname);
/*
HETATM    1  ACC CEN X   1       4.288   9.409   3.956 1.870 86.36   1.000   0.000   0.000
HETATM    2  RGC CEN X   1       5.446   9.596   2.488 1.312 63.64   0.111  -0.988  -0.107
HETATM    2  RGU CEN X   1       5.592   8.300   2.347 1.312 63.64    14 1.312
HETATM    2  RGL CEN X   1       5.301  10.892   2.629 1.312 63.64    14 1.312
HETATM    3  RGC CEN X   1       6.980  11.001  -3.368 1.417 59.09  -0.969  -0.142  -0.205
HETATM    3  RGU CEN X   1       5.608  10.800  -3.658 1.417 59.09    13 1.417
HETATM    3  RGL CEN X   1       8.353  11.202  -3.078 1.417 59.09    13 1.417
HETATM    4  RGC CEN X   1       6.011   9.512   4.252 1.393 54.55   0.077  -0.992  -0.100
HETATM    4  RGU CEN X   1       6.119   8.130   4.112 1.393 54.55    12 1.393
HETATM    4  RGL CEN X   1       5.903  10.894   4.391 1.393 54.55    12 1.393
HETATM    5  CAL CEN X   1       6.332  13.317  -1.728 1.870 50.00   1.000   0.000   0.000
HETATM    6  DON CEN X   1       7.713   8.715  -5.922 1.870 40.91   1.000   0.000   0.000
HETATM    7  ACC CEN X   1       8.761  10.607  -0.325 1.870 40.91   1.000   0.000   0.000
0123456789012345678901234567890123456789
          1         2
*/


 sprintf(mol->filename,"%s",ifname);
 end = 0;
 Nline = natom = nbond = 0;
 mol->name[0] = '\0';
 sprintf(mol->name,"%s",ifname);
 mol->Natom = mol->Nheavyatom = 0;
 mol->conmap.malloced = mol->dismap.malloced = 0;
 mol->Nmatch = 0;

 

 /*** (1) Read PDB file only to know Number of atoms ***/
 fp = fopen(ifname,"r");
 if (fp==NULL){printf("#ERROR:Can't open pdbfile '%s'\n",ifname); exit(1);}
 natom = nbond = 0;
 while ((feof(fp)==0)&&(end==0)){
   fgets(line,255,fp);
   Lline = strlen(line);
   if (line[Lline-1]=='\n'){line[Lline-1]='\0'; Lline = Lline -1;}
   if ((strncmp(line,"ATOM  ",6)==0)||(strncmp(line,"HETATM",6)==0)){
     Get_Part_Of_Line(B,line,13,15);
     if ((strcmp(B,"DON")==0)||(strcmp(B,"ACC")==0)||(strcmp(B,"MIX")==0)|| (strcmp(B,"CAL")==0)||
         (strcmp(B,"RGC")==0)||(strcmp(B,"RGU")==0)||(strcmp(B,"RGL")==0)){ natom += 1; } 
   }
 }
 fclose(fp);
  printf("#natom %d nbond %d\n",natom,nbond); 
 Malloc_MOLECULE(mol,natom,0);

/*** (2) Read PDB file ***/
 fp = fopen(ifname,"r");
 if (fp==NULL){printf("#ERROR:Can't open pdbfile '%s'\n",ifname); exit(1);}
 mol->Natom = mol->Nbond = 0;
 natom = 0;
 while ((feof(fp)==0)&&(end==0)){
   line[0] = '\0';
   fgets(line,255,fp);
   Lline = strlen(line);
   if (line[Lline-1]=='\n'){line[Lline-1]='\0'; Lline = Lline -1;}

   /** 'ATOM' or 'HETATM' lines **/
   Get_Part_Of_Line(B,line,13,15);
   if ((strncmp(line,"ATOM  ",6)==0)||(strncmp(line,"HETATM",6)==0)){
     if ((strcmp(B,"DON")==0)||(strcmp(B,"ACC")==0)||(strcmp(B,"MIX")==0)|| (strcmp(B,"CAL")==0)||
         (strcmp(B,"RGC")==0)||(strcmp(B,"RGU")==0)||(strcmp(B,"RGL")==0)){ 
         /* printf("'%s''%s' natom %d\n",line,B,natom);  */
         mol->atoms[natom].num = natom;
         if (strncmp(line,"ATOM  ",6)==0) atmhet = 'A';
         else if (strncmp(line,"HETATM",6)==0) atmhet = 'H';
         else atmhet = 'A';

         mol->atoms[natom].atmhet  = atmhet;
         Get_Part_Of_Line(B,line,6,10);  mol->atoms[natom].num_in_file = atoi(B);
         Get_Part_Of_Line(B,line,30,37); mol->atoms[natom].Pos[0] = atof(B);
         Get_Part_Of_Line(B,line,38,45); mol->atoms[natom].Pos[1] = atof(B);
         Get_Part_Of_Line(B,line,46,53); mol->atoms[natom].Pos[2] = atof(B);
         Get_Part_Of_Line(mol->atoms[natom].atomname,line,13,15);
         Get_Part_Of_Line(mol->atoms[natom].resi,line,17,19);
         Get_Part_Of_Line(mol->atoms[natom].rnum,line,22,26);
         mol->atoms[natom].chain = line[21];
         if (Lline>64){
           Get_Part_Of_Line(B,line,60,65); mol->atoms[natom].tFactor = atof(B);
         }
         else mol->atoms[natom].tFactor = -1.0;

         mol->atoms[natom].orbital = '-';
         mol->atoms[natom].DAMRL   = '-';
         Get_Part_Of_Line(B,line,13,15);
              if (strcmp(B,"DON")==0) mol->atoms[natom].DAMRL   = 'D';
         else if (strcmp(B,"ACC")==0) mol->atoms[natom].DAMRL   = 'A';
         else if (strcmp(B,"MIX")==0) mol->atoms[natom].DAMRL   = 'M';
         else if (strcmp(B,"CAL")==0) mol->atoms[natom].DAMRL   = 'L';
         else if (strcmp(B,"RGC")==0) mol->atoms[natom].DAMRL   = 'R';
         else if (strcmp(B,"RGL")==0) mol->atoms[natom].DAMRL   = '-';
         else if (strcmp(B,"RGU")==0) mol->atoms[natom].DAMRL   = '-';
         natom += 1;
     }
   }
 } /* while */
 
 fclose(fp);
 mol->Natom = natom;

} /* end of Read_DAMLring_Molecule_in_PDB() */


int Number_of_RING(HeadRingNode)
 struct RING *HeadRingNode;
{
  struct RING *rn;
  int N;

  rn = HeadRingNode;
  N = 0;
  while (rn->next != NULL){ 
    rn = rn->next;
    N += 1;
  }
  return(N);
}





void Free_RING(mol)
 struct MOLECULE *mol;
{
 struct RING *rn;

 rn = &(mol->HeadRing);
 if (rn->next != NULL){ 
   while (rn->next != NULL) rn = rn->next;
   while (rn->prev!=NULL){
     rn = rn->prev;
     if (rn->next != NULL){
       if (rn->next->num_atoms!=NULL) free(rn->next->num_atoms); 
       free(rn->next);
     }
   }
   mol->HeadRing.next = NULL;
 } 

} /* end of Free_RING() */


void Free_BLOCK(mol)
 struct MOLECULE *mol;
{
 struct RING *rn;

 rn = &(mol->HeadBlock);
 if (rn->next != NULL){ 
   while (rn->next != NULL) rn = rn->next;

   while (rn->prev!=NULL){
     rn = rn->prev;
     free(rn->next->num_atoms); 
     free(rn->next);
   }
   mol->HeadBlock.next = NULL;
 }
} /* end of Free_BLOCK() */


void Malloc_Copy_RINGs(HeadRingN,HeadRingO)
  struct RING *HeadRingN; /* new */
  struct RING *HeadRingO; /* orig */
{
/*
struct RING{
  int   num; 
  int   Natom;
  int   *num_atoms;
  char  planarity; 
  char  aromatic;
  float cpos[3];
  float nvec[3];
  float rg; 
  struct RING *prev, *next; 
};
*/

  struct RING *nr,*or;
  int i;

  printf("#Malloc_Copy_RINGs(HeadRingN,HeadRingO)\n");
  HeadRingN->next = NULL;
  HeadRingN->prev = NULL;
  HeadRingN->num  = -1;
  HeadRingN->Natom = 0;
  HeadRingN->num_atoms = NULL;
  or = HeadRingO;
  nr = HeadRingN;
  while (or->next != NULL){
    or = or->next;
    nr->next = (struct RING *)malloc(sizeof(struct RING));
    nr->next->prev = nr;
    nr = nr->next;
    nr->next = NULL;
    nr->num = or->num;
    nr->num_atoms = (int *)malloc(sizeof(int)*or->Natom); 
    nr->Natom = or->Natom;
    nr->planarity = or->planarity;
    nr->aromatic  = or->aromatic;
    nr->rg        = or->rg;
    for (i=0;i<or->Natom;++i){
      nr->num_atoms[i] = or->num_atoms[i];
    }
    for (i=0;i<3;++i){
      nr->cpos[i] = or->cpos[i];
      nr->nvec[i] = or->nvec[i];
    }
  }

} /* end of Malloc_Copy_RINGs() */




void Set_Rotatable_Bond(mol)
  struct MOLECULE *mol;
/*
 Set_Ring_Structure() and Set_Nneighbor_CHNOSP() should be peformed before execuing this function;

 mol->conmap.map[i][j] = '0' : not bonded
                       = 'R' : rotatable bond      (bond not within the ring) 
                       = 'N' : not rotatable bond  (bond within the ring)

*/
{
  int i,j,ii,jj;
  struct RING *rn;
  /* printf("#Set_Rotatable_Bond()\n");  */
  for (i=0;i<mol->Natom;++i){
    for (j=i+1;j<mol->Natom;++j){
      if (mol->conmap.map[i][j]!='0'){
         mol->conmap.map[i][j] =  mol->conmap.map[j][i] = 'R';
      }
    } 
  }

 /** Bond within the ring --> 'nonrotatble' **/
  rn = &(mol->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;
    for (i=0;i<rn->Natom;++i){
      ii = rn->num_atoms[i]; 
      for (j=i+1;j<rn->Natom;++j){
        jj = rn->num_atoms[j]; 
        if (mol->conmap.map[ii][jj]!='0') { 
          mol->conmap.map[ii][jj] = mol->conmap.map[jj][ii] = 'N';
        }
      }
    }
  }

 /*
  for (i=0;i<mol->Natom;++i){
    for (j=i+1;j<mol->Natom;++j){ 
      if (mol->conmap.map[i][j]!='0')
        printf("BOND %d %s %d %s %c\n",mol->atoms[i].num_in_file,mol->atoms[i].atomname,mol->atoms[j].num_in_file,mol->atoms[j].atomname,mol->conmap.map[i][j]);
    }
  }
  */

} /* end of Set_Rotatable_Bond() */
