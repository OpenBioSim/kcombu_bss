/*

< transform.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


functions to transform molecular conformation.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "globalvar.h"
#include "molprop.h"
#include "transform.h"
#include "PCAfit.h"
#include "energies.h"
#include "vector3.h"

/*** FUNCTIONS (GLOBAL) ***/
int Set_TRANSFORM_VARIABLE();
void Malloc_TRANSFORM_VARIABLE();
void Free_TRANSFORM_VARIABLE();
void Copy_TRANSFORM_VARIABLE();
void Rotate_TargetMolecule_by_PCA();
void Rotate_TargetMolecule_by_minimum_energy_PCA();
void Write_TRANSFORM_VARIABLE();
void Bubble_Sort_Array_of_TRANSFORM_VARIABLE();
void Set_Type_of_Rotatable_Bond();

/*** FUNCTIONS (LOCAL) ***/
static char determine_SpType();


int Set_TRANSFORM_VARIABLE(tra,mol,SetRotBond)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *mol;
  char   SetRotBond; /* 'R':set rotatable bond. otherwise not */
{
  int a,b,c,r,Nchange;
  int *toposum;
  int min_toposum,min_a;
  /**
    NOTE !: Set_Rotatable_Bond() should be executed before this function.
  ***/
  printf("#Set_TRANSFORM_VARIABLE(tra,mol,SetRotBond='%c')\n",SetRotBond);

  /** [1] Set Center Atom **/
  toposum = (int *)malloc(sizeof(int)*mol->Natom);
  min_toposum = -1;
  min_a = 0; 
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].one_char_ele != 'H'){
      toposum[a] = 0;
      for (b=0;b<mol->Natom;++b){
        if (mol->atoms[b].one_char_ele != 'H'){
          toposum[a] += mol->topodismap.map[a][b];
        }
      }

    if ((min_toposum<0) || (toposum[a]<min_toposum)) {min_toposum = toposum[a]; min_a = a;}
    /* printf("%d %s toposum %d\n",mol->atoms[a].num_in_file,mol->atoms[a].atomname,toposum[a]); */
    } 
  } 
  tra->Canum = min_a;

  /*
  printf("min_toposum %d min_a %d\n",min_toposum,mol->atoms[min_a].num_in_file);
  */


  /** [2] Count Rotatable Bond(O,Z,M) **/
  tra->Nrbond = 0; 
  if (SetRotBond=='R'){ 
   for (a=0;a<mol->Natom;++a){
      for (b=a+1;b<mol->Natom;++b){
        if ((mol->conmap.map[a][b]=='R')&&
             ( ( (PAR.HydrogenOptimize!='T')&&(mol->atoms[a].Nnei_heavy>1)&&(mol->atoms[b].Nnei_heavy>1) )  ||
               ( (PAR.HydrogenOptimize=='T')&&(mol->atoms[a].Nneighbor >1)&&(mol->atoms[b].Nneighbor >1) ) ) ){ 
          tra->Nrbond += 1;
        } 
      } 
    }
  }

  /** [3] Malloc Variables **/
  Malloc_TRANSFORM_VARIABLE(tra,tra->Nrbond,mol->Natom);

  if (SetRotBond!='R'){
    return(1);
  }

  /** [4] Set Rotatable Bond **/
  r = 0;  
  for (a=0;a<mol->Natom;++a){
    for (b=a+1;b<mol->Natom;++b){
      if ((mol->conmap.map[a][b]=='R')&&
          ( ((PAR.HydrogenOptimize!='T')&&(mol->atoms[a].Nnei_heavy>1)&&(mol->atoms[b].Nnei_heavy>1)) ||
            ((PAR.HydrogenOptimize=='T')&&(mol->atoms[a].Nneighbor >1)&&(mol->atoms[b].Nneighbor >1)) ) ){

          if (mol->topodismap.map[a][tra->Canum]<mol->topodismap.map[b][tra->Canum]){
            tra->rbZanum[r] = a; tra->rbOanum[r] = b;
          } 
          else{
            tra->rbZanum[r] = b; tra->rbOanum[r] = a;
          }
          for (c=0;c<mol->Natom;++c){
            if ((c!=a) && (c!=b) && 
                (mol->topodismap.map[c][tra->rbOanum[r]] < mol->topodismap.map[c][tra->rbZanum[r]]) ){
              tra->rbManum[r][c] = 1;
           } 
            else {tra->rbManum[r][c] = 0; }
          } 
          r += 1; 
      } 
    } 
  }

  /** [5] Bubble Sort Rotatable Bonds by the distance from the center **/
  do{ 
    Nchange = 0; 
    for (r=0;r<(tra->Nrbond-1);++r){
      if (mol->topodismap.map[tra->rbOanum[r]][tra->Canum] > mol->topodismap.map[tra->rbOanum[r+1]][tra->Canum]){
         b = tra->rbZanum[r]; 
         tra->rbZanum[r] = tra->rbZanum[r+1];
         tra->rbZanum[r+1] = b;
         b = tra->rbOanum[r]; 
         tra->rbOanum[r] = tra->rbOanum[r+1];
         tra->rbOanum[r+1] = b;
         for (a=0;a<mol->Natom;++a){
           b = tra->rbManum[r][a];
           tra->rbManum[r][a] = tra->rbManum[r+1][a];
           tra->rbManum[r+1][a] = b;
         }
         Nchange += 1;
      }
    }
  }  
  while (Nchange>0);

  printf("#CenterAtom %d %s Nrotatable_bond %d\n", mol->atoms[tra->Canum].num_in_file, mol->atoms[tra->Canum].atomname, tra->Nrbond);

/* 
  for (r=0;r<tra->Nrbond;++r){
    printf("#RBOND[%d] Z %s%d O %s%d\n",r,
mol->atoms[tra->rbZanum[r]].element, mol->atoms[tra->rbZanum[r]].num_in_file,
mol->atoms[tra->rbOanum[r]].element, mol->atoms[tra->rbOanum[r]].num_in_file);
  }
*/




  free(toposum);
  return(1);
} /* end of Set_TRANSFORM_VARIABLE() */






void Malloc_TRANSFORM_VARIABLE(tra,Nrbond,Natom)
  struct TRANSFORM_VARIABLE *tra;
  int Nrbond;
  int Natom;
{
  int a,r;

 /*
  printf("#Malloc_TRANSFORM_VARIABLE(tra,Nrbond:%d Natom:%d )\n",Nrbond,Natom);
*/
  tra->Natom         = Natom;
  tra->Nrbond        = Nrbond;
  tra->Nrbond_stamp  = 0;
  if (Nrbond>0){
    tra->rbType      = (char *)malloc(sizeof(char)*Nrbond);
    tra->rbStamp     = (char *)malloc(sizeof(char)*Nrbond);
    tra->rbPlane     = (char *)malloc(sizeof(char)*Nrbond);
    tra->rbAngle     = (float *)malloc(sizeof(float)*Nrbond);
    tra->dEd_rbAngle = (float *)malloc(sizeof(float)*Nrbond);
    tra->rbOanum     = (int   *)malloc(sizeof(int)*Nrbond);
    tra->rbZanum     = (int   *)malloc(sizeof(int)*Nrbond);

    tra->rbManum = (unsigned char **)malloc(sizeof(unsigned char*)*Nrbond);
    for (r=0;r<Nrbond;++r){
      tra->rbManum[r] = (unsigned char *)malloc(sizeof(unsigned char)*Natom);
      tra->rbAngle[r] = 0.0;
      tra->rbType[r]  = tra->rbStamp[r] = tra->rbPlane[r] = ' ';
      for (a=0;a<Natom;++a){tra->rbManum[r][a] = 0;}
    }
  }

  tra->force_atom = (float **)malloc(sizeof(float *)*Natom);
  for (a=0;a<Natom;++a){
    tra->force_atom[a] = (float *)malloc(sizeof(float)*3);
  }

  tra->atom_match.anumA = NULL;
  tra->atom_match.anumB = NULL;


} /* end of Malloc_TRANSFORM_VARIABLE() */



void Free_TRANSFORM_VARIABLE(tra)
  struct TRANSFORM_VARIABLE *tra;
{
  int r;

  for (r=0;r<tra->Nrbond;++r){
    free(tra->rbManum[r]);
  }
  free(tra->rbManum);
  free(tra->rbZanum);
  free(tra->rbOanum);
  free(tra->rbAngle);
} /* end of Free_TRANSFORM_VARIABLE() */



void Copy_TRANSFORM_VARIABLE(traN,traO)
  struct TRANSFORM_VARIABLE *traN,*traO;
{
  int a,r;

  traN->Natom        = traO->Natom;
  traN->Nrbond       = traO->Nrbond;
  traN->Nrbond_stamp = traO->Nrbond_stamp;
  traN->Canum        = traO->Canum;
  for (r=0;r<traO->Nrbond;++r){
    traN->rbStamp[r] = traO->rbStamp[r];
    traN->rbPlane[r] = traO->rbPlane[r];
    traN->rbOanum[r] = traO->rbOanum[r];
    traN->rbZanum[r] = traO->rbZanum[r];
    for (a=0;a<traO->Natom;++a) traN->rbManum[r][a] = traO->rbManum[r][a]; 
  }


  traN->Etotal      =  traO->Etotal;
  traN->Etotal_init =  traO->Etotal_init;
  traN->Eatommatch  =  traO->Eatommatch;
  traN->Eselfclash  =  traO->Eselfclash;
  traN->Eprotclash  =  traO->Eprotclash;
  traN->Evolmovlap  =  traO->Evolmovlap;
  traN->rmsd_match  =  traO->rmsd_match;
  traN->tanimoto_volume  =  traO->tanimoto_volume;
  
  traN->confnumT  =  traO->confnumT;
  traN->confnumR  =  traO->confnumR;
  traN->Ninit   =  traO->Ninit;
  traN->num_atom_match   =  traO->num_atom_match;

} /* end of Copy_TRANSFORM_VARIABLE() */





void Rotate_TargetMolecule_by_PCA(molT,molR,molP,rot_num,M)
  struct MOLECULE *molT; /* target molecule */
  struct MOLECULE *molR; /* reference molecule */
  struct MOLECULE *molP; /* receptor protein molecule */
  int    rot_num;        /* {0,1,2,3} */
  struct MATCH    *M;  
{
  char sign[4][4];
  int   a,i,j;
  float Rmat[3][3],evecB[3][3];

 /*
  printf("#Rotate_TargetMolecule_by_minimum_energy_PCA(rot_num:%d)\n",rot_num); 
 */

  Cal_Molecule_EigVec_Covariance_Matrix(molT);
  Cal_Molecule_EigVec_Covariance_Matrix(molR);

  Set_Center_of_Gravity(molT);
  Set_Center_of_Gravity(molR);

  sign[0][0] =  1; sign[0][1] =  1; sign[0][2] =  1;  /* 0-degree rotation          */
  sign[1][0] =  1; sign[1][1] = -1; sign[1][2] = -1;  /* x-axis 180-degree rotation */
  sign[2][0] = -1; sign[2][1] =  1; sign[2][2] = -1;  /* y-axis 180-degree rotation */
  sign[3][0] = -1; sign[3][1] = -1; sign[3][2] =  1;  /* z-axis 180-degree rotation */


 /** [2] Rotate molT by the rot_num rotation **/
 for (a=0;a<molT->Natom;++a){
   for (i=0;i<3;++i){
     molT->atoms[a].Pos[i] += (molR->Gcen[i]-molT->Gcen[i]);
   }
 }
 Set_Center_of_Gravity(molT);

 /* printf("#opt_k %d Eopt %f\n",opt_k,E[opt_k]); */
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){evecB[i][j] = sign[rot_num][i] * molR->eigvec[i][j];}
 }

 Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat,
   molT->eigvec[0],molT->eigvec[1],molT->eigvec[2],evecB[0],evecB[1],evecB[2]);

 Rotate_Molecule_around_Gcenter(molT,Rmat,molT->Gcen,molR->Gcen);

} /* end of Rotate_MoleculeA_by_PCA() */




void Rotate_TargetMolecule_by_minimum_energy_PCA(molT,molR,molP,M)
  struct MOLECULE *molT; /* target molecule */
  struct MOLECULE *molR; /* reference molecule */
  struct MOLECULE *molP; /* receptor protein molecule */
  struct MATCH    *M;  
{
  struct MOLECULE mA1;
  char sign[4][4];
  float E[4], Emin;
  int   opt_k,k,a,i,j;
  float Rmat[3][3],evecB[3][3];

  /*
  printf("#Rotate_TargetMolecule_by_minimum_energy_PCA()\n"); 
  */

  Cal_Molecule_EigVec_Covariance_Matrix(molT);
  Cal_Molecule_EigVec_Covariance_Matrix(molR);

  Set_Center_of_Gravity(molT);
  Set_Center_of_Gravity(molR);

  Initialize_MOLECULE(&mA1);
  Malloc_MOLECULE(&mA1,molT->Natom,molT->Nbond);
  Copy_MOLECULE(&mA1,molT);

  sign[0][0] =  1; sign[0][1] =  1; sign[0][2] =  1;  /* 0-degree rotation          */
  sign[1][0] =  1; sign[1][1] = -1; sign[1][2] = -1;  /* x-axis 180-degree rotation */
  sign[2][0] = -1; sign[2][1] =  1; sign[2][2] = -1;  /* y-axis 180-degree rotation */
  sign[3][0] = -1; sign[3][1] = -1; sign[3][2] =  1;  /* z-axis 180-degree rotation */

  Emin = 0.0;
  opt_k = 0;

  for (k=0;k<4;++k){

    Copy_MOLECULE(&mA1,molT);

    for (a=0;a<mA1.Natom;++a){
      for (i=0;i<3;++i){
        mA1.atoms[a].Pos[i] += (molR->Gcen[i]-molT->Gcen[i]);
      }
    }

    Set_Center_of_Gravity(&mA1);


   /*
     printf("#kn %d %f %f %f\n",0,mA1.atoms[0].Pos[0], mA1.atoms[0].Pos[1], mA1.atoms[0].Pos[2]);
     printf("#kn %d %f %f %f\n",1,mA1.atoms[1].Pos[0], mA1.atoms[1].Pos[1], mA1.atoms[1].Pos[2]);
    */

   for (i=0;i<3;++i){
     for (j=0;j<3;++j){ evecB[i][j] = sign[k][i] * molR->eigvec[i][j];}
   }

    Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat, molT->eigvec[0],molT->eigvec[1],molT->eigvec[2],evecB[0],evecB[1],evecB[2]);

    Rotate_Molecule_around_Gcenter(&mA1,Rmat,mA1.Gcen,molR->Gcen);

    E[k] = energy_total(&mA1,molR,molP,M);
    /* printf("#rot_num %d E %f\n",k,E[k]); */
    if ((k==0)||(E[k]<Emin)) { Emin = E[k]; opt_k = k;}
 }

 /** [2] Rotate molT by the opt_k rotation **/
 for (a=0;a<molT->Natom;++a){
   for (i=0;i<3;++i){
     molT->atoms[a].Pos[i] += (molR->Gcen[i]-molT->Gcen[i]);
   }
 }
 Set_Center_of_Gravity(molT);

 /* printf("#opt_k %d Eopt %f\n",opt_k,E[opt_k]); */
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){evecB[i][j] = sign[opt_k][i] * molR->eigvec[i][j];}
 }

 Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat,
   molT->eigvec[0],molT->eigvec[1],molT->eigvec[2],evecB[0],evecB[1],evecB[2]);

 Rotate_Molecule_around_Gcenter(molT,Rmat,molT->Gcen,molR->Gcen);
 Free_MOLECULE(&mA1);

} /* end of Rotate_TargetMolecule_by_minimum_energy_PCA() */



void Write_TRANSFORM_VARIABLE(ofname,tra,mol)
  char *ofname;
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *mol;
{
  FILE *fp;
  int r;

  printf("#Write_TRANSFORM_VARIABLE()-->'%s'\n",ofname);
  fp = fopen(ofname,"w");
  fprintf(fp,"#TRANSFORM_VARIABLE\n");
  fprintf(fp,"COMMAND %s\n",PAR.COMMAND);
  fprintf(fp,"DATE    %s\n",PAR.START_DATE);
  fprintf(fp,"Natom %d\n",tra->Natom);
  fprintf(fp,"Canum %s %d\n",mol->atoms[tra->Canum].element,mol->atoms[tra->Canum].num_in_file);
  fprintf(fp,"Nrbond       %d\n",tra->Nrbond);
  fprintf(fp,"Nrbond_stamp %d\n",tra->Nrbond_stamp);
  for (r=0;r<tra->Nrbond;++r){
    fprintf(fp,"BOND %2d [O] %s%4d [Z] %s%4d Type %c Stamp %c Plane %c\n",r,
     mol->atoms[tra->rbOanum[r]].element, mol->atoms[tra->rbOanum[r]].num_in_file,
     mol->atoms[tra->rbZanum[r]].element, mol->atoms[tra->rbZanum[r]].num_in_file,
     tra->rbType[r],tra->rbStamp[r],tra->rbPlane[r]);
  } 
  fclose(fp);

} /* end of Write_TRANSFORM_VARIABLE() */






void Bubble_Sort_Array_of_TRANSFORM_VARIABLE(N, sindex, traTtrial,SortValue)
  int N;         /* number of data */
  int *sindex;   /* sindex[N] : sorted index by increasing order of energy[] ( to be calculated ) */
  struct TRANSFORM_VARIABLE *traTtrial;
  char SortValue;  /* 'E':Etootal,'V':tanimoto_volume */
{
  int i,buff,Nchange;

  for (i=0;i<N;++i) sindex[i] = i;

  do{
    Nchange = 0;
    for (i=1;i<N;++i){
      if (((SortValue=='E')&&(traTtrial[sindex[i-1]].Etotal>traTtrial[sindex[i]].Etotal)) ||
          ((SortValue=='V')&&(traTtrial[sindex[i-1]].tanimoto_volume<traTtrial[sindex[i]].tanimoto_volume)) ){
        buff = sindex[i-1];
        sindex[i-1] = sindex[i];
        sindex[i] = buff;
        Nchange += 1;
      }
    }
  } while (Nchange>0);

} /* end of Bubble_Sort_Array_of_TRANSFORM_VARIABLE() */




void Set_Type_of_Rotatable_Bond(tra,mol,M,PlaneDecision)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *mol;
  struct MATCH *M;
  char   PlaneDecision;  /* '3'D-criteria, '2'D-criteria. (2D requires hydrogens).*/
{
  int a,r,m,z,o,nz,no,i;
  int Mmatch, nonMmatch;
  float bangZ,bangO,tole_bond_ang,tole_torsion_ang;
  char sptypeZ,sptypeO;

  printf("#Set_Type_of_Rotatable_Bond(PlaneDecision:%c)\n",PlaneDecision);

  /** [1] atom[a].mark for matched atom a **/
  for (a=0;a<mol->Natom;++a){
    mol->atoms[a].mark = 0;
  }

  for (m=0;m<M->Npair;++m){
    mol->atoms[M->anumA[m]].mark = 1;
  }


  /** [2] Count mathced atoms for rbManum for each bond **/
  for (r=0;r<tra->Nrbond;++r){
    Mmatch = nonMmatch = 0;
    for (a=0;a<mol->Natom;++a){
      if ((mol->atoms[a].mark==1) && (a!=tra->rbOanum[r])){
        if (tra->rbManum[r][a]==1){ Mmatch += 1;}
                              else{ nonMmatch += 1;}
      }
    }

    if ((Mmatch>0) && (nonMmatch>0)){ tra->rbType[r] = 'M'; }
                                else{ tra->rbType[r] = 'F'; }
/*
    printf("RBOND[%2d] Z %2d O %2d Mmatch %d nonMmatch %d type %c\n",r,
    mol->atoms[tra->rbZanum[r]].num_in_file,mol->atoms[tra->rbOanum[r]].num_in_file,Mmatch,nonMmatch,tra->rbType[r]);
  */

  }

if (PlaneDecision == '3'){

 /** [3D-plane] Set rbPlane[r] type for each rotational bond by 3D structure **/
 /**  >> Conditions for "Plane" bond << */
 /**   (1) dihedral angle  is close to 0 or 180 degree.  */
 /**   (2) Two bond angles are close to 120 degree.      */

   tole_torsion_ang  = 10.0; /* degree */
   tole_bond_ang     = 12.0; /* degree */

   for (r=0;r<tra->Nrbond;++r){
     if ((mol->atoms[tra->rbZanum[r]].Nneighbor>0)&&(mol->atoms[tra->rbZanum[r]].Nneighbor>0)){
  /*
     set 'no' and 'nz
        no        nz
          \      /
           o--- z
   */
        z = tra->rbZanum[r];
        o = tra->rbOanum[r];

        no = nz = -1;
        i = 0;
        while ((i<mol->atoms[z].Nneighbor)&&(nz<0)){
          if (o != mol->atoms[z].neighbors[i]){nz = mol->atoms[z].neighbors[i];} else {i += 1;}
        }

        i = 0;
        while ((i<mol->atoms[o].Nneighbor)&&(no<0)){
          if (z != mol->atoms[o].neighbors[i]) {no = mol->atoms[o].neighbors[i];} else {i += 1;}
        }
        tra->rbAngle[r] = 180.0*dihedral_angle(mol->atoms[nz].Pos, mol->atoms[z].Pos, mol->atoms[o].Pos, mol->atoms[no].Pos)/M_PI;
        bangZ = 180.0*bond_angle(mol->atoms[nz].Pos, mol->atoms[z].Pos, mol->atoms[o].Pos)/M_PI;
        bangO = 180.0*bond_angle(mol->atoms[z].Pos,  mol->atoms[o].Pos, mol->atoms[no].Pos)/M_PI;
        if ( ((fabs(tra->rbAngle[r]-0.0)<tole_torsion_ang)||(fabs(tra->rbAngle[r]-180.0)<tole_torsion_ang)) &&
              (fabs(bangZ-120.0)<tole_bond_ang)&&(fabs(bangO-120.0)<tole_bond_ang)){
           tra->rbPlane[r] = 'P';
        }

  /*
 *         printf("#RBOND[%2d] Z:%s%2d O:%s%2d type %c angle %+6.1f bangZ %+6.1f bangO %+6.1f Plane '%c'\n",r,
 *                    mol->atoms[tra->rbZanum[r]].element, mol->atoms[tra->rbZanum[r]].num_in_file,
 *                               mol->atoms[tra->rbOanum[r]].element, mol->atoms[tra->rbOanum[r]].num_in_file,
 *                                          tra->rbType[r],tra->rbAngle[r],
 *                                                     bangZ,bangO,tra->rbPlane[r]);
 *                                                        */

      }
   }
 }



 else if (PlaneDecision == '2'){
 /** [2D-plane] Set rbPlane[r] type for each rotational bond by 2D structure with Hydrogens **/
 /**  >> Conditions for "Plane" bond << */
 /**   (1) atom Z is not 'sp3' and azom O is not 'sp3'.  */
 /**    OR,  **/
 /**   (2) peptide bond   */

   for (r=0;r<tra->Nrbond;++r){
     z = tra->rbZanum[r];
     o = tra->rbOanum[r];
     sptypeZ = '-';
     sptypeO = '-';
     sptypeZ = determine_SpType(mol,z);
     sptypeO = determine_SpType(mol,o);
/* (1) atom Z is not 'sp3' and azom O is not 'sp3' --> planer  */
     if ((sptypeZ!=3) && (sptypeO !=3)){ tra->rbPlane[r] = 'P'; }

/* (2) For Peptide Bond --> Planar
           H
           |
           CA--C--N--CA
              ||
              O
*/
     if ((mol->atoms[z].one_char_ele=='C') && (sptypeZ=='2') && (mol->atoms[z].NneiC==1) && (mol->atoms[z].NneiO==1) &&
         (mol->atoms[o].one_char_ele=='N') && (sptypeO=='3') && (mol->atoms[z].NneiH==1) ) { tra->rbPlane[r] = 'P';}

     if ((mol->atoms[o].one_char_ele=='C') && (sptypeO=='2') && (mol->atoms[o].NneiC==1) && (mol->atoms[o].NneiO==1) &&
         (mol->atoms[z].one_char_ele=='N') && (sptypeZ=='3') && (mol->atoms[z].NneiH==1) ) { tra->rbPlane[r] = 'P';}

   }
 }

  for (r=0;r<tra->Nrbond;++r){
   printf("#RBOND[%2d] [Z] %s %s %2d %s [O] %s %s %2d %s type '%c' angle %f rbPlane '%c'\n",r,
     mol->atoms[tra->rbZanum[r]].rnum,
     mol->atoms[tra->rbZanum[r]].resi,
     mol->atoms[tra->rbZanum[r]].num_in_file,
     mol->atoms[tra->rbZanum[r]].atomname,
     mol->atoms[tra->rbOanum[r]].rnum,
     mol->atoms[tra->rbOanum[r]].resi,
     mol->atoms[tra->rbOanum[r]].num_in_file,
     mol->atoms[tra->rbOanum[r]].atomname,
     tra->rbType[r],tra->rbAngle[r],
     tra->rbPlane[r]);
 }

} /* end of Set_Type_of_Rotatable_Bond() */



char determine_SpType(mol,a)
  struct MOLECULE *mol;
  int a; /* atom_number */
{
       if ((mol->atoms[a].one_char_ele=='C')&&(mol->atoms[a].Nneighbor==4)){ return(3); }
  else if ((mol->atoms[a].one_char_ele=='C')&&(mol->atoms[a].Nneighbor==3)){ return(2); }
  else if ((mol->atoms[a].one_char_ele=='C')&&(mol->atoms[a].Nneighbor==2)){ return(1); }
  else if ((mol->atoms[a].one_char_ele=='N')&&(mol->atoms[a].Nneighbor==3)){ return(3); }
  else if ((mol->atoms[a].one_char_ele=='N')&&(mol->atoms[a].Nneighbor==2)){ return(2); }
  else if ((mol->atoms[a].one_char_ele=='N')&&(mol->atoms[a].Nneighbor==1)){ return(1); }
  else if ((mol->atoms[a].one_char_ele=='O')&&(mol->atoms[a].Nneighbor==2)){ return(3); }
  else if ((mol->atoms[a].one_char_ele=='O')&&(mol->atoms[a].Nneighbor==1)){ return(2); }
  else if ((mol->atoms[a].one_char_ele=='S')&&(mol->atoms[a].Nneighbor==2)){ return(3); }
  else if ((mol->atoms[a].one_char_ele=='S')&&(mol->atoms[a].Nneighbor==1)){ return(2); }
  return(0);
}







