/*

< qRMS.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 RMS calculation based using quaternion-rotation

 based on 
  Charles F.F. Karney "Quaternions in molecular modeling"
  E-print: arXiv:physics/0506177

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "qRMS.h"

/*** FUNCTIONS (GLOBAL) ***/
float Calculate_CRMS_MATCH_Quaternion();
float Calculate_CRMS_bwn_POS_ARRAYs();
float CRMS_MATCHed_Molecules_without_Superimpose();
float CRMS_Same_Molecules_without_Superimpose();
float Calculate_CRMS_Quaternion_BandA();
float Rotate_Molecule();
float RMSD_bwn_Rotated_and_Init();
void calculate_gcenter();
void Malloc_POS_ARRAY();
void Free_POS_ARRAY();
float Calculate_CRMS_MATCH_rotation_Z();
int  make_rot_axis_angle_from_rot_matrix();
void make_rot_matrix_from_rot_axis_angle();
void make_tvec_from_rmat_gorig_gtra();


/*** FUNCTIONS (LOCAL) ***/
static void Cal_Optimal_Rmatrix_Quaternion();
static void Cal_Rmatrix_From_Quaternion();
static double Cal_RMS();
static double Cal_RMS_by_rotation_of_pA();
static void Jacobi_Wilkinson4();
static void Set_G_to_Zero();
static void Rotation();
static void find_max_abs4();
static void Find_Minimum_Eigen_Value4();
static void equal_matrix4();
static void equal_matrix3();
static void prod_matrix4();
static void prod_matrix3();
static void Mult_Mat_Vec3();
static void print_matrix4();
static void print_matrix3();
static void make_rotmatrix_from_rot_axis_and_angle();
static void make_random_rotmatrix();
static void multiply_by_two_points_axis_random_rotmatrix();
static double normalize_vec3();


void Malloc_POS_ARRAY(PosArray,N)
  struct POS_ARRAY *PosArray;
  int N;
{
 int i;
 PosArray->Nmalloc = N;
 PosArray->pos = (double **)malloc(sizeof(double *) * N); 
 for (i=0;i<N;++i){
   PosArray->pos[i] = (double *)malloc(sizeof(double)*3);}
} /* end of Malloc_POS_ARRAY() */


void Free_POS_ARRAY(PosArray)
  struct POS_ARRAY *PosArray;
{
 int i;
 for (i=0;i<PosArray->Nmalloc;++i){
   free(PosArray->pos[i]);
 }
 free(PosArray->pos);
} /* end of Malloc_POS_ARRAY() */



float Calculate_CRMS_MATCH_Quaternion(match,mA,mB,gA,gB,Rmat,SupType)
  struct MATCH    *match; 
  struct MOLECULE *mA,*mB;  /* mA will be transformed on mB. */
  double gA[3],gB[3]; /* Center position */
  double Rmat[3][3];  /* Rotation Matrix */
  char   *SupType;    /* "BonA" or "AonB" */
{ 
  int i,j;
  float rms;
  struct POS_ARRAY mch_posT;  /* matched atoms of target (to be transformed) */
  struct POS_ARRAY mch_posR;  /* matched atoms of refererence (fixed)       */
 
  if (match->Npair==0){
    gA[0] = gA[1] = gA[2] = 0.0;
    gB[0] = gB[1] = gA[2] = 0.0;
    Rmat[0][0] = 1.0; Rmat[0][1] = 0.0; Rmat[0][2] = 0.0; 
    Rmat[1][0] = 0.0; Rmat[1][1] = 1.0; Rmat[1][2] = 0.0; 
    Rmat[2][0] = 0.0; Rmat[2][1] = 0.0; Rmat[2][2] = 1.0; 
    return(0.0);
  }

  Malloc_POS_ARRAY(&mch_posT,match->Npair);
  Malloc_POS_ARRAY(&mch_posR,match->Npair);
  mch_posT.N = mch_posR.N = match->Npair; 

  for (i=0;i<match->Npair;++i){
    for (j=0;j<3;++j){
      if (strcmp(SupType,"BonA")==0){
        mch_posT.pos[i][j] = mB->atoms[match->anumB[i]].Pos[j]; 
        mch_posR.pos[i][j] = mA->atoms[match->anumA[i]].Pos[j]; 
      }
      else{
        mch_posT.pos[i][j] = mA->atoms[match->anumA[i]].Pos[j]; 
        mch_posR.pos[i][j] = mB->atoms[match->anumB[i]].Pos[j]; 
      }
    }
  }

  if (strcmp(SupType,"BonA")==0){
    Set_G_to_Zero(&mch_posT,gB);
    Set_G_to_Zero(&mch_posR,gA);
  }
  else{
    Set_G_to_Zero(&mch_posT,gA);
    Set_G_to_Zero(&mch_posR,gB);
  }

  Cal_Optimal_Rmatrix_Quaternion(&mch_posT,&mch_posR,Rmat);
  /*
  rms = (float)Cal_RMS_by_rotation_of_pA(&mch_posT, &mch_posR,Rmat);
  printf("#rms_ori %f\n",rms);
  */

  /** Randomize for match->Npair==1 or 2 **/
  if (match->Npair==1){ 
    make_random_rotmatrix(Rmat); 
    /* print_matrix3(Rmat,"Rmat_new"); */
  } 
  else if (match->Npair==2){ 
    multiply_by_two_points_axis_random_rotmatrix(Rmat, mch_posR.pos[0], mch_posR.pos[1]);
    /* print_matrix3(Rmat,"Rmat_new"); */
  }

  Rotation(&mch_posT,Rmat);
  rms = (float)Cal_RMS(&mch_posT,&mch_posR);
  /* printf("#rms_new %f\n",rms); */
  Free_POS_ARRAY(&mch_posR);
  Free_POS_ARRAY(&mch_posT);
  return(rms);
} /* end of Calculate_CRMS_MATCH_Quaternion() */











float Calculate_CRMS_bwn_POS_ARRAYs(parrayA,parrayB)
 struct POS_ARRAY *parrayA,*parrayB; 
{ 
  int i,j;
  float rms;
  double gA[3],gB[3]; /* Center position */
  double Rmat[3][3];  /* Rotation Matrix */
  struct POS_ARRAY pAtmp,pBtmp; 

  if (parrayA->N==0){
    for (i=0;i<3;++i){
      gA[i] = gB[i] = 0.0;
      for (j=0;j<3;++j){ Rmat[i][j] = 0.0; }
      Rmat[i][i] = 1.0; 
    }
    return(0.0);
  }

 /*
  for (i=0;i<parrayA->N;++i){ 
    printf("HETATM%5d %4s%4s %s%4d    %8.3f%8.3f%8.3f\n",
       i+1,"C","MOL","A",i+1, parrayA->pos[i][0], parrayA->pos[i][1], parrayA->pos[i][2]);
  }
  for (i=0;i<parrayB->N;++i){ 
    printf("HETATM%5d %4s%4s %s%4d    %8.3f%8.3f%8.3f\n",
       i+1,"C","MOL","B",i+1, parrayB->pos[i][0], parrayB->pos[i][1], parrayB->pos[i][2]);
  }
  */
  /* [1] Copy parrayA to pAtmp, parrayB to pBtmp */
  Malloc_POS_ARRAY(&pAtmp,parrayA->N);
  pAtmp.N = parrayA->N; 
  for (i=0;i<parrayA->N;++i){ 
    for (j=0;j<3;++j){ pAtmp.pos[i][j] = parrayA->pos[i][j]; }
  }

  Malloc_POS_ARRAY(&pBtmp,parrayB->N);
  pBtmp.N = parrayB->N; 
  for (i=0;i<parrayB->N;++i){ 
    for (j=0;j<3;++j){ pBtmp.pos[i][j] = parrayB->pos[i][j]; }
  }


  /* [2] cal gA, gB and Rmat */
  Set_G_to_Zero(&pAtmp,gA);
  Set_G_to_Zero(&pBtmp,gB);
  Cal_Optimal_Rmatrix_Quaternion(&pAtmp,&pBtmp,Rmat);
  /*
  rms = (float)Cal_RMS_by_rotation_of_pA(&mch_posA, &mch_posB,Rmat);
  printf("#rms_ori %f\n",rms);
  */

  /** Randomize for match->Npair==1 or 2 **/
  if (parrayA->N==1){ 
    make_random_rotmatrix(Rmat); 
  } 
  else if (parrayA->N==2){ 
    multiply_by_two_points_axis_random_rotmatrix(Rmat, pAtmp.pos[0], pBtmp.pos[1]);
    /* print_matrix3(Rmat,"Rmat_new"); */
  }

  /* [3] calculate RMSD */
  Rotation(&pAtmp,Rmat);
  rms = (float)Cal_RMS(&pAtmp,&pBtmp);

  Free_POS_ARRAY(&pBtmp);
  Free_POS_ARRAY(&pAtmp);
  /* printf("#rms_new %f\n",rms); */
  return(rms);
} /* end of Calculate_CRMS_bwn_POS_ARRAYs() */






float CRMS_MATCHed_Molecules_without_Superimpose(match,mA,mB)
 struct MATCH    *match; 
 struct MOLECULE *mA,*mB; 
{ 
  int i,j;
  float rms,d;
 
  rms = 0.0; 
  if (match->Npair>0){
    for (i=0;i<match->Npair;++i){
      for (j=0;j<3;++j){
        d = mA->atoms[match->anumA[i]].Pos[j] - mB->atoms[match->anumB[i]].Pos[j]; 
        rms += d*d;
      }
    }
    rms = sqrt(rms/match->Npair);
  }
  return(rms);
} /* end of CRMS_MATCHed_Molecules_without_Superimpose() */




float CRMS_Same_Molecules_without_Superimpose(mA,mB)
 struct MOLECULE *mA,*mB; 
{ 
  int i,j;
  float rms,d;
 
  rms = 0.0; 
  for (i=0;i<mA->Natom;++i){
    for (j=0;j<3;++j){
      d = mA->atoms[i].Pos[j] - mB->atoms[i].Pos[j]; 
      rms += d*d;
    }
  }
  rms = sqrt(rms/mA->Natom);
  return(rms);

} /* end of CRMS_Same_Molecules_without_Superimpose(mA,mB) */





float Calculate_CRMS_Quaternion_BandA(match,mA,mB,g1,g2,Rmat)
 struct MATCH    *match; 
 struct MOLECULE *mA,*mB; 
 double g1[3],g2[3]; /* Center position */
 double Rmat[3][3];  /* Rotation Matrix */
{ 
  int i,j;
  float rms;
  struct POS_ARRAY mch_posA,mch_posB;  /* matched atoms */

  Malloc_POS_ARRAY(&mch_posA,match->Npair);
  Malloc_POS_ARRAY(&mch_posB,match->Npair);
  mch_posA.N = mch_posB.N = match->Npair; 
  for (i=0;i<match->Npair;++i){
    for (j=0;j<3;++j){
      mch_posA.pos[i][j] = mA->atoms[match->anumB[i]].Pos[j]; 
      mch_posB.pos[i][j] = mB->atoms[match->anumA[i]].Pos[j]; 
    }
  }

  Set_G_to_Zero(&mch_posA,g1);
  Set_G_to_Zero(&mch_posB,g2);
  Cal_Optimal_Rmatrix_Quaternion(&mch_posA,&mch_posB,Rmat);
  Rotation(&mch_posA,Rmat);
  rms = (float)Cal_RMS(&mch_posA,&mch_posB);
  Free_POS_ARRAY(&mch_posB);
  Free_POS_ARRAY(&mch_posA);
  return(rms);

} /* end of Calculate_CRMS_Quaternion_BandA() */





float Rotate_Molecule(mA,gA,gB,Rmat)
 struct MOLECULE *mA;
 double gA[3];
 double gB[3];
 double Rmat[3][3];
 /**
  return [rmsd between original and rotated molecules].
  **/ 
{
 int i,j;
 double pos[3],newpos[3],rpos[3];
 float rmsd;

 /*
 printf("gA %lf %lf %lf\n",gA[0],gA[1],gA[2]);
 printf("gB %lf %lf %lf\n",gB[0],gB[1],gB[2]);
 printf("R0 %lf %lf %lf\n",Rmat[0][0],Rmat[0][1],Rmat[0][2]);
 printf("R1 %lf %lf %lf\n",Rmat[1][0],Rmat[1][1],Rmat[1][2]);
 printf("R2 %lf %lf %lf\n",Rmat[2][0],Rmat[2][1],Rmat[2][2]);
 */

 rmsd = 0.0;
 /* print_matrix3(Rmat,"Rmat_rotate_molecule");  */
 for (i=0;i<mA->Natom;++i){
 /*
   printf("#%d (%f %f %f) -->" ,i,mA->atoms[i].Pos[0], mA->atoms[i].Pos[1], mA->atoms[i].Pos[2]);
 */ 
  for (j=0;j<3;++j) { pos[j] = mA->atoms[i].Pos[j] - gA[j];}

   Mult_Mat_Vec3(rpos,Rmat,pos);

   for (j=0;j<3;++j){  
     newpos[j] = rpos[j] + gB[j];
     rmsd += (newpos[j]-mA->atoms[i].Pos[j])*(newpos[j]-mA->atoms[i].Pos[j]);
     mA->atoms[i].Pos[j] = newpos[j];
   }
 /*
   printf("#(%f %f %f)\n" ,mA->atoms[i].Pos[0], mA->atoms[i].Pos[1], mA->atoms[i].Pos[2]);
 */
 }
 return(sqrt(rmsd/mA->Natom));
 
} /* end of Rotate_Molecule() */


float RMSD_bwn_Rotated_and_Init(mA,gA,gB,Rmat)
 struct MOLECULE *mA;
 double gA[3];
 double gB[3];
 double Rmat[3][3];
 /**
  return [rmsd between original and rotated molecules].
  **/ 
{
 int i,j;
 double pos[3],newpos[3],rpos[3];
 float rmsd;

 rmsd = 0.0;
 for (i=0;i<mA->Natom;++i){
   for (j=0;j<3;++j) pos[j] = mA->atoms[i].Pos[j] - gA[j];
   Mult_Mat_Vec3(rpos,Rmat,pos);
   for (j=0;j<3;++j){  
     newpos[j] = rpos[j] + gB[j];
     rmsd += (newpos[j]-mA->atoms[i].Pos[j])*(newpos[j]-mA->atoms[i].Pos[j]);
   }
 }
 return(sqrt(rmsd/mA->Natom));
 
} /* end of RMSD_bwn_Rotated_and_Init() */


void Set_G_to_Zero(parray,G)
 struct POS_ARRAY *parray;
 double G[3];
{
 int i,j;
 double sum[3];
 
 sum[0] = sum[1] = sum[2] = 0.0;
 for (i=0;i<parray->N;++i){ 
   for (j=0;j<3;++j) sum[j] += parray->pos[i][j];
  } 

 for (j=0;j<3;++j) G[j] = sum[j]/parray->N;

 for (i=0;i<parray->N;++i){
  for (j=0;j<3;++j) parray->pos[i][j] = parray->pos[i][j] - G[j];
 }
} /* end of Set_G_to_Zero() */


void Rotation(parray,mat)
 struct POS_ARRAY *parray; 
 double mat[3][3];
{
 double p[3],q[3];
 int i;

 for (i=0;i<parray->N;++i){ 
   p[0] = parray->pos[i][0];  
   p[1] = parray->pos[i][1];  
   p[2] = parray->pos[i][2];  
   Mult_Mat_Vec3(q,mat,p);
   parray->pos[i][0] = q[0]; 
   parray->pos[i][1] = q[1]; 
   parray->pos[i][2] = q[2];
 }

} /* end of Rotation() */




double Cal_RMS(pA,pB)
 struct POS_ARRAY *pA,*pB;
{ int i;
  double dx,dy,dz, RM, rms,dd;

 RM = 0.0; 
 for (i=0;i<pA->N;++i){ 
   dx = pA->pos[i][0] - pB->pos[i][0]; 
   dy = pA->pos[i][1] - pB->pos[i][1]; 
   dz = pA->pos[i][2] - pB->pos[i][2]; 
   dd = (dx*dx)+(dy*dy)+(dz*dz);
   RM += dd;
  }

 if (RM>0.0) rms = sqrt(RM/pA->N); else rms = 0.0;
 return(rms);
} /* end of Cal_RMS() */



double Cal_RMS_by_rotation_of_pA(pA, pB, R)
 struct POS_ARRAY *pA,*pB;
 double R[3][3];
{
 double posA[3],RposA[3];
 double dx,dy,dz, RM, rms,dd;
 int i,j;

 for (i=0;i<pA->N;++i){ 
   for (j=0;j<3;++j){ posA[j] = pA->pos[i][j];}  
   Mult_Mat_Vec3(RposA,R,posA);
   dx = RposA[0] - pB->pos[i][0]; 
   dy = RposA[1] - pB->pos[i][1]; 
   dz = RposA[2] - pB->pos[i][2]; 
   dd = (dx*dx)+(dy*dy)+(dz*dz);
   RM += dd;
 }
 if (RM>0.0) rms = sqrt(RM/pA->N); else rms = 0.0;
 return(rms);
} /* end of Cal_RMS_by_rotation_of_pA() */



void Cal_Optimal_Rmatrix_Quaternion(pA,pB,R)
 struct POS_ARRAY *pA,*pB; 
 double R[3][3];
{
 int i,j,k;
 double B[4][4],E[4][4],V[4][4];
 double a[3],b[3],aa[3],bb[3];
 double min_eval,min_evec[4];

/*
 printf("#int Cal_Optimal_Rmatrix_Quaternion(pA,pB,R)\n");
*/

 /**** (1) Make Matrix B  *******/ 
 /*
   a = yk+xk
   b = yk-xk
   Ak = |0.0 -b0 -b1  -b2|
        |b0  0.0 -a2   a1|
        |b1  a2  0.0  -a0|
        |b2 -a1  a0   0.0|
  
   Bk = tra[Ak] * Ak 
      = |0.0  b0  b1    b2| |0.0 -b0 -b1  -b2|
        |-b0  0.0 a2   -a1| |b0  0.0 -a2   a1|
        |-b1 -a2  0.0   a0| |b1  a2  0.0  -a0|
        |-b2  a1  -a0  0.0| |b2 -a1  a0   0.0|

      = |b0b0+b1*b1+b2*b2 a2b1-a1b2        -a2b0+a0b2       a1b0-a0b1       |
        |a2b1-a1b2        b0b0+a1*a1+b2*b2 b0b1-a0a1        b0b2-a0a2       |
        |-a2b0+a0b2       b0b1-a0a1        a0a0+b1*b1+a2*a2 b1b2-a1a2       |
        | a1b0-a0b1       b0b2-a0a2        b1b2-a1a2        a0a0+a1*a1+b2*b2|

   B = 1/N * \sum_{k}^{N} Bk

*/ 

 for (i=0;i<4;++i)
  for (j=0;j<4;++j) B[i][j] = 0.0;


 for (k=0;k<pA->N;++k){
 
  for (i=0;i<3;++i) {
   a[i] = pB->pos[k][i] + pA->pos[k][i];  
   b[i] = pB->pos[k][i] - pA->pos[k][i]; 
   aa[i] = a[i]*a[i];
   bb[i] = b[i]*b[i]; }  

  B[0][0] += bb[0]+bb[1]+bb[2];
  B[0][1] +=  a[2]*b[1] - a[1]*b[2];
  B[0][2] += -a[2]*b[0] + a[0]*b[2];
  B[0][3] +=  a[1]*b[0] - a[0]*b[1];
  
  B[1][1] += bb[0]+aa[1]+aa[2];
  B[1][2] +=  b[0]*b[1] - a[0]*a[1];
  B[1][3] +=  b[0]*b[2] - a[0]*a[2];
 
  B[2][2] += aa[0]+bb[1]+aa[2];
  B[2][3] +=  b[1]*b[2] - a[1]*a[2];
  B[3][3] += aa[0]+aa[1]+bb[2];

 } /* k */

  for (i=0;i<4;++i){
    for (j=i;j<4;++j){ 
       B[i][j] /= pA->N;
       B[j][i] = B[i][j]; 
    }
  }

 /**** (2) Calculate Eigen Value/Vector of Matrix B  *******/ 
 equal_matrix4(E,B); 
 /* print_matrix4(B,"Bmatrix"); */
 Jacobi_Wilkinson4(E,V);
 /* print_matrix4(E,"EigVal"); */
 /* print_matrix4(V,"EigVec"); */
 Find_Minimum_Eigen_Value4(E,V,&min_eval,min_evec);
 
 /**** (3) Calculate Rmat from the eigen vector quaternion ****/
 Cal_Rmatrix_From_Quaternion(R,min_evec);
 /* print_matrix3(R,"Rmat"); */


} /* end of Cal_Optimal_Rmatrix_Quaternion() */







void Find_Minimum_Eigen_Value4(E,V,min_eval,min_evec)
 double E[4][4]; /* Diagonal matrix for eigen values */
 double V[4][4]; /* Matrix for eigen vectors */
 double *min_eval;   /* minimum_evalue to be calulated */
 double min_evec[4];   /* evector for minimum_evalue to be calulated */
{
 int i,min_i;
 double ev[4];

 min_i = 0;
 for (i=0;i<4;++i){ ev[i] = E[i][i];}

      if ((ev[0]<=ev[1])&&(ev[0]<=ev[2])&&(ev[0]<=ev[3])){ min_i = 0;}
 else if ((ev[1]<=ev[0])&&(ev[1]<=ev[2])&&(ev[1]<=ev[3])){ min_i = 1;}
 else if ((ev[2]<=ev[0])&&(ev[2]<=ev[1])&&(ev[2]<=ev[3])){ min_i = 2;}
 else if ((ev[3]<=ev[0])&&(ev[3]<=ev[1])&&(ev[3]<=ev[2])){ min_i = 3;}

 *min_eval = ev[min_i];

 for (i=0;i<4;++i){min_evec[i] = V[i][min_i];}

} /* end of Find_Minimum_Eigen_Value4() */



void Cal_Rmatrix_From_Quaternion(R,q)
 double R[3][3],q[4];
{
 int i,j;
 double Q[4][4];

 for (i=0;i<4;++i) 
  for (j=i;j<4;++j) Q[i][j] = Q[j][i] = q[i]*q[j]; 
 
 R[0][0] = 2.0*(Q[0][0]+Q[1][1])-1.0;
 R[0][1] = 2.0*(Q[1][2]-Q[0][3]);
 R[0][2] = 2.0*(Q[1][3]+Q[0][2]);
 
 R[1][0] = 2.0*(Q[1][2]+Q[0][3]);
 R[1][1] = 2.0*(Q[0][0]+Q[2][2])-1.0;
 R[1][2] = 2.0*(Q[2][3]-Q[0][1]);

 R[2][0] = 2.0*(Q[1][3]-Q[0][2]);
 R[2][1] = 2.0*(Q[2][3]+Q[0][1]);
 R[2][2] = 2.0*(Q[0][0]+Q[3][3])-1.0;

 /*
 R[0][0] = 1.0-2.0*(Q[2][2]+Q[3][3]);
 R[0][1] = 2.0*(Q[1][2]+Q[0][3]);
 R[0][2] = 2.0*(Q[1][3]-Q[0][2]);
 
 R[1][0] = 2.0*(Q[1][2]-Q[0][3]);
 R[1][1] = 1.0-2.0*(Q[1][1]+Q[3][3]);
 R[1][2] = 2.0*(Q[2][3]+Q[0][1]);

 R[2][0] = 2.0*(Q[1][3]+Q[0][2]);
 R[2][1] = 2.0*(Q[2][3]-Q[0][1]);
 R[2][2] = 1.0-2.0*(Q[1][1]+Q[2][2]);
 */

} /* end of Cal_Rmatrix_From_Quaternion() */



void Jacobi_Wilkinson4(A,U)
 double A[4][4]; /* Input matrix whose eigen value/vector will be calculated. 
      After calculation, it becomes diagonal matrix of eigen values. */
 double U[4][4]; /* Output matrix corresponding to eigen vectors. 
   Eigen vector corresponds to each "column" vector of matrix U.
   For example k-th eigen vector 
   evec_k[0] =  U[0][k]; 
   evec_k[1] =  U[1][k]; 
   evec_k[2] =  U[2][k]; 
   evec_k[3] =  U[3][k]; 
   */ 
{ 
 double R[4][4],TR[4][4],BUFF[4][4];
 double max,co,si;
 double sa,r,sqrt2;
 int i,j,mi,mj,c;

 sqrt2 = sqrt(2.0);

   /* --------SETTING U[][] = I --------*/
   
   for (i=0;i<4;++i)
     for (j=0;j<4;++j) if (i==j) U[i][j] = 1.0;  else U[i][j] = 0.0;
    
    find_max_abs4(&mi,&mj,&max,A);
    c = 0;
    
    while ((max>0.00000001)&&(c<100)){
      /*--------- SET of cos sin -----------*/ 
      
      sa = (A[mi][mi] - A[mj][mj])/2;  
      r = sqrt(sa*sa + A[mi][mj]*A[mi][mj]);
      co =  sqrt(1.0+sa/r)/sqrt2;  si =  A[mi][mj]/2.0/r/co; 
      if (sa>0.0) 
       { co =  sqrt(1.0+sa/r)/sqrt2;  si =  A[mi][mj]/2.0/r/co; }
      else 
       { co =   sqrt(1.0-sa/r)/sqrt2; si = - A[mi][mj]/2.0/r/co; }
      
      /* -------- SET of Rot Matrix R and TR----- */  
      for (i=0;i<4;++i)
	for (j=0;j<4;++j) if (i==j)  R[i][j] = 1.0;  else R[i][j] =0.0;
     
       R[mi][mi] =   co; R[mi][mj] = -si;
       R[mj][mi] =   si; R[mj][mj] = co;

      for (i=0;i<4;++i)
	for (j=0;j<4;++j) TR[j][i] = R[i][j]; 

       prod_matrix4(BUFF,A,R);
       prod_matrix4(A,TR,BUFF);
       prod_matrix4(BUFF,U,R);
       equal_matrix4(U,BUFF);

       find_max_abs4(&mi,&mj,&max,A);
       ++c;
 
      } /* end of while */

}/* end of Jacobi_Wilkinson() */




void find_max_abs4(mi,mj,max,A) /* Max_Abs element A[mi][mj] = max */
 int *mi,*mj;
 double *max,A[4][4];
{ int i,j,p,q;

 *max =0.0;
  p = q = 0;
  for (i=0;i<4;++i)
    for (j=0;j<4;++j)
      if (i!=j)
       if (fabs(A[i][j])> *max)
	      { *max = fabs(A[i][j]);     p = i; q= j;}  
   *mi =p; *mj = q;
} /* end of find_max_abs4() */


void equal_matrix4(A,B)    /* A = B */
 double A[4][4],B[4][4];
{ int i,j;
  for (i=0;i<4;++i){
    for (j=0;j<4;++j){ A[i][j] = B[i][j];}
   }
}/* end of equal_matrix4() */



void equal_matrix3(A,B)    /* A = B */
 double A[3][3],B[3][3];
{ int i,j;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ A[i][j] = B[i][j];}
   }
}/* end of equal_matrix3() */

void prod_matrix4(C,A,B) /* C = A * B */
 double C[4][4],A[4][4],B[4][4];
{ int i,j,k;
  double sum;
  for (i=0;i<4;++i){
    for (j=0;j<4;++j){ 
       sum =0;
       for (k=0;k<4;++k){ sum = sum + A[i][k] * B[k][j];}
       C[i][j] = sum; 
     }
  }
} /* end of prod_matrix4() */


void prod_matrix3(C,A,B) /* C = A * B */
 double C[3][3],A[3][3],B[3][3];
{ int i,j,k;
  double sum;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
       sum =0;
       for (k=0;k<3;++k){ sum = sum + A[i][k] * B[k][j];}
       C[i][j] = sum; 
     }
  }
} /* end of prod_matrix3() */



void Mult_Mat_Vec3(y,mat,x)
 double y[3],mat[3][3],x[3];  /* This function is only for 3-D */
{ 
  y[0] = mat[0][0]*x[0] + mat[0][1]*x[1] + mat[0][2] * x[2];
  y[1] = mat[1][0]*x[0] + mat[1][1]*x[1] + mat[1][2] * x[2];
  y[2] = mat[2][0]*x[0] + mat[2][1]*x[1] + mat[2][2] * x[2]; 
}


void print_matrix4(A,comment)
 double A[4][4];
 char *comment;
{ int i,j;
  printf("\n");
  printf("#%s\n",comment);
  for (i=0;i<4;++i){
    for (j=0;j<4;++j){ 
      printf(" %8.3lf",A[i][j]);
      if (j==3){printf(";\n");} else {printf(",");} 
    }
  }
  printf("\n");
}/* end of print_matrix4() */


void print_matrix3(A,comment)
 double A[3][3];
 char *comment;
{ int i,j;
  printf("\n");
  printf("#%s\n",comment);
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
      printf(" %8.3lf",A[i][j]);
      if (j==2){printf(";\n");} else {printf(",");} 
    }
  }
  printf("\n");
}/* end of print_matrix3() */


void calculate_gcenter(mol,gcen)
  struct MOLECULE *mol;
  double gcen[3];
{
  int i,N;

  N = 0;
  gcen[0] = gcen[1] = gcen[2] = 0.0;
  for (i=0;i<mol->Natom;++i){
    if (mol->atoms[i].one_char_ele != 'H'){
       gcen[0] += mol->atoms[i].Pos[0];
       gcen[1] += mol->atoms[i].Pos[1];
       gcen[2] += mol->atoms[i].Pos[2];
       N += 1;
    }
  } 
  if (N>0){
    gcen[0] /= N; gcen[1] /= N; gcen[2] /= N;
  }
}



void make_rotmatrix_from_rot_axis_and_angle(R,axis,angle)
  double R[3][3];
  double axis[3];
  double angle; /* radian */
{
  double S[3][3];
  double SS[3][3];
  double a,b,c;
  double sint,ccost;
  int i,j;

  a = axis[0]; b = axis[1]; c = axis[2];
  S[0][0] = 0.0; S[0][1] = -c;  S[0][2] = b;
  S[1][0] = c;   S[1][1] = 0.0; S[1][2] = -a;
  S[2][0] = -b;  S[2][1] =  a;  S[2][2] = 0.0;


  SS[0][0] = -c*c-b*b; SS[0][1] = a*b;      SS[0][2] = a*c;
  SS[1][0] = a*b;      SS[1][1] = -c*c-a*a; SS[1][2] = b*c;
  SS[2][0] = a*c;      SS[2][1] = b*c;      SS[2][2] = -a*a-b*b;

  sint  = sin(angle);
  ccost = 1.0 - cos(angle);
  /* printf("angle %f sint %f ccost %f\n",angle,sint,ccost);  */
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      R[i][j] = sint*S[i][j] + ccost*SS[i][j];
      if (i==j) {R[i][j] += 1.0;}
    }
  }

} /* end of make_rotmatrix_from_rot_axis_and_angle() */



void make_random_rotmatrix(R)
  double R[3][3];
{
  double raxis[3], rangle;
  int i;
  
  for (i=0;i<3;++i){ raxis[i] = -1.0 + 2.0*((float)rand())/RAND_MAX; }
  normalize_vec3(raxis);
  rangle  = 2.0 * M_PI * ((float)rand())/RAND_MAX;
  make_rotmatrix_from_rot_axis_and_angle(R,raxis,rangle);

} /* end of make_random_rotmatrix() */




void multiply_by_two_points_axis_random_rotmatrix(R, posA, posB)
  double R[3][3]; /* (input and multiplied by random matrix) */
  double posA[3]; /* (input) */
  double posB[3]; /* (input) */
{
  double raxis[3], rangle, RposAB[3][3], RposAB_R[3][3];
  int i;

  /*
  printf("#multiply_by_two_points_axis_random_rotmatrix(posA %lf %lf %lf posB %lf %lf %lf\n",posA[0],posA[1],posA[2],posB[0],posB[1],posB[2]);
  */

  /* (1) make rotarional axis axis[] = posA[] - posB[] */
  for (i=0;i<3;++i){ raxis[i] = posA[i] - posB[i]; }
  normalize_vec3(raxis);

  /* (2) make random matrix RposAB */
  rangle  = 2.0 * M_PI * ((float)rand())/RAND_MAX;

  /* printf("#axis %lf %lf %lf rangle %lf\n",raxis[0],raxis[1],raxis[2],rangle); */
  make_rotmatrix_from_rot_axis_and_angle(RposAB, raxis, rangle);
  /* print_matrix3(RposAB,"RposAB"); */
  /* (3) R = RposAB * R */
  prod_matrix3(RposAB_R,RposAB,R);
  equal_matrix3(R,RposAB_R);
 
} /* end of multiply_two_points_axis_random_rotmatrix() */


double normalize_vec3(v)
  double v[3];
{
  double norm;
  norm = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  norm = sqrt(norm);
  if (norm<=0.0) return(0.0);
  v[0] /= norm;  v[1] /= norm;  v[2] /= norm;
  return(norm);
}


int make_rot_axis_angle_from_rot_matrix(w,R)
  double w[4];     /* rot axis [0:1:2] and rot_angle [3]  (to be calculated)*/
  double R[3][3];  /* rotation matrix (input) */
{
  double cos_t, sin_t,u[3],u_len,eps;
  int i,i0,j0,k0;

/*
  print_matrix3(R,"Rorig"); 
*/
  eps = 0.00000001;

  /* (1) calculate cos_t */
  cos_t = (R[0][0]+R[1][1]+R[2][2]-1.0)/2.0;
  if (fabs(cos_t-1) < eps){
    w[0] = 1.0; w[1] = 0.0; w[2] = 0.0; w[3] = 0.0;
    return(1);
  }

 /* (2) Decide u[i0] as positive from u[]^2 */
  for (i=0;i<3;++i){
    u[i] = (R[i][i]-cos_t)/(1.0-cos_t);
  }
  i0 = 0;
       if ((u[0]>=u[1]) && (u[0]>=u[2])) { i0 = 0;}
  else if ((u[1]>=u[0]) && (u[1]>=u[2])) { i0 = 1;}
  else if ((u[2]>=u[0]) && (u[2]>=u[1])) { i0 = 2;}

  u[i0] = sqrt(u[i0]);

  j0 = (i0+1)%3;
  k0 = (i0+2)%3;

/*
    printf("#u0 %lf u1 %lf u2 %lf\n",u[0],u[1],u[2]);
     printf("#i0 %d j0 %d k0 %d\n",i0,j0,k0); 
 */
 /* (3) Calculate sin_t from u[i0]  */
  sin_t = (R[k0][j0]-R[j0][k0])/(2.0*u[i0]);

 /* (4) Calculate u[j0] and u[k0]  */
  if (fabs(sin_t) < eps){
      u[j0] = (R[i0][j0])/(2.0*u[i0]);
      u[k0] = (R[i0][k0])/(2.0*u[i0]);
  }
  else{ 
      u[j0] = (R[i0][k0]-R[k0][i0])/(2.0*sin_t);
      u[k0] = (R[j0][i0]-R[i0][j0])/(2.0*sin_t);
  }
  /*
    printf("#cos_t %lf sin_t %lf\n",cos_t,sin_t);
    printf("#sin_t %lf\n",(R[1][0]-R[0][1])/(2.0*u[2]));
    printf("#sin_t %lf\n",(R[0][2]-R[2][0])/(2.0*u[1]));
    printf("#sin_t %lf\n",(R[2][1]-R[1][2])/(2.0*u[0]));
   */

    u_len = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    for (i=0;i<3;++i){
      w[i] = u[i]/u_len;
    }

  /* (4) Calculate rot angle from cos_t and sin_t  */
  if (cos_t>1.0){
     w[3] = 0.0; 
  }
  else if (cos_t< -1.0){
     w[3] = M_PI; 
  }
  else{
    w[3] = acos(cos_t);
  }
  if (sin_t<0.0){
    w[3] = -w[3];
  }

  /* printf("#cos_t %lf sin_t %lf w[3] %lf\n",cos_t,sin_t,w[3]);  */

  return(1);
} /* end of make_rot_axis_angle_from_rot_matrix() */


void make_rot_matrix_from_rot_axis_angle(R,w)
  double R[3][3];  /* rotation matrix (to be calculated) */
  double w[4];     /* rot axis [0:1:2] and rot_angle [3]  (input) */
{
  double Ux[3][3], Uxx[3][3],cos_t,sin_t,one_cos_t;
  int i,j;
  /* printf("#w %lf %lf %lf %lf\n",w[0],w[1],w[2],w[3]); */
  Ux[0][0] =   0.0;  Ux[0][1] = -w[2];  Ux[0][2] =  w[1];
  Ux[1][0] =  w[2];  Ux[1][1] =   0.0;  Ux[1][2] = -w[0];
  Ux[2][0] = -w[1];  Ux[2][1] =  w[0];  Ux[2][2] =  0.0;

  Uxx[0][0] =  w[0]*w[0]; Uxx[0][1] = w[0]*w[1]; Uxx[0][2] = w[0]*w[2];
  Uxx[1][0] =  w[1]*w[0]; Uxx[1][1] = w[1]*w[1]; Uxx[1][2] = w[1]*w[2];
  Uxx[2][0] =  w[2]*w[0]; Uxx[2][1] = w[2]*w[1]; Uxx[2][2] = w[2]*w[2];

 /*
  print_matrix3(Ux,"Ux");
  print_matrix3(Uxx,"Uxx");
  */
  cos_t     = cos(w[3]);
  sin_t     = sin(w[3]);
  one_cos_t = 1.0 - cos_t;
  /* printf("#cos_t %lf sin_t %lf\n",cos(w[3]), sin(w[3])); */
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      R[i][j] = sin_t * Ux[i][j] + one_cos_t*Uxx[i][j];
      if (i==j){
        R[i][j] += cos_t;
      }
      /* printf("#R[%d][%d] %lf\n",i,j,R[i][j]); */
    }
  } 
  /* print_matrix3(R,"R"); */
} /* end of make_rot_matrix_from_rot_axis_angle() */



void make_tvec_from_rmat_gorig_gtra(tvec, R,gorig,gtra)
  double  tvec[3]; /* translation vector (to be calculated) */
  double  R[3][3]; /* rotation matrix (input) */
  double  gorig[3]; /* center of original   molecule (input) */
  double  gtra[3];  /* center of transformed molecule (input) */
{
 /*
   xtra = R * (xorig - gorig) + gtra
        = R * xorig - R* gorig + gtra
        = R * xorig + tvec 
 
  Therefore, tvec =  - R* gorig + gtra
 */
  int i;
 
  Mult_Mat_Vec3(tvec,R,gorig);
  for (i=0;i<3;++i){
    tvec[i] = -tvec[i] + gtra[i];
  }



} /* end of make_tvec_from_rmat_gorig_gtra() */
