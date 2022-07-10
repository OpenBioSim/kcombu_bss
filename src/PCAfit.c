/*

< PCAfit.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  3D molecule fitting by PCA (Eigen vectors of Covariance Matrix)

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "molprop.h"

/*** FUNCTIONS (GLOBAL) ***/
void Cal_Molecule_EigVec_Covariance_Matrix();
void Write_EigVec_in_PDB();
void Cal_MATCH_for_PCAfit_Superimposed_3D_Molecule_Pair();
void Cal_Rotation_Matrix_from_Pair_of_XYZvec();
float Cal_Nmatch_for_Superimposed_3D_Molecule_Pair();
void Cal_MATCH_for_Superimposed_3D_Molecule_Pair();
void set_unit_matrix3_double();
void Rotate_Molecule_around_Gcenter();
void Rotate_Molecule_XY_flat_plane_by_PCA();

/*** FUNCTIONS (LOCAL) ***/
static float Dot_Prod3();
static void Sort_Eigen_Value3();
static void Assign_Direction_for_EigVec();
static void Jacobi_Wilkinson3();
static void find_max_abs3_double();
static void equal_matrix3_double();
static void prod_matrix3_double();
static void Mult_Mat_Vec3_float();
static void print_matrix3_float();
static void print_matrix3_double();
static void print_vector3_float();
static void make_rotational_matrix_XY_flat_plane();

void Cal_Molecule_EigVec_Covariance_Matrix(M)
 struct MOLECULE *M;
{
 int a,i,j,Natom,index[3];
 double di,dj,CovMat[3][3],EigVec[3][3];

 /** [1] Calculate Gcen **/
 for (i=0;i<3;++i) M->Gcen[i] =  0.0;
 Natom = 0; 
 for (a=0;a<M->Natom;++a){
   if (M->atoms[a].one_char_ele != 'H'){
     Natom += 1;
     for (i=0;i<3;++i) M->Gcen[i] +=  M->atoms[a].Pos[i];
   }
 }

 for (i=0;i<3;++i) M->Gcen[i] /= Natom;

 /* printf("Natom %d Gcen %f %f %f\n",Natom,M->Gcen[0], M->Gcen[1], M->Gcen[2]); */
 
 /** [2] Calculate CovMat **/
 for (i=0;i<3;++i){ 
   for (j=0;j<3;++j){  CovMat[i][j] = 0.0;}
 }
 
 for (a=0;a<M->Natom;++a){
   if (M->atoms[a].one_char_ele != 'H'){
    for (i=0;i<3;++i){ 
        di = M->atoms[a].Pos[i] - M->Gcen[i];
      for (j=i;j<3;++j){ 
        dj = M->atoms[a].Pos[j] - M->Gcen[j];
        CovMat[i][j] +=  di * dj;
      }
    }
   }
 }

 CovMat[1][0] = CovMat[0][1];
 CovMat[2][0] = CovMat[0][2];
 CovMat[2][1] = CovMat[1][2];

 /** [3] Calculate EigVal/EigVec of CovMat **/
/* print_matrix3_double(CovMat,"CovMat");  */
 Jacobi_Wilkinson3(CovMat,EigVec);
 Sort_Eigen_Value3(CovMat,index);

/* printf("#index 0: %d 1: %d 2: %d\n",index[0],index[1],index[2]);  fflush(stdout); */
 

 for (i=0;i<3;++i){ 
   M->eigval[i] = CovMat[index[i]][index[i]];
   for (j=0;j<3;++j){ 
     M->eigvec[i][j] = EigVec[j][index[i]];
   }
 }

 /*
 for (i=0;i<3;++i)
   printf("[%d] %f (%f %f %f)\n",i,M->eigval[i],M->eigvec[i][0],M->eigvec[i][1],M->eigvec[i][2]);
 */

 Assign_Direction_for_EigVec(M);
} /* end of Cal_Molecule_EigVec_Covariance_Matrix() */


void Cal_MATCH_for_PCAfit_Superimposed_3D_Molecule_Pair(m,mA,mB,ConnectGraph)
  struct MATCH *m;
  struct MOLECULE *mA; /* mA->eigval and mA->eigvec are already calculated. */
  struct MOLECULE *mB; /* mB->eigval and mB->eigvec are already calculated. */
  char   ConnectGraph;   /* 'D':disconnected, otherwise:don't care */
{
 float Rmat[3][3],evecB[3][3],maxNmatch,Nmatch;
 int i,j,k,max_k;
 struct MOLECULE mA1;
 int sign[4][3];

 /* 
 printf("#Cal_MATCH_for_PCAfit_Superimposed_3D_Molecule_Pair(m,mA,mB)\n");
 */ 
 /**[1] Find the best superimposition among the four **/

 Malloc_MOLECULE(&mA1,mA->Natom,mA->Nbond);
 Copy_MOLECULE(&mA1,mA);
 sign[0][0] =  1; sign[0][1] =  1; sign[0][2] =  1;  /* 0-degree rotation          */ 
 sign[1][0] =  1; sign[1][1] = -1; sign[1][2] = -1;  /* x-axis 180-degree rotation */
 sign[2][0] = -1; sign[2][1] =  1; sign[2][2] = -1;  /* y-axis 180-degree rotation */
 sign[3][0] = -1; sign[3][1] = -1; sign[3][2] =  1;  /* z-axis 180-degree rotation */

 maxNmatch = 0.0;
 max_k = 0;
 for (k=0;k<4;++k){
   Copy_MOLECULE(&mA1,mA);
   /*
   printf("#kn %d %f %f %f\n",0,mA1.atoms[0].Pos[0], mA1.atoms[0].Pos[1], mA1.atoms[0].Pos[2]);
   printf("#kn %d %f %f %f\n",1,mA1.atoms[1].Pos[0], mA1.atoms[1].Pos[1], mA1.atoms[1].Pos[2]);
   */
 
   for (i=0;i<3;++i){
     for (j=0;j<3;++j){ evecB[i][j] = sign[k][i] * mB->eigvec[i][j];}
   }

   Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat,
    mA->eigvec[0],mA->eigvec[1],mA->eigvec[2],evecB[0],evecB[1],evecB[2]);

   /* printf("Rmat %f %f %f\n",Rmat[0][0],Rmat[0][1],Rmat[0][2]); */

   Rotate_Molecule_around_Gcenter(&mA1,Rmat,mA1.Gcen,mB->Gcen);
  
   /* 
   printf("#km %d %f %f %f\n",0,mA1.atoms[0].Pos[0], mA1.atoms[0].Pos[1], mA1.atoms[0].Pos[2]);
   printf("#km %d %f %f %f\n",1,mA1.atoms[1].Pos[0], mA1.atoms[1].Pos[1], mA1.atoms[1].Pos[2]);
    */
 
   Nmatch = Cal_Nmatch_for_Superimposed_3D_Molecule_Pair(&mA1,mB);
   if (Nmatch > maxNmatch){maxNmatch = Nmatch; max_k = k;}
 }

 /* printf("#max_k %d maxMatch %d\n",max_k,maxNmatch); */
 /** [2] Rotate mA1 by the max_k rotation **/
 Copy_MOLECULE(&mA1,mA);
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){ evecB[i][j] = sign[max_k][i] * mB->eigvec[i][j];}
 }
 
 Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat,
   mA->eigvec[0],mA->eigvec[1],mA->eigvec[2],evecB[0],evecB[1],evecB[2]);
 
 Rotate_Molecule_around_Gcenter(&mA1,Rmat,mA1.Gcen,mB->Gcen);

 /** [3] Find  atomic MATCH  **/
 Cal_MATCH_for_Superimposed_3D_Molecule_Pair(m,&mA1,mB,ConnectGraph);

 Free_MOLECULE(&mA1); 

} /* end of Cal_MATCH_for_PCAfit_Superimposed_3D_Molecule_Pair() */







void Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat,Io,Jo,Ko,In,Jn,Kn)
  float Rmat[3][3];
  float Io[3],Jo[3],Ko[3];   /* original XYZ axis */
  float In[3],Jn[3],Kn[3];   /* new XYZ axis */
{
/*
      |Io0 Jo0 Ko0|
 Xo = |Io1 Jo1 Ko1|
      |Io2 Jo2 Ko2|

      |In0 Jn0 Kn0|
 Xn = |In1 Jn1 Kn1|
      |In2 Jn2 Kn2|

 We want to calculate the matrix R, satisfying R*Xo = Xn.
 Therefore,  R = Xn * inv[Xo] = Xn * tr[Xo].

                   |In0 Jn0 Kn0| |Io0 Io1 Io2|
 R = Xn * tr[Xo] = |In1 Jn1 Kn1|*|Jo0 Jo1 Jo2|
                   |In2 Jn2 Kn2| |Ko0 Ko1 Ko2|

*/

 Rmat[0][0] = In[0]*Io[0] + Jn[0]*Jo[0] + Kn[0]*Ko[0];
 Rmat[0][1] = In[0]*Io[1] + Jn[0]*Jo[1] + Kn[0]*Ko[1];
 Rmat[0][2] = In[0]*Io[2] + Jn[0]*Jo[2] + Kn[0]*Ko[2];

 Rmat[1][0] = In[1]*Io[0] + Jn[1]*Jo[0] + Kn[1]*Ko[0];
 Rmat[1][1] = In[1]*Io[1] + Jn[1]*Jo[1] + Kn[1]*Ko[1];
 Rmat[1][2] = In[1]*Io[2] + Jn[1]*Jo[2] + Kn[1]*Ko[2];

 Rmat[2][0] = In[2]*Io[0] + Jn[2]*Jo[0] + Kn[2]*Ko[0];
 Rmat[2][1] = In[2]*Io[1] + Jn[2]*Jo[1] + Kn[2]*Ko[1];
 Rmat[2][2] = In[2]*Io[2] + Jn[2]*Jo[2] + Kn[2]*Ko[2];

} /* end of Cal_Rotation_Matrix_from_Pair_of_XYZvec() */


void Rotate_Molecule_around_Gcenter(mol,Rmat,Gcen_orig,Gcen_new)
  struct MOLECULE *mol;
  float Rmat[3][3],Gcen_orig[3], Gcen_new[3];
{
  int a,i;
  float PosG[3];

  for (a=0;a<mol->Natom;++a){
    for (i=0;i<3;++i) PosG[i] = mol->atoms[a].Pos[i] - Gcen_orig[i];
    Mult_Mat_Vec3_float(mol->atoms[a].Pos,Rmat,PosG);
    for (i=0;i<3;++i) mol->atoms[a].Pos[i] += Gcen_new[i];
  }

} /* end of Rotate_Molecule_around_Gcenter() */



float Dot_Prod3(a,b)
 float a[3],b[3];
{
 return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void Write_EigVec_in_PDB(ofname,M)
 char *ofname;
 struct MOLECULE *M;
{
 FILE *fpo;
 float s;
 printf("#Write_EigVec_in_PDB()-->'%s'\n",ofname);
 fpo = fopen(ofname,"w");
 
 s = 4.0;
 fprintf(fpo,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
  1," CA ","CEN",'V',1, M->Gcen[0], M->Gcen[1],M->Gcen[2]);
 fprintf(fpo,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
  2," CA ","XXX",'V',2,M->Gcen[0]+s*M->eigvec[0][0], M->Gcen[1]+s*M->eigvec[0][1],M->Gcen[2]+s*M->eigvec[0][2]);
 fprintf(fpo,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
  3," CA ","YYY",'V',3,M->Gcen[0]+s*M->eigvec[1][0], M->Gcen[1]+s*M->eigvec[1][1],M->Gcen[2]+s*M->eigvec[1][2]);
 fprintf(fpo,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
  4," CA ","ZZZ",'V',4,M->Gcen[0]+s*M->eigvec[2][0], M->Gcen[1]+s*M->eigvec[2][1],M->Gcen[2]+s*M->eigvec[2][2]);
 fclose(fpo);

} /* end of Write_EigVec_in_PDB() */


void Assign_Direction_for_EigVec(M)
  struct MOLECULE *M;
 /*
  Count number of atoms for both side of EigVec.
  Assign "plus" direction for the side having more atoms.
 
 */ 
{
 int i,j,a,Npos,Nneg;
 float gvec[3];
 

 for (i=0;i<3;++i){
   Npos = Nneg = 0;
   for (a=0;a<M->Natom;++a){
     if (M->atoms[a].one_char_ele != 'H'){
       for (j=0;j<3;++j)  gvec[j] = M->atoms[a].Pos[j] - M->Gcen[j];
       if (Dot_Prod3(gvec,M->eigvec[i])>0.0) Npos += 1; else Nneg += 1;
     }
   }
 /* printf("axis [%d]  Npos %d Nneg %d\n",i,Npos,Nneg); */
  if (Npos<Nneg){ 
       for (j=0;j<3;++j)  M->eigvec[i][j] = -M->eigvec[i][j];
  }
 }

} /* end of Assign_Direction_for_EigVec() */




void Sort_Eigen_Value3(E,index)
 double E[3][3]; /* Diagonal matrix for eigen values (E[0][0],E[1][1],E[2][2])*/
 int    index[3]; /* Sort index by increaising order of E[0][0],E[1][1],E[2][2] */
{
 int i,min_i,max_i;
 double ev0,ev1,ev2;

 ev0 = E[0][0]; ev1 = E[1][1]; ev2 = E[2][2];

 /* printf("#ev0 %lf ev1 %lf ev2 %lf\n",ev0,ev1,ev2); */
 
 min_i = max_i = 0;

      if ((ev0<=ev1)&&(ev0<=ev2)) min_i = 0;
 else if ((ev1<=ev0)&&(ev1<=ev2)) min_i = 1;
 else if ((ev2<=ev0)&&(ev2<=ev1)) min_i = 2;

      if ((ev0>=ev1)&&(ev0>=ev2)&&(min_i!=0)) max_i = 0;
 else if ((ev1>=ev0)&&(ev1>=ev2)&&(min_i!=1)) max_i = 1;
 else if ((ev2>=ev0)&&(ev2>=ev1)&&(min_i!=2)) max_i = 2;

 index[0] = max_i;
 index[2] = min_i;
 for (i=0;i<3;++i)
   if ((i!=min_i)&&(i!=max_i)) index[1] = i; 

} /* end of Sort_Eigen_Value3() */




void Jacobi_Wilkinson3(A,U)
 double A[3][3]; /* Input matrix whose eigen value/vector will be calculated. 
      After calculation, diagonal elements A[0][0],A[1][1],A[2][2] become eigen values. */
 double U[3][3]; /* Output matrix corresponding to eigen vectors. 
   Eigen vector corresponds to each "column" vector of matrix U.
   For example k-th eigen vector 
   evec_k[0] =  U[0][k]; 
   evec_k[1] =  U[1][k]; 
   evec_k[2] =  U[2][k]; 
   */ 
{ 
 double R[3][3],TR[3][3],BUFF[3][3];
 double max,co,si;
 double wa,sa,r,sqrt2;
 int i,j,mi,mj,c;

 sqrt2 = sqrt(2.0);

   /* --------SETTING U[][] = I --------*/
   
   for (i=0;i<3;++i)
     for (j=0;j<3;++j) if (i==j) U[i][j] = 1.0;  else U[i][j] = 0.0;
    
    find_max_abs3_double(&mi,&mj,&max,A);
    c = 0;
    
    while ((max>0.00000001)&&(c<100)){
      /*--------- SET of cos sin -----------*/ 
      
      wa = (A[mi][mi] + A[mj][mj])/2;  
      sa = (A[mi][mi] - A[mj][mj])/2;  
      r = sqrt(sa*sa + A[mi][mj]*A[mi][mj]);
      co =  sqrt(1.0+sa/r)/sqrt2;  si =  A[mi][mj]/2.0/r/co; 
      if (sa>0.0) 
       { co =  sqrt(1.0+sa/r)/sqrt2;  si =  A[mi][mj]/2.0/r/co; }
      else 
       { co =   sqrt(1.0-sa/r)/sqrt2; si = - A[mi][mj]/2.0/r/co; }
      
      /* -------- SET of Rot Matrix R and TR----- */  
      for (i=0;i<3;++i)
	for (j=0;j<3;++j) if (i==j)  R[i][j] = 1.0;  else R[i][j] =0.0;
     
       R[mi][mi] =   co; R[mi][mj] = -si;
       R[mj][mi] =   si; R[mj][mj] = co;

      for (i=0;i<3;++i)
	for (j=0;j<3;++j) TR[j][i] = R[i][j]; 

       prod_matrix3_double(BUFF,A,R);
       prod_matrix3_double(A,TR,BUFF);
       prod_matrix3_double(BUFF,U,R);
       equal_matrix3_double(U,BUFF);

       find_max_abs3_double(&mi,&mj,&max,A);
       ++c;
 
      } /* end of while */

}/* end of Jacobi_Wilkinson3() */




void find_max_abs3_double(mi,mj,max,A) /* Max_Abs element A[mi][mj] = max */
 int *mi,*mj;
 double *max,A[3][3];
{ int i,j,p,q;

  p = q = 0;
  *max =0.0;
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      if (i!=j)
       if (fabs(A[i][j])> *max)
	      { *max = fabs(A[i][j]);     p = i; q= j;}  
   *mi =p; *mj = q;
} /* end of find_max_abs3() */


void equal_matrix3_double(A,B)    /* A = B */
 double A[3][3],B[3][3];
{ int i,j;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j)  A[i][j] = B[i][j];
   }
}/* end of equal_matrix3() */


void prod_matrix3_double(C,A,B) /* C = A * B */
 double C[3][3],A[3][3],B[3][3];
{ int i,j,k;
  double sum;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
       sum =0;
       for (k=0;k<3;++k)  sum = sum + A[i][k] * B[k][j];
       C[i][j] = sum; }
   }
 } /* end of prod_matrix3() */


void Mult_Mat_Vec3_float(y,mat,x)
 float y[3],mat[3][3],x[3];  /* y = mat * x (This function is only for 3-D) */
{ y[0] = mat[0][0]*x[0] + mat[0][1]*x[1] + mat[0][2] * x[2];
  y[1] = mat[1][0]*x[0] + mat[1][1]*x[1] + mat[1][2] * x[2];
  y[2] = mat[2][0]*x[0] + mat[2][1]*x[1] + mat[2][2] * x[2]; }


void print_matrix3_float(A,comment)
 float A[3][3];
 char *comment;
{ int i,j;
  printf("\n");
  printf("#%s\n",comment);
  printf("[");
  for (i=0;i<3;++i){ 
    printf("[");
    for (j=0;j<3;++j){ 
      printf(" %.3f",A[i][j]); 
      if (j==2) printf("]"); else printf(",");
     }
    if (i==2) printf("]"); else printf(",");
  }
  printf("\n");
}/* end of print_matrix3_float() */


void print_matrix3_double(A,comment)
 double A[3][3];
 char *comment;
{ int i,j;
  printf("\n");
  printf("#%s\n",comment);
  printf("[");
  for (i=0;i<3;++i){ 
    printf("[");
    for (j=0;j<3;++j){ 
      printf(" %.3lf",A[i][j]); 
      if (j==2) printf("]"); else printf(",");
     }
    if (i==2) printf("]"); else printf(",");
  }
  printf("\n");
}/* end of print_matrix3_double() */


void print_vector3_float(a,comment)
 float a[3];
 char *comment;
{ 
  printf("#%s [[%.3f],[%.3f],[%.3f]]\n",comment,a[0],a[1],a[2]);
}/* end of print_vector3() */


void set_unit_matrix3_double(A)
   double A[3][3];
{
  int i,j;
  for (i=0;i<3;++i){ 
   for (j=0;j<3;++j) if (i==j) A[i][j] = 1.0; else A[i][j] = 0.0;
 }
}



float Cal_Nmatch_for_Superimposed_3D_Molecule_Pair(mA,mB)
  struct MOLECULE *mA, *mB;
{
 int a,b,i,NmatchA, NmatchB;
 float del,DD,DDthre,Nmatch;

 DDthre = PAR.Tole_distan * PAR.Tole_distan;
 /* (1) Initialize */
 for (a=0;a<mA->Natom;++a) mA->atoms[a].mark = 0;
 for (b=0;b<mB->Natom;++b) mB->atoms[b].mark = 0;

 /* (2) Check atomtype and distance */
 for (a=0;a<mA->Natom;++a) { 
   if (mA->atoms[a].one_char_ele != 'H'){
     for (b=0;b<mB->Natom;++b){
       if (mB->atoms[b].one_char_ele != 'H'){
         if (strncmp(mA->atoms[a].atomtype,mB->atoms[b].atomtype,4)==0){
            DD = 0.0; 
            for (i=0;i<3;++i) { del = mA->atoms[a].Pos[i] - mB->atoms[b].Pos[i]; DD += del*del;}
            if (DD<=DDthre){                  
              mA->atoms[a].mark = 1;
              mB->atoms[b].mark = 1;
            }
         }
       }
     }
   }
 }

 /* (2) Check atomtype and distance */
 NmatchA = NmatchB = 0;
 for (a=0;a<mA->Natom;++a) NmatchA += mA->atoms[a].mark;
 for (b=0;b<mB->Natom;++b) NmatchB += mB->atoms[b].mark;
 Nmatch = (float)(NmatchA + NmatchB)/2.0; 
 /* 
 printf("#NmatchA %d NmatchB %d Nmatch %f\n",NmatchA,NmatchB,Nmatch);
 */ 
 return(Nmatch);

} /* end of Cal_Nmatch_for_Superimposed_3D_Molecule_Pair() */




void Cal_MATCH_for_Superimposed_3D_Molecule_Pair(m,mA,mB,ConnectGraph)
  struct MATCH    *m;
  struct MOLECULE *mA, *mB;
  char   ConnectGraph;   /* 'D':disconnected, otherwise:don't care */
{
 int a,b,aa,bb,min_aa,min_bb,i,round;
 float del,DD,DDthre,minD;
 float **Dmap; 
 char  hit, connect_ok; 

/*
 printf("#Cal_MATCH_for_Superimposed_3D_Molecule_Pair(NheavyA %d B %d)\n",mA->Nheavyatom, mB->Nheavyatom);
 */
 if (mA->Nheavyatom>mB->Nheavyatom){
   Malloc_MATCH(m,mA->Nheavyatom);
 }
 else{
   Malloc_MATCH(m,mB->Nheavyatom);
 }
 DDthre = PAR.Tole_distan * PAR.Tole_distan;

 /** [1] malloc and initialize Dmap[a][b] = -1.0 **/
 Dmap = (float **)malloc(sizeof(float *)*mA->Nheavyatom);
 for (aa=0;aa<mA->Nheavyatom;++aa){
   Dmap[aa] = (float *)malloc(sizeof(float)*mB->Nheavyatom);
 }

 for (aa=0;aa<mA->Nheavyatom;++aa) {
   for (bb=0;bb<mB->Nheavyatom;++bb) Dmap[aa][bb] = -1.0;
 }
 

 /** [2] calculate dmap only for distance(a,b)<=Dthre **/
 for (aa=0;aa<mA->Nheavyatom;++aa) { 
   for (bb=0;bb<mB->Nheavyatom;++bb){
       a = mA->num_frm_heavy[aa];
       b = mB->num_frm_heavy[bb];
       if (strncmp(mA->atoms[a].atomtype,mB->atoms[b].atomtype,4)==0){
         DD = 0.0; 
         for (i=0;i<3;++i) { del = mA->atoms[a].Pos[i] - mB->atoms[b].Pos[i]; DD += del*del;}
         if (DD <= DDthre) Dmap[aa][bb] = sqrt(DD);
      }
  }
 }

 /** [3] Find atom pairs by repeating to chose the closest pair **/
 m->Npair = 0;
 round = 0; 
 do{
   hit = 0; 
   min_aa = min_bb = 0;
   minD = -1.0;
   for (aa=0;aa<mA->Nheavyatom;++aa) { 
      for (bb=0;bb<mB->Nheavyatom;++bb){
          /* printf("aa %d bb %d Dmap %f\n",aa,bb,Dmap[aa][bb]); */
          if (Dmap[aa][bb]>=0.0){
            if ((minD<0.0) || (Dmap[aa][bb]<minD)){
              connect_ok = 1;
              /* connection check */ 
              if (ConnectGraph=='D'){
                a = mA->num_frm_heavy[aa];
                b = mB->num_frm_heavy[bb];
                for (i=0;i<m->Npair;++i){
                  if (((mA->conmap.map[m->anumA[i]][a]=='0')&&(mB->conmap.map[m->anumB[i]][b]!='0')) || 
                      ((mA->conmap.map[m->anumA[i]][a]!='0')&&(mB->conmap.map[m->anumB[i]][b]=='0')) ) {connect_ok = 0;}
                }
                /* printf("connection check %d %d --> connect_ok %d\n",a,b,connect_ok); */
              }

              if (connect_ok==1){ 
                minD = Dmap[aa][bb]; min_aa = aa; min_bb = bb; hit = 1;
              }
            }
          }
      }
   }


   if (hit==1){

     a = mA->num_frm_heavy[min_aa];
     b = mB->num_frm_heavy[min_bb];
     m->anumA[m->Npair] =  a;
     m->anumB[m->Npair] =  b;

    /*
     printf("#pair %d/%d (%d %d '%s')-(%d %d '%s') minD %f\n",
       m->Npair,m->Npair_malloc,min_aa,a,mA->atoms[a].atomname,min_bb,b,mB->atoms[b].atomname,minD); 

     printf("#(%f %f %f) (%f %f %f)\n",
       mA->atoms[a].Pos[0], mA->atoms[a].Pos[1], mA->atoms[a].Pos[2],
       mB->atoms[b].Pos[0], mB->atoms[b].Pos[1], mB->atoms[b].Pos[2]);
    */

     for (aa=0;aa<mA->Nheavyatom;++aa){Dmap[aa][min_bb] = -1.0;}
     for (bb=0;bb<mB->Nheavyatom;++bb){Dmap[min_aa][bb] = -1.0;}
     m->Npair += 1;
   }
   round += 1;

  } while (hit==1);


 /** [4] Free Dmap **/
 for (aa=0;aa<mA->Nheavyatom;++aa) free(Dmap[aa]);
 free(Dmap);

} /* end of Cal_MATCH_for_Superimposed_3D_Molecule_Pair() */




void Rotate_Molecule_XY_flat_plane_by_PCA(mol)
  struct MOLECULE *mol;
{
  float Rmat[3][3],Gcen_new[3];
  /*
  printf("#Rotate_Molecule_XY_flat_plane_by_PCA('%s')\n",mol->name); fflush(stdout);
  */
  Gcen_new[0] = Gcen_new[1] = Gcen_new[2] = 0.0;
  Cal_Molecule_EigVec_Covariance_Matrix(mol);
  make_rotational_matrix_XY_flat_plane(Rmat,mol->eigval,mol->eigvec);
  Rotate_Molecule_around_Gcenter(mol,Rmat,mol->Gcen,Gcen_new);
  mol->Gcen[0] = mol->Gcen[1] = mol->Gcen[2] = 0.0;
} /* end of Rotate_molecule_XY_flat_plane_by_PCA() */





void make_rotational_matrix_XY_flat_plane(Rmat,eval,evec)
  float Rmat[3][3];
  float eval[3];      /* eval[0] > eval[1] > eval[2] */
  float evec[3][3];   /* eval[0]--> evec[0][], eval[1]-->evec[1][], eval[2]-->evec[2][] */           
{
  /*
   x-axis : the largest  width  (largest  eigen value) eval0,evec0
   y-axis : the middle   width  (middle   eigen value) eval1,evec1
   z-axis : the smallest width  (smallest eigen value) eval2,evec2
 */
  int i,j;

  for (i=0;i<3;++i){
     for (j=0;j<3;++j){
      Rmat[i][j] = evec[i][j];
    }
  }
} /* end of  make_rotational_matrix_XY_flat_plane() */



