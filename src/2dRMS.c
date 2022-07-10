/*

< 2dRMS.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 RMS calculation for 2D-MOLECULES


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "qRMS.h"

/*** FUNCTIONS (GLOBAL) ***/
float Calculate_CRMS_MATCH_rotation_Z();


/*** FUNCTIONS (LOCAL) ***/
static void setup_identity_matrix();
static void setup_Zrot_Rmat();
static void multi_mat_vec_double();
static void multi_mat3x3_double();




float Calculate_CRMS_MATCH_rotation_Z(match,mA,mB,gA,gB,Rmat)
 struct MATCH    *match; 
 struct MOLECULE *mA,*mB; 
 double gA[3],gB[3]; /* Center position */
 double Rmat[3][3];  /* Rotation Matrix */
{ 
  int i,j,x;
  float rms;
  double posA[3],posB[3],rposA[3];
  double A,B,sin_phi, cos_phi,norm,phi,theta,sqA,sqB,rms_ang; 
  double Rx[2][3][3]; /* Rx[0][][]:just identity matrix, Rx[1][][]:180-degree rot.around X-axis */  
  double Rz[2][3][3];  /* theta-rotation around Z-axis for two x values Rz[x=0][][] and Rz[x=1][][]. */
  double rms_x[2];
 
  printf("#Calculate_CRMS_MATCH_rotation_Z(match,mA,mB,gA,gB,Rmat)\n");
  /** [0] Initialize gA[],gB[], Rmat[][] **/
  gA[0] = gA[1] = gA[2] = gB[0] = gB[1] = gB[2] = 0.0;
  setup_identity_matrix(Rmat);
  setup_identity_matrix(Rx[0]);
  Rx[1][0][0] = 1.0; Rx[1][0][1] =  0.0; Rx[1][0][2] =  0.0; 
  Rx[1][1][0] = 0.0; Rx[1][1][1] = -1.0; Rx[1][1][2] =  0.0; 
  Rx[1][2][0] = 0.0; Rx[1][2][1] =  0.0; Rx[1][2][2] = -1.0; 

  if (match->Npair==0){
    return(0.0);
  }

  /** [1] Calculate gA[] and gB[] **/
  for (i=0;i<match->Npair;++i){
    for (j=0;j<3;++j){
      gA[j] += mA->atoms[match->anumA[i]].Pos[j];
      gB[j] += mB->atoms[match->anumB[i]].Pos[j];
    }
  }

  for (j=0;j<3;++j){
     gA[j] /= match->Npair;
     gB[j] /= match->Npair;
  }

  for (x=0;x<2;++x){
    /** [x-2] Calculate A and B **/
    A = 0.0;  B = 0.0; 
    posA[2] = posB[2] = 0.0;  
    sqA = sqB = 0.0;
    for (i=0;i<match->Npair;++i){
      for (j=0;j<2;++j){
        posA[j] = mA->atoms[match->anumA[i]].Pos[j] - gA[j];
        posB[j] = mB->atoms[match->anumB[i]].Pos[j] - gB[j];
      }
      multi_mat_vec_double(rposA,Rx[x],posA);
      A   += posB[0]*rposA[0] + posB[1]*rposA[1]; 
      B   += posB[1]*rposA[0] - posB[0]*rposA[1];  
      sqA += rposA[0]*rposA[0] + rposA[1]*rposA[1];  
      sqB += posB[0]*posB[0] + posB[1]*posB[1];  
    }

    /** [x-3] Calculate sin_p and cos_p **/
    /* printf("#sqA %lf sqB %lf A %lf B %lf\n",sqA,sqB,A,B); */
    norm = sqrt(A*A + B*B);
    if (norm == 0.0){ return(0);}
    sin_phi = A/norm; 
    cos_phi = B/norm;

    /** [x-4] Calculate phi and theta **/
    phi = acos(cos_phi);
    if (sin_phi < 0.0){
      phi = -phi;
    } 
     /*
    printf("#sin_phi %lf(%lf) cos_phi %lf(%lf) phi %lf\n",sin_phi,sin(phi), cos_phi, cos(phi),phi);
     */
    theta = M_PI/2.0 - phi;
    setup_Zrot_Rmat(Rz[x],theta);
    /* print_matrix3(Rz[x],"Rz");  */

    /** [x-5] Calculate RMSD **/
    multi_mat3x3_double(Rmat,Rz[x],Rx[x]);
    rms_x[x] = 0.0; 
    posA[2] = posB[2] = 0.0;  
    for (i=0;i<match->Npair;++i){
      for (j=0;j<2;++j){
        posA[j] = mA->atoms[match->anumA[i]].Pos[j] - gA[j];
        posB[j] = mB->atoms[match->anumB[i]].Pos[j] - gB[j];
      }
      multi_mat_vec_double(rposA,Rmat,posA);
      for (j=0;j<3;++j){
        rms_x[x]  += (rposA[j] - posB[j]) * (rposA[j] - posB[j]);
      }
    }
    rms_x[x] = sqrt(rms_x[x]/match->Npair); 
    rms_ang  = sqrt( (sqA+sqB - 2.0*norm*sin(theta+phi))/match->Npair); 
    printf("#x %d rms_x %lf phi %lf theta %lf phi+theta %lf rms_phi_theta %lf\n",
      x,rms_x[x],phi,theta,phi+theta,rms_ang);
 } 
  /** [6] Take the smallest RMSD **/
  if (rms_x[0] < rms_x[1]){
    multi_mat3x3_double(Rmat,Rz[0],Rx[0]);
    rms = rms_x[0];
  }
  else{
    multi_mat3x3_double(Rmat,Rz[1],Rx[1]);
    rms = rms_x[1];
  }

/*
  printf("#gA %lf %lf %lf gB %lf %lf %lf\n",gA[0],gA[1],gA[2],gB[0],gB[1],gB[2]);
  print_matrix3(Rmat,"Rmat");  
 */
  return(rms);
} /* end of Calculate_CRMS_MATCH_rotation_Z() */


void setup_identity_matrix(R)
  double R[3][3];
{
  R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0; 
  R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0; 
  R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0; 
} /* end of setup_identity_matrix() */


void setup_Zrot_Rmat(R,theta)
  double R[3][3];
  double theta;
{
  double sin_th, cos_th;
  sin_th = sin(theta);
  cos_th = cos(theta);
  R[0][0] = cos_th; R[0][1] = -sin_th; R[0][2] = 0.0; 
  R[1][0] = sin_th; R[1][1] =  cos_th; R[1][2] = 0.0; 
  R[2][0] = 0.0;    R[2][1] = 0.0;     R[2][2] = 1.0; 
} /* end of setup_Zrot_Rmat() */

void multi_mat_vec_double(y,A,x) /* y = A x */
  double y[3],A[3][3],x[3];
{
  int i,j;
  for (i=0;i<3;++i){
    y[i] = 0.0;
    for (j=0;j<3;++j){
      y[i] += A[i][j] * x[j];
    }
  }
} /* end of multi_mat_vec_double() */


void multi_mat3x3_double(C,A,B) /* C = A B */
  double C[3][3],A[3][3],B[3][3];
{
  int i,j,k;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      C[i][j] = 0.0;
      for (k=0;k<3;++k){
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
} /* end of multi_mat3x3_double() */

