/*

< vector3.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

functions to for handling 3D vector


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalvar.h"

/** FUNCTION (GLOBAL) **/
void sub_vec3();
void copy_vec3();
float normalize_vec3();
float length_vec3();
void cross_vec3();
float dot_prod3();
float distance_vec3();
float bond_angle();
float dihedral_angle();
float angle_bwn_vec();
void make_rotmatrix_from_rot_axis_and_angle();
void rotate_vec_around_center();
void transform_vec_by_Rmat_oldO_newO();
void multi_mat_vec3();
void multi_mat3x3();
void make_transpose_matrix3x3();
void show_matrix3x3();
void set_identity_matrix3x3();




void sub_vec3(c,a,b)   /* c = a - b */
 float c[3],a[3],b[3];
{ c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2]; }


void copy_vec3(a,b)   /* a = b */
 float a[3],b[3];
{ a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2]; }


float normalize_vec3(v)
  float v[3];
{
  float norm;
  norm = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  norm = sqrt(norm);
  if (norm<=0.0) return(0.0);
  v[0] /= norm;  v[1] /= norm;  v[2] /= norm;  
  return(norm);
} 


float length_vec3(v)
  float v[3];
{
  float norm;
  norm = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  norm = sqrt(norm);
  if (norm<=0.0) return(0.0);
  return(norm);
} 


void cross_vec3(c,a,b)   /* c = a x b */
 float c[3],a[3],b[3];
{ c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0]; }



float dot_prod3(a,b)
 float a[3],b[3];
{ return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]); }



float distance_vec3(a,b)
 float a[3],b[3];
{ 
 float d;
 d =   (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]);
 return(sqrt(d));
}


float bond_angle(A,O,B)
 float A[3],O[3],B[3]; /* position of atoms A--O--B */
{
  float OA[3],OB[3];

  sub_vec3(OA,A,O);
  normalize_vec3(OA);
  sub_vec3(OB,B,O);
  normalize_vec3(OB);
  return(acos(dot_prod3(OA,OB)));
}





float dihedral_angle(A,B,C,D)
 float A[3],B[3],C[3],D[3]; /* position of atoms A,B,C,D */
 /* return dihedral angle in radian */
{
 float CB[3],CD[3],BC[3],BA[3],x[3],y[3],z[3];
 float prod_bc,prodCBz,dih;

 /*  A\   /D
       B-C    

 A  D
 | /
 BC
 
 x[]: normal vec for the plane A-B-C.
 y[]: normal vec for the plane B-C-D.

 x  y
 | /
 BC
 
 z = x x y

*/

 sub_vec3(CB,B,C);
 sub_vec3(CD,D,C);
 sub_vec3(BA,A,B);
 sub_vec3(BC,C,B);
 cross_vec3(x,BA,BC); normalize_vec3(x);
 cross_vec3(y,CB,CD); normalize_vec3(y);
 prod_bc = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
 if (prod_bc> 1.0) prod_bc =  1.0; 
 if (prod_bc<-1.0) prod_bc = -1.0;
 dih = acos(prod_bc);
 
 /*
 Because acos() cannot distinguigh theta and -theta,
 The dihedral angle may be -acos(prod_bc). 
 By the definition of dihedral angle, if the angle is positive,
 z = x x y should be parallel to the vector CB. 
 */ 
 cross_vec3(z,x,y);
 prodCBz = CB[0]*z[0] + CB[1]*z[1] + CB[2]*z[2];
 if (prodCBz>0.0){dih = -dih + 2.0*M_PI;}
 if (dih>=M_PI)  {dih -= 2.0*M_PI;}
 /*
 printf("\n");
 printf("#A %f %f %f\n",A[0],A[1],A[2]);
 printf("#B %f %f %f\n",B[0],B[1],B[2]);
 printf("#C %f %f %f\n",C[0],C[1],C[2]);
 printf("#D %f %f %f --> dif %f\n",D[0],D[1],D[2],dih);
 */ 
 return(dih);
} /* end of dihedral_angle() */



float angle_bwn_vec(A,B)
 float A[3],B[3];
{ float dot_prod,normA,normB,cosAB;
  normA    = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  normB    = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
  dot_prod = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  cosAB = dot_prod/sqrt(normA)/sqrt(normB);
  if (cosAB < 0.0) cosAB *= -1.0;
  return(acos(cosAB));
}




void make_rotmatrix_from_rot_axis_and_angle(R,axis,angle)
  float R[3][3];
  float axis[3];
  float angle; /* radian */
{
  float S[3][3];
  float SS[3][3];
  float a,b,c;
  float sint,ccost;
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
 /*
  printf("#R0 %f %f %f\n",R[0][0],R[0][1],R[0][2]);
  printf("#R1 %f %f %f\n",R[1][0],R[1][1],R[1][2]);
  printf("#R2 %f %f %f\n",R[2][0],R[2][1],R[2][2]);
 */

} /* end of make_rotmatrix_from_rot_axis_and_angle() */



void rotate_vec_around_center(x,R,cen)
  float x[3];
  float R[3][3];
  float cen[3]; 
{
  int i,j;
  float xx[3];
  
  for (i=0;i<3;++i) xx[i] = x[i];
  for (i=0;i<3;++i){
    x[i] = cen[i];
    for (j=0;j<3;++j){ 
      x[i] += R[i][j] * (xx[j]-cen[j]); 
    }
  }
}



void transform_vec_by_Rmat_oldO_newO(x,R,oldO,newO)
  float x[3];
  float R[3][3];
  float oldO[3],newO[3]; 
{
  int i,j;
  float xx[3];
  for (i=0;i<3;++i){ xx[i] = x[i];}

  for (i=0;i<3;++i){
    x[i] = newO[i];
    for (j=0;j<3;++j){ 
      x[i] += R[i][j] * (xx[j]-oldO[j]); 
    }
  }
}


void multi_mat_vec3(y,A,x) /* y = A x */
  float y[3],A[3][3],x[3];
{
  int i,j;
  for (i=0;i<3;++i){
    y[i] = 0.0;
    for (j=0;j<3;++j){ 
      y[i] += A[i][j] * x[j]; 
    }
  }
}


void multi_mat3x3(C,A,B) /* C = A B */
  float C[3][3],A[3][3],B[3][3];
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
}

void make_transpose_matrix3x3(TA,A)
  float TA[3][3],A[3][3];
{
  int i,j;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){ 
      TA[i][j] = A[j][i];
    }
  }
}



void show_matrix3x3(A, comment)
  float A[3][3];
  char *comment;
{
  printf("#>%s\n",comment);
  printf("#%8.3f %8.3f %8.3f\n",A[0][0],A[0][1],A[0][2]);
  printf("#%8.3f %8.3f %8.3f\n",A[1][0],A[1][1],A[1][2]);
  printf("#%8.3f %8.3f %8.3f\n",A[2][0],A[2][1],A[2][2]);
}

void set_identity_matrix3x3(A)
  float A[3][3];
{
  int i,j;
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      if (i==j) A[i][j] = 1.0;
           else A[i][j] = 0.0;
    }
  }
}
