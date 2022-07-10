/*
 
 <disgeo2D.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 Functions for making 2D(XY)-coordinates for the UNIT
 with complicated ring structure.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "molprop.h"
#include "ioLINE.h"
#include "ioSDF.h"
#include "moltopodis.h"
#include "gen2D.h"


/*** FUNCTIONS (GLOBAL) ***/
int set_XY_to_UNIT_by_distance_geometry();
float restraint_energy_with_derivative();
float restraint_energy();

/*** FUNCTIONS (LOCAL) ***/
static int find_two_atoms_in_regular_polygons();
static void Jacobi_Wilkinson();
static void find_max_abs_nondiagonal();
static void equal_matrix();
static void prod_matrix();
static void find_max_diagonal();
static void bubble_sort_diagonal();



int set_XY_to_UNIT_by_distance_geometry(mol,unit,HeadRing)
  struct MOLECULE *mol;
  struct UNIT *unit;
  struct UNIT *HeadRing;  /* Head of pseudo-SSSR for the "unit" */
{
  struct FLOAT2DMAP Dmap,Bmap,Umap;
  int i,j,ii,jj,ix,iy,iz,iw,*index;
  float *Aix,*Axj,Axx,lambdaX,lambdaY,lambdaZ,lambdaW;
  float *dEx,*dEy,E0,E1,StepSize,noise_magni;
  int Nrepeat,Nrepeat_max; 
  char end;

  printf("#set_XY_to_UNIT_by_distance_geometry(unit->num %d Natom %d)[",unit->num,unit->Natom);
  for (i=0;i<unit->Natom;++i) {
   printf(" %d%s", mol->atoms[unit->num_atoms[i]].num_in_file,mol->atoms[unit->num_atoms[i]].atomname);
  }  
  printf("]\n");
 

  Malloc_FLOAT2DMAP(&Dmap,unit->Natom); 
  Malloc_FLOAT2DMAP(&Bmap,unit->Natom); 
  Malloc_FLOAT2DMAP(&Umap,unit->Natom); 
  Aix = (float *)malloc(sizeof(float)*unit->Natom);
  Axj = (float *)malloc(sizeof(float)*unit->Natom);
  index = (int *)malloc(sizeof(int)*unit->Natom);
  dEx = (float *)malloc(sizeof(float)*unit->Natom);
  dEy = (float *)malloc(sizeof(float)*unit->Natom);

 /** (1) Set up Dmap by topological distance map **/
  for (i=0;i<unit->Natom;++i){
    ii = unit->num_atoms[i];
    for (j=0;j<unit->Natom;++j){
      jj = unit->num_atoms[j];
      Dmap.map[i][j] = PAR.LengthOfBond2D * (float)mol->topodismap.map[ii][jj];
/*
      printf("(%d %s)--(%d %s) : %f\n",
mol->atoms[ii].num_in_file, mol->atoms[ii].atomname,
mol->atoms[jj].num_in_file, mol->atoms[jj].atomname,Dmap.map[i][j]);
 */
    }
  }

  /** (2) Prepare Bmap by Young-Householder transformation **/
  Axx = 0.0;
 
  for (i=0;i<unit->Natom;++i){
    Aix[i] = 0.0;
    Axj[i] = 0.0;
    for (j=0;j<unit->Natom;++j){
      Aix[i] += Dmap.map[i][j] * Dmap.map[i][j];
      Axj[i] += Dmap.map[j][i] * Dmap.map[j][i];
      Axx    += Dmap.map[i][j] * Dmap.map[i][j];
    } 
  }

  Axx = Axx/unit->Natom/unit->Natom;
  for (i=0;i<unit->Natom;++i){
    Aix[i] = Aix[i]/unit->Natom;
    Axj[i] = Axj[i]/unit->Natom;
  }

  for (i=0;i<unit->Natom;++i){
    for (j=0;j<unit->Natom;++j){
      Bmap.map[i][j] = -0.5*(Dmap.map[i][j]*Dmap.map[i][j] - Aix[i] - Axj[j] - Axx);
    }
  }

  /** (3) Caluculte Eigen values/vectors **/
 /*
  for (i=0;i<Bmap.N;++i){
    for (j=0;j<Bmap.N;++j){
      printf("B(before) %d %d %f\n",i,j,Bmap.map[i][j]);
    }
  } 
 */

  Jacobi_Wilkinson(&Bmap,&Umap);
 /*
  for (i=0;i<Bmap.N;++i){
    for (j=0;j<Bmap.N;++j){
      printf("B(after) %d %d %f\n",i,j,Bmap.map[i][j]);
    }
  }
  for (i=0;i<Umap.N;++i){
    for (j=0;j<Umap.N;++j){
      printf("U %d %d %f\n",i,j,Umap.map[i][j]);
    }
  }
  */

  bubble_sort_diagonal(index,&Bmap);
  find_max_diagonal(&ix,&iy,&lambdaX,&lambdaY,&Bmap);
 /* 
  printf("ix %d lambdaX %f %f\n",ix,lambdaX,Bmap.map[ix][ix]);
  printf("iy %d lambdaY %f %f\n",iy,lambdaY,Bmap.map[iy][iy]);
  for (i=0;i<Bmap.N;++i){
    printf("[%d] lambda %f\n",index[i],Bmap.map[index[i]][index[i]]);
  } 
 */

  ix = index[1]; iy = index[2]; iz = index[3]; iw = index[4];

  /** (4) Set XY from eigen vectors with the first and second largest eigen value */
  lambdaX = sqrt(Bmap.map[ix][ix]); 
  lambdaY = sqrt(Bmap.map[iy][iy]); 
  lambdaZ = sqrt(Bmap.map[iz][iz]); 
  lambdaW = sqrt(Bmap.map[iw][iw]);  
  noise_magni = 0.05;

  for (i=0;i<unit->Natom;++i){
    unit->pos[i][0] = lambdaX * Umap.map[i][ix];
    /*  
    unit->pos[i][1] = lambdaY * Umap.map[i][iy];  
    unit->pos[i][1] = 0.9*lambdaY * Umap.map[i][iy] + 0.1*lambdaZ * Umap.map[i][iz];  
    */
    unit->pos[i][1] = 0.9*lambdaY * Umap.map[i][iy] + 0.1*lambdaZ * Umap.map[i][iz] 
                      + noise_magni * (-1.0 + 2.0*rand()/RAND_MAX);  

    /* printf("[%d] X %f %f Y %f %f\n",i,unit->pos[i][0],Umap.map[i][ix],unit->pos[i][1],Umap.map[i][iy]); */
  }

  printf("#restraint_energy %f\n",restraint_energy(mol,unit,HeadRing));

 /* Write_SDF_UNIT("dg_init.sdf",unit,mol,"");  */
  
  /*** (5) Steepest descent optimization ***/
  Nrepeat_max = 10 * unit->Natom;
  Nrepeat = 0;
  StepSize = 0.01 * PAR.LengthOfBond2D;  
  end = 0; 
  while ((Nrepeat < Nrepeat_max)&&(end==0)){
    restraint_energy_with_derivative(&E0,dEx,dEy,mol,unit,HeadRing,StepSize);
    for (i=0;i<unit->Natom;++i){
      unit->pos[i][0] -= dEx[i] * StepSize * unit->Natom;
      unit->pos[i][1] -= dEy[i] * StepSize * unit->Natom;
    }
    E1 = restraint_energy(mol,unit,HeadRing);
    printf("#[%d] E0 %f E1 %f\n",Nrepeat,E0,E1);
    if (E1>E0){
      end = 1;
      for (i=0;i<unit->Natom;++i){
        unit->pos[i][0] += dEx[i] * StepSize * unit->Natom;
        unit->pos[i][1] += dEy[i] * StepSize * unit->Natom;
      }
    } 
    Nrepeat += 1;
  }

  unit->set_localXY = 'T';
  
  /* Write_SDF_UNIT("dg_fin.sdf",unit,mol,"");  */
  
  free(dEx);
  free(dEy);
  free(Axj);
  free(Aix);
  Free_FLOAT2DMAP(&Dmap); 
  Free_FLOAT2DMAP(&Bmap); 
  Free_FLOAT2DMAP(&Umap); 
  free(index);
  return(1);
} /* end of set_XY_to_UNIT_by_distance_geometry() */



float restraint_energy(mol,unit,HeadRing)
  struct MOLECULE *mol;
  struct UNIT *unit;
  struct UNIT *HeadRing;
{
  int i,j,ii,jj;
  float E,dx,dy,dd,ddtar,LL,weight;
  int Natom_poly, topodiff;
  E = 0.0;
  LL = PAR.LengthOfBond2D * PAR.LengthOfBond2D;
 
  for (i=0;i<unit->Natom;++i){
    ii = unit->num_atoms[i];
    for (j=i+1;j<unit->Natom;++j){
      jj = unit->num_atoms[j];
      dx = unit->pos[i][0]-unit->pos[j][0]; 
      dy = unit->pos[i][1]-unit->pos[j][1]; 
      dd = dx*dx + dy*dy;
      ddtar = -1.0;
      if (mol->conmap.map[ii][jj]!='0'){ ddtar = LL; weight = 3.0; }
      else{
        find_two_atoms_in_regular_polygons(&Natom_poly,&topodiff,ii,jj,HeadRing);
        if (Natom_poly>=3){
/*
           printf("%d %s - %d %s Natom_poly %d topodiff %d\n",
             mol->atoms[ii].num_in_file, mol->atoms[ii].atomname,
             mol->atoms[jj].num_in_file, mol->atoms[jj].atomname,Natom_poly, topodiff);
*/
          if (topodiff==2){
            ddtar = 2.0 * LL * (1 + cos(2*M_PI/Natom_poly));
            weight = 1.0;
          }
           else if ((topodiff==3)&&(Natom_poly==6)){ 
            ddtar = 4.0 * LL;
            weight = 1.0;
          }
        }
      }

      if (ddtar > 0.0){
        E += weight * (dd - ddtar) * (dd - ddtar);  
      }
    }
  }

  return(E);
} /* end of restraint_energy() */



float restraint_energy_with_derivative(E,dEx,dEy,mol,unit,HeadRing)
  float  *E;   /* energy ( to be calculated) */
  float  *dEx; /* derivative by x [unit->Natom] (to be calculated) */
  float  *dEy; /* derivative by y [unit->Natom] (to be calculated) */
  struct MOLECULE *mol;
  struct UNIT *unit;
  struct UNIT *HeadRing;
{
  int i,j,ii,jj;
  float dx,dy,dd,ddtar,LL,weight,norm;
  int Natom_poly, topodiff;
 
  for (i=0;i<unit->Natom;++i) dEx[i] = dEy[i] = 0.0;  
 
  *E = 0.0;
  LL = PAR.LengthOfBond2D * PAR.LengthOfBond2D;
 
  for (i=0;i<unit->Natom;++i){
    ii = unit->num_atoms[i];
    for (j=0;j<unit->Natom;++j){
      if (i!=j){
        jj = unit->num_atoms[j];
        dx = unit->pos[i][0]-unit->pos[j][0]; 
        dy = unit->pos[i][1]-unit->pos[j][1]; 
        dd = dx*dx + dy*dy;
        ddtar = -1.0;
        if (mol->conmap.map[ii][jj]!='0'){ ddtar = LL; weight = 3.0;}
        else{
          find_two_atoms_in_regular_polygons(&Natom_poly,&topodiff,ii,jj,HeadRing);
          if (Natom_poly>=3){
            if (topodiff==2){
              ddtar = 2.0 * LL * (1 + cos(2*M_PI/Natom_poly));
              weight = 1.0;
            }
            else if ((topodiff==3)&&(Natom_poly==6)){ 
              ddtar = 4.0 * LL;
              weight = 1.0;
            }
          }
        }

        if (ddtar > 0.0){
          if (i<j) *E += weight * (dd - ddtar) * (dd - ddtar);  
          dEx[i] += 4.0*weight*(dd-ddtar)*dx;    
          dEy[i] += 4.0*weight*(dd-ddtar)*dy;    
        }
      }
    }
  }

  /** Normalization of [dEx,dEy] **/
  norm = 0.0;
  for (i=0;i<unit->Natom;++i){
    norm += dEx[i]*dEx[i] + dEy[i]*dEy[i];
  }
  norm = sqrt(norm);
  for (i=0;i<unit->Natom;++i){
     dEx[i] = dEx[i]/norm;
     dEy[i] = dEy[i]/norm;
  }

  
  return(*E);

} /* end of restraint_energy_with_derivative() */







int find_two_atoms_in_regular_polygons(Natom_poly, topodiff,anum1,anum2,HeadRing)
  int *Natom_poly;   /* Natom for the ring including anum1 and anum2 (to be calculated) */
  int *topodiff;    /* topological difference bwn anum1 and anum2 (1,2,3..)  */
  int anum1, anum2;  /* atom number 1 and 2 (given input) */
  struct UNIT *HeadRing; /* (given input) */
  /* return (Natom of the regular polygon including anum1 and anum2) */
{
  struct UNIT *rn;
  int i, i1, i2, t1,t2;

  *Natom_poly = *topodiff = 0;
  rn = HeadRing;
  
  while (rn->next != NULL){
    rn = rn->next;
    i1 = i2 = -1;
    for (i=0;i<rn->Natom;++i){
      if (rn->num_atoms[i]==anum1) i1 = i;
      if (rn->num_atoms[i]==anum2) i2 = i;
    }
    if ((i1>=0)&&(i2>=0)){
      *Natom_poly = rn->Natom;
      if (i1>=i2){ 
        t1 = abs(i1-i2)%rn->Natom;
        t2 = abs(i2+rn->Natom-i1)%rn->Natom;}
      else if (i1<i2){ 
        t1 = abs(i2-i1)%rn->Natom;
        t2 = abs(i1+rn->Natom-i2)%rn->Natom; }
      
       if (t1<t2) *topodiff = t1; else *topodiff = t2; 
 
       if ((i1==0)&&(i2==(rn->Natom-1))) *topodiff = 1;
       if ((i2==0)&&(i1==(rn->Natom-1))) *topodiff = 1;
      return(rn->Natom);
    }
  } 

  return(0);
} /* end of find_two_atoms_in_regular_polygons() */



void Jacobi_Wilkinson(A,U)
 struct FLOAT2DMAP *A; /* Input matrix whose eigen value/vector will be calculated. 
      After calculation, it becomes diagonal matrix of eigen values. */
 struct FLOAT2DMAP *U; /* Output matrix corresponding to eigen vectors. 
   Eigen vector corresponds to each "column" vector of matrix U.
   For example k-th eigen vector 
   evec_k[0] =  U[0][k]; 
   evec_k[1] =  U[1][k]; 
   evec_k[2] =  U[2][k]; 
   evec_k[3] =  U[3][k]; 
   */
{
 float max,co,si;
 float wa,sa,r,sqrt2,eps;
 int i,j,mi,mj,c,maxc;
 struct FLOAT2DMAP R,TR,BUFF;

 printf("#Jacobi_Wilkinson(N %d)\n",A->N);
 sqrt2 = sqrt(2.0);
 Malloc_FLOAT2DMAP(&R,A->N); 
 Malloc_FLOAT2DMAP(&TR,A->N); 
 Malloc_FLOAT2DMAP(&BUFF,A->N); 

   /* --------SETTING U[][] = I --------*/

 for (i=0;i<U->N;++i){
   for (j=0;j<U->N;++j){if (i==j) U->map[i][j] = 1.0;  else U->map[i][j] = 0.0;}
 }

 maxc = 1000;
 eps  = 0.00001;
 find_max_abs_nondiagonal(&mi,&mj,&max,A);
 c = 0;

 while ((max>eps)&&(c<maxc)){
   /*--------- SET of cos sin -----------*/

   wa = (A->map[mi][mi] + A->map[mj][mj])/2.0;
   sa = (A->map[mi][mi] - A->map[mj][mj])/2.0;
   r = sqrt(sa*sa + A->map[mi][mj]*A->map[mi][mj]);
   co =  sqrt(1.0+sa/r)/sqrt2;  si =  A->map[mi][mj]/2.0/r/co;
   if (sa>0.0)
    { co =  sqrt(1.0+sa/r)/sqrt2;  si =  A->map[mi][mj]/2.0/r/co; }
    else
    { co =   sqrt(1.0-sa/r)/sqrt2; si = - A->map[mi][mj]/2.0/r/co; }

   /* -------- SET of Rot Matrix R and TR----- */
   for (i=0;i<A->N;++i){
     for (j=0;j<A->N;++j){if (i==j)  R.map[i][j] = 1.0;  else R.map[i][j] =0.0;}
   }
    R.map[mi][mi] =   co; R.map[mi][mj] = -si;
    R.map[mj][mi] =   si; R.map[mj][mj] = co;

   for (i=0;i<A->N;++i){
     for (j=0;j<A->N;++j) TR.map[j][i] = R.map[i][j];
   }

   prod_matrix(&BUFF,A,&R);
   prod_matrix(A,&TR,&BUFF);
   prod_matrix(&BUFF,U,&R);
   equal_matrix(U,&BUFF);

   find_max_abs_nondiagonal(&mi,&mj,&max,A);
   /* printf("#c %d max %f\n",c,max); */
   ++c;

  } /* end of while */

 Free_FLOAT2DMAP(&R);
 Free_FLOAT2DMAP(&TR);
 Free_FLOAT2DMAP(&BUFF);

}/* end of Jacobi_Wilkinson() */


void find_max_abs_nondiagonal(mi,mj,max,A) /* Max_Abs element A[mi][mj] = max */
 int *mi,*mj;
 float *max;
 struct FLOAT2DMAP *A;
{ int i,j,p,q;
 *max =0.0;
  p = q = 0;
  for (i=0;i<A->N;++i){
    for (j=0;j<A->N;++j){
      if (i!=j){
       if (fabs(A->map[i][j])> *max)
              { *max = fabs(A->map[i][j]);     p = i; q= j;}
      }
    }
  }
 *mi =p; *mj = q;
} /* end of find_max_abs_nondiagonal() */



void equal_matrix(A,B)    /* A = B */
  struct FLOAT2DMAP *A, *B;
{ int i,j;
  for (i=0;i<A->N;++i){
    for (j=0;j<A->N;++j)  A->map[i][j] = B->map[i][j];
   }
}/* end of equal_matrix() */


void prod_matrix(C,A,B) /* C = A * B */
  struct FLOAT2DMAP  *C, *A, *B;
{ int i,j,k;
  float sum;
  for (i=0;i<A->N;++i){
    for (j=0;j<A->N;++j){
       sum =0;
       for (k=0;k<A->N;++k)  sum = sum + A->map[i][k] * B->map[k][j];
       C->map[i][j] = sum; }
   }
} /* end of prod_matrix() */



void find_max_diagonal(i1,i2,max1,max2,A) /* find maximum A[i][i] */
 int *i1,*i2;
 float *max1, *max2;
 struct FLOAT2DMAP *A;
{ int i;
  if (A->map[0][0]>A->map[1][1]){ 
   *max1 = A->map[0][0];
   *max2 = A->map[1][1];
   *i1 = 0; *i2 = 1;
  }
  else{
   *max1 = A->map[1][1];
   *max2 = A->map[0][0];
   *i1 = 1; *i2 = 0;
  }

  for (i=2;i<A->N;++i){
    if (A->map[i][i] > *max1) { 
      *max2 = *max1;        
      *i2 = *i1;
      *max1 = A->map[i][i];
      *i1 = i;
    }
    else if (A->map[i][i] > *max2){ 
      *max2 = A->map[i][i];
      *i2 = i;
    }
    /* printf("%d %f max1 %f max2 %f\n",i,A->map[i][i],*max1, *max2); */
  }
} /* end of find_max_diagonal() */



void bubble_sort_diagonal(index,A)
 int     *index;
 struct FLOAT2DMAP *A;
{
  int i,Nchange,b;
  for (i=0;i<A->N;++i) index[i] = i;
 
  Nchange = 0;
  do{
    Nchange = 0;
    for (i=1;i<A->N;++i){
      if (A->map[index[i-1]][index[i-1]] < A->map[index[i]][index[i]]){
        b = index[i];
        index[i] = index[i-1];
        index[i-1] = b;
        Nchange += 1;
      }
    } 
  }while (Nchange>0);
} /* end of bubble_sort_diagonal() */
