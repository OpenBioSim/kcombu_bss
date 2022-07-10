/*

< stamp_transf.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


functions to transform molecular conformation
by stamping conformations

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string.h>
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "globalvar.h"
#include "molprop.h"
#include "vector3.h"
#include "transform.h"
#include "energies.h"
#include "PCAfit.h"
#include "buildup.h"
#include "moltopodis.h"
#include "ioSDF.h"
#include "ioPDB.h"

/*** FUNCTIONS (GLOBAL) ***/
int Stamp_RotAngle_for_Reference_Dihedral_Angle();
int Stamp_sp3_Chirality();
int Check_3D_sp3_Chirality();
int Modify_MATCH_For_3D_sp3_Chirality_with_Permutation();
int Stamp_Five_or_Six_Atom_Ring_Geometry();
int Check_Five_or_Six_Atom_Ring_Geometry();
void Superimpose_Just_Two_Gcenters();
int Random_Change_sp3_Chirality();
int Random_Change_Flip_Of_Fragments();
int Change_Flip_Of_Fragment();
int Check_Bond_Angle_Difference();

/*** FUNCTIONS (LOCAL) ***/
static void set_four_neighbor_atom_geometry();
static float chirality_three_neighbor_atoms();
static void chose_rotatable_bonded_atoms_among_4neighbors();
static void random_chose_rotatable_bonded_atoms_among_4neighbors();
static void exchange_positions_of_two_atom_groups();
static int  check_ring_overlap();
static int check_ring_atom_overlap();
static void set_ring_num_atoms_in_rotating_order();
static void make_rotmatrix_by_atom_triplet_matching();
static int  check_equal_ring_chirality();
static void set_RingConnectAtoms_array();
static int check_ring_planarity();
static int check_neighbor_atoms_planarity();
static void reflection();
static void rot_matrix_of_axisAO_for_atomBtoC();
static void make_b_on_a_from_MATCH();
static int overlap_permuted_index_heavy_and_marked_atoms();
static void make_permu_index_from_permu_index_heavy();
static void show_atoms_in_two_rings();
static void write_pos6_in_pdb();




int Stamp_RotAngle_for_Reference_Dihedral_Angle(rnum,traP,M,molP,molQ)
  int    rnum;              /* number of rotational bond (0..tra->Nrbond-1) */
  struct TRANSFORM_VARIABLE *traP;
  struct MATCH *M; 
  struct MOLECULE *molP; /* target    molecule */
  struct MOLECULE *molQ; /* reference molecule */
{
  int pz,po,pzz,poo;
  int qz,qo,qzz,qoo;
  int i,j,m;
  float dihP,dihQ;
  char hit;
  /*
   molP(target):

     pzz          poo
        \         /
         \       /
          pz==>po

   molQ (reference): 

     qzz          qoo
        \         /
         \       /
          qz==>qo

 
    pz,qz : start atoms of rotating bond
    po,qo : end atoms   of rotating bond (origin of rotation)
  */

  /*** [1] set-up (po,pz) and (qo,qz) ***/

  traP->rbAngle[rnum] = 0.0;

  qz = qo = qzz = qoo = pzz = poo =  -1;
  pz = traP->rbZanum[rnum];
  po = traP->rbOanum[rnum];
  for (m=0;m<M->Npair;++m){
    if (M->anumA[m]==pz){ qz = M->anumB[m]; }
    if (M->anumA[m]==po){ qo = M->anumB[m]; }
  } 

  if ((qz<0) || (qo<0)){ return(0);}

  /* if bond qz-qo is not rotatble, exit. */
  if (molQ->conmap.map[qz][qo]!='R') return(0);
 
 
  /*** [2] find pzz and qzz (the pzz- and qzz-connected atom pairs, first appeared in M->anumA, anumB ) ***/
  hit = 0;
  i = 0; 
  while ((i<molP->atoms[pz].Nnei_heavy) && (hit==0)){
    pzz = molP->atoms[pz].neighbors[i];  
    j = 0;
    while ((j<molQ->atoms[qz].Nnei_heavy) && (hit==0)){
      qzz = molQ->atoms[qz].neighbors[j];
      for (m=0;m<M->Npair;++m){
        if ((M->anumA[m]==pzz) && (M->anumB[m]==qzz)
            &&(pzz!=pz)&&(pzz!=po)&&(qzz!=qz)&&(qzz!=qo)){
        hit = 1; m = M->Npair + 1;}
      }
      if (hit==0) j += 1;
    }    
    i += 1; 
  }

  if (hit==0) return(0);



  /*** [3] find poo and qoo  (the poo- and qoo-connected atom pairs, first appeared in M->anumA, anumB ) ***/
  hit = 0;
  i = 0; 
  while ((i<molP->atoms[po].Nnei_heavy) && (hit==0)){
    poo = molP->atoms[po].neighbors[i];  
    j = 0;
    while ((j<molQ->atoms[qo].Nnei_heavy) && (hit==0)){
      qoo = molQ->atoms[qo].neighbors[j];
      for (m=0;m<M->Npair;++m){
        if ((M->anumA[m]==poo) && (M->anumB[m]==qoo)
            &&(poo!=pz)&&(poo!=po)&&(qoo!=qz)&&(qoo!=qo))
            {hit = 1; m = M->Npair + 1;}
      }
      if (hit==0) j += 1;
    }    
    i += 1; 
  }
 
  if (hit==0) return(0);


  /* [4] Stamp dihedral angle with dihQ - dihP !!*/ 
  if ((pzz>=0) && (poo>=0) && (qzz>=0) && (qoo>=0)){
  
    dihP = dihedral_angle(molP->atoms[pzz].Pos, molP->atoms[pz ].Pos,
                          molP->atoms[po ].Pos, molP->atoms[poo].Pos);
    dihQ = dihedral_angle(molQ->atoms[qzz].Pos, molQ->atoms[qz ].Pos,
                          molQ->atoms[qo ].Pos, molQ->atoms[qoo].Pos);

    printf("#STAMP_RBOND[%2d] P(%4d-%4d-%4d-%4d) Q(%4d-%4d-%4d-%4d)",rnum,
       molP->atoms[pzz].num_in_file, molP->atoms[pz].num_in_file,
       molP->atoms[po].num_in_file, molP->atoms[poo].num_in_file,
       molQ->atoms[qzz].num_in_file, molQ->atoms[qz].num_in_file,
       molQ->atoms[qo].num_in_file, molQ->atoms[qoo].num_in_file);
    printf("dihP %+6.1f dihQ %+6.1f dihQ-P:rotating angle %+6.1f)\n",180.0*dihP/M_PI, 180.0*dihQ/M_PI,180.0*(dihQ-dihP)/M_PI);

    traP->rbAngle[rnum] = dihQ - dihP; 
    traP->rbStamp[rnum] = 'S';

  } 

 return(1);

} /* end of  Stamp_RotAngle_for_Reference_Dihedral_Angle() */







int Stamp_sp3_Chirality(molA, molB,M)
  struct MOLECULE *molA; /* target    */
  struct MOLECULE *molB; /* reference */
  struct MATCH *M;
{
 /*

  >> CAUTION <<
   Set_Topological_Distance() and Set_Rotatable_Bond()
   should be performed before this function.
 
 */ 
  int a,b,m,i,aa,k;
  int *b_on_a;
  int  Nmatch_nei, Nstamp;
  int  match_neiA[4],match_neiB[4];  /* atom number of four neighboring atoms.  If no neighobor match_neiA[n] = -1. */
  int  ch_i,ch_j; /* {0,1,2,3} */
  float rA[4][3],chiralA,rB[4][3],chiralB;
  float max_chiral_plane; /* angstrom */

  max_chiral_plane = 0.5;
  ch_i = ch_j = -1;

  printf("#Stamp_sp3_Chirality(Npair %d);\n",M->Npair);
  if (M->Npair==0) return(0);

 /** (1) prepare b_on_a[anum_in_molA] := anum_in_molB for matched atom pairs **/
  b_on_a = (int *)malloc(sizeof(int)*molA->Natom);
  make_b_on_a_from_MATCH(b_on_a, M, molA);

  /** (2) Check each atom pair, and if the pair satisfies the following condition 
          change the chirality of molecule A 
          (i)   a.Nneighbor>=3 and b.Nneighbor>=3.
          (ii)  At least three of the neighobors are matched.
          (iii) focused atoms and neighbor atoms are not on the same plane (|chirality|>chiral_tole angstrom)
          (iv)  Signs of chirality are different.
 **/

  Nstamp = 0; 
  for (m=0;m<M->Npair;++m){

    a = M->anumA[m];
    b = M->anumB[m];

    if ((molA->atoms[a].Nneighbor>=3)&&(molB->atoms[b].Nneighbor>=3)){
      Nmatch_nei = 0;
      for (k=0;k<4;++k) match_neiA[k] = match_neiB[k] = -1;

      for (i=0;i<molA->atoms[a].Nneighbor;++i){
        aa = molA->atoms[a].neighbors[i];
        if (b_on_a[aa]>=0){
          if (Nmatch_nei<4){
            match_neiA[Nmatch_nei] = aa;
            match_neiB[Nmatch_nei] = b_on_a[aa];
          }
          Nmatch_nei += 1;
        }
      }

     if ((Nmatch_nei==3)||(Nmatch_nei==4)){
        set_four_neighbor_atom_geometry(molA,a,match_neiA[0],match_neiA[1],match_neiA[2],match_neiA[3],rA,&chiralA);
        set_four_neighbor_atom_geometry(molB,b,match_neiB[0],match_neiB[1],match_neiB[2],match_neiB[3],rB,&chiralB);

        /** Do chiratity change  **/
        if ((chiralA*chiralB<0.0) && (fabs(chiralA)>max_chiral_plane)&& (fabs(chiralB)>max_chiral_plane)){

         chose_rotatable_bonded_atoms_among_4neighbors(molA,a,Nmatch_nei,match_neiA,&ch_i,&ch_j);

         if ((ch_i>=0)&&(ch_j>=0)){

           printf("#STAMP_SP3_CHIRAL A (%2s %4d %4d %4d %4d) Nnei %d B (%2s %4d %4d %4d %4d) Nnei %d Nmatch %d chiralA %f B %f exchange atoms %d and %d.\n",
            molA->atoms[a].element,
            molA->atoms[a].num_in_file,
            molA->atoms[match_neiA[0]].num_in_file, molA->atoms[match_neiA[1]].num_in_file, molA->atoms[match_neiA[2]].num_in_file,
            molA->atoms[a].Nneighbor,
            molA->atoms[b].element,
            molB->atoms[b].num_in_file,
            molB->atoms[match_neiB[0]].num_in_file, molB->atoms[match_neiB[1]].num_in_file, molB->atoms[match_neiB[2]].num_in_file, 
            molB->atoms[b].Nneighbor,
            Nmatch_nei, chiralA,chiralB, molA->atoms[match_neiA[ch_i]].num_in_file, molA->atoms[match_neiA[ch_j]].num_in_file);

            Nstamp += 1;
            exchange_positions_of_two_atom_groups(molA,a,Nmatch_nei,match_neiA,rA,ch_i, ch_j);

         } /* ch_i and ch_j */
        } /* chirality conditions */
      } /* Nmatch_nei >=3 */
    } /* Nneighbor >= 3 */
  }  /* m */

 free(b_on_a);
 return(Nstamp);
} /* end of Stamp_sp3_Chirality() */
















void set_four_neighbor_atom_geometry(mol,numP,num0,num1,num2,num3,rvec,chiral)
  struct MOLECULE *mol;          /* molecule    (given) */
  int numP,num0,num1,num2,num3;  /* atom number (given)*/
   /* numP: for focused atom, num0,num1,num2:for neighbor atoms of the focused atom */
  float rvec[4][3];     /* four bond vectors (to be calulated)*/
  float *chiral;        /* chirality (to be calculated). 
                         If chiral <0.0 then stereo_parity = 1, otherwise stereo_parity = 2. */  
{
  float nvec01[3];
  int k;
/*
  p: focused atom.
  x0, x1, x2, x3 : neighbor atoms
  r0, r1, r2, r3 : directional vector from the focused atoms to the neighbor atom
  
  ri = (xi-p)/|xi-p|  

  [position]   [edge vector] 
      x3          r3
      |            |
  x0--p--x2    r0--p--r2
      |            |
      x1          r1

  nvec01 = normalized(r0 cross r1)

  [chirality] = (r0 cross r1) dot r2 

  *The chiratity can be regarded as a distance between the atom n2 and the plane with normal vec b0 x b1  !!

  If the atom x2 on    the plane r0 x r1, the chirality has positive value.
  If the atom x2 under the plane r0 x r1, the chirality has negative value.

  If the signs of chirality are different, chirarities are different.

  ** The atom x3 is not necessary for deciding chirarity.

*/

  if (num0>=0) sub_vec3(rvec[0],mol->atoms[num0].Pos, mol->atoms[numP].Pos);
  if (num1>=0) sub_vec3(rvec[1],mol->atoms[num1].Pos, mol->atoms[numP].Pos);
  if (num2>=0) sub_vec3(rvec[2],mol->atoms[num2].Pos, mol->atoms[numP].Pos);
  if (num3>=0) sub_vec3(rvec[3],mol->atoms[num2].Pos, mol->atoms[numP].Pos);

  if ((num0>=0) && (num1>=0) && (num2>=0)){
    cross_vec3(nvec01,rvec[0],rvec[1]);
    normalize_vec3(nvec01);
    *chiral = dot_prod3(nvec01,rvec[2]);

    normalize_vec3(rvec[0]);
    normalize_vec3(rvec[1]);
    normalize_vec3(rvec[2]);

    if (num3>=0) sub_vec3(rvec[3],mol->atoms[num3].Pos, mol->atoms[numP].Pos);
   /* If num3 is not assigned(a hydrogen atom may be missing), make a DUMMY atom */
    else{
      for (k=0;k<3;++k) {rvec[3][k] = -(rvec[0][k] + rvec[1][k] + rvec[2][k]);}
    }
    normalize_vec3(rvec[3]);
  }


} /* end of set_four_neighbor_atom_geometry() */




float chirality_three_neighbor_atoms(mol,numP,num0,num1,num2)
  struct MOLECULE *mol;      /* molecule    (given) */
  int numP,num0,num1,num2;  /* atom number (given)*/
   /* numP: for focused atom, num0,num1,num2:for neighbor atoms of the focused atom */
{
  float nvec01[3];
  float rvec[3][3];     /* four bond vectors (to be calulated)*/
  float chiral;
/*
  p: focused atom.
  x0, x1, x2, (x3) : neighbor atoms
  r0, r1, r2, (r3) : directional vector from the focused atoms to the neighbor atom
  
  ri = (xi-p)/|xi-p|  

  [position]   [edge vector] 
      x0          r0
      |            |
  x2--p--x1    r2--p--r1
      |            |
     (x3)         (r3)

  nvec01 = normalized(r0 cross r1)

  [chirality] = nvec01 dot r2 

  *The chiratity can be regarded as a distance between the atom x2 and the plane with normal vec r0 x r1  !!

  If the atom x2 on    the plane r0 x r1, the chirality has positive value.
  If the atom x2 under the plane r0 x r1, the chirality has negative value.

  If the signs of chirality are different, chirarities are different.

  ** The atom x3 is not necessary for deciding chirarity.

*/

  sub_vec3(rvec[0],mol->atoms[num0].Pos, mol->atoms[numP].Pos);
  sub_vec3(rvec[1],mol->atoms[num1].Pos, mol->atoms[numP].Pos);
  sub_vec3(rvec[2],mol->atoms[num2].Pos, mol->atoms[numP].Pos);

  cross_vec3(nvec01,rvec[0],rvec[1]);
  normalize_vec3(nvec01);
  normalize_vec3(rvec[2]);
  chiral = dot_prod3(nvec01,rvec[2]);
  return(chiral);
} /* end of chirality_three_neighbor_atoms() */











void chose_rotatable_bonded_atoms_among_4neighbors(mol,o,Nnei_match,anum_nei,ch_i,ch_j)
  struct MOLECULE *mol; /* (given) */
  int o;                /* atom number for the focused atom (given) */ 
  int Nnei_match;       /* 3 or 4. (given). */
  int anum_nei[4];      /* atom numbers for neighbors */ 
  int *ch_i, *ch_j;     /* change i and j among {0,1,2,3} (to be calculated) */
  /*
    Bonds [o]-[ch_i], and [o]-[ch_j] should be 'rotatable'. 
  */ 
{
  int i,j;

  *ch_i = *ch_j = -1;

  if (Nnei_match==3){
/*
    for (i=0;i<3;++i){ 
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;} 
    }
    if (*ch_i>=0) *ch_j = 3;
*/
    for (i=0;i<3;++i){ 
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;} /* if o-i is rotatable */
    }
    for (j=0;j<3;++j){ 
      if ((j!=*ch_i)&&(mol->conmap.map[o][anum_nei[j]]=='R')) {*ch_j = j;} /* if o-j is rotatable */
     }

  }

  if (Nnei_match==4){
    for (i=0;i<4;++i){ 
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;} /* if o-i is rotatable */
    }
    for (j=0;j<4;++j){ 
      if ((j!=*ch_i)&&(mol->conmap.map[o][anum_nei[j]]=='R')) {*ch_j = j;} /* if o-j is rotatable */
     }
  } 

} /* end of chose_rotatable_bonded_atoms_among_4neighbors() */



void random_chose_rotatable_bonded_atoms_among_4neighbors(mol,o,Nnei,anum_nei,ch_i,ch_j)
  struct MOLECULE *mol; /* (given) */
  int o;                /* atom number for the focused atom */ 
  int Nnei;             /* 3 or 4. (given). */
  int anum_nei[4];      /* atom numbers for neighbors */ 
  int *ch_i, *ch_j;     /* change i and j among {0,1,2,3} (to be calculated) */
{
  int i,j,n,N;

  N = 10;

  if (Nnei==3){
    *ch_i = -1; 
    *ch_j = 3;  /* The last dummy atom is always chosen as 'ch_j' */
    n = 0;
    while ((n<N)&&((*ch_i)==-1)){
      i = rand()%3;
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;}  /* if o-i is rotatble */
      n += 1; 
    }
  }

  if (Nnei==4){
    *ch_i = *ch_j = -1;
    n = 0;
    while ((n<N)&&((*ch_i)==-1)){
      i = rand()%4;
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;} /* if a-i is rotatable */
      n += 1; 
    }
    n = 0;
    while ((n<N)&&((*ch_j)==-1)){
      j = rand()%4;
      if ((j!=*ch_i)&&(mol->conmap.map[o][anum_nei[j]]=='R')) {*ch_j = j;} /* if a-j is rotatable */
      n += 1; 
    } 
  }

} /* end of random_chose_rotatable_bonded_atoms_among_4neighbors() */








void exchange_positions_of_two_atom_groups(mol,o,Nnei,anum_nei,r,ch_i, ch_j)
  struct MOLECULE *mol;
  int o;            /* atom number for the focused atom */ 
  int Nnei;         /* 3 or 4. (given). */
  int anum_nei[4];  /* atom numbers for neighbors (given) */ 
  float r[4][3];    /* four directional vectors (given) */ 
  int ch_i, ch_j;   /* change i and j among {0,1,2,3} (given) */
{
/*
      4                   2
      |                   |
      j--3                i--1 
     /      2    =>      /      4 
    /      /            /      /
---o------i---1     ---o------j---3

 R : angle : rotation from i to j 
 RI:-angle : rotation from j to i

*/
  float rot_axis[3],angle,R[3][3],RI[3][3];
  int k;

  cross_vec3(rot_axis,r[ch_i],r[ch_j]);
  normalize_vec3(rot_axis);
  angle = acos(dot_prod3(r[ch_i],r[ch_j]));
  make_rotmatrix_from_rot_axis_and_angle(R, rot_axis, angle);
  make_rotmatrix_from_rot_axis_and_angle(RI,rot_axis,-angle);

  printf("#exchange_positions_of_two_atom_groups(o %d ch_i %d ch_j %d angle %f)\n",
    mol->atoms[o].num_in_file,mol->atoms[anum_nei[ch_i]].num_in_file,mol->atoms[anum_nei[ch_j]].num_in_file,angle); 

  for (k=0;k<mol->Natom;++k){
    if (k!=o){
    
 /* 
       printf("#ATOM %s %d o %d ch_i %d ch_j %d topodismap %d_vs_%d %d %d_vs_%d %d %d_vs_%d %d\n",
         mol->atoms[k].element,
         mol->atoms[k].num_in_file, 
         mol->atoms[o].num_in_file, 
         mol->atoms[anum_nei[ch_i]].num_in_file, 
         mol->atoms[anum_nei[ch_j]].num_in_file, 
         mol->atoms[k].num_in_file, mol->atoms[o].num_in_file, 
         mol->topodismap.map[k][o],
         mol->atoms[k].num_in_file, mol->atoms[anum_nei[ch_i]].num_in_file, 
         mol->topodismap.map[k][anum_nei[ch_i]], 
         mol->atoms[k].num_in_file, mol->atoms[anum_nei[ch_j]].num_in_file, 
         mol->topodismap.map[k][anum_nei[ch_j]]
        );
 */

       if (mol->topodismap.map[k][anum_nei[ch_i]] < mol->topodismap.map[k][o]){
          rotate_vec_around_center(mol->atoms[k].Pos,R,mol->atoms[o].Pos);
          /* printf("#atom %s%d ROT %f\n",mol->atoms[k].element,mol->atoms[k].num_in_file,angle);  */
       }
       else if ((mol->topodismap.map[k][anum_nei[ch_j]] < mol->topodismap.map[k][o])){
         rotate_vec_around_center(mol->atoms[k].Pos,RI,mol->atoms[o].Pos);
         /* printf("#atom %s%d IROT %f\n",mol->atoms[k].element,mol->atoms[k].num_in_file,-angle); */
       }

/*
       else if ((Nnei==4)&&(mol->topodismap.map[k][anum_nei[ch_j]] < mol->topodismap.map[k][o])){
         rotate_vec_around_center(mol->atoms[k].Pos,RI,mol->atoms[o].Pos);
          printf("#atom %s%d IROT %f\n",mol->atoms[k].element,mol->atoms[k].num_in_file,-angle);  
       }

       else if ((Nnei==3)&&(mol->conmap.map[o][k] !='0') && (mol->atoms[k].one_char_ele=='H') ){
         rotate_vec_around_center(mol->atoms[k].Pos,RI,mol->atoms[o].Pos);
          printf("#atom %s%d ROT %f\n",mol->atoms[k].element,mol->atoms[k].num_in_file,-angle);  
       }
*/
    }
  } 

} /* end of exchange_positions_of_two_atom_groups() */









int Random_Change_sp3_Chirality(mol,Pexchange)
  struct MOLECULE *mol;
  float  Pexchange;     /* Probability for exchange */
{
 /*

  >> CAUTION <<
   Set_Topological_Distance() and Set_Rotatable_Bond()
   should be performed before this function.
 
 */ 
  int a,aa,i,k,ch_i,ch_j, Nchange;
  float r[4][3],chiral,Prand;
  int Nnei,anum_nei[4];
  float max_chiral_plane; /* angstrom */
  max_chiral_plane = 0.5;

  Nchange = 0;
  ch_i = ch_j = -1;
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].Nneighbor>=3){
      Nnei = 0;
      for (k=0;k<4;++k) anum_nei[k] = -1;

      for (i=0;i<mol->atoms[a].Nneighbor;++i){
        aa = mol->atoms[a].neighbors[i];
        anum_nei[Nnei] = aa;
        Nnei += 1;
      }

      set_four_neighbor_atom_geometry(mol,a,anum_nei[0],anum_nei[1],anum_nei[2],anum_nei[3],r,&chiral);


      Prand = rand()/(float)RAND_MAX;

      if ((Prand <= Pexchange) && (chiral > max_chiral_plane)){ 
        printf("#CHIRATIY_CHANGE %d\n",mol->atoms[a].num_in_file);
        random_chose_rotatable_bonded_atoms_among_4neighbors(mol,a,Nnei,anum_nei,&ch_i,&ch_j);
        if ((ch_i>=0) && (ch_j>=0)){
          exchange_positions_of_two_atom_groups(mol,a,Nnei,anum_nei,r,ch_i, ch_j);
          Nchange += 1;
        }
      }
   }
  } /* a */

  printf("#Random_Change_sp3_Chirality(mol,Pexchange, Nchange %d)\n",Nchange);
  return(Nchange);

} /* end of Random_Change_sp3_Chirality() */





int Random_Change_Flip_Of_Fragments(mol,Pexchange)
  struct MOLECULE *mol;
  float  Pexchange;     /* Probability for exchange */
{
 /*

  >> CAUTION <<
   Set_Topological_Distance() and Set_Rotatable_Bond()
   should be performed before this function.
 
 */ 
  int Nchange,i,j; 
  int a,b,c,d,m,l;  /* atom number of atom A,B,C,D,E,F,G,H,M,L */
  struct RING  *rn;
  float CA[3],DB[3],angCA_DB,Prand; 
  int pathAB[100],Natom_pathAB;

  /*
  <Atom of fragment> 

   C-F1-D     : fixed
   A-M-F2-L-B : to be flipped

       E
    C--A--M
  /    F   \
  |         | 
 F1         F2
  |         | 
  \    G   / 
    D--B--L
       H
 
   X:mid_point of C-D.
   Y:mid_point of M-L.

  */

  Nchange = 0;
  rn = &(mol->HeadRing); 
  while (rn->next != NULL){
    rn = rn->next;
 /*
    printf("RING[%d] Natom %d plane %d {",rn->num,rn->Natom,check_ring_planarity(mol,rn->Natom,rn->num_atoms));
    for (i=0;i<rn->Natom;++i){
      printf("%s%d ",mol->atoms[rn->num_atoms[i]].element,mol->atoms[rn->num_atoms[i]].num_in_file);
    }
    printf("}\n");
 */

    if (check_ring_planarity(mol,rn->Natom,rn->num_atoms)==0){
    for (i=0;i<rn->Natom;++i){
      a = rn->num_atoms[i];
      c = rn->num_atoms[(i-1+rn->Natom)%rn->Natom];
      m = rn->num_atoms[(i+1+rn->Natom)%rn->Natom];
      for (j=i+2;j<rn->Natom;++j){
        b = rn->num_atoms[j];
        d = rn->num_atoms[(j+1+rn->Natom)%rn->Natom];
        l = rn->num_atoms[(j-1+rn->Natom)%rn->Natom];
        if ((mol->topodismap.map[a][b]>=2)){
          set_shortest_path_between_two_atoms(mol,a,b,pathAB,&Natom_pathAB);
          /*
          for (s=0;s<Natom_pathAB;++s){
            printf("#PATH %s%d\n",mol->atoms[pathAB[s]].element,mol->atoms[pathAB[s]].num_in_file); 
          }
          */
          m = pathAB[1];
          l = pathAB[Natom_pathAB-2];
         
          if (rn->num_atoms[(i-1+rn->Natom)%rn->Natom]==m) c =  rn->num_atoms[(i+1+rn->Natom)%rn->Natom];
                                                     else  c =  rn->num_atoms[(i-1+rn->Natom)%rn->Natom];
          
          if (rn->num_atoms[(j-1+rn->Natom)%rn->Natom]==l) d =  rn->num_atoms[(j+1+rn->Natom)%rn->Natom];
                                                     else  d =  rn->num_atoms[(j-1+rn->Natom)%rn->Natom];

          sub_vec3(CA,mol->atoms[a].Pos,mol->atoms[c].Pos); normalize_vec3(CA);
          sub_vec3(DB,mol->atoms[b].Pos,mol->atoms[d].Pos); normalize_vec3(DB);
          angCA_DB = acos(dot_prod3(CA,DB))*180.0/M_PI;
         
         /* 
          printf("#FLIP_CANDIDATE A:%s%d B:%s%d C:%s%d D:%s%d M:%s%d L:%s%d angCA_DB %f\n",
              mol->atoms[a].element, mol->atoms[a].num_in_file, mol->atoms[b].element, mol->atoms[b].num_in_file,
              mol->atoms[c].element, mol->atoms[c].num_in_file, mol->atoms[d].element, mol->atoms[d].num_in_file,
              mol->atoms[m].element, mol->atoms[m].num_in_file, mol->atoms[l].element, mol->atoms[l].num_in_file, angCA_DB);
           */
          Prand = (float)(rand())/RAND_MAX;

          if ( (Prand<=Pexchange)&&(angCA_DB <=PAR.minAngleFlipOfFragment)
               &&(check_neighbor_atoms_planarity(mol,a)==0)&&(check_neighbor_atoms_planarity(mol,b)==0)
               && (check_ring_atom_overlap(mol,rn,a)==0) && (check_ring_atom_overlap(mol,rn,b)==0)
               &&(mol->atoms[a].one_char_ele == 'C') && (mol->atoms[b].one_char_ele == 'C')){
            
              Change_Flip_Of_Fragment(mol,a,b,c,d,m,l);
              Nchange += 1; 
           }
        } 
      } /* j */ 
    } /* i */
   } /* plane ring */

  } /* rn */

  printf("#Random_Change_Flip_Of_Fragments(Pexchange %f Nchange %d)\n",Pexchange,Nchange);
  return(Nchange);

} /* end of Random_Change_Flip_Of_Fragments() */







int Change_Flip_Of_Fragment(mol,a,b,c,d,m,l)
  struct MOLECULE *mol;
  int a,b,c,d,m,l;  /* atom number of atom A,B,C,D,M,L */
{
 /*

  >> CAUTION <<
   Set_Topological_Distance() and Set_Rotatable_Bond()
   should be performed before this function.
 
 */ 
  int e,f,g,h;
  int k,s; 
  float CA[3],DB[3],angCA_DB; 
  float X[3],XA[3],XB[3],normABX[3]; 
  float Y[3],YA[3],YB[3],normABY[3]; 
  float PosMorig[3],PosLorig[3],Rmat[3][3];
  float bangCAMorig,bangCAMflip,bangDBLorig,bangDBLflip;
  /*
 
 Based on 
   
   Mekenyan,O;Pavlov,T;Grancharov,V;Todorov,M;Schmieder,P;Veith G. 
   2D-3D Migration of Large Chemical Inventories with Conformational Multiplication.
   J.Chem.Inf.Model., 2005, 283-292.   

   Vanio,M,J;Johnson, M.S. Generating Conformer Ensembles Using a Multiobjective Genetic Algorithm.
   J.Chem.Inf.Model., 2007, 47, 2462-2474.

  <Atom of fragment> 

   C-F1-D     : fixed
   A-M-F2-L-B : to be flipped

       E
    C--A--M
  /    F   \
  |         | 
 F1         F2
  |         | 
  \    G   / 
    D--B--L
       H
 
   X:mid_point of C-D.
   Y:mid_point of M-L.

  */


          
  sub_vec3(CA,mol->atoms[a].Pos,mol->atoms[c].Pos);
  normalize_vec3(CA);
  sub_vec3(DB,mol->atoms[b].Pos,mol->atoms[d].Pos);
  normalize_vec3(DB);
  angCA_DB = acos(dot_prod3(CA,DB))*180.0/M_PI;
          

   for (k=0;k<3;++k){ 
     PosMorig[k] = mol->atoms[m].Pos[k];
     PosLorig[k] = mol->atoms[l].Pos[k];
   }

    /* (1) set up protein e,f,g,h */
   e = f = g = h = -1; 
   for (s=0;s<mol->atoms[a].Nneighbor;++s){
      if ((mol->atoms[a].neighbors[s]!=c) && (mol->atoms[a].neighbors[s]!=m)){
        if (e<0) e = mol->atoms[a].neighbors[s];
            else f = mol->atoms[a].neighbors[s];
     }
   }
   for (s=0;s<mol->atoms[b].Nneighbor;++s){
     if ((mol->atoms[b].neighbors[s]!=d) && (mol->atoms[b].neighbors[s]!=l)){
       if (g<0) g = mol->atoms[b].neighbors[s];
           else h = mol->atoms[b].neighbors[s];
     }
   }


    bangCAMorig = bond_angle(mol->atoms[c].Pos, mol->atoms[a].Pos, mol->atoms[m].Pos)*180.0/M_PI;
    bangDBLorig = bond_angle(mol->atoms[d].Pos, mol->atoms[b].Pos, mol->atoms[l].Pos)*180.0/M_PI;

    /* (2) mark flipping atoms */
    for (s=0;s<mol->Natom;++s){
       mol->atoms[s].mark = 0;
       if ((s!=a)&&(s!=b)&&(s!=c)&&(s!=d)){

        if  ((mol->topodismap.map[s][m]<mol->topodismap.map[s][a])&&(mol->topodismap.map[s][l]<mol->topodismap.map[s][b])){
           mol->atoms[s].mark = 1; /* connected neighbors of atoms M and L */
         }

        if ((s!=m)&&(s!=l)&&(mol->topodismap.map[s][a]<mol->topodismap.map[s][c])&&(mol->topodismap.map[s][a]<mol->topodismap.map[s][m])){
          mol->atoms[s].mark = 2; /* connected neighbors of atom A */
         }

         if ((s!=m)&&(s!=l)&&(mol->topodismap.map[s][b]<mol->topodismap.map[s][d])&&(mol->topodismap.map[s][b]<mol->topodismap.map[s][l])){
           mol->atoms[s].mark = 3; /* connected neighbors of atom B */
         }
         /*
         if (mol->atoms[s].mark != 0){
           printf("  FLIPPING_ATOM %s%d mark %d\n",mol->atoms[s].element,mol->atoms[s].num_in_file,mol->atoms[s].mark);
         }
         */
       }
    } 

    /* (3) the first reflection using AB-X */
    for (k=0;k<3;++k) X[k] = (mol->atoms[c].Pos[k] + mol->atoms[d].Pos[k])/2.0;         
    sub_vec3(XA,mol->atoms[a].Pos,X); normalize_vec3(XA);
    sub_vec3(XB,mol->atoms[b].Pos,X); normalize_vec3(XB);
    cross_vec3(normABX,XA,XB); normalize_vec3(normABX);    
 
    for (s=0;s<mol->Natom;++s){
      if (mol->atoms[s].mark==1){ 
        reflection(mol->atoms[s].Pos,normABX,X);
      }
    } 
    /*
    Write_SDF_Molecule("flip1.sdf",mol,"1st flip");
     */         
    /* (4) the second reflection using AB-Y */
    for (k=0;k<3;++k) Y[k] = (mol->atoms[m].Pos[k] + mol->atoms[l].Pos[k])/2.0;
    sub_vec3(YA,mol->atoms[a].Pos,Y); normalize_vec3(YA);
    sub_vec3(YB,mol->atoms[b].Pos,Y); normalize_vec3(YB);
    cross_vec3(normABY,YA,YB);  normalize_vec3(normABY);   
    for (s=0;s<mol->Natom;++s){
      if (mol->atoms[s].mark==1){ 
        reflection(mol->atoms[s].Pos,normABY,Y);
      }
    } 
    /* (5) rotate atoms around the atom A */

    /*
    printf("#M %f %f %f -> %f %f %f\n",PosMorig[0],PosMorig[1],PosMorig[2],mol->atoms[m].Pos[0],mol->atoms[m].Pos[1],mol->atoms[m].Pos[2]);
    */ 
    rot_matrix_of_axisAO_for_atomBtoC(Rmat,mol->atoms[c].Pos, mol->atoms[a].Pos,PosMorig,mol->atoms[m].Pos);
    for (s=0;s<mol->Natom;++s){
      if (mol->atoms[s].mark==2){
        rotate_vec_around_center(mol->atoms[s].Pos,Rmat,mol->atoms[a].Pos);
      }
    } 

   /* (6) rotate atoms around the atom B */
/*
   printf("#L %f %f %f -> %f %f %f\n",PosLorig[0],PosLorig[1],PosLorig[2],mol->atoms[l].Pos[0],mol->atoms[l].Pos[1],mol->atoms[l].Pos[2]);
 */ 
    rot_matrix_of_axisAO_for_atomBtoC(Rmat,mol->atoms[d].Pos, mol->atoms[b].Pos,PosLorig,mol->atoms[l].Pos);
   for (s=0;s<mol->Natom;++s){
   if (mol->atoms[s].mark==3){
        rotate_vec_around_center(mol->atoms[s].Pos,Rmat,mol->atoms[b].Pos);
      }
   } 

  /*
   Write_SDF_Molecule("flip2.sdf",mol,"2nd flip");
   */
   bangCAMflip = bond_angle(mol->atoms[c].Pos, mol->atoms[a].Pos, mol->atoms[m].Pos)*180.0/M_PI;
   bangDBLflip = bond_angle(mol->atoms[d].Pos, mol->atoms[b].Pos, mol->atoms[l].Pos)*180.0/M_PI;
   printf("#Change_Flip_Of_Fragment A:%s%2d B:%s%2d C:%s%2d D:%s%2d M:%s%2d L:%s%2d angCA_DB %5.1f bangCAM %5.1f->%5.1f  DBL %5.1f->%5.1f\n",
              mol->atoms[a].element, mol->atoms[a].num_in_file,
              mol->atoms[b].element, mol->atoms[b].num_in_file,
              mol->atoms[c].element, mol->atoms[c].num_in_file,
              mol->atoms[d].element, mol->atoms[d].num_in_file,
              mol->atoms[m].element, mol->atoms[m].num_in_file, 
              mol->atoms[l].element, mol->atoms[l].num_in_file, angCA_DB,bangCAMorig,bangCAMflip,bangDBLorig,bangDBLflip);
   return(1);

} /* end of Change_Flip_Of_Fragment() */













void reflection(pos,norm,orig)
  float pos[3];  /* position to be reflected */
  float norm[3]; /* norm vector     of the reflection plane.*/
  float orig[3]; /* origin position of the reflection plane.*/
{
  float p_o[3],dprod;
  int k;
 /*
  p' = p - 2 [norm dot (pos-orig)] norm 
 */

  sub_vec3(p_o,pos,orig);
  dprod = dot_prod3(norm,p_o); 
  for (k=0;k<3;++k) pos[k] = pos[k] - 2.0 * dprod * norm[k];

} /* end of reflection() */






void rot_matrix_of_axisAO_for_atomBtoC(Rmat,A,O,B,C)
  float Rmat[3][3];
  float A[3],O[3],B[3],C[3];
{
  /*
        C 
       / 
   A==O   /|\ rangle(BtoC)
       \
        B
  */ 
  float OA[3],OB[3],OC[3],normAOB[3],normAOC[3],rot_axis[3],rangle,AOBxAOC[3];
  
  sub_vec3(OA,A,O); normalize_vec3(OA);
  sub_vec3(OB,B,O); normalize_vec3(OB);
  sub_vec3(OC,C,O); normalize_vec3(OC);
  cross_vec3(normAOB,OA,OB);  normalize_vec3(normAOB);
  cross_vec3(normAOC,OA,OC);  normalize_vec3(normAOC);
  cross_vec3(AOBxAOC,normAOB,normAOC); normalize_vec3(AOBxAOC);
  rangle = acos(dot_prod3(normAOB,normAOC));
  sub_vec3(rot_axis,O,A); normalize_vec3(rot_axis);
  /*
  printf("rangle %f\n",180.0*rangle/M_PI);
  printf("dprod raxis AOBxAOC %f\n",dot_prod3(rot_axis,AOBxAOC));
  */ 
  if (dot_prod3(rot_axis,AOBxAOC) <0.0) rangle = -rangle;
  make_rotmatrix_from_rot_axis_and_angle(Rmat, rot_axis, rangle);

} /* end of rot_matrix_of_axisAO_for_atomBtoC() */


int Check_Five_or_Six_Atom_Ring_Geometry(M,molA,molB)
  struct MATCH *M; 
  struct MOLECULE *molA; /* target    molecule */
  struct MOLECULE *molB; /* reference molecule */
{
  struct RING *rn;
  int *b_on_a;
  int i,Nmatch,Nring_stamp,Nring_diffcon;
  int Natom_ring;  /* 5 or 6 */
  int   num_atomsB[6]; 
 
  printf("#Check_Five_or_Six_Atom_Ring_Geometry()\n");
  Nring_diffcon = 0;

 /** [1] prepare b_on_a[anum_in_molA] := anum_in_molB for matched atom pairs **/
  b_on_a = (int *)malloc(sizeof(int)*molA->Natom);
  make_b_on_a_from_MATCH(b_on_a, M, molA);

 /** [2] Scan the rings **/
  Nring_stamp = 0;
  rn = &(molA->HeadRing); 
  while (rn->next != NULL){
    rn = rn->next;
    Nmatch = 0;
    Natom_ring = rn->Natom;
    if (((Natom_ring ==5)||(Natom_ring ==6))&&(check_ring_overlap(molA,rn)==0)&&
        (check_ring_planarity(molA,rn->Natom,rn->num_atoms)==0)){

        for (i=0;i<rn->Natom;++i){
          if (b_on_a[rn->num_atoms[i]]>=0){ 
            num_atomsB[Nmatch] = b_on_a[rn->num_atoms[i]];
            Nmatch += 1;
          }
        }

        if ((Nmatch == rn->Natom)&&(check_ring_planarity(molB,rn->Natom,num_atomsB)==0)){
          set_ring_num_atoms_in_rotating_order(molA,rn);
          if (check_equal_ring_chirality(molA,molB,rn->num_atoms,Natom_ring,b_on_a)==0){
            Nring_diffcon += 1;
          }
      }  /* Nmatch == rn->Natom) */
   } /* rn->Natom = 5 or 6 */
 } /* rn */

 printf("#Nring_stamp %d\n",Nring_stamp);
 free(b_on_a);

 return(Nring_diffcon);

} /* end of Check_Five_or_Six_Atom_Ring_Geometry() */






int Stamp_Five_or_Six_Atom_Ring_Geometry(M,molA,molB)
  struct MATCH *M; 
  struct MOLECULE *molA; /* target    molecule */
  struct MOLECULE *molB; /* reference molecule */
{
  struct RING *rn;
  int *b_on_a;
  int a,i,Nmatch,pre_i,next_i,Nring_stamp;
  int Natom_ring;  /* 5 or 6 */
  float R[3][3],oldO[3],newO[3];
  char  *RngConAtoms;
  int   num_atomsB[6]; 
  float rposA[6][3],rposB[6][3];

  printf("#Stamp_Five_or_Six_Atom_Ring_Geometry()\n");
 /** [1] prepare b_on_a[anum_in_molA] := anum_in_molB for matched atom pairs **/
  b_on_a = (int *)malloc(sizeof(int)*molA->Natom);
  make_b_on_a_from_MATCH(b_on_a, M, molA);

  RngConAtoms = (char *)malloc(sizeof(char)*molA->Natom);
  

 /** [2] Scan the rings **/
  Nring_stamp = 0;
  rn = &(molA->HeadRing); 
  while (rn->next != NULL){
    rn = rn->next;
    Nmatch = 0;
    Natom_ring = rn->Natom;
    if (((Natom_ring ==5)||(Natom_ring ==6))&&(check_ring_overlap(molA,rn)==0)&&
        (check_ring_planarity(molA,rn->Natom,rn->num_atoms)==0)){

        for (i=0;i<rn->Natom;++i){
          if (b_on_a[rn->num_atoms[i]]>=0){ 
            num_atomsB[Nmatch] = b_on_a[rn->num_atoms[i]];
            Nmatch += 1;
          }
        }

        if ((Nmatch == rn->Natom)&&(check_ring_planarity(molB,rn->Natom,num_atomsB)==0)){
          set_ring_num_atoms_in_rotating_order(molA,rn);

          if (check_equal_ring_chirality(molA,molB,rn->num_atoms,Natom_ring,b_on_a)==0){
            show_atoms_in_two_rings("TRANSFORM_RING_GEOMETRY", molA, molB,rn->num_atoms,Natom_ring,b_on_a);

            for (i=0;i<Natom_ring;++i) {
              copy_vec3(rposA[i],molA->atoms[rn->num_atoms[i]].Pos);
              copy_vec3(rposB[i],molB->atoms[b_on_a[rn->num_atoms[i]]].Pos);
            }

           /** (1) Transform molB by Tr[b(012)->a(012)] **/ 
           make_rotmatrix_by_atom_triplet_matching(R,oldO,newO,
               rposB[0],  rposB[1],  rposB[2],  rposA[0],  rposA[1],  rposA[2] ); 

           for (i=0;i<Natom_ring;++i){ transform_vec_by_Rmat_oldO_newO(rposB[i],R,oldO,newO); }
           /*
           write_pos6_in_pdb("init_rposA.pdb",rposA,Natom_ring,1);
           write_pos6_in_pdb("init_rposB.pdb",rposB,Natom_ring,2);
           */

           /** (2) Transform molA by Tr[a(i-1,i,i+1) --> b(i-1,i,i+1)] **/ 
           for (i=0;i<Natom_ring;++i) {
             if (i!=1){
               set_RingConnectAtoms_array(molA,i,rn->num_atoms,Natom_ring, RngConAtoms);
               if (i==0)              pre_i  = Natom_ring-1; else pre_i  = i - 1;
               if (i==(Natom_ring-1)) next_i = 0; else next_i = i + 1;
               make_rotmatrix_by_atom_triplet_matching(R,oldO,newO, 
                   rposA[pre_i], rposA[i], rposA[next_i], rposB[pre_i], rposB[i], rposB[next_i]);

               printf("#tranform  (%d-%d-%d) into (%d %d %d)\n",
                 molA->atoms[rn->num_atoms[pre_i]].num_in_file, molA->atoms[rn->num_atoms[i]].num_in_file, molA->atoms[rn->num_atoms[next_i]].num_in_file,
                 molB->atoms[b_on_a[rn->num_atoms[pre_i]]].num_in_file, molB->atoms[b_on_a[rn->num_atoms[i]]].num_in_file, molB->atoms[b_on_a[rn->num_atoms[next_i]]].num_in_file);

               for (a=0;a<molA->Natom;++a){
                 if (RngConAtoms[a]==1){
                   transform_vec_by_Rmat_oldO_newO(molA->atoms[a].Pos,R,oldO,newO);
                 }
               }
             }
             /*
             sprintf(ofname,"%dA.pdb",molA->atoms[rn->num_atoms[i]].num_in_file);
             Write_PDB_Molecule(ofname,molA,'w','A',1,-1,'F','F',"");
             */
           }

          Nring_stamp += 1;

        } /* chirarities differ */
      }  /* Nmatch == rn->Natom) */
   } /* rn->Natom = 5 or 6 */
 } /* rn */

 /*
 sprintf(ofname,"finA.pdb");
 Write_PDB_Molecule(ofname,molA,'w','A',1,-1,'F','F',"");
 */

 printf("#Nring_stamp %d\n",Nring_stamp);
 free(RngConAtoms);
 free(b_on_a);
 return(Nring_stamp);

} /* end of Stamp_Five_or_Six_Atom_Ring_Geometry() */









void Superimpose_Just_Two_Gcenters(molA,molB)
  struct MOLECULE *molA,*molB;
{
  float gA[3],gB[3];
  int a,b,k;

  gA[0] = gA[1] = gA[2] = 0.0;
  gB[0] = gB[1] = gB[2] = 0.0;

  for (a=0;a<molA->Natom;++a){
    for (k=0;k<3;++k) gA[k] += molA->atoms[a].Pos[k]; 
  }
  for (k=0;k<3;++k) gA[k] /= molA->Natom;

  for (b=0;b<molB->Natom;++b){
    for (k=0;k<3;++k) gB[k] += molB->atoms[b].Pos[k]; 
  }
  for (k=0;k<3;++k) gB[k] /= molB->Natom;


  for (a=0;a<molA->Natom;++a){
    for (k=0;k<3;++k)  molA->atoms[a].Pos[k] += (gB[k]-gA[k]); 
  }
  
} /* end of Superimpose_Just_Two_Gcenters() */





int check_ring_overlap(mol,rn)
  struct MOLECULE *mol;
  struct RING *rn;
{
  struct RING *sn;
  int i,j;

  sn = &(mol->HeadRing); 

  while (sn->next != NULL){
    sn = sn->next;
    if (sn->num != rn->num){
      for (i=0;i<rn->Natom;++i){
        for (j=0;j<sn->Natom;++j){
          if (rn->num_atoms[i]==sn->num_atoms[j]){return(1); /* overlap !! */}
        } 
      }
    } 
  }
 return(0); /* non-overlap */
} /* end of check_ring_overlap() */


int check_ring_atom_overlap(mol,rn,ao)
  struct MOLECULE *mol;
  struct RING *rn;
  int    ao;  /* atom number of focused atom belong to rn */
{
  struct RING *sn;
  int j;

  sn = &(mol->HeadRing); 

  while (sn->next != NULL){
    sn = sn->next;
    if (sn->num != rn->num){
      for (j=0;j<sn->Natom;++j){
        if (sn->num_atoms[j]==ao) { return(1); /* overlap !! */}
      } 
    } 
  }
 return(0); /* non-overlap */
} /* end of check_ring_atom_overlap() */





void set_ring_num_atoms_in_rotating_order(mol,rn)
  struct MOLECULE *mol;
  struct RING *rn;
{
  int *num_atoms_orig,i,j;
  char *used_atoms,hit;
 
  num_atoms_orig = (int *)malloc(sizeof(int) * rn->Natom);
  used_atoms     = (char *)malloc(sizeof(char) * rn->Natom);

  for (i=0;i<rn->Natom;++i) {
    num_atoms_orig[i] = rn->num_atoms[i];
    rn->num_atoms[i] = -1;
    used_atoms[i]     = 0;
  }

  rn->num_atoms[0] = num_atoms_orig[0]; 
  used_atoms[0] = 1;
 
  for (i=1;i<rn->Natom;++i){
    hit = 0;
    j = 0;
    while ((j<rn->Natom) && (hit==0)){
      if ((used_atoms[j]==0) && (mol->conmap.map[rn->num_atoms[i-1]][num_atoms_orig[j]]!='0')){
        rn->num_atoms[i] = num_atoms_orig[j];
        used_atoms[j] = 1; 
        hit = 1;
      }
      else {j += 1;} 
    } 
    if (hit==0){printf("#ERROR(set_ring_num_atoms_in_rotating_number):i %d\n",i); exit(1);}
  }

  free(used_atoms);
  free(num_atoms_orig);
} /* end of set_ring_num_atoms_in_rotating_order() */



int check_equal_ring_chirality(molA, molB,anumAarray,Natom_ring,b_on_a)
  struct MOLECULE *molA, *molB;
  int *anumAarray;  /* int[Natom_ring] */
  int Natom_ring;
  int *b_on_a;
{
  float nvA[3],nvB[3],vg0[3],vg1[3],vi1[3];
  float chiralA,chiralB,chiral_tole,gposA[3],gposB[3];
  int i,k,equal;
 /*
     4--5
    /    \
   0  g   3
    \    /
     1--2 

   Normal vec of the ring plane is generated by vg0 x vg1.

 */

  chiral_tole = 0.5;

  for (k=0;k<3;++k){gposA[k]=0.0; gposB[k]=0.0;} 
  for (i=0;i<Natom_ring;++i){
    for (k=0;k<3;++k){
      gposA[k] += molA->atoms[anumAarray[i]].Pos[k]; 
      gposB[k] += molB->atoms[b_on_a[anumAarray[i]]].Pos[k]; 
    } 
  }

  for (k=0;k<3;++k){
    gposA[k] /= Natom_ring; 
    gposB[k] /= Natom_ring;
  } 
  
  sub_vec3(vg0,molA->atoms[anumAarray[0]].Pos, gposA);
  sub_vec3(vg1,molA->atoms[anumAarray[1]].Pos, gposA);
  cross_vec3(nvA,vg0,vg1); normalize_vec3(nvA);

  sub_vec3(vg0,molB->atoms[b_on_a[anumAarray[0]]].Pos, gposB);
  sub_vec3(vg1,molB->atoms[b_on_a[anumAarray[1]]].Pos, gposB);
  cross_vec3(nvB,vg0,vg1); normalize_vec3(nvB);

  equal = 1;
  for (i=3;i<Natom_ring;++i){
    sub_vec3(vi1,molA->atoms[anumAarray[i]].Pos, gposA);
    chiralA = dot_prod3(nvA,vi1);
    sub_vec3(vi1,molB->atoms[b_on_a[anumAarray[i]]].Pos, gposB);
    chiralB = dot_prod3(nvB,vi1);

    if   (fabs(chiralA-chiralB)>chiral_tole){ 
       equal = 0;
       printf("#RING_CONF_DIFFERENCE %s %d chiral %f %s %d chiral %f\n",
          molA->atoms[anumAarray[i]].element,molA->atoms[anumAarray[i]].num_in_file,chiralA,
          molB->atoms[b_on_a[anumAarray[i]]].element,molB->atoms[b_on_a[anumAarray[i]]].num_in_file,chiralB);
    }
 }
 return(equal);
} /* end of check_equal_ring_chirality() */



void show_atoms_in_two_rings(title, molA, molB,anumAarray,Natom_ring,b_on_a)
  char  *title;
  struct MOLECULE *molA, *molB;
  int *anumAarray;  /* int[Natom_ring] */
  int Natom_ring;
  int *b_on_a;
{
  int i;

  printf("#%s ",title);
  printf("(");
  for (i=0;i<Natom_ring;++i){
     printf(" %s %d", molA->atoms[anumAarray[i]].element,molA->atoms[anumAarray[i]].num_in_file);
  }
  printf(")");
  printf("(");
  for (i=0;i<Natom_ring;++i){
     printf(" %s %d", molB->atoms[b_on_a[anumAarray[i]]].element,molB->atoms[b_on_a[anumAarray[i]]].num_in_file);
  }
  printf(")\n");
} /* end of show_atoms_in_two_rings() */








void make_rotmatrix_by_atom_triplet_matching(R, oldO, newO,posAx,posAo,posAy,posBx,posBo,posBy)
  float  R[3][3];   /* rotational matrix */ 
  float  oldO[3],newO[3];  /* oldO:origin for molA, newO:origin of molB */
  float posAx[3],posAo[3],posAy[3];   /* Atom Position of molecule A (to be rotated)   */
  float posBx[3],posBo[3],posBy[3];   /* Atom position of molecule B (fixed reference) */
{
  float xA[3],yA[3],zA[3]; 
  float xB[3],yB[3],zB[3]; 
  float XA[3][3],XB[3][3]; 
  /*
  XA, XB : matrix for x,y,z-axis vectors   
  
   XA[0],XB[0] : x-axis vector
   XA[1],XB[1] : y-axis vector
   XA[2],XB[2] : z-axis vector

   R = XB x tr[XA] 
  */ 
  float yyA[3],yyB[3];  
  float dprod;
  int i,j,k;

  /* (1) Calculate XA[0],XA[1],XA[2] **/
  for (k=0;k<3;++k) oldO[k] = posAo[k];
  sub_vec3(xA,posAx,oldO);
  normalize_vec3(xA);
  sub_vec3(yyA,posAy,oldO);
  dprod = dot_prod3(yyA,xA);
  for (k=0;k<3;++k) yA[k] = yyA[k] - dprod * xA[k];
  normalize_vec3(yA);
  cross_vec3(zA,xA,yA);
  
  XA[0][0] = xA[0]; XA[0][1] = yA[0]; XA[0][2] = zA[0];
  XA[1][0] = xA[1]; XA[1][1] = yA[1]; XA[1][2] = zA[1];
  XA[2][0] = xA[2]; XA[2][1] = yA[2]; XA[2][2] = zA[2];


  /* (2) Calculate XB[0],XB[1],XB[2] **/
  for (k=0;k<3;++k) newO[k] = posBo[k];
  sub_vec3(xB,posBx,newO);
  normalize_vec3(xB);
  sub_vec3(yyB,posBy,newO);
  dprod = dot_prod3(yyB,xB);
  for (k=0;k<3;++k) yB[k] = yyB[k] - dprod * xB[k];
  normalize_vec3(yB);
  cross_vec3(zB,xB,yB);
  XB[0][0] = xB[0]; XB[0][1] = yB[0]; XB[0][2] = zB[0];
  XB[1][0] = xB[1]; XB[1][1] = yB[1]; XB[1][2] = zB[1];
  XB[2][0] = xB[2]; XB[2][1] = yB[2]; XB[2][2] = zB[2];

  /* (3)  R = XB x tr[XA]  */
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      R[i][j] = 0.0;
      for (k=0;k<3;++k){
        R[i][j] = R[i][j] + XB[i][k] * XA[j][k];
     }
    }
  } 

} /* end of make_rotmatrix_by_atom_triplet_matching() */









void set_RingConnectAtoms_array(mol,focused_i,ring_anums,Natom_ring,ring_con_atoms)
  struct MOLECULE *mol;   /* (input) */
  int  focused_i;         /* focused_number(0,1,2,3,4,5) in the ring_anums (input) */
  int  *ring_anums;       /* [Natom_ring] other atoms in the rings (input) */
  int  Natom_ring;        /* number of atoms in the ring */
  char *ring_con_atoms;   /* [Natom] ring contact atoms array (to be calculated) */
{
  int a,i,Tfocus,minT;

  for (a=0;a<mol->Natom;++a){
    Tfocus = mol->topodismap.map[a][ring_anums[focused_i]]; 
    minT = mol->Natom;
    for (i=0;i<Natom_ring;++i){
      if ((i != focused_i)&&(ring_anums[i] != a)){
        if (mol->topodismap.map[a][ring_anums[i]]<minT){
          minT = mol->topodismap.map[a][ring_anums[i]];
        }
      } 
    } 
    if (Tfocus<minT){ring_con_atoms[a] = 1;} else {ring_con_atoms[a] = 0;} 
  } 

 /*
  printf("#[%d]",mol->atoms[ring_anums[focused_i]].num_in_file); 
  for (a=0;a<mol->Natom;++a){
    if (ring_con_atoms[a]==1){ printf(" %d",mol->atoms[a].num_in_file); }
  } 
  printf("\n");
  */
} /* end of set_RingConnectAtoms_array() */



int check_ring_planarity(mol,Natom_ring,num_atoms)
  struct MOLECULE *mol;
  int    Natom_ring;
  int    *num_atoms;
/*
  return(1) : if the ring is plane.
  return(0) : if the ring is not plane 
*/
{
  int i,o,a,b,c;
  float angle,angle_tole,OA[3],OB[3],OC[3],Nvec[3];
  int plane;

  angle_tole = 8.0; /** degree **/

  if (Natom_ring < 4){ return(1);}

  o = num_atoms[0];
  a = num_atoms[1];
  b = num_atoms[2];
  sub_vec3(OA,mol->atoms[a].Pos, mol->atoms[o].Pos);
  sub_vec3(OB,mol->atoms[b].Pos, mol->atoms[o].Pos);
  cross_vec3(Nvec,OA,OB);

  plane = 1;

  for (i=3;i<Natom_ring;++i){
    c = num_atoms[i];
    sub_vec3(OC,mol->atoms[c].Pos, mol->atoms[o].Pos);
    angle =  angle_bwn_vec(Nvec,OC)*180.0/M_PI;
    if (fabs(angle-90.0) > angle_tole) plane = 0;
    /* printf("angle %f\n",angle);  */
  }
  return(plane);
} /* end of check_ring_planarity() */



int check_neighbor_atoms_planarity(mol,o)
  struct MOLECULE *mol;
  int    o;  /* atom number of the focused atom */
/*
  return(1) : if atom o and neighbors are on the plane.
  return(0) :                         are not on the plane 
*/
{
  float angle,angle_tole,OA[3],OB[3],OC[3],OD[3],Nvec[3];
  int i,plane,a,b,c,d;
 /*
     c 
     |
  a--o--b
     |
     d
  ** Atom a and b should be heavy atoms. 
 */
 
  angle_tole = 8.0; /** degree **/
  a = b = c = d = -1;
  plane = 1;
 
  for (i=0;i<mol->atoms[o].Nneighbor;++i){
    if (mol->atoms[mol->atoms[o].neighbors[i]].one_char_ele != 'H'){
           if (a<0) a = mol->atoms[o].neighbors[i];
      else if (b<0) b = mol->atoms[o].neighbors[i];
      else if (c<0) c = mol->atoms[o].neighbors[i];
      else if (d<0) d = mol->atoms[o].neighbors[i];
    }
  } 

  sub_vec3(OA,mol->atoms[a].Pos, mol->atoms[o].Pos); normalize_vec3(OA); 
  sub_vec3(OB,mol->atoms[b].Pos, mol->atoms[o].Pos); normalize_vec3(OB);
  cross_vec3(Nvec,OA,OB); normalize_vec3(Nvec);

  if (c>=0){
    sub_vec3(OC,mol->atoms[c].Pos, mol->atoms[o].Pos);
    angle =  angle_bwn_vec(Nvec,OC)*180.0/M_PI;
    if (fabs(angle-90.0) > angle_tole) plane = 0;
  }

  if (d>=0){
    sub_vec3(OD,mol->atoms[d].Pos, mol->atoms[o].Pos);
    angle =  angle_bwn_vec(Nvec,OD)*180.0/M_PI;
    if (fabs(angle-90.0) > angle_tole) plane = 0;
  }

  return(plane);
} /* end of check_neighbor_atoms_planarity() */


  
void make_b_on_a_from_MATCH(b_on_a, M, molA)
  int *b_on_a;  /* int[mol->Natom] (should be already malloced) **/
  struct MATCH *M;
  struct MOLECULE *molA;
{
  int a,m;

  for (a=0;a<molA->Natom;++a){ 
    b_on_a[a] = -1; 
  }
  for (m=0;m<M->Npair;++m){
    if (M->anumA[m]<molA->Natom){
      b_on_a[M->anumA[m]] = M->anumB[m];
    }
  }

} /* end of make_b_on_a_from_MATCH() */





int Check_3D_sp3_Chirality(molA, molB,M, OutStdOut)
  struct MOLECULE *molA; /* target    */
  struct MOLECULE *molB; /* reference */
  struct MATCH *M;
  char   OutStdOut;     /* if 'T', output in stdout */
/*
 OUTPUT: 
  return (Dchiral) (Number of atoms with different 3D chirality).

  molA->atoms[a].mark :=1, molB->atoms[a].mark :=1 for the atom pairs (a,b) 
  if chiralities of (a,b) differ, or (a,b) are neighbors of atoms with different chiralities.

*/
{
  int a,b,m,i,aa,k,Nmatch_nei, Dchiral;
  int *b_on_a;
  int  match_neiA[4],match_neiB[4];  /* atom number of four neighboring atoms.  If no neighobor match_neiA[n] = -1. */
  float chiralA,chiralB;
  float max_chiral_plane;            /* angstrom */

  max_chiral_plane = 0.5;

  printf("#Check_3D_sp3_Chirality(Npair %d);\n",M->Npair); 
  if (M->Npair==0) return(0);

 /** (1) prepare b_on_a[anum_in_molA] := anum_in_molB for matched atom pairs **/
  b_on_a = (int *)malloc(sizeof(int)*molA->Natom);
  make_b_on_a_from_MATCH(b_on_a, M, molA);
  for (a=0;a<molA->Natom;++a){ molA->atoms[a].mark = 0; }
  for (b=0;b<molB->Natom;++b){ molB->atoms[b].mark = 0; }


  /** (2) Check each atom pair, and if the pair satisfies the following condition 
          change the chirality of molecule A 
          (i)   a.Nnei_heavy>=3 and b.Nnei_heavy>=3.
          (ii)  At least three of the neighobors are matched.
          (iii) focused atoms and neighbor atoms are not on the same plane (|chirality|>chiral_tole angstrom)
          (iv)  Signs of chirality are different.
 **/

  Dchiral = 0; 
  for (m=0;m<M->Npair;++m){
    a = M->anumA[m];
    b = M->anumB[m];

    if ((molA->atoms[a].Nnei_heavy>=3)&&(molB->atoms[b].Nnei_heavy>=3)){
      Nmatch_nei = 0;
      for (k=0;k<4;++k) match_neiA[k] = match_neiB[k] = -1;

      for (i=0;i<molA->atoms[a].Nneighbor;++i){
        aa = molA->atoms[a].neighbors[i];
        if (molA->atoms[a].one_char_ele != 'H'){
          if (b_on_a[aa]>=0){
            if (Nmatch_nei<4){
              match_neiA[Nmatch_nei] = aa;
              match_neiB[Nmatch_nei] = b_on_a[aa];
            }
            Nmatch_nei += 1;
          }
        }
      }

     if ((Nmatch_nei==3)||(Nmatch_nei==4)){
        chiralA = chirality_three_neighbor_atoms(molA,a,match_neiA[0],match_neiA[1],match_neiA[2]);
        chiralB = chirality_three_neighbor_atoms(molB,b,match_neiB[0],match_neiB[1],match_neiB[2]);

        /** Chiratities are different  **/
        if ((chiralA*chiralB<0.0) && (fabs(chiralA)>max_chiral_plane)&& (fabs(chiralB)>max_chiral_plane)){
           if (OutStdOut=='T'){
             printf("#DIFFENRENT_3D_CHIRAL %2s %4d(%4d %4d %4d) Nnei %d chiral %f A --- %2s %4d (%4d %4d %4d) Nnei %d chiral %f A.\n",
              molA->atoms[a].element, molA->atoms[a].num_in_file,
              molA->atoms[match_neiA[0]].num_in_file, molA->atoms[match_neiA[1]].num_in_file, molA->atoms[match_neiA[2]].num_in_file,
              molA->atoms[a].Nnei_heavy, chiralA, 
              molA->atoms[b].element, molB->atoms[b].num_in_file,
              molB->atoms[match_neiB[0]].num_in_file, molB->atoms[match_neiB[1]].num_in_file, molB->atoms[match_neiB[2]].num_in_file, 
              molB->atoms[b].Nnei_heavy, chiralB);
          }
            molA->atoms[a].mark = 1;
            molB->atoms[b].mark = 1;
            for (k=0;k<3;++k){
              molA->atoms[match_neiA[k]].mark = 1;
              molB->atoms[match_neiB[k]].mark = 1;
            } 

            Dchiral += 1;
        } /* chirality conditions */
      } /* Nmatch_nei >=3 */
    } /* Nnei_heavy >= 3 */
  }  /* m */

 if ((OutStdOut=='T')&&(Dchiral==0)){
    printf("#IDENTICAL_3D_CHIRAL\n");
 }

 free(b_on_a);
 return(Dchiral);
} /* end of Check_3D_sp3_Chirality() */





int Modify_MATCH_For_3D_sp3_Chirality_with_Permutation(molA, molB,M)
  struct MOLECULE *molA; /* target    */
  struct MOLECULE *molB; /* reference */
  struct MATCH *M;
{
  int m,Dchiral0,Dchiral_permu,Dchiral_permu_min;
  struct PERMUTATION *pA,*pB;
  char chiral_ok;
  struct MATCH  newM;
  int *indexA,*indexB; /* full atom index */

  printf("#Modify_MCS_For_3D_sp3_Chirality_with_Permutation(Npair %d);\n",M->Npair);
  if (M->Npair==0) return(0);

  Dchiral0 = Check_3D_sp3_Chirality(molA, molB,M, 'T');
  if (Dchiral0==0){ return(0);}

  Malloc_MATCH(&newM,M->Npair);
  Copy_MATCH(&newM,M);
  indexA = (int *)malloc(sizeof(int)*(molA->Natom));
  indexB = (int *)malloc(sizeof(int)*(molB->Natom));

  Dchiral_permu_min = Dchiral0;
  chiral_ok = 0;
  pA = &(molA->permuhead);
  while ((pA->next != NULL) && (chiral_ok==0)){
    pA = pA->next;
    pB = &(molB->permuhead);
    while ((pB->next != NULL) && (chiral_ok==0)){
      pB = pB->next;
      if ((overlap_permuted_index_heavy_and_marked_atoms(molA,pA)>0) || (overlap_permuted_index_heavy_and_marked_atoms(molB,pB)>0)){
         make_permu_index_from_permu_index_heavy(indexA,molA,pA->index_heavy);
         make_permu_index_from_permu_index_heavy(indexB,molB,pB->index_heavy);
         for (m=0;m<M->Npair;++m){
           newM.anumA[m] = indexA[M->anumA[m]];
           newM.anumB[m] = indexB[M->anumB[m]];
         } 
         Dchiral_permu = Check_3D_sp3_Chirality(molA,molB, &newM,'T');
         if (Dchiral_permu<Dchiral_permu_min){
           Dchiral_permu_min = Dchiral_permu;
           printf("#MODIFY ATOM MATCH FOR 3D CHIRALITY !!\n");
           Copy_MATCH(M,&newM);
         }
         if (Dchiral_permu==0){ chiral_ok = 1; }
      }
   } /* pB */
 } /* pA */

 free(indexB);
 free(indexA);
 Free_MATCH(&newM);
 Dchiral_permu_min = Check_3D_sp3_Chirality(molA, molB,M,'T');
 if (Dchiral_permu_min>0){
   printf("#CANNOT MODIFY ATOM MATCH FOR THE IDENTIAL 3D CHIRALITY.\n");
 }
 return(Dchiral_permu_min);

} /* end of Modify_MATCH_For_3D_sp3_Chirality_with_Permutation() */


int overlap_permuted_index_heavy_and_marked_atoms(mol,per)
  struct MOLECULE *mol;
  struct PERMUTATION *per;
{
  int h,a,Npermu_mark;

  Npermu_mark = 0;
  for (h=0;h<per->N;++h){
    if (per->index_heavy[h]!=h){
      a = mol->num_frm_heavy[h];
      if (mol->atoms[a].mark==1){ Npermu_mark +=1; }
    }
  }
  return(Npermu_mark);
} /* end of overlap_permuted_index_heavy_and_marked_atoms() */


void make_permu_index_from_permu_index_heavy(index,mol,index_heavy)
  int *index;           /* [0,1,...,molA->Natom-1]      (for all the atoms) */
  struct MOLECULE *mol;
  int *index_heavy;     /* [0,1,...,molA->Nheavyatom-1] (only for the heavy atoms) */
{
  int a,h,aa;

  for (a=0;a<mol->Natom;++a){index[a] = a;}
  for (h=0;h<mol->Nheavyatom;++h){
      a  = mol->num_frm_heavy[h];
      aa = mol->num_frm_heavy[index_heavy[h]];
      index[a] = aa; 
  }
} /* end of make_permu_index_from_permu_index_heavy() */


int Check_Bond_Angle_Difference(molA, molB,M, OutStdOut)
  struct MOLECULE *molA; /* target    */
  struct MOLECULE *molB; /* reference */
  struct MATCH *M;
{
  int i,j,k;
  int ai,aj,ak,bi,bj,bk;
  float bangA,bangB,diff_bang;
  /*
   ai    ak   bi    bk
     \  /       \  /
      aj         bj
  */

  printf("#Check_Bond_Angle_Difference(molA, molB,M->Npair %d)\n",M->Npair);
  for (i=0;i<M->Npair;++i){
    ai = M->anumA[i];
    bi = M->anumB[i];
    for (j=0;j<M->Npair;++j){
      aj = M->anumA[j];
      bj = M->anumB[j];
      if ((i!=j)&&(molA->conmap.map[ai][aj] != '0') && (molB->conmap.map[bi][bj] != '0')){
        for (k=0;k<M->Npair;++k){
          ak = M->anumA[k];
          bk = M->anumB[k];
          if ((i!=k)&&(j!=k)
              && (molA->conmap.map[aj][ak] != '0') && (molB->conmap.map[bj][bk] != '0') 
              && (molA->conmap.map[ai][ak] == '0') && (molB->conmap.map[bi][bk] == '0')){
            bangA = bond_angle(molA->atoms[ai].Pos, molA->atoms[aj].Pos, molA->atoms[ak].Pos)*180.0/M_PI;
            bangB = bond_angle(molB->atoms[bi].Pos, molB->atoms[bj].Pos, molB->atoms[bk].Pos)*180.0/M_PI;
            diff_bang = fabs(bangA-bangB);
            printf("#BOND_ANGLE %2s %4d %2s %4d %2s %4d -- %2s %4d %2s %4d %2s %4d bangA %6.2f bangB %6.2f diff %6.2f (degree)\n",
               molA->atoms[ai].element, molA->atoms[ai].num_in_file,
               molA->atoms[aj].element, molA->atoms[aj].num_in_file,
               molA->atoms[ak].element, molA->atoms[ak].num_in_file,
               molB->atoms[bi].element, molB->atoms[bi].num_in_file,
               molB->atoms[bj].element, molB->atoms[bj].num_in_file,
               molB->atoms[bk].element, molB->atoms[bk].num_in_file,
               bangA,bangB,diff_bang);
          }
        } /* k */
      }
    } /* j */
  } /* i */

 return(0);
} /* end of Check_Bond_Angle_Difference() */
 

void write_pos6_in_pdb(ofname,pos,R,model_num)
  char  *ofname;
  float pos[6][3];
  int   R;
  int   model_num;
{
  FILE *fpo; 
  int i;

  printf("#write_pos6_in_pdb()-->'%s'\n",ofname);
  fpo = fopen(ofname,"w");
  fprintf(fpo, "REMARK  write_pos_6_in_pdb() -->'%s'\n",ofname);
  fprintf(fpo, "REMARK  Natom_in_ring %d model_num %d\n",R,model_num);
  if (model_num>0) fprintf(fpo,"MODEL    %5d\n",model_num);
  for (i=0;i<R;++i){
    fprintf(fpo,"HETATM%5d %4s%4s %d%5s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",i+1,"C","MOL",model_num,"1", pos[i][0], pos[i][1],pos[i][2],0.0,0.0,"C");
  }
  fprintf(fpo, "ENDMDL\n");
  fclose(fpo);

} /* end of write_pos6_in_pdb() */
