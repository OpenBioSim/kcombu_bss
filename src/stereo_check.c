/*
 
 < stereo_check.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
  functions for checking stereo parities.
 
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
#include "vector3.h"
#include "transform.h"
#include "energies.h"
#include "PCAfit.h"
#include "buildup.h"
#include "moltopodis.h"
#include "ioSDF.h"
#include "ioLINE.h"
#include "ioSMILES.h"
#include "qRMS.h"



/** functions (global) **/
void Set_stereo_parity_of_3D_molecule();
void check_3D_stereo_parity_vs_assigned_parity();
int  check_stereo_parity_for_MATCH();
int neighbors_permutation_number();

/** functions (local) **/
static int check_neighbor_atoms_planarity();
static float chirality_from_four_neighbors();
static float chirality_from_four_neighbors_RMSD();
static void normalize_tetra_POS_ARRAY();




void Set_stereo_parity_of_3D_molecule(mol)
  struct MOLECULE *mol;
{
  int i,j,k,Nsame_rank_pair;
  char chiral_center;
  float chiral_value;

/*
 ### Condition of "chiral_center": ###
    (1) Nneighbor = 3 or 4
    (2) Any neighbor atom pairs have not the same "rank". 
    (3) If Nneighbor==3, these 3 neighbor atoms are not planar.
        It is because a dummy hydrogen atom cannot be generated for the 3-planar atoms. 

  To confirm the "rank" condition (3),  Set_Unique_Extended_Connectivity(mol,-1,0) has to be exececuted.

*/
  printf("#Set_stereo_parity_of_3D_molecule('%s')\n",mol->filename);

  for (i=0;i<mol->Natom;++i){ mol->atoms[i].constr_num = 0;}

  Set_Unique_Extended_Connectivity(mol,-1,0);

  for (i=0;i<mol->Natom;++i){
    chiral_center = 0;
    chiral_value  = 0.0;
    Nsame_rank_pair = 0;
    for (j=0;j<mol->atoms[i].Nneighbor;++j){
      for (k=j+1;k<mol->atoms[i].Nneighbor;++k){
        if (mol->atoms[mol->atoms[i].neighbors[j]].rank == mol->atoms[mol->atoms[i].neighbors[k]].rank) { 
          Nsame_rank_pair += 1;
        } 
      }
    } 

  /* This condition is modified by T.Kawabata (2018/11/01) **/
   if (Nsame_rank_pair==0){
     if ((mol->atoms[i].Nneighbor==4) 
        || ((mol->atoms[i].Nneighbor==3) && (check_neighbor_atoms_planarity(mol,i)==0))){
        chiral_center = 1; 
     }
   }
/*  >> OLD CONDITION<<
   if ((mol->atoms[i].Nneighbor>=3) && (check_neighbor_atoms_planarity(mol,i)==0) && (Nsame_rank_pair==0)){
     chiral_center = 1; 
   }
*/

  if (chiral_center==1){
    /*
    chiral_value =  chirality_from_four_neighbors(mol,i);
    */
    chiral_value =  chirality_from_four_neighbors_RMSD(mol,i);
  }

        if (chiral_value>0.0){ mol->atoms[i].stereo_parity = '1'; }
   else if (chiral_value<0.0){ mol->atoms[i].stereo_parity = '2'; }
   else { mol->atoms[i].stereo_parity = '0'; }

/*
   printf("#%s %d Nneighbor %d Nsame_rank_pair %d neighbor_planar %d chiral_center %d chiral_value %f stereo_parity %c\n",
       mol->atoms[i].atomname,mol->atoms[i].num_in_file,mol->atoms[i].Nneighbor,
      Nsame_rank_pair,
      check_neighbor_atoms_planarity(mol,i),
      chiral_center,chiral_value, mol->atoms[i].stereo_parity);

   for (j=0;j<mol->atoms[i].Nneighbor;++j){
      k = mol->atoms[i].neighbors[j];
      printf("#neighbors_for_%s %d[%d] %s %d (%f %f %f)\n", 
        mol->atoms[i].atomname,mol->atoms[i].num_in_file,
        j,mol->atoms[k].atomname, mol->atoms[k].num_in_file, mol->atoms[k].Pos[0],mol->atoms[k].Pos[1],mol->atoms[k].Pos[2]); 
   }
 */
  /*
   printf("#%-3d %2s rk %2d Nnei %d Nsame_rk_pair %d ch_cen %d ch_val %+5.2f parity %c\n",
     mol->atoms[i].num_in_file,mol->atoms[i].element,
     mol->atoms[i].rank, mol->atoms[i].Nneighbor,Nsame_rank_pair,chiral_center,chiral_value,mol->atoms[i].stereo_parity); 
  */ 
 }

} /* end of Set_stereo_parity_of_3D_molecule() */





void check_3D_stereo_parity_vs_assigned_parity(mol)
  struct MOLECULE *mol;
{
  int i, Nchiral_atoms, Nchiral_atoms_agree, Nchiral_atoms_disagree;
  char stereo_parity_3D;
  float chiral_value;

  printf("#check_3D_stereo_parity_vs_assigned_parity()\n");

  Nchiral_atoms = Nchiral_atoms_agree =  Nchiral_atoms_disagree = 0;

  for (i=0;i<mol->Natom;++i){
    if ((mol->atoms[i].stereo_parity == '1') || (mol->atoms[i].stereo_parity == '2')){
      Nchiral_atoms += 1;
    /*
      chiral_value =  chirality_from_four_neighbors(mol,i);
    */
      chiral_value =  chirality_from_four_neighbors_RMSD(mol,i);

          if (chiral_value>0.0){ stereo_parity_3D = '1'; }
     else if (chiral_value<0.0){ stereo_parity_3D = '2'; }
     else { stereo_parity_3D = '0'; }

     printf("#STEREO_PARITY %d %s stereo_parity %c parity_3D  %c ",
         mol->atoms[i].num_in_file, mol->atoms[i].element, mol->atoms[i].stereo_parity, stereo_parity_3D);  
     if (mol->atoms[i].stereo_parity == stereo_parity_3D){
       printf(" AGREE\n");
       Nchiral_atoms_agree += 1;
     }
     else{
       printf(" DISAGREE\n");
       Nchiral_atoms_disagree += 1;
     }
    } 
 }
 printf("#N_CHIRAL_ATOMS %d\n",Nchiral_atoms);
 printf("#N_CHIRAL_ATOMS_AGREE %d\n",Nchiral_atoms_agree);
 printf("#N_CHIRAL_ATOMS_DISAGREE %d\n",Nchiral_atoms_disagree);

} /* end of check_3D_stereo_parity_vs_assigned_parity(mol) */






int check_stereo_parity_for_MATCH(m,molA,molB,type,header_string)
  struct MATCH *m;
  struct MOLECULE *molA,*molB;
  char   type;   /* 'C':just checking, 'D':delete unmatch stereo atom pairs */
  char   *header_string;  /* For example, "INIT_" or "FINAL_".   */
/*
    (molA->atoms[a].stereo_parity ) and (molB->atoms[b].stereo_parity ) 
    have to be assigned before executing this function.

  OUTPUT:
    molA->atoms[a].mark ==  0: otherwise 
    molA->atoms[a].mark ==  1: chirarity-agree
    molA->atoms[a].mark ==  2: chirarity-agree

*/

{
  int i,j,k,a,b, missA, missB, agree, disagree;
  int Natom_chiralA,Natom_chiralB, Natom_chiral_agree, Natom_chiral_disagree, Natom_chiral_missA, Natom_chiral_missB;
  int *b_on_a,Npermu, *anumA0,*anumB0;
  int hvyneiA[10],NhvyneiA,hvyneiB[10],NhvyneiB;  /* heavy atoms neighbors checking for matching permutations */
  char *delmark;

  printf("#check_stereo_parity_for_MATCH(type %c header_string '%s')\n",type,header_string);
  b_on_a = (int *)malloc(sizeof(int)*molA->Natom);
  for (a=0;a<molA->Natom;++a){ 
    b_on_a[a] = -1;
    molA->atoms[a].mark = 0;
  }
  for (i=0;i<m->Npair;++i){ b_on_a[m->anumA[i]] = m->anumB[i]; }


  anumA0 = anumB0  = NULL;
  delmark = NULL;

  if (type=='D'){
    anumA0  = (int *)malloc(sizeof(int)*m->Npair);
    anumB0  = (int *)malloc(sizeof(int)*m->Npair);
    delmark = (char *)malloc(sizeof(char)*m->Npair);
    for (i=0;i<m->Npair;++i){ 
      anumA0[i] = m->anumA[i];
      anumB0[i] = m->anumB[i];
      delmark[i] = 0;
    }
  }


  Natom_chiralA = Natom_chiralB =  Natom_chiral_agree = Natom_chiral_disagree = 0;
  Natom_chiral_missA = Natom_chiral_missB = 0;

  for (i=0;i<m->Npair;++i){
    a = m->anumA[i];
    b = m->anumB[i];
    /** count NhvyneiA and NhvyneiB **/
    NhvyneiA = NhvyneiB = 0;
    for (j=0;j<molA->atoms[a].Nneighbor;++j){
      k = molA->atoms[a].neighbors[j];
      if (molA->atoms[k].one_char_ele != 'H'){
        hvyneiA[NhvyneiA] = k;
        hvyneiB[NhvyneiB] = b_on_a[k];
        if (hvyneiB[NhvyneiB]>=0){ NhvyneiB += 1;}
        NhvyneiA +=1;
      }
   }
   disagree = agree = missA = missB =0;
   Npermu = 0;
 /* 
  printf("#molA %s %d stereo %c molB %s %d stereo %c\n",
     molA->atoms[a].atomname,molA->atoms[a].num_in_file,molA->atoms[a].stereo_parity,
     molB->atoms[b].atomname,molB->atoms[b].num_in_file,molB->atoms[b].stereo_parity);
 */ 
    if ((NhvyneiA == NhvyneiB) && ((molA->atoms[a].Nneighbor==4)||(molB->atoms[b].Nneighbor==4)) 
        && (  (molA->atoms[a].stereo_parity == '1')||(molA->atoms[a].stereo_parity=='2') 
            ||(molB->atoms[b].stereo_parity == '1')||(molB->atoms[b].stereo_parity=='2'))) { 

      if ((molA->atoms[a].stereo_parity == '1')||(molA->atoms[a].stereo_parity=='2')){ 
        Natom_chiralA += 1; 
        if ((molB->atoms[b].stereo_parity != '1')&&(molB->atoms[b].stereo_parity!='2')){ 
          missB = 1; 
          Natom_chiral_missB += 1;
        }
      }
 
      if ((molB->atoms[b].stereo_parity == '1')||(molB->atoms[b].stereo_parity=='2')){ 
        Natom_chiralB += 1; 
        if ((molA->atoms[a].stereo_parity != '1')&&(molA->atoms[a].stereo_parity!='2')){ 
          missA = 1; 
          Natom_chiral_missA += 1;
        }
      }

 
      if ( ((molA->atoms[a].stereo_parity == '1')||(molA->atoms[a].stereo_parity=='2')) &&
           ((molB->atoms[b].stereo_parity == '1')||(molB->atoms[b].stereo_parity=='2')) ){

        Npermu = neighbors_permutation_number(NhvyneiA, hvyneiA, hvyneiB);
        if ((((Npermu%2)==0)&&(molA->atoms[a].stereo_parity != molB->atoms[b].stereo_parity))||
            (((Npermu%2)==1)&&(molA->atoms[a].stereo_parity == molB->atoms[b].stereo_parity))){
          disagree = 1;
        } 
        else{
          agree = 1;
        }
        molA->atoms[a].mark = 1;
        if (disagree==1){
          molA->atoms[a].mark = 2;
          Natom_chiral_disagree += 1;
        } 
        else{
          Natom_chiral_agree += 1;
        }
 
        if ((type=='D') && (disagree==1)){
           printf("#delmark [%d] m->Npair %d\n",i,m->Npair); fflush(stdout);
           delmark[i] = 1;
        }
     }

    /*** OUTPUT ***/
        printf("#STEREO_CHECK (%s%-2d %s%-2d)",
          molA->atoms[a].element,molA->atoms[a].num_in_file,
          molB->atoms[b].element,molB->atoms[b].num_in_file);

        printf(" parity (%c %c) Npermu %d",
           molA->atoms[a].stereo_parity, molB->atoms[b].stereo_parity,Npermu);
        if (disagree==1){ 
          printf(" DISAGREE");
        } 
        if (agree==1){ 
          printf(" AGREE   ");
        }
        if (missA==1){ 
          printf(" MISS_A   ");
        }
        if (missB==1){ 
          printf(" MISS_B   ");
        }

        printf(" (");
        for (j=0;j<NhvyneiA;++j){ printf(" %s%-2d",
          molA->atoms[hvyneiA[j]].element,molA->atoms[hvyneiA[j]].num_in_file);
         }
        printf(")(");
        for (j=0;j<NhvyneiB;++j){ printf(" %s%-2d",
          molB->atoms[hvyneiB[j]].element,molB->atoms[hvyneiB[j]].num_in_file);
        }
        printf(")\n");

    }
  }

 free(b_on_a);
 printf("#%sN_CHIRAL_ATOMS_A %d\n",header_string,Natom_chiralA);
 printf("#%sN_CHIRAL_ATOMS_B %d\n",header_string,Natom_chiralB);
 printf("#%sN_CHIRAL_ATOMS_MISS_A %d\n",header_string,Natom_chiral_missA);
 printf("#%sN_CHIRAL_ATOMS_MISS_B %d\n",header_string,Natom_chiral_missB);
 printf("#%sN_CHIRAL_ATOMS_AGREE    %d\n",header_string,Natom_chiral_agree);
 printf("#%sN_CHIRAL_ATOMS_DISAGREE %d\n",header_string,Natom_chiral_disagree);

 /* type=='D':delete unmatch stereo atom pairs */
 if (type=='D'){
   j = 0;
   for (i=0;i<m->Npair;++i){
     if (delmark[i]==0){
       m->anumA[j] = anumA0[i];
       m->anumB[j] = anumB0[i];
       j += 1;
     } 
   }
   m->Npair = j;

 }
 free(delmark);
 free(anumB0);
 free(anumA0);
 return(Natom_chiral_disagree);

} /* end of check_stereo_parity_for_MATCH() */






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

  Nvec = oa X ob 

  if (angle bwn Nvec and oc) = 90 degree, (a,b,c) are on the plane.
  if (angle bwn Nvec and od) = 90 degree, (a,b,d) are on the plane.

 */
  angle_tole = 8.0; /** degree **/
  a = b = c = d = -1;
  plane = 1;

/*
  printf("#%s %d\n",mol->atoms[o].atomname,mol->atoms[o].num_in_file);
*/
  for (i=0;i<mol->atoms[o].Nneighbor;++i){
    if (mol->atoms[mol->atoms[o].neighbors[i]].one_char_ele != 'H'){
           if (a<0){ a = mol->atoms[o].neighbors[i];}
      else if (b<0){ b = mol->atoms[o].neighbors[i];}
      else if (c<0){ c = mol->atoms[o].neighbors[i];}
      else if (d<0){ d = mol->atoms[o].neighbors[i];}
    }
  }

  if ((a>=0) && (b>=0)){
    sub_vec3(OA,mol->atoms[a].Pos, mol->atoms[o].Pos); normalize_vec3(OA);
    sub_vec3(OB,mol->atoms[b].Pos, mol->atoms[o].Pos); normalize_vec3(OB);
    cross_vec3(Nvec,OA,OB); normalize_vec3(Nvec);
  }


  if (c>=0){
    sub_vec3(OC, mol->atoms[c].Pos, mol->atoms[o].Pos);
    angle =  angle_bwn_vec(Nvec,OC)*180.0/M_PI;
    if (fabs(angle-90.0) > angle_tole){ plane = 0;}
    /* printf("#angle_c %f\n",angle); */
  }

  else if (d>=0){
    sub_vec3(OD, mol->atoms[d].Pos, mol->atoms[o].Pos);
    angle =  angle_bwn_vec(Nvec,OD)*180.0/M_PI;
    if (fabs(angle-90.0) > angle_tole){ plane = 0;}
    /* printf("#angle_d %f\n",angle); */
  }

  return(plane);
} /* end of check_neighbor_atoms_planarity() */





float chirality_from_four_neighbors(mol,o)
  struct MOLECULE *mol;          /* molecule    (given) */
  int    o;  /* atom number of the focused atom */
{
/*
  o: focused atom.
  0, 1, 2, 3 : neighbor atoms (order of atomic number)
      
[position] 
       0      
       |
    3--o--1
       |
       2 

   n = 01 cross 12
   n = n /||n||
  [chirality]  = n dot (3-o) 

 The chiratity can be regarded as a distance bwn the atom 4 and the plane with normal vec 12 x 23.
 
 If atom 3 above the plane 01 x 12, the chirality >0, 012 is clockwise.        (stereo_parity = 1)
 If atom 3 below the plane 01 x 12, the chirality <0, 012 is counter-clockwise (stereo_parity = 2).
 If the signs of chirality are different, chirarities are different.
 
                                    
 ** If the atom 3 is not available, (This is modified in 2018/10/24) 
     (0-o)+(1-o)+(2-o)+(3-o) = 0+1+2+3-4o = zero;
      Then 3 = 4o - (0+1+2)
 */

  float ev01[3],ev12[3],nvec[3],evo3[3],pos3[3],chirality;
  int i,j,k;
  int Nneighbor,neighbors[10];

  Nneighbor = 0;
  for (i=0;i<mol->atoms[o].Nneighbor;++i){
    j = mol->atoms[o].neighbors[i];
    neighbors[Nneighbor] = j;
    Nneighbor += 1;

/* This part was modified in 2018/11/02 by T.Kawabata 
   Hydrogen atoms should be always included to determine chilarity,
   not regarding "includeHinMCS" flag. 
 */
 
/*  [Old program] 
 
    if ((PAR.includeHinMCS=='T') || (mol->atoms[j].one_char_ele != 'H')){
      neighbors[Nneighbor] = j;
      Nneighbor += 1;
    }
*/

  }
 
  if (Nneighbor<3){ return(0.0);}

  sub_vec3(ev01,mol->atoms[neighbors[1]].Pos, mol->atoms[neighbors[0]].Pos); normalize_vec3(ev01);
  sub_vec3(ev12,mol->atoms[neighbors[2]].Pos, mol->atoms[neighbors[1]].Pos); normalize_vec3(ev12);
  cross_vec3(nvec,ev01,ev12); normalize_vec3(nvec);
  if (Nneighbor>=4){
     sub_vec3(evo3,mol->atoms[neighbors[3]].Pos, mol->atoms[o].Pos);
  }
  else{
    for (k=0;k<3;++k) {
/*
      pos3[k] =-(mol->atoms[neighbors[0]].Pos[k] + mol->atoms[neighbors[1]].Pos[k] + mol->atoms[neighbors[2]].Pos[k])
               + 2.0*mol->atoms[o].Pos[k];
    This equation is modified as follows (2018/10/24)
*/
      pos3[k] = 4.0*mol->atoms[o].Pos[k] -(mol->atoms[neighbors[0]].Pos[k] + mol->atoms[neighbors[1]].Pos[k] + mol->atoms[neighbors[2]].Pos[k]);
    }
    sub_vec3(evo3, pos3, mol->atoms[o].Pos);
  }
  normalize_vec3(evo3);
  chirality = dot_prod3(nvec,evo3);

 /*
  printf("#%s %2d nvec %+6.2f %+6.2f %+6.2f evo3 %+6.2f %+6.2f %+6.2f chiral %f",
      mol->atoms[o].element,mol->atoms[o].num_in_file,nvec[0],nvec[1],nvec[2],evo3[0],evo3[1],evo3[2],chirality);
  printf(" Nneighbor %d",mol->atoms[o].Nneighbor);
  printf("\n");
 */

  return(chirality);

} /* end of chirality_from_four_neighbors() */




float chirality_from_four_neighbors_RMSD(mol,o)
  struct MOLECULE *mol;          /* molecule    (given) */
  int    o;  /* atom number of the focused atom */
{
/*
  o: focused atom.
  0, 1, 2, 3 : neighbor atoms (order of atomic number)
      
[position] 
       2        3      
       |        |
    3--o--0  4--0--1
       |        |
       1        2 

  Prepare clock-wise Tetra (clckT) and anti-clock-wise Tetera (antiT) 
      2                 2   
       \                 \
        \                 \
3-down- 0-----1   4-down- 0-----1
       /                 /
      up                up
     /                 /
    4                 3 
  >> clckT <<      >> antiT <<


      3                3      
      |                |
      |                |
 4 -- 0-down-1    4 -- 0-down-1
      \                \ 
       up               up
         \               \
          2              2
  >> clckT <<      >> antiT <<


  chiralitly = RMSD(targT, antiT) - RMSD(targT, clckT)

  return(chirality)
  
  if (chilarity > 0.0):
     RMSD(targT, antiT) > RMSD(targT, clckT)   ===> targT is more similar to clckT.
  else:
     RMSD(targT, antiT) < RMSD(targT, clckT)   ===> targT is more similar to antiT.


  This algorithm is more robust than the sign-volume algorithm, for skewed tetrahedral atoms.

 */
  int i,j,k;
  int Nneighbor,neighbors[10];
  struct POS_ARRAY clckT, antiT; /* ideal tetrahedral  5 atoms. clckT:clock-wise, antiT:anti-clock-wise */
  struct POS_ARRAY targT;        /* target tetrahedral 5 atoms */
  float clckRMSD, antiRMSD,chirality;
  double sqrt2,sqrt3;

  /** [1] Setup clckT and antiT ***/ 
  Malloc_POS_ARRAY(&clckT,5);
  Malloc_POS_ARRAY(&antiT,5);
  clckT.N = 5;
  sqrt2 = sqrt(2.0);
  sqrt3 = sqrt(3.0);
  clckT.pos[0][0] = 0.0;      clckT.pos[0][1] = 0.0;       clckT.pos[0][2] = 0.0; 
  clckT.pos[1][0] = 1.0;      clckT.pos[1][1] = 0.0;       clckT.pos[1][2] = 0.0; 
  clckT.pos[2][0] = -1.0/3.0; clckT.pos[2][1] = 2*sqrt2/3; clckT.pos[2][2] = 0.0; 
  clckT.pos[3][0] = -1.0/3.0; clckT.pos[3][1] =  -sqrt2/3; clckT.pos[3][2] = -sqrt2/sqrt3; 
  clckT.pos[4][0] = -1.0/3.0; clckT.pos[4][1] =  -sqrt2/3; clckT.pos[4][2] =  sqrt2/sqrt3; 

  antiT.N = 5;
  for (i=0;i<5;++i){
    for (k=0;k<3;++k){
      antiT.pos[i][k] = clckT.pos[i][k];
    }
  } 
  antiT.pos[3][2] = clckT.pos[4][2];
  antiT.pos[4][2] = clckT.pos[3][2];

  /*** [2] Setup targT ***/
  Malloc_POS_ARRAY(&targT,5);
  targT.N = 5;
  for (k=0;k<3;++k){
    targT.pos[0][k] = mol->atoms[o].Pos[k]; 
  }

  Nneighbor = 0;
  for (i=0;i<mol->atoms[o].Nneighbor;++i){
    j = mol->atoms[o].neighbors[i];
    neighbors[Nneighbor] = j;
    Nneighbor += 1;
    for (k=0;k<3;++k){
      targT.pos[Nneighbor][k] = mol->atoms[j].Pos[k];
    }
  }
 
  if (Nneighbor<3){ return(0.0);}

  if (Nneighbor==3){
    for (k=0;k<3;++k) {
      targT.pos[4][k] = 4.0*mol->atoms[o].Pos[k] 
                          - (mol->atoms[neighbors[0]].Pos[k] + mol->atoms[neighbors[1]].Pos[k] + mol->atoms[neighbors[2]].Pos[k]);
    }
  }

  normalize_tetra_POS_ARRAY(&targT);

  clckRMSD = Calculate_CRMS_bwn_POS_ARRAYs(&targT, &clckT);
  antiRMSD = Calculate_CRMS_bwn_POS_ARRAYs(&targT, &antiT);
  chirality = antiRMSD - clckRMSD; 
 
/* 
  printf("#%s %2d clckRMSD %lf antiRMSD %lf chiral %f\n",
      mol->atoms[o].element,mol->atoms[o].num_in_file,clckRMSD, antiRMSD,chirality);
 */
 /*
  printf("#%s %2d nvec %+6.2f %+6.2f %+6.2f evo3 %+6.2f %+6.2f %+6.2f chiral %f",
      mol->atoms[o].element,mol->atoms[o].num_in_file,nvec[0],nvec[1],nvec[2],evo3[0],evo3[1],evo3[2],chirality);
  printf(" Nneighbor %d",mol->atoms[o].Nneighbor);
  printf("\n");
 */

  Free_POS_ARRAY(&targT);
  Free_POS_ARRAY(&antiT);
  Free_POS_ARRAY(&clckT);

  return(chirality);

} /* end of chirality_from_four_neighbors_RMSD() */


void normalize_tetra_POS_ARRAY(T)
  struct POS_ARRAY *T;
{
  int i,k;
  double norm;

  if (T->N == 5){
    for (i=1;i<5;++i){
      norm = 0.0;
      for (k=0;k<3;++k){
        T->pos[i][k] =  T->pos[i][k] - T->pos[0][k];
        norm += T->pos[i][k] * T->pos[i][k];
      }
      norm = sqrt(norm);
      for (k=0;k<3;++k){
        T->pos[i][k] =  T->pos[i][k]/norm;
      }
    }

    T->pos[0][0] =  T->pos[0][1] = T->pos[0][2] = 0.0;

  }
  else{
   printf("#normalize_tetra_POS_ARRAY(T):T->N %d is not 5.\n",T->N);
  }
} /* end of normalize_tetra_POS_ARRAY() */






int neighbors_permutation_number(N,A,B)
  int N;
  int *A; /* array of neighboring atoms for molecule A */
  int *B; /* array of neighboring atoms for molecule B */
/*

 return (ntrans)
 
 ntrans: number of transpose to get equivalent order B[] for A[].

ex1) molA       molB 
      1          2
      |          |
      a..4       b..5       A = [1,2,3,4], B = [2,3,4,5] => ntrans = 0
    //  \      //  \ 
   3     2    4     3 
   clock       clock     ==> agree


ex2) molA       molB 
      1          2
      |          |
      a..4       b..5       A = [1,2,3,4], B = [2,4,3,5] => ntrans = 1
    //  \      //  \ 
   3     2    3     4 
   clock      counter-clk     ==> agree


ex3) molA       molB 
      1          2
      |          |
      a..4       b..3       A = [1,2,3,4], B = [2,4,5,3] => ntrans = 2
    //  \      //  \ 
   3     2    5     4 
   clock      clock   ==> agree

ex4) molA       molB 
      1          2
      |          |
      a..4       b==5       A = [1,2,3,4], B = [2,3,4,5] => ntrans = 0
    //  \      .  \ 
   3     2    4     3 
    clock     counter-clock   ==> disagree


ex5) molA       molB 
      1          2
      |          |
      a..4       b==3       A = [1,2,3,4], B = [2,4,5,3] => ntrans = 2
    //  \      .  \ 
   3     2    5    4 
    clock      counter-clk   ==> disagree


*/
{
  int *B0;
  int i,Ntranspose,ntrans,buff;

  B0 = (int *)malloc(sizeof(int)*N);
  for (i=0;i<N;++i){
    B0[i] = B[i];
  }
  Ntranspose = 0;
  do{
    ntrans = 0;
    for (i=1;i<N;++i){
      /* if orders of A[] and B[] are different, exchange B[i-1] and B[i] */
      if (((A[i-1]<A[i]) && (B0[i-1]>B0[i])) ||
          ((A[i-1]>A[i]) && (B0[i-1]<B0[i])) ){
        buff = B0[i-1]; 
        B0[i-1] = B0[i];
        B0[i] = buff; 
        ntrans += 1;
        Ntranspose += 1;
      }
    } 
  } while (ntrans>0);

  free(B0);
  return(Ntranspose);

} /* end of neighbors_permutation_number() */

