/*

< stereo_adjust.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 functions for adjusting stereo parities.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
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
#include "ioLINE.h"
#include "stereo_check.h"

/*** FUNCTIONS (GLOBAL) ***/
int Adjust_3D_Stereo_for_2D_Stereo_Reference();
int Change_Flip_Of_Fragment();

/*** FUNCTIONS (LOCAL) ***/
static void Reflection_of_All_Atoms();
static int Exchange_sp3_Chirality_for_Stereo_Parity();
static int Pyramid_Inversion_sp3_Chirality_for_Stereo_Parity();
static int Fold_Six_Atom_Ring_for_sp3_Chirality_and_Stereo_Parity();
static void chose_rotatable_bonded_atoms_among_4neighbors();
static void set_reference_matched_neighbor_atoms();
static void exchange_positions_of_two_atom_groups();
static int  inversion_of_pyramid();
static void chose_pybase_from_match_nei();
static void normal_vector_of_ring();
static void make_r_on_t_from_MATCH();
static int check_in_the_same_ring_bwn_two_atom_groups();
static int check_in_ring();



int Adjust_3D_Stereo_for_2D_Stereo_Reference(molT, molR, m_model, PyramidInversion,AllReflectChiral,FoldHexaRing)
  struct MOLECULE *molT;
  struct MOLECULE *molR;
  struct MATCH *m_model;
  char  PyramidInversion;
  char  AllReflectChiral;
  char  FoldHexaRing;
{
   int  round,Nmodified,Ndisagree,Ndisagree_orig, Ndisagree_fin, Nround_max;
   char Transform[20];
   int  Nexchange_pyramid,Nreflect, Nfoldsix; 

  printf("#Adjust_3D_Stereo_for_2D_Stereo_Reference(PyramidInversion %c AllReflectChiral %c FoldHexaRing %c\n",
          PyramidInversion, AllReflectChiral,FoldHexaRing);

   Nround_max = 12;
   Transform[ 0] = 'E';  /* Exchange */
   Transform[ 1] = 'p';  /* pyramid  */
   Transform[ 2] = 'P';  /* Pyramid enforce */
   Transform[ 3] = 'R';  /* Reflection */
   Transform[ 4] = 'E';  /* Exchange */
   Transform[ 5] = 'p';  /* pyramid */
   Transform[ 6] = 'P';  /* Pyramid enforce */
   Transform[ 7] = 'S';  /* FoldSix */
   Transform[ 9] = 'E';  /* Exchange */
   Transform[10] = 'p';  /* pyramid */
   Transform[11] = 'P';  /* Pyramid enforce */


   Set_stereo_parity_of_3D_molecule(molT);
   Ndisagree = check_stereo_parity_for_MATCH(m_model, molT, molR,'C',"INIT_");
   printf("#INIT_Eselfclash %f\n",energy_selfclash(molT));
   printf("#INIT_Nclash_atom_pair %d\n",count_selfclash(molT));
   
   if (Ndisagree==0){
     printf("#NO ATOM REQUIRED FOR STEREO PARITY CORRECTION.\n");
     return(0);
   }

   Nexchange_pyramid = Ndisagree_orig =   Nreflect = Nfoldsix = 0;

   round = 0;
   while ((round<Nround_max)&&(Ndisagree>0)){
     printf("#>ROUND[%d] transform %c\n",round,Transform[round]);
     /* Copy_MOLECULE_XYZ(molT,&molTo); */

     if ((AllReflectChiral=='T')&&(Transform[round]=='R')){
        Reflection_of_All_Atoms(molT);
        Nreflect += 1;
     }
     else if ((FoldHexaRing=='T') && (Transform[round]=='S')){
       Nmodified = Fold_Six_Atom_Ring_for_sp3_Chirality_and_Stereo_Parity(molT, molR,m_model);
       if (Nmodified>0){
         Nfoldsix += 1;
       }
     }
     else if ((PyramidInversion=='T') && (Transform[round]=='p')){
       Nmodified = Pyramid_Inversion_sp3_Chirality_for_Stereo_Parity(molT, molR, m_model,'C');
     }
     else if ((PyramidInversion=='T') && (Transform[round]=='P')){
       Nmodified = Pyramid_Inversion_sp3_Chirality_for_Stereo_Parity(molT, molR, m_model,'-');
     }
     else if (Transform[round]=='E'){
       Nmodified = Exchange_sp3_Chirality_for_Stereo_Parity(molT, molR, m_model,'-');
     }

     Set_stereo_parity_of_3D_molecule(molT);
     Ndisagree = check_stereo_parity_for_MATCH(m_model, molT, molR,'C',"");
     round += 1;
   }
  Set_stereo_parity_of_3D_molecule(molT);
  Ndisagree_fin = check_stereo_parity_for_MATCH(m_model, molT, molR,'C',"FINAL_");
  printf("#N_ROUND   %d\n",round);
  printf("#N_REFLECT %d\n",Nreflect);
  printf("#N_FOLDSIX %d\n",Nfoldsix);
  printf("#FINAL_Eselfclash %f\n",energy_selfclash(molT));
  printf("#FINAL_Nclash_atom_pair %d\n",count_selfclash(molT));
  return(Ndisagree_orig -Ndisagree_fin);

} /* end of Adjust_3D_Stereo_for_2D_Stereo_Reference() */











int Exchange_sp3_Chirality_for_Stereo_Parity(molT, molR, M)
  struct MOLECULE *molT;    /* target    */
  struct MOLECULE *molR;    /* reference */
  struct MATCH *M;
{
  int a,b,m;
  int *r_on_t;
  int  Ndisagree, Ndisagree_exchange,Ndisagree_rest;
  int  Nmatch_nei, Nmodified, Nmodified_exchange;
  int  match_neiT[4],match_neiR[4];  /* atom number of four neighboring atoms.  If no neighobor match_neiT[n] = -1. */
  int  ch_i,ch_j;
  
  /*
 
  Before execute this function
    check_stereo_parity_for_MATCH(m,molA,molB,'C',"") should be perfomed 
   to obrain molT->atoms[].mark. 
 
    molT->atoms[a].mark ==  1: chirarity-agree
    molT->atoms[a].mark ==  2: chirarity-agree
    molT->atoms[a].mark ==  0: otherwise 
  */
  printf("#Exchange_sp3_Chirality_for_Stereo_Parity()\n");

  if (M->Npair==0){ return(0);}

  /** [1] prepare r_on_t[anum_in_molT] := anum_in_molR for matched atom pairs **/
  r_on_t = (int *)malloc(sizeof(int)*molT->Natom);
  make_r_on_t_from_MATCH(r_on_t, M, molT);


  /** [2] scan the match ***/
  Ndisagree = Ndisagree_exchange = Ndisagree_rest = 0; 
  Nmodified = Nmodified_exchange = 0;

  for (m=0;m<M->Npair;++m){

    a = M->anumA[m];
    b = M->anumB[m];

    if ((molT->atoms[a].Nneighbor>=3)&&(molR->atoms[b].Nneighbor>=3)&&
        ((molR->atoms[b].stereo_parity =='1')||(molR->atoms[b].stereo_parity=='2'))&&
         (molT->atoms[a].mark == 2)
        ){

       Ndisagree += 1; 

       set_reference_matched_neighbor_atoms(molT, a, r_on_t, &Nmatch_nei,match_neiT,match_neiR);
       chose_rotatable_bonded_atoms_among_4neighbors(molT,a, Nmatch_nei, match_neiT,&ch_i,&ch_j);

         /** (1) Exchange positions ch_i and ch_j **/
       if ((ch_i>=0)&&(ch_j>=0)){
          Ndisagree_exchange += 1;
          exchange_positions_of_two_atom_groups(molT,a,Nmatch_nei,match_neiT,ch_i, ch_j);
          Nmodified += 1;
          Nmodified_exchange += 1;
       } 
      /*
       show_chiralities_molT_molR(molT, molR,a,b,match_neiT,match_neiR,0.0,'-');
      */
    } 

  } 

 free(r_on_t);
 return(Nmodified);
} /* end of Exchange_sp3_Chirality_for_Stereo_Parity() */







int Pyramid_Inversion_sp3_Chirality_for_Stereo_Parity(molT, molR, M, CheckType)
  struct MOLECULE *molT;    /* target    */
  struct MOLECULE *molR;    /* reference */
  struct MATCH *M;
  char CheckType;   /* if 'C', checking angle bwn 0-1-2 and o-included rings is < norm_angle_max */
{
  int a,b,m;
  int *r_on_t;
  int  Ndisagree, Ndisagree_pyramid,Ndisagree_rest;
  int  Nmatch_nei, Nmodified, Nmodified_pyramid, Ninversion;
  int  match_neiT[4],match_neiR[4];  /* atom number of four neighboring atoms.  If no neighobor match_neiT[n] = -1. */
  int  ch_i,ch_j;
  int pybaseT[3];  /* atom numbers of three pyramid-base atoms */
  
  /*
 
  Before execute this function
    check_stereo_parity_for_MATCH(m,molA,molB,'C',"") should be perfomed 
   to obrain molT->atoms[].mark. 
 
    molT->atoms[a].mark ==  1: chirarity-agree
    molT->atoms[a].mark ==  2: chirarity-agree
    molT->atoms[a].mark ==  0: otherwise 
  */
  printf("#Pyramid_Inversion_sp3_Chirality_for_Stereo_Parity()\n");

  if (M->Npair==0){ return(0);}
  /** [1] prepare r_on_t[anum_in_molT] := anum_in_molR for matched atom pairs **/
  r_on_t = (int *)malloc(sizeof(int)*molT->Natom);
  make_r_on_t_from_MATCH(r_on_t, M, molT);


  /** [2] scan the match ***/
  Ndisagree = Ndisagree_pyramid = Ndisagree_rest = 0; 
  Nmodified = Nmodified_pyramid = 0;

  for (m=0;m<M->Npair;++m){
    a = M->anumA[m];
    b = M->anumB[m];

    if ((molT->atoms[a].Nneighbor>=3)&&(molR->atoms[b].Nneighbor>=3)&&
        ((molR->atoms[b].stereo_parity =='1')||(molR->atoms[b].stereo_parity=='2'))&&
         (molT->atoms[a].mark == 2)
        ){

       Ndisagree += 1; 

       set_reference_matched_neighbor_atoms(molT, a, r_on_t, &Nmatch_nei,match_neiT,match_neiR);
       chose_rotatable_bonded_atoms_among_4neighbors(molT,a, Nmatch_nei, match_neiT,&ch_i,&ch_j);

       /** (1) Inversion of Pyramids **/
       if ((ch_i<0)||(ch_j<0)){
         chose_pybase_from_match_nei(pybaseT,ch_i,ch_j,Nmatch_nei,match_neiT);
         Ninversion = inversion_of_pyramid(molT,a,pybaseT,CheckType);
         if (Ninversion==1){
           Ndisagree_pyramid +=  1;
           Nmodified += 1;
           Nmodified_pyramid += 1;
         }
         else { 
           Ndisagree_rest += 1;
         }
       }
    }
  }

 free(r_on_t);
 return(Nmodified);
} /* end of Pyramid_Inversion_sp3_Chirality_for_Stereo_Parity() */











int subnum_one_for_exchange_possible_atoms(molT, molR, M)
  struct MOLECULE *molT;    /* target    */
  struct MOLECULE *molR;    /* reference */
  struct MATCH *M;
{
  int a,b,m;
  int *r_on_t;
  int  Nmatch_nei, ch_i,ch_j;
  int  match_neiT[4],match_neiR[4];  /* atom number of four neighboring atoms.  If no neighobor match_neiT[n] = -1. */
  /*
 
  Before executing this function
    check_stereo_parity_for_MATCH(m,molA,molB,'C',"") should be perfomed 
   to obrain molT->atoms[].mark. 
 
    molT->atoms[a].mark ==  1: chirarity-agree
    molT->atoms[a].mark ==  2: chirarity-disagree
    molT->atoms[a].mark ==  0: otherwise 

  Output  
    molT->atoms[a].subnum ==  1:  chirarity-disagree and it will be modified by exchange.
 
  */
  if (M->Npair==0){ return(0);}
  for (a=0;a<molT->Natom;++a){
    molT->atoms[a].subnum = 0;
  } 
  /** [1] prepare r_on_t[anum_in_molT] := anum_in_molR for matched atom pairs **/
  r_on_t = (int *)malloc(sizeof(int)*molT->Natom);
  make_r_on_t_from_MATCH(r_on_t, M, molT);
  
  /** [2] scan atom pairs ***/
  for (m=0;m<M->Npair;++m){

    a = M->anumA[m];
    b = M->anumB[m];

    if ((molT->atoms[a].Nneighbor>=3)&&(molR->atoms[b].Nneighbor>=3)&&
        ((molR->atoms[b].stereo_parity =='1')||(molR->atoms[b].stereo_parity=='2'))&&
         (molT->atoms[a].mark == 2)
        ){

       set_reference_matched_neighbor_atoms(molT, a, r_on_t, &Nmatch_nei,match_neiT,match_neiR);
       chose_rotatable_bonded_atoms_among_4neighbors(molT,a, Nmatch_nei, match_neiT,&ch_i,&ch_j);

         /** (1) Exchange positions ch_i and ch_j **/
       if ((ch_i>=0)&&(ch_j>=0)){
          molT->atoms[a].subnum = 1;
       } 
    } 
  }  
  return(1);
} /* end of subnum_one_for_exchange_possible_atoms() */










int Fold_Six_Atom_Ring_for_sp3_Chirality_and_Stereo_Parity(molT, molR, M)
  struct MOLECULE *molT;    /* target    */
  struct MOLECULE *molR;    /* reference */
  struct MATCH *M;
{
  int s,t,tt,r,rr;
  int *r_on_t;
  struct RING *rn;
  int NmatchB,Ndisagree,Nmodify;
  int a,b,c,d,e,f,g,h,i,j; /* atom number for molT */
  float AD[3],AG[3],BE[3],BH[3],ADxAG[3],angleAD_BE,angleAG_BH;
  float tole_angle_parallel; /* degree */
  float rotaxisBA[3],rangle,rangleDAG,rangleEBH,Rmat[3][3],RImat[3][3];
  int   minDabc,minDdef,Da,Db,Dc,Dd,De,Df,Dg,Dh;
  char buffstr[32],out_string[256];
  int numA[3],numB[3],insamering;
 
  printf("#Fold_Six_Atom_Ring_for_sp3_Chirality_and_Stereo_Parity()\n");
  /*
 
  This process is mainly for modifying "morphine" (KEGG_DRUG: D08233) structure.

  Before executing this function, two functions, 
    check_stereo_parity_for_MATCH(m,molA,molB,'C',"")
    subnum_one_for_exchange_possible_atoms(molT, molR, M)
 should be perfomed to obrain molT->atoms[].mark. and molT->atoms[].subnum 
 
    molT->atoms[a].mark ==  1: chirarity-agree
    molT->atoms[a].mark ==  2: chirarity-disagree
    molT->atoms[a].mark ==  0: otherwise 
    molT->atoms[a].subnum ==  1:  chirarity-disagree and it will be modified by exchange.


 << rotation scheme >> 

  *Chirarities of atoms A and B should be modified.
      I                    I
      |                    |
    --A----D             - A--G
   /  |     \          /   |     
  /   G      \        /    |      
 C           F  ==> C      D      
  \         /        \     |        
   \       /          \    | 
 J- B-----E          J-B-H | 
    |                  |   F
    H                  |  /  
                       E-/
  */

  tole_angle_parallel = 30.0*M_PI/180.0;

  r_on_t = (int *)malloc(sizeof(int)*molT->Natom);
  make_r_on_t_from_MATCH(r_on_t, M, molT);
  Nmodify = 0;

  rn = &(molT->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;

    if (rn->Natom==6){
      NmatchB = Ndisagree = 0;
      for (s=0;s<rn->Natom;++s){
        if (r_on_t[rn->num_atoms[s]]>=0){NmatchB += 1;}
        if (molT->atoms[rn->num_atoms[s]].mark!=0){Ndisagree += 1;}
      }

      if ((NmatchB==rn->Natom)&&(Ndisagree>=1)){
        for (s=0;s<6;++s){
          a = b = c = d = e = f = g = h = i = j = -1;
          angleAD_BE = angleAG_BH = M_PI;

          a = rn->num_atoms[s];
          b = rn->num_atoms[(s+2)%rn->Natom];
          c = rn->num_atoms[(s+1+rn->Natom)%rn->Natom];
          d = rn->num_atoms[(s-1+rn->Natom)%rn->Natom];
          f = rn->num_atoms[(s-2+rn->Natom)%rn->Natom];
          e = rn->num_atoms[(s+3)%rn->Natom];

/*
          if (((molT->atoms[a].mark==-1)||(molT->atoms[a].mark==2)||
               (molT->atoms[b].mark==-1)||(molT->atoms[b].mark==2)) 
             &&(molT->atoms[c].mark!=-1)&&(molT->atoms[c].mark!=2) ){
*/
          if ((molT->atoms[a].mark == 2)&&(molT->atoms[b].mark==2)&&(molT->atoms[c].mark!=2)){

            sub_vec3(AD,molT->atoms[d].Pos, molT->atoms[a].Pos); normalize_vec3(AD);
            sub_vec3(BE,molT->atoms[e].Pos, molT->atoms[b].Pos); normalize_vec3(BE);
            angleAD_BE = acos(dot_prod3(AD,BE));

            for (t=0;t<molT->atoms[a].Nneighbor;++t){
              tt = molT->atoms[a].neighbors[t];
              if ((tt!=c)&&(tt!=d)&&(check_in_ring(molT,tt)==0)){
                for (r=0;r<molT->atoms[b].Nneighbor;++r){
                  rr = molT->atoms[b].neighbors[r];
                  if ((rr!=c)&&(rr!=e)&&(check_in_ring(molT,rr)==0)){
                    sub_vec3(AG,molT->atoms[tt].Pos, molT->atoms[a].Pos); normalize_vec3(AG);
                    sub_vec3(BH,molT->atoms[rr].Pos, molT->atoms[b].Pos); normalize_vec3(BH);
                    angleAG_BH = acos(dot_prod3(AG,BH));
                    if (fabs(angleAG_BH)<=tole_angle_parallel){ g = tt; h = rr; }
                  }
                }
              } 
            }

            if ((g>=0) && (h>=0)){
              for (t=0;t<molT->atoms[a].Nneighbor;++t){
                tt = molT->atoms[a].neighbors[t];
                if ((tt!=c)&&(tt!=d)&&(tt!=g)) i = tt;
              }

              for (t=0;t<molT->atoms[b].Nneighbor;++t){
                tt = molT->atoms[b].neighbors[t];
                if ((tt!=c)&&(tt!=e)&&(tt!=h)) j = tt;
              }
            }

            printf("#CANDIDATE [A]%s %2d [B]%s %2d [C]%s %2d [D]%s %2d [E]%s %2d [F]%s %2d",
              molT->atoms[a].element, molT->atoms[a].num_in_file,
              molT->atoms[b].element, molT->atoms[b].num_in_file,
              molT->atoms[c].element, molT->atoms[c].num_in_file,
              molT->atoms[d].element, molT->atoms[d].num_in_file,
              molT->atoms[e].element, molT->atoms[e].num_in_file,
              molT->atoms[f].element, molT->atoms[f].num_in_file);
            if (g>=0) printf("[G]%s %2d ", molT->atoms[g].element, molT->atoms[g].num_in_file);
            if (h>=0) printf("[H]%s %2d ", molT->atoms[h].element, molT->atoms[h].num_in_file);
            if (i>=0) printf("[I]%s %2d ", molT->atoms[i].element, molT->atoms[i].num_in_file);
            if (j>=0) printf("[J]%s %2d ", molT->atoms[j].element, molT->atoms[j].num_in_file);
            numA[0] = d;   numA[1] = f;   numA[2] = e;   
            numB[0] = c;   numB[1] = i;   numB[2] = j;   
            insamering = check_in_the_same_ring_bwn_two_atom_groups(molT,rn,3,numA,3,numB);
            printf(" angAD_BE %.1f angAG_BH %.1f samering %d\n",angleAD_BE*180/M_PI,angleAG_BH*180/M_PI,insamering);
 
 
            if ((a>=0)&&(b>=0)&&(c>=0)&&(d>=0)&&(e>=0)&&(f>=0)&&(g>=0)&&
               (fabs(angleAD_BE)<=tole_angle_parallel)&&(fabs(angleAG_BH)<=tole_angle_parallel)&&(insamering==0)){
              
              /*** DO ROTATION ***/
              Nmodify += 1;
              molT->atoms[a].mark = 3;
              molT->atoms[b].mark = 3;
              sub_vec3(rotaxisBA,molT->atoms[b].Pos, molT->atoms[a].Pos); normalize_vec3(rotaxisBA);
              rangleDAG = acos(dot_prod3(AD,AG));
              rangleEBH = acos(dot_prod3(BE,BH));
              rangle = (rangleDAG + rangleEBH)/2.0;
              cross_vec3(ADxAG,AD,AG); normalize_vec3(ADxAG);

              printf("#DO_MODIFY_RING_CONFORMATION [A]%s %d [B]%s %d\n",
                molT->atoms[a].element, molT->atoms[a].num_in_file,
                molT->atoms[b].element, molT->atoms[b].num_in_file);
              printf("#DO MODIFICATION rangleCAF %f DBG %f ACxAF_dot_AB %f\n",rangleDAG,rangleEBH,dot_prod3(ADxAG,rotaxisBA));
              if (dot_prod3(ADxAG,rotaxisBA)<0.0) rangle = -rangle; 
              make_rotmatrix_from_rot_axis_and_angle(Rmat, rotaxisBA,  rangle);
              make_rotmatrix_from_rot_axis_and_angle(RImat,rotaxisBA, -rangle);

              for (t=0;t<molT->Natom;++t){
                 Da = molT->topodismap.map[a][t];
                 Db = molT->topodismap.map[b][t];
                 Dc = molT->topodismap.map[c][t];
                 Dd = molT->topodismap.map[d][t];
                 De = molT->topodismap.map[e][t];
                 Df = molT->topodismap.map[f][t];
                 Dg = molT->topodismap.map[g][t];
                 Dh = molT->topodismap.map[h][t];


                 minDabc = Da;
                 minDdef = Dd;
                      if ((Da<=Db)&&(Da<=Dc)) minDabc = Da;
                 else if ((Db<=Da)&&(Db<=Dc)) minDabc = Db;
                 else if ((Dc<=Da)&&(Dc<=Db)) minDabc = Dc;
               
                      if ((Dd<=De)&&(Dd<=Df)) minDdef = Dd;
                 else if ((De<=Dd)&&(De<=Df)) minDdef = De;
                 else if ((Df<=Dd)&&(Df<=De)) minDdef = Df;

                 /** (1) rotate C, E,and D with 'angle' **/ 
                 if ((t==d)||(t==e)||(t==f) || (minDdef < minDabc)){
                  rotate_vec_around_center(molT->atoms[t].Pos,Rmat,molT->atoms[a].Pos);
                 }

                /** (1) rotate F and G with '-angle' **/ 
                if ((t==g)||(Dg<Da)||(t==h)||(Dh<Db)){
                  rotate_vec_around_center(molT->atoms[t].Pos,RImat,molT->atoms[a].Pos);
                }


              }

             /*
             sprintf(buffstr,"hex%d.sdf",Nmodify); 
             Write_SDF_Molecule(buffstr,molT,"");
             */ 
            } /* do rotation */

          } /* a.mark = -1 and b.mark = -1 */

        }  /* i */

      } /* matched */
    } /* 6 */ 

  }  /* rn */
  
  out_string[0] = '\0';
  for (s=0;s<molT->Natom;++s){
     if (molT->atoms[s].mark==3){
       sprintf(buffstr," %d",molT->atoms[s].num_in_file);
       strcat(out_string,buffstr);
     }
  }

  if (out_string[0] != '\0') Add_key_value_string_to_KEYVALUEs("MODIFY_FOLD_SIX_RING",out_string,&(molT->HeadProperty));
  return(Nmodify);
} /* end of Fold_Six_Atom_Ring_for_sp3_Chirality_and_Stereo_Parity() */





int check_in_the_same_ring_bwn_two_atom_groups(mol,exclude_ring,NatomA,num_atomsA,NatomB,num_atomsB)
  struct MOLECULE *mol;
  struct RING     *exclude_ring;
  int NatomA;
  int *num_atomsA;
  int NatomB;
  int *num_atomsB;
{
  struct RING *rn;
  int i,j;
  char hitA,hitB;

  rn = &(mol->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;
    if (rn->num != exclude_ring->num){
      hitA = hitB = 0;
      for (i=0;i<rn->Natom;++i){
        for (j=0;j<NatomA;++j){
          if (rn->num_atoms[i]==num_atomsA[j]) hitA = 1;
        } 
        for (j=0;j<NatomB;++j){
          if (rn->num_atoms[i]==num_atomsB[j]) hitB = 1;
        } 
      }
      if ((hitA==1)&&(hitB==1)) return(1);
    }
  }
  return(0);
} /* end of check_in_the_same_ring_bwn_two_atom_groups() */






int check_in_ring(mol,o)
  struct MOLECULE *mol;
  int o; /* atom number of focused atom */
{
  struct RING *rn;
  int i;

  rn = &(mol->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;
    for (i=0;i<rn->Natom;++i){
      if (rn->num_atoms[i]==o) return(1);
    } 
  }
  return(0);
} /* end of check_in_ring() */







void chose_rotatable_bonded_atoms_among_4neighbors(mol,o,Nnei,anum_nei,ch_i,ch_j)
  struct MOLECULE *mol; /* (given) */
  int o;                /* atom number for the focused atom */ 
  int Nnei;             /* 3 or 4. (given). */
  int anum_nei[4];      /* atom numbers for neighbors */ 
  int *ch_i, *ch_j;     /* change i and j among {0,1,2,3} (to be calculated) */
  /*
    Bonds [o]-[ch_i], and [o]-[ch_j] should be 'rotatable'. 
  */ 
{
  int i,j;

  *ch_i = *ch_j = -1;

  if (Nnei==3){
    for (i=0;i<3;++i){ 
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;}  /* if o-i is rotatable */
    }
    if (*ch_i>=0) *ch_j = 3;  /* The last dummy atom is always chosen as 'ch_j' */
  }

  if (Nnei==4){
    for (i=0;i<4;++i){ 
      if (mol->conmap.map[o][anum_nei[i]]=='R') {*ch_i = i;} /* if a-i is rotatable */
    }
    for (j=0;j<4;++j){ 
      if ((j!=*ch_i)&&(mol->conmap.map[o][anum_nei[j]]=='R')) {*ch_j = j;} /* if a-j is rotatable */
     }
  } 

} /* end of chose_rotatable_bonded_atoms_among_4neighbors() */








void set_reference_matched_neighbor_atoms(molT,o, r_on_t, Nmatch_nei, match_neiT, match_neiR)
  struct MOLECULE *molT;  /* target molecule (input) */
  int o;              /* atom number of focused atom of molT (input) */
  int *r_on_t;        /* (input) */
  int *Nmatch_nei;    /* (to be calculaed) */
  int match_neiT[4];  /* target  (to be calculated) */  
  int match_neiR[4];  /* reference  (to be calculated) */
{
  int i,k,tt,change,buff;

  /** (1) make match_neiT[] and match_neiR[] **/
  *Nmatch_nei = 0;
  for (k=0;k<4;++k){ 
    match_neiT[k] = match_neiR[k] =  -1;
  }

  for (i=0;i<molT->atoms[o].Nneighbor;++i){
    tt = molT->atoms[o].neighbors[i];
    if (r_on_t[tt]>=0){
      match_neiT[*Nmatch_nei] = tt;
      match_neiR[*Nmatch_nei] = r_on_t[tt];
      *Nmatch_nei += 1;
    }
 }
 
  /** (2) bubble sort match_neiT[] by match_neiR[] **/
  change = 0;
  do{
    change = 0;
    for (i=1;i<(*Nmatch_nei);++i){
      if (match_neiR[i-1]>match_neiR[i]){
        buff = match_neiT[i-1];
        match_neiT[i-1] = match_neiT[i];
        match_neiT[i] = buff;
        buff = match_neiR[i-1];
        match_neiR[i-1] = match_neiR[i];
        match_neiR[i] = buff;
        change += 1;
      }
    }
  } while (change>0);

 /*** (3) If extra hydrogen atom exists, add it to match_neiT[] */
 /** this part is added in 2018/11/02 by T.Kawabata **/ 
  for (i=0;i<molT->atoms[o].Nneighbor;++i){
    tt = molT->atoms[o].neighbors[i];
    if ((r_on_t[tt]<0) && (molT->atoms[tt].one_char_ele=='H')){
      match_neiT[*Nmatch_nei] = tt;
      match_neiR[*Nmatch_nei] = -1;
      *Nmatch_nei += 1;
    }
  }

} /* end of set_reference_matched_neighbor_atoms() */



void exchange_positions_of_two_atom_groups(mol,o,Nnei,anum_nei,ch_i, ch_j)
  struct MOLECULE *mol;
  int o;                /* atom number for the focused atom */ 
  int Nnei;             /* 3 or 4. (given). */
  int anum_nei[4];      /* atom numbers for neighbors (given) */ 
  int ch_i, ch_j;     /* change i and j among {0,1,2,3} (given) */
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
  float ri[3],rj[3];
  float g[3];  /* gravity center of neighbors */
  int i,k,N;


  printf("#exchange_positions_of_two_atom_groups(o:%s %d Nnei:%d ch_i:%s %d ch_j:%s %d)\n",
     mol->atoms[o].atomname,mol->atoms[o].num_in_file, Nnei,
     mol->atoms[ch_i].atomname,mol->atoms[ch_i].num_in_file,
     mol->atoms[ch_j].atomname,mol->atoms[ch_j].num_in_file);


  g[0] = g[1] = g[2] = 0;  
  N = 0;
  for (i=0;i<4;++i){ 
    if (anum_nei[i]>=0){ 
      N += 1;
      for (k=0;k<3;++k) g[k] += mol->atoms[anum_nei[i]].Pos[k];
    }
  }
  for (k=0;k<3;++k) g[k] = g[k]/N;
  /* 
  printf("N %d\n",N);
  printf("g %f %f %f\n",g[0],g[1],g[2]);
   */ 
  if (anum_nei[ch_i]>=0){ sub_vec3(ri,mol->atoms[anum_nei[ch_i]].Pos, mol->atoms[o].Pos); }
                  else  { sub_vec3(ri,mol->atoms[o].Pos, g);}
  if (anum_nei[ch_j]>=0){ sub_vec3(rj,mol->atoms[anum_nei[ch_j]].Pos, mol->atoms[o].Pos); }
                  else  { sub_vec3(rj,mol->atoms[o].Pos, g);}

  normalize_vec3(ri);
  normalize_vec3(rj);
  /* 
  printf("ri %f %f %f rj %f %f %f\n",ri[0],ri[1],ri[2],rj[0],rj[1],rj[2]);  
  */
 
  cross_vec3(rot_axis,ri,rj);
  normalize_vec3(rot_axis);
  angle = acos(dot_prod3(ri,rj));
  make_rotmatrix_from_rot_axis_and_angle(R, rot_axis, angle);
  make_rotmatrix_from_rot_axis_and_angle(RI,rot_axis,-angle);
  for (k=0;k<mol->Natom;++k){
    if (k!=o){
      if (mol->topodismap.map[k][anum_nei[ch_i]] < mol->topodismap.map[k][o]){
        rotate_vec_around_center(mol->atoms[k].Pos,R,mol->atoms[o].Pos);
       }
       else if ((Nnei==4)&&(mol->topodismap.map[k][anum_nei[ch_j]] < mol->topodismap.map[k][o])){
         rotate_vec_around_center(mol->atoms[k].Pos,RI,mol->atoms[o].Pos);
       }
       else if ((Nnei==3)&&(mol->conmap.map[o][k] !='0') && (mol->atoms[k].one_char_ele=='H') ){
         rotate_vec_around_center(mol->atoms[k].Pos,RI,mol->atoms[o].Pos);
        }
      }
    } 

} /* end of exchange_positions_of_two_atom_groups() */



int inversion_of_pyramid(mol,o,pybase,CheckType)
  struct MOLECULE *mol;
  int o;            /* atom number for the focused atom */ 
  int pybase[3];    /* atom numbers for pyramid base. 0,1,2:base, 3:vertex */ 
  char CheckType;   /* if 'C', checking angle bwn 0-1-2 and o-included rings is < norm_angle_max */
/*
          x
          | 
          o
         /|\
        / | \
pybase 0  1  2  ==>  0  1  2 
                     \  |  /
                      \ | /
                        o
                        |
                        x
  If CheckType=='C',
   if angle bwn 0-1-2 and o-included rings is > norm_angle_max, do not perform the inversion. 
   In other words, if atom 'o' collides with another atom in the ring of the atom 'o', do not perform the inversion. 

*/ 
{
 int i,k,hit,accept;
 float posG[3],rot_axis[3],Rmat[3][3],nvec_pybase[3],nvec_ring[3],angle,nvec[3],vA[3],vB[3];
 float norm_angle_max;  /* degree */
 struct RING *rn; 
 
 /*
 printf(">> o %d pybase %d %d %d",
  mol->atoms[o].num_in_file,
  mol->atoms[pybase[0]].num_in_file, mol->atoms[pybase[1]].num_in_file, mol->atoms[pybase[2]].num_in_file);
 if (pybase[3]<0) printf(" ---\n"); else printf(" %d\n",mol->atoms[pybase[3]].num_in_file);
 */

 if (pybase[0]<0) {return(0);}

 norm_angle_max = 70.0;
 /** [1] caluclate posG[] **/
 for (k=0;k<3;++k){
   posG[k] = (mol->atoms[pybase[0]].Pos[k] + mol->atoms[pybase[1]].Pos[k] + mol->atoms[pybase[2]].Pos[k])/3.0;
 }

 /** [2] calculate Rmat[][] **/
 sub_vec3(rot_axis,posG,mol->atoms[pybase[0]].Pos); 
 normalize_vec3(rot_axis);
 make_rotmatrix_from_rot_axis_and_angle(Rmat, rot_axis, M_PI);

 /** [3] calculate nvec_pybase[] **/
 nvec_pybase[0] = nvec_pybase[1] = nvec_pybase[2] = 0.0;
 for (i=0;i<3;++i){
   sub_vec3(vA,mol->atoms[pybase[i]].Pos,posG);
   sub_vec3(vB,mol->atoms[pybase[(i+1)%3]].Pos,posG);
   cross_vec3(nvec,vA,vB); normalize_vec3(nvec);
   for (k=0;k<3;++k) nvec_pybase[k] += nvec[k];
 }
 normalize_vec3(nvec_pybase);


 /** [4] Checking angle between pyramid base plane and plane of ring through the atom o  **/
 accept = 1;
 if (CheckType == 'C'){
   rn = &(mol->HeadRing);
   while (rn->next != NULL){
     rn = rn->next;
     hit = 0;
     for (i=0;i<rn->Natom;++i){
       if (rn->num_atoms[i]==o) hit = 1;
     }
     if (hit == 1){
       normal_vector_of_ring(nvec_ring, mol, rn);
       angle = acos(dot_prod3(nvec_pybase,nvec_ring))*180.0/M_PI;
       if ((180.0-angle)<angle){  angle = 180.0 - angle;}
       /* printf("#angle %f\n",angle);  */
       if (angle>norm_angle_max){ 
         accept = 0;
       } 
     }
   }
 }

 /** [5] Do Inversion **/
 printf("#inversion_of_pyramid(): CheckType %c accept %d o (%s %d) pybase (%s %d|%s %d|%s %d)\n",
    CheckType, accept,
    mol->atoms[o].atomname, mol->atoms[o].num_in_file,
    mol->atoms[pybase[0]].atomname, mol->atoms[pybase[0]].num_in_file,
    mol->atoms[pybase[1]].atomname, mol->atoms[pybase[1]].num_in_file,
    mol->atoms[pybase[2]].atomname, mol->atoms[pybase[2]].num_in_file);
 
 if (accept == 1){
   for (i=0;i<mol->Natom;++i){
     if ((mol->topodismap.map[i][o] < mol->topodismap.map[i][pybase[0]]) &&
         (mol->topodismap.map[i][o] < mol->topodismap.map[i][pybase[1]]) &&
         (mol->topodismap.map[i][o] < mol->topodismap.map[i][pybase[2]]) ){
       rotate_vec_around_center(mol->atoms[i].Pos,Rmat,posG);
    }
   }
   return(1); 
 }
 else{
   return(0);
 }
} /* end of inversion_of_pyramid() */



void normal_vector_of_ring(normalvec,mol,rn)
  float normalvec[3];
  struct MOLECULE *mol;
  struct RING *rn;
{
  int i,k;
  float posG[3],vecA[3],vecB[3],norm[3];

 /*
  printf("\nRING [");
  for (i=0;i<rn->Natom;++i){
    printf(" %d",mol->atoms[rn->num_atoms[i]].num_in_file);
  }
  printf("] ");
 */
 
  for (k=0;k<3;++k) {posG[k] = normalvec[k] = 0.0; }
  for (i=0;i<rn->Natom;++i){
    for (k=0;k<3;++k) {posG[k] += mol->atoms[rn->num_atoms[i]].Pos[k];} 
  }
  for (k=0;k<3;++k) {posG[k] /= rn->Natom;}

  for (i=0;i<rn->Natom;++i){
    sub_vec3(vecA,mol->atoms[rn->num_atoms[i]].Pos,posG);
    sub_vec3(vecB,mol->atoms[rn->num_atoms[(i+1)%rn->Natom]].Pos,posG);
    cross_vec3(norm,vecA,vecB);
    normalize_vec3(norm);
    for (k=0;k<3;++k) {normalvec[k] += norm[k];}
  } 

  normalize_vec3(normalvec);

} /* end of normal_vector_of_ring() */







void Reflection_of_All_Atoms(mol)
  struct MOLECULE *mol;
{
  int i;


  printf("#DO_ALL_ATOM_REFLECTION (minus X-coordinate) !!!\n");
  /*
   THIS IS VERY SIMPLE AND CRUDE WAY TO MODIFY STEREO PARITY !! 
  */
  for (i=0;i<mol->Natom;++i){
    mol->atoms[i].Pos[0] = -mol->atoms[i].Pos[0];
  }
  
  Add_key_value_string_to_KEYVALUEs("MODIFY_REFLECT","Refletion_of_All_Atoms",&(mol->HeadProperty));

} /* end of Reflection_of_All_Atoms() */

 


 
void make_r_on_t_from_MATCH(r_on_t, M, molT)
  int *r_on_t;  /* int[mol->Natom] (should be already malloced) **/
  struct MATCH *M;
  struct MOLECULE *molT;
{
  int a,m;

  for (a=0;a<molT->Natom;++a){ 
    r_on_t[a] = -1; 
  }
  for (m=0;m<M->Npair;++m){
    r_on_t[M->anumA[m]] = M->anumB[m];
  }

} /* end of make_r_on_t_from_MATCH() */


void chose_pybase_from_match_nei(pybase,ch_i,ch_j,Nmatch_nei,match_nei)
  int pybase[3]; /* to be calculated */
  int ch_i, ch_j;
  int Nmatch_nei;
  int match_nei[4];
{
  int i,Npybase;

  pybase[0] = pybase[1] = pybase[2] = -1;
  if (((ch_i<0)||(ch_j<0))&&(Nmatch_nei==3)){ 
   for (i=0;i<3;++i) pybase[i] =  match_nei[i];
  } 

  if ((ch_i>=0)&&(ch_j<0)&&(Nmatch_nei==4)) {
      Npybase = 0;       
      for (i=0;i<4;++i){
      if (i!=ch_i) {pybase[Npybase] = match_nei[i]; Npybase += 1;}
    }
  }
} /* end of chose_pybase_from_match_nei() */
