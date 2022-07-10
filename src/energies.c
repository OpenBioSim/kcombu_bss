/*

< energies.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

functions of energies to transform molecular conformation


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
#include "transform.h"
#include "vector3.h"
#include "stamp_transf.h"

/*** FUNCTIONS (GLOBAL) ***/
int  Transform_MOLECULE_by_RotBond();
void Transform_MOLECULE_by_Translation_Rotation();
float energy_total();
void  show_energy_total();
void cal_energy_for_transform();
float energy_atommatch();
float energy_selfclash();
float energy_protclash();
int   mark_ligand_neighbor_receptor_atoms_for_Eprotclash();
float energy_protatrct();
float energy_volmovlap();
float energy_tpdisrest();
float sphere_volume_overlap();
void  show_selfclash();
int   count_selfclash();
void  show_protclash();
float fixed_rmsd_atommatch();
struct MATCH *chose_MATCH_with_min_rmsd();
void  cal_Force_and_Torque();
void  cal_dEd_rbAngle();
void Random_Transform_MOLECULE_by_RotBond();

/*** FUNCTIONS (LOCAL) ***/



int Transform_MOLECULE_by_RotBond(rnum,rbangle,tra,mol)
  int rnum; 
  float rbangle;  /* radian */
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *mol;
{
  int z,o,k,a,natom_tra;
  float axis[3],R[3][3];
 /*
  printf("#void Transform_MOLECULE_by_RotBond(rnum:%d ,rbangle:%f,tra,mol)\n",rnum,rbangle);
 */
  z = tra->rbZanum[rnum];  
  o = tra->rbOanum[rnum];  
  for (k=0;k<3;++k){ 
   axis[k] = mol->atoms[o].Pos[k] - mol->atoms[z].Pos[k];  
  } 
  normalize_vec3(axis);
  /* printf("axis %f %f %f\n",axis[0],axis[1],axis[2]); */
  make_rotmatrix_from_rot_axis_and_angle(R,axis,rbangle);
  natom_tra = 0;
  for (a=0;a<mol->Natom;++a){
    /* printf("tra->rbMrnum[%d][%d]:%d\n",rnum,a,tra->rbManum[rnum][a]); */
    if (tra->rbManum[rnum][a]==1){ 
      rotate_vec_around_center(mol->atoms[a].Pos, R, mol->atoms[o].Pos);
      natom_tra += 1;
    }
  }

 return(natom_tra);
} /* end of Transform_MOLECULE_by_RotBond() */





void Transform_MOLECULE_by_Translation_Rotation(dT,rotaxis,rotangle,tra,mol)
  float  dT[3],rotaxis[3];
  float rotangle; /* radian */
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *mol;
/*
  Rotation is done around the center atom of the molecule (xc).
    xi'  = xi + dT
    xi'' = R[dW](xi' - xc') + xc'
     =>
    xi'' = R[dW](xi + dT - (xc+dT)) + xc + dT 
         = R[dW](xi-xc) + xc + dT 
*/ 
{
  int k,a;
  float R[3][3];

  /** (1) Translation by dT **/
  for (a=0;a<mol->Natom;++a){
    for (k=0;k<3;++k){
      mol->atoms[a].Pos[k] += dT[k];
    } 
  }
 /** (2) Rotation by dW **/
 if (rotangle>0.0){
  make_rotmatrix_from_rot_axis_and_angle(R,rotaxis,rotangle);
  for (a=0;a<mol->Natom;++a){
    rotate_vec_around_center(mol->atoms[a].Pos,R,mol->atoms[tra->Canum].Pos);
  }
 }

} /* end of Transform_MOLECULE_by_Translation_Rotation() */










float energy_total(molT,molR,molP,M)
  struct MOLECULE *molT,*molR,*molP;
  struct MATCH *M;
{
  float Eatommatch, Eselfclash,Eprotclash, Eprotatrct, Evolmovlap, Etpdisrest,E;
     
  Eatommatch = Eselfclash = Eprotclash = Eprotatrct = Evolmovlap = Etpdisrest = 0.0;

  if (PAR.WEatommatch > 0.0){ Eatommatch = energy_atommatch(molT,molR,M); }

  if (PAR.WEselfclash > 0.0){ Eselfclash = energy_selfclash(molT); }

  if (PAR.WEprotclash > 0.0){ Eprotclash = energy_protclash(molT,molP); }

  if (PAR.WEprotatrct > 0.0){ Eprotatrct = energy_protatrct(molT,molP); }
  
  if (PAR.WEvolmovlap > 0.0){ Evolmovlap = energy_volmovlap(molT,molR); }

  if (PAR.WEtpdisrest > 0.0){ Etpdisrest = energy_tpdisrest(molT,molP); }
 
  E = 0.0; 
  E += PAR.WEatommatch * Eatommatch;
  E += PAR.WEselfclash * Eselfclash;
  E += PAR.WEprotclash * Eprotclash;
  E += PAR.WEprotatrct * Eprotatrct;
  E += PAR.WEvolmovlap * Evolmovlap;
  E += PAR.WEtpdisrest * Etpdisrest;
  return(E);

} /* end of energy_total() */



void show_energy_total(molT,molR,molP,M)
  struct MOLECULE *molT,*molR,*molP;
  struct MATCH *M;
{
  float Eatommatch, Eselfclash,Eprotclash, Eprotatrct, Evolmovlap, Etpdisrest,E;
     
  Eatommatch = Eselfclash = Eprotclash = Eprotatrct = Evolmovlap = Etpdisrest = 0.0;

  Eatommatch = energy_atommatch(molT,molR,M); 
  Eselfclash = energy_selfclash(molT); 
  Eprotclash = energy_protclash(molT,molP); 
  Eprotatrct = energy_protatrct(molT,molP); 
  Evolmovlap = energy_volmovlap(molT,molR); 
  Etpdisrest = energy_tpdisrest(molT,molP); 
 
  E = 0.0; 
  E += PAR.WEatommatch * Eatommatch;
  E += PAR.WEselfclash * Eselfclash;
  E += PAR.WEprotclash * Eprotclash;
  E += PAR.WEprotatrct * Eprotatrct;
  E += PAR.WEvolmovlap * Evolmovlap;
  E += PAR.WEtpdisrest * Etpdisrest;
  printf("#Eatommatch %f\n",Eatommatch);
  printf("#Eselfclash %f\n",Eselfclash);
  printf("#Eprotclash %f\n",Eprotclash);
  printf("#Eprotatrct %f\n",Eprotatrct);
  printf("#Evolmovlap %f\n",Evolmovlap);
  printf("#Etpdisrest %f\n",Etpdisrest);
  printf("#Natommatch %d\n",M->Npair);
} /* end of show_energy_total() */




void cal_energy_for_transform(tra,molT,molR,molP,M)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *molT; /* 'T'arget */
  struct MOLECULE *molR; /* 'R'eference */
  struct MOLECULE *molP; /* recepttor 'P'rotein */
  struct MATCH *M;
{
  tra->Eatommatch = tra->Eselfclash = tra->Eprotclash = tra->Eprotatrct = tra->Evolmovlap = tra->Etpdisrest = 0.0;

  tra->Eatommatch =  energy_atommatch(molT,molR,M); 

  tra->Eselfclash =  energy_selfclash(molT); 

  tra->Eprotclash =  energy_protclash(molT,molP);
  
  tra->Eprotatrct =  energy_protatrct(molT,molP);
  
  tra->Evolmovlap =  energy_volmovlap(molT,molR);
 
  tra->Etpdisrest =  energy_tpdisrest(molT,molP);

  tra->Nselfclash =  count_selfclash(molT);

  tra->Etotal =   PAR.WEatommatch * tra->Eatommatch + PAR.WEselfclash * tra->Eselfclash 
                + PAR.WEprotclash * tra->Eprotclash + PAR.WEprotatrct * tra->Eprotatrct 
                + PAR.WEvolmovlap * tra->Evolmovlap + PAR.WEtpdisrest * tra->Etpdisrest;

  tra->tanimoto_volume = -energy_volmovlap(molT,molR)/(-energy_volmovlap(molT,molT)-energy_volmovlap(molR,molR) + energy_volmovlap(molT,molR));

} /* end of cal_energy_for_transform() */









float energy_atommatch(molT,molR,M)
  struct MOLECULE *molT,*molR;
  struct MATCH *M;
{
  int p,a,b,k;
  float d[3],sum; 

  sum = 0.0;

  for (p=0;p<M->Npair;++p){
    a = M->anumA[p];
    b = M->anumB[p];
    for (k=0;k<3;++k){
      d[k] = molT->atoms[a].Pos[k] - molR->atoms[b].Pos[k];
    } 
    sum += d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; 
  }
  /*
  printf("ene %f RMSD %f\n",0.5*sum,sqrt(sum/M->Npair));
  */ 
  return(0.5*sum);
} /* end of energy_atommatch() */




float energy_selfclash(mol)
  struct MOLECULE *mol;
{
  int a,b,k;
  float d,dis,disrad,E,DDthre; 
/*
 Eselfclash = 1 - 1/(1+exp[-alpha(d-Ra-Rb+tole)]
*/
  DDthre = PAR.Dmax_clash * PAR.Dmax_clash;
  E = 0.0;
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].one_char_ele != 'H'){
      for (b=a+1;b<mol->Natom;++b){
        if (mol->atoms[b].one_char_ele != 'H'){
          if ((mol->conmap.map[a][b]=='0')&&(mol->topodismap.map[a][b]>2)){
            dis = 0.0;
            for (k=0;k<3;++k){
              d = mol->atoms[a].Pos[k] - mol->atoms[b].Pos[k];
              dis += d*d;
            }
            if (dis < DDthre){
              dis = sqrt(dis); 
              disrad = dis - mol->atoms[a].radius - mol->atoms[b].radius + PAR.tole_dis_clash;
              E += 1.0 - 1.0/(1+exp(-PAR.alpha_sigmoid_clash *disrad)); 
            }
          }
        }
      }
    }
  }
  return(E);
} /* end of energy_selfclash() */




float energy_protclash(molT, molP)
  struct MOLECULE *molT; /* target molecule  */
  struct MOLECULE *molP; /* receptor molecule */
{
  int a,p,k;
  float d,dis,disrad,E,DDthre; 

/*
 Eprotclash = 1 - 1/(1+exp[-alpha(d-Ra-Rr+tole)]

  NOTE :: energy_proclash only consider receptor atoms[r] with molP->atoms[r].mark=1 !!


 */
  E = 0.0;
  DDthre = PAR.Dmax_clash*PAR.Dmax_clash;

  for (p=0;p<molP->Natom;++p){
    if (molP->atoms[p].one_char_ele != 'H'){
      if (molP->atoms[p].mark==1){   /* <-- only considers mark==1 molP atoms !! */
        for (a=0;a<molT->Natom;++a){
          if (molT->atoms[a].one_char_ele != 'H'){
            dis = 0.0;
            for (k=0;k<3;++k){
              d = molT->atoms[a].Pos[k] - molP->atoms[p].Pos[k];
              dis += d*d;
            }
            if (dis<DDthre){
              dis = sqrt(dis); 
              disrad = dis - molT->atoms[a].radius - molP->atoms[p].radius + PAR.tole_dis_clash;
              E += 1.0 - 1.0/(1+exp(-PAR.alpha_sigmoid_clash *disrad)); 
            }
          }
        }
      }
    }
  }

  return(E);
} /* end of energy_protclash() */





int mark_ligand_neighbor_receptor_atoms_for_Eprotclash(molR,molL,Dthre)
  struct MOLECULE *molR; /* Receptor molecule */
  struct MOLECULE *molL; /* Ligand   molecule */
  float  Dthre;          /* Distance threshold for neighbor */
/*
  Note: both Eprotclash and Eprotatrct require mark[] = 1.

*/
{
  int r,i,k,Nrepatm_nei;
  float DDthre,d,DD;

  DDthre = Dthre * Dthre;
  Nrepatm_nei = 0;
  for (r=0;r<molR->Natom;++r){
    molR->atoms[r].mark = 0;
    i = 0;
    while ((i<molL->Natom) && (molR->atoms[r].mark==0)){
      DD = 0.0;
      for (k=0;k<3;++k){
        d = molR->atoms[r].Pos[k] - molL->atoms[i].Pos[k];
        DD += d*d;
      }
      if (DD<=DDthre){
         molR->atoms[r].mark = 1;
      }
      i += 1;
    }
    if (molR->atoms[r].mark == 1){
      Nrepatm_nei += 1;
    }
  }

 return(Nrepatm_nei);

} /* end of mark_ligand_neighbor_receptor_atoms_for_Eprotclash() */



/*
float energy_protatrct(molT, molP)
  struct MOLECULE *molT;  target molecule  
  struct MOLECULE *molP;  receptor (protin) molecule 
{
  int a,p,k;
  float d,dis,disrad,E,DDthre; 

 Eprotatrct =   1 - 1/(1+exp[-alpha(d-Ra-Rp+tole)]
              -(1 - 1/(1+exp[-alpha(d-Ra-Rp-tole)])
            =  -1/(1+exp[-alpha(d-Ra-Rp+tole)] +  1/(1+exp[-alpha(d-Ra-Rp-tole)]

( Only consider receptor atoms[r] with molP->atoms[p].mark=1.)

  E = 0.0;
  DDthre = PAR.Dmax_clash*PAR.Dmax_clash;

  for (p=0;p<molP->Natom;++p){
    if (molP->atoms[p].mark==1){
      for (a=0;a<molT->Natom;++a){
        dis = 0.0;
        for (k=0;k<3;++k){
          d = molT->atoms[a].Pos[k] - molP->atoms[p].Pos[k];
          dis += d*d;
        }
        if (dis<DDthre){
          dis = sqrt(dis); 
          disrad = dis - molT->atoms[a].radius - molP->atoms[p].radius;
          E += - 1.0/(1+exp(-PAR.alpha_sigmoid_clash*(disrad + PAR.tole_dis_clash)));
          E +=   1.0/(1+exp(-PAR.alpha_sigmoid_clash*(disrad - PAR.tole_dis_clash))); 
        }
      }
    }
  }
  return(E);
} 
*/

float energy_protatrct(molT, molP)
  struct MOLECULE *molT;  /* target molecule  */
  struct MOLECULE *molP;  /* receptor (protein) molecule */
{
  int a,p,k;
  float d,dis,disrad,E,DDthre;
 
/*
 Eprotatrct =   1 - 1/(1+exp[-alpha(d-Ra-Rp)]
              -(1 - 1/(1+exp[-alpha(d-Ra-Rp-tole)])
            =  -1/(1+exp[-alpha(d-Ra-Rp)] +  1/(1+exp[-alpha(d-Ra-Rp-tole)]

( Only consider receptor atoms[r] with molP->atoms[p].mark=1.)
*/
  E = 0.0;
  DDthre = PAR.Dmax_clash*PAR.Dmax_clash;

  for (p=0;p<molP->Natom;++p){
    if (molP->atoms[p].one_char_ele != 'H'){
      if (molP->atoms[p].mark==1){
        for (a=0;a<molT->Natom;++a){
          if (molT->atoms[a].one_char_ele != 'H'){
            dis = 0.0;
            for (k=0;k<3;++k){
              d = molT->atoms[a].Pos[k] - molP->atoms[p].Pos[k];
              dis += d*d;
            }
            if (dis<DDthre){
              dis = sqrt(dis); 
              disrad = dis - molT->atoms[a].radius - molP->atoms[p].radius;
              E += - 1.0/(1+exp(-PAR.alpha_sigmoid_clash*(disrad)));
              E +=   1.0/(1+exp(-PAR.alpha_sigmoid_clash*(disrad - PAR.tole_dis_atrct))); 
            }
          }
        }
      }
    }
  }
/*
  printf("#energy_protatrct tole_dis_atrct %f %f %e\n",PAR.tole_dis_atrct,E,E);
 */
  return(E);
} /* end of energy_protatrct() */






float energy_volmovlap(molA,molB)
  struct MOLECULE *molA,*molB;
{
/*

 If density of the i-th atom is represented as:

  fi(r) = p * exp[-pi[3p/4pi]^(2/3) (r-rj)^2/Ri^2]  

 The overlap of i-th and j-th atoms is represented as:
  
  ov(fi,fj) = 
   4/3*pi*p*[Ri^2*Rj^2/(Ri^2+Rj^2)]^3/2 exp[-pi*(3p/4pi)**(2/3) (ri-rj)^2/(Ri^2+Rj^2)]  

  Evolovlap = -ov(fi,fj)

*/
  int a,b,k;
  float d,dd,ene,rrA,rrB; 
  float Cvo,alpha_vo;  /* Evo = -Cvo * exp(-alpha_vo*dd) */ 
  float c1vo;          /* constant :: 4/3*pi*p           */
  float c2vo;          /* constant :: pi*(3p/4pi)**(2/3) */
  float Watomtype;     /* Weight for atomtypes */

  c1vo = 4.0/3.0 * M_PI * PAR.p_volmovlap;
  c2vo = M_PI * pow(3.0*PAR.p_volmovlap/4.0/M_PI,2.0/3.0);

  ene = 0.0;
  for (a=0;a<molA->Natom;++a){
    if (molA->atoms[a].one_char_ele != 'H'){
      rrA = molA->atoms[a].radius * molA->atoms[a].radius;
      for (b=0;b<molB->Natom;++b){
        if (molB->atoms[b].one_char_ele != 'H'){ 
            if (PAR.atomtype_volmovlap!='C') Watomtype=1.0;
          else Watomtype = Score_bwn_atomtypes(molA->atoms[a].atomtype,molB->atoms[b].atomtype);

          rrB = molB->atoms[b].radius * molB->atoms[b].radius;
          dd = 0.0; 
          for (k=0;k<3;++k){
             d = molA->atoms[a].Pos[k] - molB->atoms[b].Pos[k];
             dd += d*d;
           }
          Cvo      = c1vo * pow(rrA*rrB/(rrA+rrB),3.0/2.0);
          alpha_vo = c2vo/(rrA+rrB);
          ene -= Watomtype * Cvo*exp(-alpha_vo*dd); 
          /* printf("a %d b %d %lf\n",molA->atoms[a].num_in_file,b,C*exp(-alpha*dd));  */
         }
      }
    }
  }
  return(ene);
} /* end of energy_volmovlap() */


float energy_tpdisrest(molT, molP)
  struct MOLECULE *molT; /* target molecule  */
  struct MOLECULE *molP; /* receptor molecule */
{
  int i,t,cal_rest;
  float d,dis,E; 

  E = 0.0;
  for (t=0;t<3;++t){
    cal_rest = 1;
    if (PAR.string_tpdisrest[t][0]=='\0'){cal_rest = 0;}
    if ((PAR.Tatom_num_tpdisrest[t]<0) || (PAR.Tatom_num_tpdisrest[t]>=molT->Natom)){cal_rest = 0;}
    if ((PAR.Patom_num_tpdisrest[t]<0) || (PAR.Patom_num_tpdisrest[t]>=molP->Natom)){cal_rest = 0;}

    if (cal_rest == 1){
      dis = 0.0;
      for (i=0;i<3;++i){
        d = molT->atoms[PAR.Tatom_num_tpdisrest[t]].Pos[i] - molP->atoms[PAR.Patom_num_tpdisrest[t]].Pos[i];
        dis += d*d;
      }
  
      if (dis > PAR.Dupper_tpdisrest[t]*PAR.Dupper_tpdisrest[t]){
        dis = sqrt(dis);
        E += 0.5 * (dis - PAR.Dupper_tpdisrest[t]) * (dis - PAR.Dupper_tpdisrest[t]); 
      }
   }
 }
 return(E); 

} /* end of energy_tpdisrest() */



float sphere_volume_overlap(molA,molB)
  struct MOLECULE *molA,*molB;
{
  int a,b,k;
  float d,dAB,vol,rA,rB,rrA,rrB,rrrA,rrrB,eps; 

  vol  = 0.0;
  eps = 0.00000001;
  for (a=0;a<molA->Natom;++a){
    if (molA->atoms[a].one_char_ele != 'H'){
      rA = molA->atoms[a].radius;
      rrA  = rA * rA;
      rrrA = rA * rrA;
      for (b=0;b<molB->Natom;++b){
        if (molB->atoms[b].one_char_ele != 'H'){
          dAB = 0.0; 
          rB = molB->atoms[b].radius;
          rrB = rB * rB;
          rrrB = rB * rrB;
          for (k=0;k<3;++k){
             d = molA->atoms[a].Pos[k] - molB->atoms[b].Pos[k];
             dAB += d;
           }
           dAB = sqrt(dAB);
           if ((dAB < (rA+rB)) && (dAB>eps)){
            /*
            vol += M_PI/12*(rA+rB-dAB)*(rA+rB-dAB)*(dAB + 2*(rA+rB)-3/dAB * (rA-rB)*(rA-rB)); 
            printf("dAB %lf vol %lf\n",dAB,vol);
            */
             vol += M_PI/12.0*(8.0*(rrrA+rrrB)-3.0*(rrA-rrB)*(rrA-rrB)/dAB - 6.0*dAB*(rrA+rrB) + dAB*dAB*dAB); 
           
           }
           if (dAB < eps){
             if (rA<rB) vol += 4.0/3.0*M_PI*rrrA;
                  else  vol += 4.0/3.0*M_PI*rrrB;
           }

         }
      }
    }
  }
  return(vol);
} /* end of sphere_volume_overlap() */



void show_selfclash(mol,header_string)
  struct MOLECULE *mol;
  char   *header_string;
{
  int a,b,k;
  float d,dis,disrad,E,DDthre; 
/*
 Eselfclash = 1 - 1/(1+exp[-alpha(d-Ra-Rb+tole)]
*/

  DDthre = PAR.Dmax_clash*PAR.Dmax_clash;
  E = 0.0;


  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].one_char_ele != 'H'){
      for (b=a+1;b<mol->Natom;++b){
        if (mol->atoms[b].one_char_ele != 'H'){
          if ((mol->conmap.map[a][b]=='0')&&(mol->topodismap.map[a][b]>2)){
            dis = 0.0;
            for (k=0;k<3;++k){
              d = mol->atoms[a].Pos[k] - mol->atoms[b].Pos[k];
              dis += d*d;
            }
            if (dis < DDthre){
              dis = sqrt(dis); 
              disrad = dis - mol->atoms[a].radius - mol->atoms[b].radius + PAR.tole_dis_clash;
              E += 1.0 - 1.0/(1+exp(-PAR.alpha_sigmoid_clash *disrad)); 
              if (disrad < 0.0){
                printf("#SELFCRASH%s %4d %s R %f vs %4d %s R %f tole %f dis %f disrad %f\n",
                  header_string,
                  mol->atoms[a].num_in_file,mol->atoms[a].atomname,mol->atoms[a].radius,
                  mol->atoms[b].num_in_file,mol->atoms[b].atomname,mol->atoms[b].radius,
                  PAR.tole_dis_clash,dis,disrad);
               }
            }
          }
        }
      }
    }
  }

} /* end of show_selfclash() */





int count_selfclash(mol)
  struct MOLECULE *mol;
{
  int a,b,k,Npair_clash;
  float d,dis,disrad,DDthre; 
/*
 Eselfclash = 1 - 1/(1+exp[-alpha(d-Ra-Rb+tole)]
*/
  Npair_clash = 0;
  DDthre = PAR.Dmax_clash*PAR.Dmax_clash;
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].one_char_ele != 'H'){
      for (b=a+1;b<mol->Natom;++b){
        if (mol->atoms[b].one_char_ele != 'H'){
          if ((mol->conmap.map[a][b]=='0')&&(mol->topodismap.map[a][b]>2)){
            dis = 0.0;
            for (k=0;k<3;++k){
              d = mol->atoms[a].Pos[k] - mol->atoms[b].Pos[k];
              dis += d*d;
            }
            if (dis < DDthre){
              dis = sqrt(dis); 
              disrad = dis - mol->atoms[a].radius - mol->atoms[b].radius + PAR.tole_dis_clash;
              if (disrad < 0.0){ Npair_clash += 1; }
            }
          }
        }
      }
    }
  }
  return(Npair_clash);
} /* end of count_selfclash() */






void show_protclash(molT, molP)
  struct MOLECULE *molT; /* target molecule  */
  struct MOLECULE *molP; /* receptor molecule */
{
  int a,r,k;
  float d,dis,disrad,e,E,DDthre; 
  E = 0.0;
  printf("#show_protclash(Dmax_clash %f tole_dis_clash %f)\n",PAR.Dmax_clash, PAR.tole_dis_clash);
  DDthre = PAR.Dmax_clash*PAR.Dmax_clash;

  for (r=0;r<molP->Natom;++r){
    if (molP->atoms[r].mark==1){
      for (a=0;a<molT->Natom;++a){
        dis = 0.0;
        for (k=0;k<3;++k){
          d = molT->atoms[a].Pos[k] - molP->atoms[r].Pos[k];
          dis += d*d;
        }
        if (dis<DDthre){
          dis = sqrt(dis); 
          disrad = dis - molT->atoms[a].radius - molP->atoms[r].radius + PAR.tole_dis_clash;
          e = 1.0 - 1.0/(1+exp(-PAR.alpha_sigmoid_clash *disrad)); 
          E += e;

        if (disrad < 0.0){
           printf("#RCPTCRASH %4d %s R %f vs %4d %s R %f tole %f dis %f disrad %f eprotclash %f\n",
             molT->atoms[a].num_in_file,molT->atoms[a].atomname, molT->atoms[a].radius,
             molP->atoms[r].num_in_file,molP->atoms[r].atomname, molP->atoms[r].radius,
             PAR.tole_dis_clash,dis,disrad,e);
        }
       }
     }
    }
  }
  printf("#Eprotclash %f %f\n",E, energy_protclash(molT,molP)); 

} /* end of show_protclash() */



float fixed_rmsd_atommatch(molA,molB,M)
  struct MOLECULE *molA,*molB;
  struct MATCH *M;
{
  int p,a,b,k;
  float d[3],sum,rmsd; 

  if (M->Npair==0) return(0.0);

  sum = 0.0;
  for (p=0;p<M->Npair;++p){
    a = M->anumA[p];
    b = M->anumB[p];
    for (k=0;k<3;++k){
      d[k] = molA->atoms[a].Pos[k] - molB->atoms[b].Pos[k];
    } 
    sum += d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; 
  }
  rmsd = sqrt(sum/M->Npair);
  return(rmsd);

} /* end of fixed_rmsd_atommatch() */


struct MATCH  *chose_MATCH_with_min_rmsd(molA,molB,Mlist)
  struct MOLECULE *molA,*molB;
  struct MATCH *Mlist;
{
  float rmsd,rmsd_min;
  struct MATCH *M,*optM;
  int N;

  printf("#float chose_MATCH_with_min_rmsd(molA,molB,Mlist)\n");
  N = 0;
  M = Mlist;
  rmsd_min = -1.0;
  optM = NULL;

  while (M->next!=NULL){
    M = M->next;
    rmsd =  fixed_rmsd_atommatch(molA,molB,M);
    N += 1;
    if ((rmsd_min < 0.0) || (rmsd < rmsd_min)){ rmsd_min = rmsd; optM = M;}
  }
  return(optM);
}



void cal_Force_and_Torque(tra, molT, molR,molP, M)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *molT,*molR, *molP;
  struct MATCH *M;
{
  int i,j,t,r,k,m,p;
  float x_r[3],x_c[3],tor[3],x_p[3],Rij,Rijrad,Rxp,Rxprad,expR,dSdR,DDthre_clash;
  float rrT,rrR,RRtr,t_r[3],pi_pj[3],force_nume[3],delta,E0,E1;
  char  check_numerical_force;
  float Cvo,alpha_vo;  /* Evo = -Cvo * exp(-alpha_vo*dd) */ 
  float c1vo;          /* constant :: 4/3*pi*p           */
  float c2vo;          /* constant :: pi*(3p/4pi)**(2/3) */
  float Watomtype;     /* Weight for atomtypes */

  DDthre_clash = PAR.Dmax_clash * PAR.Dmax_clash;

  /*** [0] initialize **/ 
  for (k=0;k<3;++k){ 
    tra->Force[k] = tra->Torque[k] = 0.0;
    for (m=0;m<M->Npair;++m){
    }
  }
  for (t=0;t<molT->Natom;++t){ 
    for (k=0;k<3;++k){ tra->force_atom[t][k] = 0.0; }
  }

 /*** [1] Calculate force_atom[][] for Eatommatch **/ 
 if (PAR.WEatommatch>0.0){
   for (m=0;m<M->Npair;++m){
     t = M->anumA[m];
     r = M->anumB[m];
     if ((t>=0) && (r>=0)){
       sub_vec3(x_r,molT->atoms[t].Pos,molR->atoms[r].Pos); 
       for (k=0;k<3;++k){tra->force_atom[t][k] -= PAR.WEatommatch * x_r[k];}
     }
   }
 } 

 /*** [2] Calculate force_atom[][] for Eselfclash **/ 
 if (PAR.WEselfclash>0.0){
   for (i=0;i<molT->Natom;++i){
     if (molT->atoms[i].one_char_ele != 'H'){
       for (j=i+1;j<molT->Natom;++j){
         if (molT->atoms[j].one_char_ele != 'H'){
           if ((i!=j)&&(molT->conmap.map[i][j]=='0')&&(molT->topodismap.map[i][j]>2)){
             sub_vec3(pi_pj, molT->atoms[i].Pos, molT->atoms[j].Pos);
             Rij = pi_pj[0]*pi_pj[0] + pi_pj[1]*pi_pj[1] + pi_pj[2]*pi_pj[2]; 
             if (Rij<DDthre_clash){
               Rij = sqrt(Rij); 
               Rijrad = Rij - molT->atoms[i].radius - molT->atoms[j].radius + PAR.tole_dis_clash;
               expR = exp(-PAR.alpha_sigmoid_clash*Rijrad);
               dSdR = - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
               for (k=0;k<3;++k){tra->force_atom[i][k] -= PAR.WEselfclash * dSdR * pi_pj[k] /Rij;}
               for (k=0;k<3;++k){tra->force_atom[j][k] -= PAR.WEselfclash * dSdR * pi_pj[k] /Rij;}
            } 
           } 
         } 
        }
      }
   }
 }

 /*** [3] Calculate force_atom[][] for Eprotclash **/ 
 if (PAR.WEprotclash > 0.0){
   for (t=0;t<molT->Natom;++t){
     if (molT->atoms[t].one_char_ele != 'H'){
       for (p=0;p<molP->Natom;++p){
         if (molP->atoms[p].one_char_ele != 'H'){
           if (molP->atoms[p].mark==1){
             sub_vec3(x_p,molT->atoms[t].Pos,molP->atoms[p].Pos); 
             Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
             if (Rxp<DDthre_clash){
               Rxp = sqrt(Rxp);
               Rxprad = Rxp - molT->atoms[t].radius - molP->atoms[p].radius + PAR.tole_dis_clash;
               expR = exp(-PAR.alpha_sigmoid_clash*Rxprad);
               dSdR = - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
               for (k=0;k<3;++k){
                 tra->force_atom[t][k] -= PAR.WEprotclash * dSdR * x_p[k] /Rxp;
               }
             }
           }
         }
       }
     } 
   }
 }

 /*** [4] Calculate force_atom[][] for Eprotatrct **/ 
 if (PAR.WEprotatrct > 0.0){
   for (t=0;t<molT->Natom;++t){
     if (molT->atoms[t].one_char_ele != 'H'){
       for (p=0;p<molP->Natom;++p){
         if (molP->atoms[p].mark==1){
           if (molP->atoms[p].one_char_ele != 'H'){
             sub_vec3(x_p,molT->atoms[t].Pos,molP->atoms[p].Pos); 
             Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
             if (Rxp<DDthre_clash){
               Rxp = sqrt(Rxp);
               Rxprad = Rxp - molT->atoms[t].radius - molP->atoms[p].radius;
               dSdR = 0.0;
               expR = exp(-PAR.alpha_sigmoid_clash*(Rxprad));
               dSdR += - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
               expR = exp(-PAR.alpha_sigmoid_clash*(Rxprad - PAR.tole_dis_atrct));
               dSdR +=   PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
               for (k=0;k<3;++k){
                 tra->force_atom[t][k] -= PAR.WEprotatrct * dSdR * x_p[k] /Rxp;
               }
             }
           }
         }
       }
     } 
   }
 }
/*
 if (PAR.WEprotatrct > 0.0){
   for (t=0;t<molT->Natom;++t){
     for (p=0;p<molP->Natom;++p){
       if (molP->atoms[p].mark==1){
         sub_vec3(x_p,molT->atoms[t].Pos,molP->atoms[p].Pos); 
         Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
         if (Rxp<DDthre_clash){
           Rxp = sqrt(Rxp);
           Rxprad = Rxp - molT->atoms[t].radius - molP->atoms[p].radius;
           dSdR = 0.0;
           expR = exp(-PAR.alpha_sigmoid_clash*(Rxprad + PAR.tole_dis_clash));
           dSdR += - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
           expR = exp(-PAR.alpha_sigmoid_clash*(Rxprad - PAR.tole_dis_clash));
           dSdR +=   PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
           for (k=0;k<3;++k){
             tra->force_atom[t][k] -= PAR.WEprotatrct * dSdR * x_p[k] /Rxp;
           }
         }
       }
     } 
   }
 }
*/


 
 /*** [5] Calculate force_atom[][] for Evolmovlap **/ 
 if (PAR.WEvolmovlap>0.0){
   c1vo = 4.0/3.0 * M_PI * PAR.p_volmovlap;
   c2vo = M_PI * pow(3.0*PAR.p_volmovlap/4.0/M_PI,2.0/3.0);
   for (t=0;t<molT->Natom;++t){
     if (molT->atoms[t].one_char_ele != 'H'){
       rrT = molT->atoms[t].radius * molT->atoms[t].radius;
       for (r=0;r<molR->Natom;++r){
         if (molR->atoms[r].one_char_ele != 'H'){
           if (PAR.atomtype_volmovlap!='C') Watomtype=1.0;
           else Watomtype = Score_bwn_atomtypes(molT->atoms[t].atomtype,molR->atoms[r].atomtype);
           rrR = molR->atoms[r].radius * molR->atoms[r].radius;
           sub_vec3(t_r,molT->atoms[t].Pos,molR->atoms[r].Pos); 
           RRtr = t_r[0]*t_r[0] + t_r[1]*t_r[1] + t_r[2]*t_r[2]; 
           Cvo      = c1vo * pow(rrT*rrR/(rrT+rrR),3.0/2.0);
           alpha_vo = c2vo/(rrT+rrR);
           for (k=0;k<3;++k){
             tra->force_atom[t][k] -= PAR.WEvolmovlap*2.0*Watomtype*Cvo*alpha_vo*t_r[k]*exp(-alpha_vo*RRtr);
           }
         }
       }
     }
   }
 } 


 /*** [6] Calculate force_atom[][] for Etpdisrest **/ 
 if (PAR.WEtpdisrest>0.0){
   for (t=0;t<3;++t){
     if (PAR.string_tpdisrest[t][0] != '\0'){
       sub_vec3(x_p,molT->atoms[PAR.Tatom_num_tpdisrest[t]].Pos,molP->atoms[PAR.Patom_num_tpdisrest[t]].Pos); 
       Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
       if ((Rxp > PAR.Dupper_tpdisrest[t] * PAR.Dupper_tpdisrest[t]) && (Rxp>0.0)){
         Rxp = sqrt(Rxp);
         for (k=0;k<3;++k){ 
           tra->force_atom[PAR.Tatom_num_tpdisrest[t]][k] -=  PAR.WEtpdisrest * (Rxp - PAR.Dupper_tpdisrest[t])/Rxp * x_p[k];
         }
       }
     }
   }
 }


 check_numerical_force = 'F';

 if (check_numerical_force=='T'){ 
   printf("\n");
   delta = 0.001;
   for (t=0;t<molT->Natom;++t){
     for (k=0;k<3;++k){
       E0 = energy_total(molT,molR,molP,M);
       molT->atoms[t].Pos[k] += delta;
       E1 = energy_total(molT,molR,molP,M);
       molT->atoms[t].Pos[k] -= delta;
       force_nume[k] = -(E1-E0)/delta;
      }
     printf("ATOM %d %s force_atom %f %f %f numerical_grad %f %f %f\n",
         molT->atoms[t].num_in_file,molT->atoms[t].atomname,tra->force_atom[t][0], tra->force_atom[t][1], tra->force_atom[t][2],
         force_nume[0], force_nume[1], force_nume[2]);
   } 
 }


 
 /*** [7] Calulate tra->Force[] and tra->Torque[] ***/
 for (t=0;t<molT->Natom;++t){
   sub_vec3(x_c,molT->atoms[t].Pos, molT->atoms[tra->Canum].Pos);
   cross_vec3(tor, x_c,tra->force_atom[t]);
   for (k=0;k<3;++k){
     tra->Force[k]  += tra->force_atom[t][k];
     tra->Torque[k] += tor[k];
   }
 }
  
  normalize_vec3(tra->Force);
  normalize_vec3(tra->Torque);

  /*
  printf("#cal_Force_and_Torque() F %f %f %f T %f %f %f\n",
  tra->Force[0], tra->Force[1], tra->Force[2],
  tra->Torque[0], tra->Torque[1], tra->Torque[2]);
  */
} /* end of cal_Force_and_Torque() */





void cal_dEd_rbAngle(tra, molT, molR, molP, M)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *molT; /* target   */
  struct MOLECULE *molR; /* reference */
  struct MOLECULE *molP; /* receptor protein*/
  struct MATCH *M;
{
 int m,r,i,j,k,n;
 float x_r[3],o_z[3],x_o[3],dx_dphi[3],xi_xj[3],xi_o[3],x_p[3],norm;
 float E0,E1,dAng,Rij,Rxp,RRxr,Rxprad,Rijrad,dSdR,expR,DDthre_clash;
 char  check_numerical_dEdAngle;
 float rrT,rrR;
 float Cvo,alpha_vo;  /* Evo = -Cvo * exp(-alpha_vo*dd) */ 
 float c1vo;          /* constant :: 4/3*pi*p           */
 float c2vo;          /* constant :: pi*(3p/4pi)**(2/3) */
 float Watomtype;     /* Weight for atomtypes */
 
/*
 for (i=0;i<molT->Natom;++i){
   printf("cal_dE_dRbAngle molT %d %s %f\n",molT->atoms[i].num_in_file,
    molT->atoms[i].atomname,
    molT->atoms[i].radius);
 }

 for (i=0;i<molR->Natom;++i){
   printf("cal_dE_dRbAngle molR %d %s %f\n",molR->atoms[i].num_in_file,
    molR->atoms[i].atomname,
    molR->atoms[i].radius);
 }

 printf("#cal_dEd_rbAngle(tra, molT, molR, M)\n");
 */
  DDthre_clash = PAR.Dmax_clash*PAR.Dmax_clash; 

  for (r=0;r<tra->Nrbond;++r){ tra->dEd_rbAngle[r] = 0.0; }
 /* (1) Calculate rbAngle-derivative of Eatommatch */
 /*   dE/dphi = (xi-ri) dot (o_z x (xi-ok) */

 if (PAR.WEatommatch>0.0){
   for (k=0;k<tra->Nrbond;++k){
     sub_vec3(o_z, molT->atoms[tra->rbOanum[k]].Pos, molT->atoms[tra->rbZanum[k]].Pos);
     normalize_vec3(o_z);
     for (m=0;m<M->Npair;++m){
       i = M->anumA[m];
       j = M->anumB[m];
       if ((i>=0) && (j>=0) && (tra->rbManum[k][i]==1)){
         sub_vec3(x_r,molT->atoms[i].Pos, molR->atoms[j].Pos);
         sub_vec3(x_o,molT->atoms[i].Pos, molT->atoms[tra->rbOanum[k]].Pos);
         cross_vec3(dx_dphi,o_z,x_o);
         tra->dEd_rbAngle[k] += PAR.WEatommatch*dot_prod3(x_r,dx_dphi);
       }
     }
   }
 }
 /* (2) Calculate rbAngle-derivative of Eselfclash */
 /*   dE/dphi = dS/dR * (pi-qi)/R dot (o_z x (pi-ok)) */
 /*   dS/dR = -A*exp[-A(R-tau)]/(1+exp[-A(R-tau)])**2 */

 if (PAR.WEselfclash>0.0){
   for (r=0;r<tra->Nrbond;++r){
     sub_vec3(o_z, molT->atoms[tra->rbOanum[r]].Pos, molT->atoms[tra->rbZanum[r]].Pos);
     normalize_vec3(o_z);
     for (i=0;i<molT->Natom;++i){
       for (j=i+1;j<molT->Natom;++j){
         if ((molT->conmap.map[i][j]=='0')&&(molT->topodismap.map[i][j]>2)){ 
           if ( ((tra->rbManum[r][i]==1)&&(tra->rbManum[r][j]!=1))||
                ((tra->rbManum[r][i]!=1)&&(tra->rbManum[r][j]==1))   ) {
             sub_vec3(xi_xj, molT->atoms[i].Pos, molT->atoms[j].Pos);
             Rij = xi_xj[0]*xi_xj[0] + xi_xj[1]*xi_xj[1] + xi_xj[2]*xi_xj[2]; 
             if (Rij<DDthre_clash){
               Rij = sqrt(Rij); 
               Rijrad = Rij - molT->atoms[i].radius - molT->atoms[j].radius + PAR.tole_dis_clash;
               expR = exp(-PAR.alpha_sigmoid_clash*Rijrad);
               dSdR = - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
               sub_vec3(xi_o,molT->atoms[i].Pos, molT->atoms[tra->rbOanum[r]].Pos);
               cross_vec3(dx_dphi,o_z,xi_o);
/*
        printf("%f %f Rij %f Rijrad %f expR %f dSdR %f\n",molT->atoms[i].radius,molT->atoms[j].radius,Rij,Rijrad,expR,dSdR); 
*/

               if ((tra->rbManum[r][i]==1)&&(tra->rbManum[r][j]!=1)) 
                      tra->dEd_rbAngle[r] += PAR.WEselfclash * dSdR * dot_prod3(xi_xj,dx_dphi)/Rij;
                 else tra->dEd_rbAngle[r] -= PAR.WEselfclash * dSdR * dot_prod3(xi_xj,dx_dphi)/Rij;
             }
           } 
         }
       }
     }
   }
 }

 /* (3) Calculate rbAngle-derivative of Eprotclash */
 /*   dE/dphi = dS/dR * (pi-sj)/R dot (o_z x (pi-ok)) */
 /*   dS/dR = -A*exp[-A(R-tau)]/(1+exp[-A(R-tau)])**2 */
 if (PAR.WEprotclash>0.0){
   for (r=0;r<tra->Nrbond;++r){
     sub_vec3(o_z, molT->atoms[tra->rbOanum[r]].Pos, molT->atoms[tra->rbZanum[r]].Pos);
     normalize_vec3(o_z);
     for (i=0;i<molT->Natom;++i){
       if (tra->rbManum[r][i]==1){
         sub_vec3(xi_o,molT->atoms[i].Pos, molT->atoms[tra->rbOanum[r]].Pos);
         cross_vec3(dx_dphi,o_z,xi_o);

         for (j=0;j<molP->Natom;++j){
           if (molP->atoms[j].mark==1){
             sub_vec3(x_p, molT->atoms[i].Pos, molP->atoms[j].Pos);
             Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
             if (Rxp<DDthre_clash){ 
               Rxp = sqrt(Rxp);
               Rxprad = Rxp - molT->atoms[i].radius - molP->atoms[j].radius + PAR.tole_dis_clash;
               expR = exp(-PAR.alpha_sigmoid_clash*Rxprad);
               dSdR = - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));
               tra->dEd_rbAngle[r] += PAR.WEprotclash * dSdR * dot_prod3(x_p,dx_dphi)/Rxp;
             } 
            }
         }
       }
     }
   }
 }

 /* (4) Calculate rbAngle-derivative of Eprotatrct */
 /*   dE/dphi = dS/dR * (pi-sj)/R dot (o_z x (pi-ok)) */
 /*   dS/dR = -A*exp[-A(R+tau)]/(1+exp[-A(R+tau)])**2  + A*exp[-A(R-tau)]/(1+exp[-A(R-tau)])**2 */
 if (PAR.WEprotatrct>0.0){
   for (r=0;r<tra->Nrbond;++r){
     sub_vec3(o_z, molT->atoms[tra->rbOanum[r]].Pos, molT->atoms[tra->rbZanum[r]].Pos);
     normalize_vec3(o_z);
     for (i=0;i<molT->Natom;++i){
       if (tra->rbManum[r][i]==1){
         sub_vec3(xi_o,molT->atoms[i].Pos, molT->atoms[tra->rbOanum[r]].Pos);
         cross_vec3(dx_dphi,o_z,xi_o);

         for (j=0;j<molP->Natom;++j){
           if (molP->atoms[j].mark==1){
             sub_vec3(x_p, molT->atoms[i].Pos, molP->atoms[j].Pos);
             Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
             if (Rxp<DDthre_clash){ 
               Rxp = sqrt(Rxp);
               Rxprad = Rxp - molT->atoms[i].radius - molP->atoms[j].radius;
               dSdR = 0.0; 

               expR = exp(-PAR.alpha_sigmoid_clash*(Rxprad+PAR.tole_dis_clash));
               dSdR += - PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));

               expR = exp(-PAR.alpha_sigmoid_clash*(Rxprad-PAR.tole_dis_clash));
               dSdR +=   PAR.alpha_sigmoid_clash * expR /((1.0+expR)*(1.0+expR));

               tra->dEd_rbAngle[r] += PAR.WEprotatrct * dSdR * dot_prod3(x_p,dx_dphi)/Rxp;
             } 
            }
         }
       }
     }
   }
 }





 /* (5) Calculate rbAngle-derivative of Evolmovlap */
 /*   dE/dphi = 2Cexp[-alpha*(pi-qj)] (pi-qj) dot (o_z x (pi-ok) */

 if (PAR.WEvolmovlap>0.0){
   c1vo = 4.0/3.0 * M_PI * PAR.p_volmovlap;
   c2vo = M_PI * pow(3.0*PAR.p_volmovlap/4.0/M_PI,2.0/3.0);

   for (r=0;r<tra->Nrbond;++r){
     sub_vec3(o_z, molT->atoms[tra->rbOanum[r]].Pos, molT->atoms[tra->rbZanum[r]].Pos);
     normalize_vec3(o_z);
      n = 0;
     for (i=0;i<molT->Natom;++i){
       if ((molT->atoms[i].one_char_ele != 'H')&&(tra->rbManum[r][i]==1)){ 
         sub_vec3(x_o,molT->atoms[i].Pos, molT->atoms[tra->rbOanum[r]].Pos);
         rrT = molT->atoms[i].radius * molT->atoms[i].radius;
         for (j=0;j<molR->Natom;++j){
           if (molR->atoms[j].one_char_ele != 'H'){ 
             if (PAR.atomtype_volmovlap!='C') Watomtype=1.0;
             else Watomtype = Score_bwn_atomtypes(molT->atoms[i].atomtype,molR->atoms[j].atomtype);
             rrR = molR->atoms[j].radius * molR->atoms[j].radius;
             sub_vec3(x_r,molT->atoms[i].Pos, molR->atoms[j].Pos);
             cross_vec3(dx_dphi,o_z,x_o);
             RRxr = x_r[0]*x_r[0] + x_r[1]*x_r[1] + x_r[2]*x_r[2];
             Cvo      = c1vo * pow(rrT*rrR/(rrT+rrR),3.0/2.0);
             alpha_vo = c2vo/(rrT+rrR);
             tra->dEd_rbAngle[r] += PAR.WEvolmovlap*2.0*Watomtype*Cvo*alpha_vo*exp(-alpha_vo*RRxr)*dot_prod3(x_r,dx_dphi); 
             n += 1;
           }
         }
       }
     }
   }
 }


 /* (6) Calculate rbAngle-derivative of Etpdisrest */
 /*   dE/dphi = -(Dxp-Dupp)/Dxp (x-p) dot (o_z x (x-ok)) */

 if (PAR.WEtpdisrest>0.0){
/*
   for (r=0;r<tra->Nrbond;++r){
     for (i=0;i<molT->Natom;++i){
        printf("#BOND[%d] Z %d O %d rbManum[%d][%d]:%d\n",r,tra->rbZanum[r],tra->rbOanum[r],r,i,tra->rbManum[r][i]);
     }
   }
*/
   for (k=0;k<3;++k){   
     if (PAR.string_tpdisrest[k][0] != '\0'){
       sub_vec3(x_p, molT->atoms[PAR.Tatom_num_tpdisrest[k]].Pos, molP->atoms[PAR.Patom_num_tpdisrest[k]].Pos);
        Rxp = x_p[0]*x_p[0] + x_p[1]*x_p[1] + x_p[2]*x_p[2]; 
        /* printf("#Rxp %f Dupper %f\n",sqrt(Rxp),PAR.Dupper_tpdisrest); */
        if ((Rxp>PAR.Dupper_tpdisrest[k]) && (Rxp>0.0)){
        Rxp = sqrt(Rxp);
        for (r=0;r<tra->Nrbond;++r){
           /* printf("#bond[%d] rbManum[%d][%d]:%d\n",r,r,PAR.Tatom_num_tpdisrest,tra->rbManum[r][PAR.Tatom_num_tpdisrest]); */
           if (tra->rbManum[r][PAR.Tatom_num_tpdisrest[k]]==1){
             /* printf("#IN bond[%d]\n",r); */
             sub_vec3(o_z, molT->atoms[tra->rbOanum[r]].Pos, molT->atoms[tra->rbZanum[r]].Pos);
             normalize_vec3(o_z);
             sub_vec3(xi_o,molT->atoms[PAR.Tatom_num_tpdisrest[k]].Pos, molT->atoms[tra->rbOanum[r]].Pos);
             cross_vec3(dx_dphi,o_z,xi_o);
             tra->dEd_rbAngle[r] += PAR.WEtpdisrest*(Rxp - PAR.Dupper_tpdisrest[k])/Rxp * dot_prod3(x_p,dx_dphi);
           }
         }
       }
     }
   } 
 }

 check_numerical_dEdAngle = 'F';

 if (check_numerical_dEdAngle=='T'){ 
   printf("\n");
   /* dAng = 0.01; */
   dAng = 0.001; 
   for (r=0;r<tra->Nrbond;++r){
     E0 = energy_total(molT,molR,molP,M);
     Transform_MOLECULE_by_RotBond(r,dAng,tra,molT);
     E1 = energy_total(molT,molR,molP,M);
     n = Transform_MOLECULE_by_RotBond(r,-dAng,tra,molT);
     printf("RBOND[%2d] CenAtom %d Z %2d O %2d dEdAngle %+13.10f (E1-E0)/dang %+13.10f n %d\n",r,
       molT->atoms[tra->Canum].num_in_file,
       molT->atoms[tra->rbZanum[r]].num_in_file,molT->atoms[tra->rbOanum[r]].num_in_file,tra->dEd_rbAngle[r],(E1-E0)/dAng,n);
   }
 }


 /** (7) Normalize ||(dEdAngle[0],dEdAngle[1],...,dEdangle[K-1])|| = 1; **/
 norm = 0.0; 
 for (r=0;r<tra->Nrbond;++r){
   norm += tra->dEd_rbAngle[r]*tra->dEd_rbAngle[r];
 }
 if (norm>0.0){
   norm = sqrt(norm/tra->Nrbond);
   for (r=0;r<tra->Nrbond;++r){
      tra->dEd_rbAngle[r] = tra->dEd_rbAngle[r]/norm;
   }
 }

} /* end of cal_dEd_rbAngle() */



















void Random_Transform_MOLECULE_by_RotBond(tra,molA,maxRbAngle,FixStamp,FixPlane)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *molA;/* molecule to be randomly transformed */
  float  maxRbAngle;    /* maximum deformed angle (degree) */
  char   FixStamp;      /* 'F':fixing, 'R'otating stamped angles */
  char   FixPlane;      /* 'F':fixing, 'R'otating plane angles. 'I':allow inversing (180 degree rotation) */
{
  int r;
  float rand_ang; /* radian */
  char rot;

  printf("#Random_Transform_MOLECULE_by_RotBond(tra,molA,maxRbAngle:%f FixStamp '%c' FixPlane '%c')\n",maxRbAngle,FixStamp,FixPlane);

  for (r=0;r<tra->Nrbond;++r){
    rot = 0;
    if ( ((FixStamp!='F')||(tra->rbStamp[r]!='S'))  &&
         ((FixPlane!='F')||(tra->rbPlane[r]!='P'))) { rot = 1;}

/*
    printf("#rot_bond %d rbStamp %c rbPlane %c atomO->Z %d->%d rot %d\n",r,
      tra->rbStamp[r],tra->rbPlane[r], molA->atoms[tra->rbOanum[r]].num_in_file,molA->atoms[tra->rbZanum[r]].num_in_file,rot);
 */
    if (rot==1){
      rand_ang = (-1.0 + 2.0*(float)rand()/(float)RAND_MAX) * M_PI * maxRbAngle /180.0;
      if ((FixPlane=='I')&&(tra->rbPlane[r]=='P')){
        if ((rand()%2)==0){rand_ang = 0.0;} else {rand_ang = M_PI;}
        /* printf("#FixPlane FLIP %f\n",rand_ang); */
      }
      /* printf("#rand_ang %f\n",rand_ang); */
      printf("#rot_bond %d rbStamp %c rbPlane %c atomO->Z %d->%d rot %d rand_ang %.2f (degree)\n",r,
        tra->rbStamp[r],tra->rbPlane[r], molA->atoms[tra->rbOanum[r]].num_in_file,molA->atoms[tra->rbZanum[r]].num_in_file,rot,rand_ang*180/M_PI);
      Transform_MOLECULE_by_RotBond(r,rand_ang,tra,molA);
    }
  }
} /* end of Random_Transform_MOLECULE_by_RotBond() */


