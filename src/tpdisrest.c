/*

< tpdisrest.c >

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


functions for tpdisrest (target-protein distance restraint)


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
#include "transform.h"
#include "vector3.h"
#include "energies.h"
#include "stamp_transf.h"
#include "PCAfit.h"
#include "qRMS.h"
#include "io_match.h"

/*** FUNCTIONS (GLOBAL) ***/
void setup_variables_from_string_tpdisrest();
void mk_string_status_tpdisrest();
void mk_short_string_status_tpdisrest();
void Transform_MOLECULE_by_RandomChosenAtoms_Translation_Random_Rotation();
void Transform_MOLECULE_by_Random_3PointPairs_with_TPDISREST_restraint();


/*** FUNCTIONS (LOCAL) ***/
static void random_select_target_atom();
static void random_select_reference_atom_among_tpdisrest_satisfied_atoms();
static void make_XYZ_from_Po_Px_Py();


void setup_variables_from_string_tpdisrest(molT,molP)
  struct MOLECULE *molT;
  struct MOLECULE *molP;
{
  int t,Nword,Wsta[100],Wend[100];
  char buff[100],proper_rest;
  struct ATOM *atmT,*atmP;

/*

 FORMAT OF STRING_TPDISREST: 
   string_tpdisrest[] = (Tarom_num):(Patom_num):(Dupper):('S'elect)

*/
  printf("#setup_variables_from_string_tpdisrest(molT,molP)\n");

  proper_rest = 1;
  if ((PAR.string_tpdisrest[0][0]== '\0') && ((PAR.string_tpdisrest[1][0] != '\0') || (PAR.string_tpdisrest[2][0]!= '\0'))){
    proper_rest = 0;
  }

  if ((PAR.string_tpdisrest[0][0] != '\0') && (PAR.string_tpdisrest[1][0] == '\0') && (PAR.string_tpdisrest[2][0]!= '\0')){
    proper_rest = 0;
  }
 
  if (proper_rest==0){
    printf("#ERROR:Improper tpdisrest.");
    for (t=0;t<3;++t){ printf(" -tppair%d '%s' ",t,PAR.string_tpdisrest[t]); }
    printf("\n");
    exit(1);
  }

  for (t=0;t<3;++t){
     if (PAR.string_tpdisrest[t][0] != '\0'){
       Split_to_Words(PAR.string_tpdisrest[t],':',&Nword,Wsta,Wend,100);
       Get_Part_Of_Line(buff,PAR.string_tpdisrest[t],Wsta[0],Wend[0]);  
       PAR.Tatom_num_tpdisrest[t] = atom_num_from_num_in_file(molT,atoi(buff));
    
       Get_Part_Of_Line(buff,PAR.string_tpdisrest[t],Wsta[1],Wend[1]);  
       PAR.Patom_num_tpdisrest[t] = atom_num_from_num_in_file(molP,atoi(buff));
 
       Get_Part_Of_Line(buff,PAR.string_tpdisrest[t],Wsta[2],Wend[2]);  
       PAR.Dupper_tpdisrest[t] = atof(buff);
 
      if ((PAR.Tatom_num_tpdisrest[t]<0) || (PAR.Patom_num_tpdisrest[t]<0)){
        printf("#ERROR:string_tpdisrest[%d] '%s' => Tatom_num:%d Patom_num:%d Dupper:%f\n",t,PAR.string_tpdisrest[t],PAR.Tatom_num_tpdisrest[t],PAR.Patom_num_tpdisrest[t],PAR.Dupper_tpdisrest[t]);
        exit(1);
       }

      atmT = &(molT->atoms[PAR.Tatom_num_tpdisrest[t]]);
      atmP = &(molP->atoms[PAR.Patom_num_tpdisrest[t]]);
      
      printf("#TPDISREST[%d] atomT %5d %4s atomP %5d %4s %3s %3s Dupper %8.3f\n",t,
        atmT->num_in_file, atmT->atomname, atmP->num_in_file, atmP->atomname, atmP->resi, atmP->rnum, PAR.Dupper_tpdisrest[t]);
     }
  }
} /* end of setup_variables_from_string_tpdisrest() */


void mk_string_status_tpdisrest(retstr,molT,molP)
  char   *retstr;
  struct MOLECULE *molT;
  struct MOLECULE *molP;
{
  int t,i;
  char cal_rest,status[32],buff[1024];
  float d,dis,ene;
  struct ATOM *atmT,*atmP;

  retstr[0] = '\0';
  for (t=0;t<3;++t){
    cal_rest = 1;
    status[0] = '\0'; 
    if (PAR.string_tpdisrest[t][0]=='\0'){cal_rest = 0;}
    if ((PAR.Tatom_num_tpdisrest[t]<0) || (PAR.Tatom_num_tpdisrest[t]>=molT->Natom)){cal_rest = 0;}
    if ((PAR.Patom_num_tpdisrest[t]<0) || (PAR.Patom_num_tpdisrest[t]>=molP->Natom)){cal_rest = 0;}

    if (cal_rest == 1){
      dis = 0.0;
      atmT = &(molT->atoms[PAR.Tatom_num_tpdisrest[t]]);
      atmP = &(molP->atoms[PAR.Patom_num_tpdisrest[t]]);
      for (i=0;i<3;++i){
         d = atmT->Pos[i] - atmP->Pos[i];
        dis += d*d;
      }
      dis = sqrt(dis);
      ene = 0.0;
 /*
      printf("#atmT %f %f %f atmP %f %f %f dis %f\n",
  atmT->Pos[0], atmT->Pos[1], atmT->Pos[2],
  atmP->Pos[0], atmP->Pos[1], atmP->Pos[2],dis);
 */
      if (dis <= PAR.Dupper_tpdisrest[t]){ sprintf(status,"SATISFIED"); }
      else {
        sprintf(status,"NOT-SATISFIED");
        ene = 0.5*(dis - PAR.Dupper_tpdisrest[t])*(dis - PAR.Dupper_tpdisrest[t]);
      }
     /*
      printf("#TPDISREST[%d] atomT %d %s atomP %d %s %s %s Dupper %f dis %f ene %f %s\n",t,
        atmT->num_in_file, atmT->atomname, atmP->num_in_file, 
        atmP->atomname, atmP->resi, atmP->rnum, PAR.Dupper_tpdisrest[t], dis,ene,status);
      */
      sprintf(buff,"#TPDISREST[%d] atomT %d %s atomP %d %s %s %s Dupper %f dis %f ene %f %s\n",t,
        atmT->num_in_file, atmT->atomname, atmP->num_in_file, 
        atmP->atomname, atmP->resi, atmP->rnum, PAR.Dupper_tpdisrest[t], dis,ene,status);
      strcat(retstr,buff);
   }
 }

} /* end of mk_string_status_tpdisrest() */




void mk_short_string_status_tpdisrest(retstr,t,molT,molP)
  char   *retstr;  /* return string */ 
  int    t;  /* restraint number (0,1,2) */
  struct MOLECULE *molT;
  struct MOLECULE *molP;
{
  int i;
  char cal_rest,status;
  float d,dis;
  struct ATOM *atmT,*atmP;

  retstr[0] = '\0';
  cal_rest = 1;
  status = '-'; 
  if (PAR.string_tpdisrest[t][0]=='\0'){cal_rest = 0;}
  if ((PAR.Tatom_num_tpdisrest[t]<0) || (PAR.Tatom_num_tpdisrest[t]>=molT->Natom)){cal_rest = 0;}
  if ((PAR.Patom_num_tpdisrest[t]<0) || (PAR.Patom_num_tpdisrest[t]>=molP->Natom)){cal_rest = 0;}

  if (cal_rest == 1){
    dis = 0.0;
    atmT = &(molT->atoms[PAR.Tatom_num_tpdisrest[t]]);
    atmP = &(molP->atoms[PAR.Patom_num_tpdisrest[t]]);
    for (i=0;i<3;++i){
       d = atmT->Pos[i] - atmP->Pos[i];
      dis += d*d;
    }
    dis = sqrt(dis);
    if (dis <= PAR.Dupper_tpdisrest[t]){ status = 'S'; }
    else {
      status = 'B';
      /* ene = 0.5*(dis - PAR.Dupper_tpdisrest[t])*(dis - PAR.Dupper_tpdisrest[t]); */
    }
    sprintf(retstr,"tp%d_%c_%.2f",t,status,dis);
  }

} /* end of mk_short_string_status_tpdisrest() */



void Transform_MOLECULE_by_Random_3PointPairs_with_TPDISREST_restraint(tra,molT,molR,molP)
  struct TRANSFORM_VARIABLE *tra;
  struct MOLECULE *molT;/* target molecule to be randomly transformed */
  struct MOLECULE *molR;/* reference molecule for the template */
  struct MOLECULE *molP;/* protein molecule for the receptor */
{
  float dd,min_dd,Dt01,Dt02,Dt12,Dr01,Dr02,Dr12;
  int r,t,min_aR;
  int anumT[3],anumR[3];
  float Xt[3],Yt[3],Zt[3],Xr[3],Yr[3],Zr[3],Rmat[3][3],Ot[3],Or[3];

  printf("#Transform_MOLECULE_by_Random_3PointPairs_with_TPDISREST_restraint(restraint '%s' '%s' '%s')\n",
       PAR.string_tpdisrest[0], PAR.string_tpdisrest[1], PAR.string_tpdisrest[2]);
  /** [0-th pair] decide (anumT[0], anumR[0])  **/
  if (PAR.string_tpdisrest[0][0] != '\0'){
    anumT[0] = PAR.Tatom_num_tpdisrest[0];
    random_select_reference_atom_among_tpdisrest_satisfied_atoms(anumR,0,molR, molP);
  } 
  else{
    anumT[0] = rand()%molT->Natom;
    anumR[0] = rand()%molR->Natom;
  }

  /** [1st pair] decide (anumT[1], anumR[1])  **/
  if (PAR.string_tpdisrest[1][0] != '\0'){
    anumT[1] = PAR.Tatom_num_tpdisrest[1];
    random_select_reference_atom_among_tpdisrest_satisfied_atoms(anumR,1,molR, molP);
  } 
  else{ 
    random_select_target_atom(anumT,1,molT);
    Dt01 = distance_vec3(molT->atoms[anumT[0]].Pos,molT->atoms[anumT[1]].Pos); 

    min_dd = -1.0;  
    min_aR = 0;
    for (r=0;r<molR->Natom;++r){
        if (r != anumR[0]){
      
        Dr01 = distance_vec3(molR->atoms[anumR[0]].Pos,molR->atoms[r].Pos);
        dd = (Dr01-Dt01)*(Dr01-Dt01);
        if ((min_dd<0.0) || (dd<min_dd)){
          min_dd = dd;
          min_aR = r;
        } 
      }
    }
    anumR[1] = min_aR;
  }

  /** [2nd pair] decide (anumT[2], anumR[2])  **/
  if (PAR.string_tpdisrest[2][0] != '\0'){
    anumT[2] = PAR.Tatom_num_tpdisrest[2];
    random_select_reference_atom_among_tpdisrest_satisfied_atoms(anumR,2,molR, molP);
  }
  else{ 
    random_select_target_atom(anumT,2,molT);
    Dt02 = distance_vec3(molT->atoms[anumT[0]].Pos,molT->atoms[anumT[2]].Pos); 
    Dt12 = distance_vec3(molT->atoms[anumT[1]].Pos,molT->atoms[anumT[2]].Pos); 
    min_dd = -1.0;  
    min_aR = 0;
    for (r=0;r<molR->Natom;++r){
      if ((r != anumR[0])&&(r != anumR[1])){
        Dr02 = distance_vec3(molR->atoms[anumR[0]].Pos,molR->atoms[r].Pos);
        Dr12 = distance_vec3(molR->atoms[anumR[1]].Pos,molR->atoms[r].Pos);
        dd = (Dr02-Dt02)*(Dr02-Dt02) + (Dr12-Dt12)*(Dr12-Dt12);
        if ((min_dd<0.0) || (dd<min_dd)){
          min_dd = dd;
          min_aR = r;
        } 
      }
    }
    anumR[2] = min_aR;
  }

  Dt01 = distance_vec3(molT->atoms[anumT[0]].Pos,molT->atoms[anumT[1]].Pos); 
  Dt02 = distance_vec3(molT->atoms[anumT[0]].Pos,molT->atoms[anumT[2]].Pos); 
  Dt12 = distance_vec3(molT->atoms[anumT[1]].Pos,molT->atoms[anumT[2]].Pos); 
  Dr01 = distance_vec3(molR->atoms[anumR[0]].Pos,molR->atoms[anumR[1]].Pos); 
  Dr02 = distance_vec3(molR->atoms[anumR[0]].Pos,molR->atoms[anumR[2]].Pos); 
  Dr12 = distance_vec3(molR->atoms[anumR[1]].Pos,molR->atoms[anumR[2]].Pos); 

  printf("#T(%d %d %d) (%f %f %f) R(%d %d %d) (%f %f %f)\n",anumT[0],anumT[1],anumT[2],Dt01,Dt02,Dt12,anumR[0],anumR[1],anumR[2],Dr01,Dr02,Dr12);
 
  /* (6) make local coodinate (Xt,Yt,Zt) from anumT[0,1,2], and make (Xr,Yr,Zr) from anumR[0,1,2].  */
 
  make_XYZ_from_Po_Px_Py(Xt,Yt,Zt,molT->atoms[anumT[0]].Pos, molT->atoms[anumT[1]].Pos, molT->atoms[anumT[2]].Pos);
  make_XYZ_from_Po_Px_Py(Xr,Yr,Zr,molR->atoms[anumR[0]].Pos, molR->atoms[anumR[1]].Pos, molR->atoms[anumR[2]].Pos);


  /** (7) Make Rmat, transforming (Xt,Yt,Zt) into (Xr,Yr,Zr). */
  Cal_Rotation_Matrix_from_Pair_of_XYZvec(Rmat,Xt,Yt,Zt,Xr,Yr,Zr);
  /* 
  make_transpose_matrix3x3(TRmat,Rmat);
  multi_mat3x3(Amat,Rmat,TRmat);
  show_matrix3x3(Rmat,"Rmat"); 
  show_matrix3x3(TRmat,"TRmat"); 
  show_matrix3x3(Amat,"Rmat*TRmat"); 
  */
  
  /** (8) Transoform all the atoms in molT */
  copy_vec3(Ot,molT->atoms[anumT[0]].Pos);
  copy_vec3(Or,molR->atoms[anumR[0]].Pos);
  
  for (t=0;t<molT->Natom;++t){
    transform_vec_by_Rmat_oldO_newO(molT->atoms[t].Pos,Rmat,Ot,Or);
  }


} /* end of Transform_MOLECULE_by_Random_3PointPairs_with_TPDISREST_restraint() */



void random_select_reference_atom_among_tpdisrest_satisfied_atoms(anumR,rest_num,molR, molP)
  int    anumR[3];
  int    rest_num;  /* restraint_number (0,1,2) */
  struct MOLECULE *molR; /* reference molecule (pocket) */
  struct MOLECULE *molP; /* protein receptor molecule */
{
  float r_p[3];
  float dd,min_dd;
  int r,min_aR,Nratom_withinDupper,*ratom_withinDupper;

  ratom_withinDupper = (int *)malloc(sizeof(int)*molR->Natom);

  for (r=0;r<molR->Natom;++r){ratom_withinDupper[r] = -1;}

  Nratom_withinDupper = 0;
  min_dd = -1.0;
  min_aR = -1;

  for (r=0;r<molR->Natom;++r){
    if (  (rest_num==0) 
       || ((rest_num==1) && (r != anumR[0]))
       || ((rest_num==2) && (r != anumR[0]) && (r != anumR[1])) ){

      dd = 0.0;
      sub_vec3(r_p, molR->atoms[r].Pos, molP->atoms[PAR.Patom_num_tpdisrest[rest_num]].Pos);
      dd = r_p[0]*r_p[0] + r_p[1]*r_p[1] + r_p[2]*r_p[2];
      if (dd <= PAR.Dupper_tpdisrest[rest_num] * PAR.Dupper_tpdisrest[rest_num]){
        ratom_withinDupper[Nratom_withinDupper] = r;
        Nratom_withinDupper += 1;
      }
      if ((min_aR<0) || (dd < min_dd)){
        min_aR  = r;
        min_dd = dd;
      }
    }
  }

  if (Nratom_withinDupper > 0){
    anumR[rest_num] = ratom_withinDupper[rand()%Nratom_withinDupper];
    printf("#aR %d is chosen to satisfy tpdisrest among the Nratom_withinDupper %d.\n",anumR[0],Nratom_withinDupper);
  }
  else{
    anumR[rest_num] = min_aR;
    printf("#aR %d is chosen, but not satisfy tpdisrest. Dr_p:%f\n",anumR[0],sqrt(min_dd));
  }

  free(ratom_withinDupper);

} /* end of random_select_reference_atom_among_tpdisrest_satisfied_atoms() */



void random_select_target_atom(anumT,rest_num,molT)
  int    anumT[3];  /* atom numbers of molT  */
  int    rest_num;  /* restraint_number (0,1,2) */
  struct MOLECULE *molT;
{
  int i,j,rnum;
  rnum = rand()%(molT->Natom-rest_num);

  anumT[rest_num] = -1;
  i = 0; j = 0;
  while ((anumT[rest_num]==-1)&&(i < molT->Natom) && ( j < (molT->Natom-rest_num))){
    if ((rest_num==0) ||
        ((rest_num==1)&&(i != anumT[0])) ||
        ((rest_num==2)&&(i != anumT[0])&&(i != anumT[1]))){
      if (rnum==j){ anumT[rest_num] = i; }
      ++j; 
    }
    ++i; 
  }

} /* end of void random_select_target_atom() */



void make_XYZ_from_Po_Px_Py(X,Y,Z,Po,Px,Py)
  float X[3],Y[3],Z[3];   /* The three axes for XYZ coodinates */
  float Po[3];            /* Origin */
  float Px[3];            /* X[] ==>  Px-Po */
  float Py[3];            /* Y[] ==> (Py-Po) - [X dot (Py-Po)]X */
{
  float yo[3],dprod_yo_X;
  int k;

  sub_vec3(X,Px,Po);
  normalize_vec3(X);
  sub_vec3(yo,Py,Po);
  dprod_yo_X = dot_prod3(X,yo);
  for (k=0;k<3;++k){
    Y[k] = yo[k] - dprod_yo_X * X[k];
  }  
  normalize_vec3(Y);
  cross_vec3(Z,X,Y);

} /* end of make_XYZ_from_Po_Px,Py() */

