/*
 
  <gen2D.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  for generating two-dimentional (x,y) coordinate
  only from connection table.
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "molprop.h"
#include "ioLINE.h"
#include "ioSDF.h"
#include "moltopodis.h"
#include "PCAfit.h"
#include "gen2D.h"
#include "disgeo2D.h"
#include "vector3.h"

/** LOCAL STRUCTURE for UNIT OF ATOMS **/

/*** FUNCTIONS (GLOBAL) ***/
int  Generate_2D_Coordinates();
int  set_hydrogen_XY();
void Write_SDF_UNIT();
void Free_UNITs();
void Rescaling_Corrdinates_of_Molecule();

/*** FUNCTIONS (LOCAL) ***/
static int Assemble_UNITs_to_make_global_XYcordinates();
static int combine_new_UNIT_into_seed_atoms();
static void Decompose_Atoms_into_UNITs();
static int set_XY_to_ring_UNIT();
static void set_localXY_in_regular_polygon();
static int combine_new_ring_into_seed_ring();
static void make_ring_to_unit_index();
static int find_shared_atoms_bwn_two_rings();
static int set_XY_to_chain_UNIT();
static int mark_neighbors_for_smallest_cycle_in_UNIT();
static void add_new_UNIT_to_tail();
static void setup_neighbors_for_UNITs();
static void show_UNITs();
static void show_UNIT();
static int Number_of_UNITs();
static void normalize_2Dvec();
static void congestion_score();
static float point_congestion_score();
static float congestion_in_ring_for_dir();
static float congestion_in_marked_atoms_for_dir();
static void cal_mid_point();
static void cal_angular_mid_point();
static void cal_angular_three_division_point();
static int  Collision_Resolving();
static void copy_posXY_from_MOLECULE();
static void copy_posXY_to_MOLECULE();
static int  inversion();
static int check_bond_distance();
static float myacos();


int Generate_2D_Coordinates(mol)
  struct MOLECULE *mol;
{
  struct UNIT HeadUnit;
  struct UNIT *un;
  int i, Natmpair_clash;
  float score;
  char fname[256];
  struct RING *rn;

  printf("#Generate_2D_Coordinates('%s' Natom %d)\n",mol->name,mol->Natom);
  for (i=0;i<mol->Natom;++i){
    mol->atoms[i].Pos[0] = mol->atoms[i].Pos[1] = mol->atoms[i].Pos[2] = 0.0; 
    mol->atoms[i].ring = ' ';  
 }

  if (mol->Nheavyatom==1){
    set_hydrogen_XY(mol);
    return(1);
  }

  /** [0] Set mol->atoms[].ring : 'R':in large (>6)ring, 'r':in small (<=6) ring **/  
  rn = &(mol->HeadRing);
  while (rn->next != NULL){
    rn = rn->next;
    /* printf("#RING[%d] Natom %d\n",rn->num,rn->Natom);  */
    for (i=0;i<rn->Natom;++i){
      if (rn->Natom<=6){ mol->atoms[rn->num_atoms[i]].ring = 'r';}
      else if (mol->atoms[rn->num_atoms[i]].ring == '0'){ mol->atoms[rn->num_atoms[i]].ring = 'R';}
    }
  }
 /*
  for (i=0;i<mol->Natom;++i){
    printf("#%s%d ring '%c'\n",mol->atoms[i].atomname,mol->atoms[i].num_in_file,mol->atoms[i].ring);
  }
 */
  /** [1] Decomposing atoms into UNITs **/
  Decompose_Atoms_into_UNITs(mol,&HeadUnit);
  /*
  show_UNITs(mol,&HeadUnit,"list of UNIT");
  */
  un = &HeadUnit; 

  while (un->next != NULL){
    un = un->next;
    if (un->type=='R') set_XY_to_ring_UNIT(mol,un);
    if (un->type=='C') set_XY_to_chain_UNIT(mol,un);
    sprintf(fname,"%d%c.sdf",un->num,un->type);
    /* Write_SDF_UNIT(fname,un,mol); */
  }

  /*
  show_UNITs(mol,&HeadUnit,"list of UNIT");
  */

  /** [2] Assembel UNITs  **/
  Assemble_UNITs_to_make_global_XYcordinates(mol,&HeadUnit);
  
  /* Write_SDF_Molecule("before_resolv.sdf",mol,""); */
  
  /** [3]  Collision Resolving  **/
  congestion_score(mol,&score,&Natmpair_clash);
  printf("#before_resolve %s Natmpair_clash %d score %f\n",mol->filename,Natmpair_clash,score);
  if (PAR.CollisionResolve2D == 'T') { 
    Collision_Resolving(mol,&HeadUnit);
    congestion_score(mol,&score,&Natmpair_clash);
  }
  set_hydrogen_XY(mol); 
  printf("#final_score %s Natmpair_clash %d score %f\n",mol->filename,Natmpair_clash,score);

  Rotate_Molecule_XY_flat_plane_by_PCA(mol);
  Free_UNITs(&HeadUnit);

  return(mol->Nheavyatom);

} /* end of Generate_2D_coordinates() */



int Assemble_UNITs_to_make_global_XYcordinates(mol,HeadUnit)
  struct MOLECULE *mol;
  struct UNIT *HeadUnit;
{
  struct UNIT *xn,*cn,*on, *sn, *an; 
  int maxNatom_chain,maxNatom_other,Nunit,Nunit_rest;
  int i,j,ii,match,ishare_shareA,anum_share; 
  char ofname[128],combok; 
  float score;
  int Natmpair_clash; 
  printf("#Assemble_UNITs_to_make_global_XYcordinates(mol,HeadUnit)\n");

  ofname[0] = '\0';
  for (i=0;i<mol->Natom;++i){ mol->atoms[i].mark = 0;}

  /** [1] find seed UNITs ***/
  /*** (1st priority) : longest chain unit (Natom>=3) **/
  /*** (2nd priority) : largest unit of all type */

  Nunit = Number_of_UNITs(HeadUnit);
  printf("#Nunit %d\n", Nunit);

  xn = HeadUnit;
  cn = on = NULL;
  maxNatom_chain = maxNatom_other = -1;

  while (xn->next != NULL){
    xn = xn->next;
 /*
    printf("#unit[%d] type %c Natom %d\n",xn->num,xn->type,xn->Natom);
  */
    if ((xn->type=='C')&&(xn->Natom>=3)){
      if (xn->Natom > maxNatom_chain){
         maxNatom_chain = xn->Natom;
         cn = xn;
      }
    }
    else{
      if (xn->Natom > maxNatom_other){
         maxNatom_other = xn->Natom;
         on = xn;
      }
    }

  }  

       if (cn!=NULL){ sn = cn;}
  else if (on != NULL) {sn = on;}
  else {printf("#ERROR(Assemble_UNIT('%s')):no proper seed UNIT.\n",mol->name); return(0);}

  sn->set_globalXY = 'T';
  /*
  show_UNIT(mol,sn,"seed unit ");
  */

  for (i=0;i<sn->Natom;++i){
    ii = sn->num_atoms[i];
    mol->atoms[ii].Pos[0] = sn->pos[i][0];
    mol->atoms[ii].Pos[1] = sn->pos[i][1];
    mol->atoms[ii].Pos[2] = 0.0;
    mol->atoms[ii].mark = 1;
  }

  Nunit_rest = Nunit -1;
  sprintf(ofname,"tmp/%d.sdf",Nunit_rest);
  /* Write_SDF_Molecule(ofname,mol,"");  */
 
  /*** [2] combine other units into seed, one by one ***/
  while (Nunit_rest>0){ 
    match = 0; ishare_shareA = anum_share = -1;
    an = HeadUnit;
    while ((an->next!=NULL)&&(match==0)){
      an = an->next;
      if (an->set_globalXY != 'T'){
        for (i=0;i<an->Natom;++i){
          for (j=0;j<mol->Natom;++j){
            if ((mol->atoms[j].mark==1) && (an->num_atoms[i]==j)){
              ishare_shareA = i; anum_share = j; match = 1;
            }
          }
       }
      }
    }
 
    if (match==1){
    /*
      printf("#>Nunit_rest %d shared atom %s%d\n",Nunit_rest,mol->atoms[anum_share].atomname, mol->atoms[anum_share].num_in_file);
      show_UNIT(mol,an,"add unit ");
  */
      combok = combine_new_UNIT_into_seed_atoms(mol,an,ishare_shareA,anum_share);
      if (combok==0){ return(0);}
      congestion_score(mol,&score,&Natmpair_clash);
      an->set_globalXY = 'T';
      Nunit_rest -= 1;
      sprintf(ofname,"tmp/%d.sdf",Nunit_rest);
      /* Write_SDF_Molecule(ofname,mol,""); */  
    }
  }

 return(1);

} /* end of Assemble_UNITs_to_make_global_XYcordinates() */



int combine_new_UNIT_into_seed_atoms(mol,addU,ishare_a,anum_share)
  struct MOLECULE *mol;
  struct UNIT *addU;
  int    ishare_a;    /* index of the shared atom for addU */
  int    anum_share;  /* atom number of shared atom    */
{
  float Os[2],Oa[2]; /* shared atom XY coordinates */
  float Bs[2],Ba[2]; /* shared bond vector coordinates */
  float dpos[2];
  float max_angle,angle,bdir[2],bdirs[4][2]; 
  int i,ii,j, ina[3];
  int Nneighbor_share, anum_nei_seed[10]; 
  float th,c_t,s_t;
  char typeS[8];  /* type for seed.   ('T'erminal or 'M'id-point) */
  char typeA[8];  /* type adding unit.('T'erminal or 'M'id-point) */
  float  ConSc,minConSc;
  char seedRing; 
  /*
  printf("#combine_new_UNIT_into_seed_atoms(addU:%d %c ishare_a:%d anum_share %d)\n",
      addU->num,addU->type,ishare_a,anum_share);
  */
  max_angle = -1.0;
  /*********************************/
  /** (1) Setup Oa and Ba for addU */
  /*********************************/
  
  Oa[0] = addU->pos[ishare_a][0];
  Oa[1] = addU->pos[ishare_a][1];
  ina[0] = addU->neighbors[ishare_a][0];
  ina[1] = addU->neighbors[ishare_a][1];
  ina[2] = addU->neighbors[ishare_a][2];
 
  /* 
  printf("addU share %d neighbor %d", mol->atoms[addU->num_atoms[ishare_a]].num_in_file, mol->atoms[addU->num_atoms[ina0]].num_in_file);
  if (addU->Nneighbor[ishare_a]>=2) printf(" %d",mol->atoms[addU->num_atoms[ina1]].num_in_file);
  printf("\n");
  */
  sprintf(typeA,"-");
  if (addU->Nneighbor[ishare_a]==1){
    Ba[0] = Oa[0] - addU->pos[ina[0]][0];
    Ba[1] = Oa[1] - addU->pos[ina[0]][1];
    normalize_2Dvec(Ba);
    sprintf(typeA,"T");
  }
  else if ((addU->Nneighbor[ishare_a]==2)&&(addU->type!='R')){
    Ba[0] = -(addU->pos[ina[0]][0] + addU->pos[ina[1]][0])/2.0 + Oa[0];
    Ba[1] = -(addU->pos[ina[0]][1] + addU->pos[ina[1]][1])/2.0 + Oa[1];
    normalize_2Dvec(Ba);
    sprintf(typeA,"2C");
  }
  else if ((addU->Nneighbor[ishare_a]==2)&&(addU->type=='R')&&(mol->atoms[anum_share].Nnei_heavy==3)){
    Ba[0] =   addU->pos[ina[0]][1] - addU->pos[ina[1]][1];
    Ba[1] = -(addU->pos[ina[0]][0] - addU->pos[ina[1]][0]);
    normalize_2Dvec(Ba);
    if (congestion_in_ring_for_dir(addU,ishare_a,Ba[0],Ba[1])>congestion_in_ring_for_dir(addU,ishare_a,-Ba[0],-Ba[1])){
      Ba[0] = -Ba[0];
      Ba[1] = -Ba[1];
    }
    sprintf(typeA,"2R");
  }
  else if ((addU->Nneighbor[ishare_a]==2)&&(addU->type=='R')&&(mol->atoms[anum_share].Nnei_heavy==4)){
    cal_angular_three_division_point(bdirs,Oa,addU->pos[ina[0]],addU->pos[ina[1]]);
    minConSc = -1.0;
    for (j=0;j<4;++j){
      ConSc  = congestion_in_ring_for_dir(addU,ishare_a,bdirs[j][0],bdirs[j][1]);
      if ((minConSc<0.0)||( ConSc < minConSc)){ minConSc = ConSc; Ba[0] = bdirs[j][0]; Ba[1] = bdirs[j][1]; }
    }
    sprintf(typeA,"2Rfour");
  }
  else if ((addU->Nneighbor[ishare_a]==3)&&(addU->type!='R')){
   max_angle = -1.0;
   for (j=0;j<3;++j){
     cal_angular_mid_point(&angle,bdir,Oa,addU->pos[ina[j]],addU->pos[ina[(j+1)%3]],addU->pos[ina[(j+2)%3]]);
     if ((max_angle<0.0)||(angle > max_angle)){
       max_angle = angle;
       Ba[0] = bdir[0];
       Ba[1] = bdir[1];
     }
   }
   sprintf(typeA,"3C");
  }
  else if ((addU->Nneighbor[ishare_a]==3)&&(addU->type=='R')){
   minConSc = -1.0;
   for (j=0;j<3;++j){
     cal_mid_point(bdir,Oa,addU->pos[ina[j]],addU->pos[ina[(j+1)%3]]);

     ConSc  = congestion_in_ring_for_dir(addU,ishare_a,bdir[0],bdir[1]);
     if ((minConSc<0.0)||( ConSc < minConSc)){ minConSc = ConSc; Ba[0] = bdir[0]; Ba[1] = bdir[1]; }

     ConSc  = congestion_in_ring_for_dir(addU,ishare_a,-bdir[0],-bdir[1]);
     if ((minConSc<0.0)||( ConSc < minConSc)){ minConSc = ConSc; Ba[0] = -bdir[0]; Ba[1] = -bdir[1]; }

     bdir[0] = Oa[0] - addU->pos[ina[j]][0];
     bdir[1] = Oa[1] - addU->pos[ina[j]][1];
     ConSc  = congestion_in_ring_for_dir(addU,ishare_a,bdir[0],bdir[1]);
     if ((minConSc<0.0)||( ConSc < minConSc)){ minConSc = ConSc; Ba[0] = bdir[0]; Ba[1] = bdir[1]; }


   }
   sprintf(typeA,"3R");
  }

  else {printf("#ERROR(combine):for addU[%d] Nneighbor = %d \n",ishare_a,addU->Nneighbor[ishare_a]); return(0);}
  /* printf("#Oa %f %f Ba %f %f\n",Oa[0],Oa[1],Ba[0],Ba[1]); */

  /*********************************/
  /** (2) Setup Os and Bs for seed */
  /*********************************/

  Os[0] = mol->atoms[anum_share].Pos[0];
  Os[1] = mol->atoms[anum_share].Pos[1];
  if ((mol->atoms[anum_share].ring == '-') || (mol->atoms[anum_share].ring==' ')) seedRing = 'F';
            else seedRing = 'T';

  Nneighbor_share = 0;
  for (i=0;i<mol->Natom;++i){ 
    if ((mol->atoms[i].mark==1) && (mol->conmap.map[anum_share][i]!='0')){
      anum_nei_seed[Nneighbor_share] = i;
      Nneighbor_share += 1;
    }
  }

  
  if (Nneighbor_share==1){
    Bs[0] =  Os[0] - mol->atoms[anum_nei_seed[0]].Pos[0];
    Bs[1] =  Os[1] - mol->atoms[anum_nei_seed[0]].Pos[1];
    normalize_2Dvec(Bs);
    sprintf(typeS,"1");
  }
  else if ((Nneighbor_share ==2)&&(seedRing!='T')){
    Bs[0] = -(mol->atoms[anum_nei_seed[0]].Pos[0]  + mol->atoms[anum_nei_seed[1]].Pos[0])/2.0 + Os[0];
    Bs[1] = -(mol->atoms[anum_nei_seed[0]].Pos[1]  + mol->atoms[anum_nei_seed[1]].Pos[1])/2.0 + Os[1];
    normalize_2Dvec(Bs);
    sprintf(typeS,"2C");
  }
  else if ((Nneighbor_share ==2)&&(seedRing=='T')&&(mol->atoms[anum_share].Nnei_heavy==3)&&(mol->atoms[anum_share].ring=='r')){
    Bs[0] =  (mol->atoms[anum_nei_seed[0]].Pos[1]  - mol->atoms[anum_nei_seed[1]].Pos[1]);
    Bs[1] = -(mol->atoms[anum_nei_seed[0]].Pos[0]  - mol->atoms[anum_nei_seed[1]].Pos[0]);
    normalize_2Dvec(Bs);
    if (congestion_in_marked_atoms_for_dir(mol,anum_share,Bs[0],Bs[1],'T') > congestion_in_marked_atoms_for_dir(mol,anum_share,-Bs[0],-Bs[1],'T')){
      Bs[0] = -Bs[0];
      Bs[1] = -Bs[1];
    }
    sprintf(typeS,"2r");
  }
  else if ((Nneighbor_share ==2)&&(seedRing=='T')&&(mol->atoms[anum_share].Nnei_heavy==3)){
    Bs[0] =  (mol->atoms[anum_nei_seed[0]].Pos[1]  - mol->atoms[anum_nei_seed[1]].Pos[1]);
    Bs[1] = -(mol->atoms[anum_nei_seed[0]].Pos[0]  - mol->atoms[anum_nei_seed[1]].Pos[0]);
    normalize_2Dvec(Bs);
    if (congestion_in_marked_atoms_for_dir(mol,anum_share,Bs[0],Bs[1],'-') > congestion_in_marked_atoms_for_dir(mol,anum_share,-Bs[0],-Bs[1],'-')){
      Bs[0] = -Bs[0];
      Bs[1] = -Bs[1];
    }
    sprintf(typeS,"2R");
  }
  else if ((Nneighbor_share ==2)&&(seedRing=='T')&&(mol->atoms[anum_share].Nnei_heavy==4)){
    cal_angular_three_division_point(bdirs,Os,mol->atoms[anum_nei_seed[0]].Pos,mol->atoms[anum_nei_seed[1]].Pos);
    minConSc = -1.0;
    for (j=0;j<4;++j){
      ConSc =  congestion_in_marked_atoms_for_dir(mol,anum_share,bdirs[j][0],bdirs[j][1],'-');
      if ((minConSc<0)||(ConSc < minConSc)){ minConSc = ConSc; Bs[0] = bdirs[j][0]; Bs[1] = bdirs[j][1]; }
    } 
    sprintf(typeS,"2Rfour");
  }
  else if ((Nneighbor_share==3)&&(seedRing!='T')){ 
   max_angle = -1.0;
   for (j=0;j<3;++j){
     printf("#%d %d %d\n",
         mol->atoms[anum_nei_seed[j]].num_in_file, mol->atoms[anum_nei_seed[(j+1)%3]].num_in_file, mol->atoms[anum_nei_seed[(j+2)%3]].num_in_file);
     cal_angular_mid_point(&angle,bdir,Os,
            mol->atoms[anum_nei_seed[j]].Pos, mol->atoms[anum_nei_seed[(j+1)%3]].Pos, mol->atoms[anum_nei_seed[(j+2)%3]].Pos);
       if ((max_angle<0.0)||(angle > max_angle)){
         max_angle = angle;
         Bs[0] = bdir[0];
         Bs[1] = bdir[1];
       }
   }
   sprintf(typeS,"3C");
  }
  else if ((Nneighbor_share==3)&&(seedRing=='T')){ 
   minConSc = -1.0;
   for (j=0;j<3;++j){
     cal_mid_point(bdir,Os, mol->atoms[anum_nei_seed[j]].Pos, mol->atoms[anum_nei_seed[(j+1)%3]].Pos);
     ConSc =  congestion_in_marked_atoms_for_dir(mol,anum_share,bdir[0],bdir[1],'-');
     if ((minConSc<0)||(ConSc < minConSc)){ minConSc = ConSc; Bs[0] = bdir[0]; Bs[1] = bdir[1]; }
     ConSc =  congestion_in_marked_atoms_for_dir(mol,anum_share,-bdir[0],-bdir[1],'-');
    if ((minConSc<0)||(ConSc < minConSc)){ minConSc = ConSc; Bs[0] = -bdir[0]; Bs[1] = -bdir[1]; }
     bdir[0] = Os[0] - mol->atoms[anum_nei_seed[j]].Pos[0];
     bdir[1] = Os[1] - mol->atoms[anum_nei_seed[j]].Pos[1];
     normalize_2Dvec(bdir);
     ConSc =  congestion_in_marked_atoms_for_dir(mol,anum_share,bdir[0],bdir[1],'-');
     if ((minConSc<0)||(ConSc < minConSc)){ minConSc = ConSc; Bs[0] = bdir[0]; Bs[1] = bdir[1]; }
   }
   sprintf(typeS,"3R");
  }

  else {printf("#ERROR(combine):for seed atom [%d %s Nneighbor = %d] \n",
    mol->atoms[anum_share].num_in_file,mol->atoms[anum_share].atomname,Nneighbor_share); return(0);}

  if (PAR.OutCalProcess=='T'){
    printf("#add_unit %d %c typeA '%s' seed typeS '%s'\n",addU->num,addU->type,typeA,typeS);
    printf("#Oa %f %f Ba %f %f Os %f %f Bs %f %f\n",Oa[0],Oa[1], Ba[0],Ba[1], Os[0],Os[1],Bs[0],Bs[1]); 
  } 
 /*** (3) Calculated rotated addU XY, and copy into mol->atoms[] **/
  /***  add_pos_new = R(theta) * (add_pos_old - Oa) + Os **/

  th = myacos(-Ba[0]*Bs[0] -Ba[1]*Bs[1]);
  if ((Ba[0]*Bs[1]-Ba[1]*Bs[0])>0.0) th = -th; 
  c_t = cos(th);
  s_t = sin(th);

  /* printf("th %f c_t %f s_t %f\n",th,c_t,s_t);  */

  for (i=0;i<addU->Natom;++i){
    dpos[0] = addU->pos[i][0] - Oa[0];
    dpos[1] = addU->pos[i][1] - Oa[1];
    addU->pos[i][0] = (c_t * dpos[0] - s_t * dpos[1]) + Os[0];
    addU->pos[i][1] = (s_t * dpos[0] + c_t * dpos[1]) + Os[1];
    ii = addU->num_atoms[i];
    mol->atoms[ii].Pos[0] = addU->pos[i][0];
    mol->atoms[ii].Pos[1] = addU->pos[i][1];
    mol->atoms[ii].mark = 1;
  }

  return(1);

} /* end of combine_new_UNIT_into_seed_atoms() */




void Decompose_Atoms_into_UNITs(mol,HeadUnit)
  struct MOLECULE *mol;
  struct UNIT *HeadUnit;
{
  struct RING *rn;
  int a,b,i,j,ii,jj,minD,min_a,min_b;
  struct UNSIGNED_CHAR2DMAP  used_bond;
  int *path, Natom_in_path,Natom_nonused;
  char used, *term_atom;

/*
  printf("#Decompose_Atoms_into_UNITs(mol,HeadUnit)\n");
*/
  Malloc_UNSIGNED_CHAR2DMAP(&used_bond,mol->Natom);
  path = (int *)malloc(sizeof(int)*mol->Natom);
  term_atom = (char *)malloc(sizeof(char)*mol->Natom);

  HeadUnit->next = HeadUnit->prev = NULL;
  HeadUnit->num = -1;
  for (a=0;a<mol->Natom;++a){
    mol->atoms[a].mark = 0;
    mol->atoms[a].ring = ' ';
    for (b=0;b<mol->Natom;++b){ 
       used_bond.map[a][b] = used_bond.map[b][a] = 0; 
    }
  }
  
  /*** [1] copy RingBlocks into the HeadUnit ***/
  rn = &(mol->HeadBlock);
  while (rn->next != NULL){
    rn = rn->next;
    add_new_UNIT_to_tail(HeadUnit, rn->Natom, rn->num_atoms,'R');
    for (i=0;i<rn->Natom;++i){
      ii = rn->num_atoms[i];
      mol->atoms[ii].mark = 1;
      mol->atoms[ii].ring = 'r';
      /* printf("#i %d ii %d mark %d\n",i,ii,mol->atoms[ii].mark); */
      for (j=i+1;j<rn->Natom;++j){
        jj = rn->num_atoms[j];
        if (mol->conmap.map[ii][jj] != '0'){ 
          used_bond.map[ii][jj] = used_bond.map[jj][ii] = 1; 
        }
      }
    }
  }

  /*** [2] Take the longest non-used chain UNIT, one by one. ***/

  do{
    /* (2-1) mark used_atom[a] = used_atom[b] = 0 for in non-marked bonds with used_bond[a][b]==0 */
    for (a=0;a<mol->Natom;++a){ term_atom[a] = 0; }

    for (a=0;a<mol->Natom;++a){
      if (mol->atoms[a].one_char_ele != 'H'){
        for (b=0;b<mol->Natom;++b){
          if (mol->atoms[b].one_char_ele != 'H'){
            if ((mol->conmap.map[a][b]!='0')&&(used_bond.map[a][b] ==0)){ 
              term_atom[a] = term_atom[b] = 1;
            }
          }
         }
       }
    }
    Natom_nonused = 0;
    for (a=0;a<mol->Natom;++a){ 
       if ((term_atom[a]==1)&&(mol->atoms[a].one_char_ele != 'H')) {Natom_nonused += 1;}
       /* printf("[%d] %d %s used %d\n",a,mol->atoms[a].num_in_file,mol->atoms[a].atomname,used_atom[a]); */
    } 
    /* printf("#Natom_nonused %d\n",Natom_nonused); */

    /* (2-2) Find the longest chain consisting of non-used bonds */
    if (Natom_nonused>0){
      minD = -1;
      min_a = min_b = 0; 
      for (a=0;a<mol->Natom;++a){
        if ((term_atom[a]==1)&&(mol->atoms[a].one_char_ele != 'H')) {
          for (b=a+1;b<mol->Natom;++b){
            if ((term_atom[b]==1)&&(mol->atoms[b].one_char_ele != 'H')) {
              set_shortest_path_between_two_atoms(mol,a,b,path,&Natom_in_path);
              used = 0;
              /* check used bond on the path*/
              for (i=1;i<Natom_in_path;++i){
                if (used_bond.map[path[i-1]][path[i]] ==1){ used = 1;}
              }
 
              /* check used atom on the middle (non-terminal) of path */
              for (i=1;i<(Natom_in_path-1);++i){ 
                 if (mol->atoms[path[i]].mark ==1) used = 1; 
              }
              if ((used==0)&&((minD<0)||(mol->topodismap.map[a][b]>minD))){
          /*
              if ((used==0)&&
                  ((minD<0)
                   ||(mol->topodismap.map[a][b]>minD)
                   ||((mol->topodismap.map[a][b]==minD)&&
                      ( (Number_of_Element(mol->atoms[a].element)+Number_of_Element(mol->atoms[b].element))
                       <(Number_of_Element(mol->atoms[min_a].element)+Number_of_Element(mol->atoms[min_b].element))))
                  )){
          */
                minD = mol->topodismap.map[a][b];
                min_a = a;
                min_b = b;
              }      
            } 
          }
        }
      }
  
    /* (2-3) Add the longest chain into HeadUnit */
      /* printf("#minD %d\n",minD); */
      set_shortest_path_between_two_atoms(mol,min_a,min_b,path,&Natom_in_path);
      add_new_UNIT_to_tail(HeadUnit,Natom_in_path,path,'C');
      /* show_UNITs(mol, HeadUnit, "Decomposing now"); */

      for (i=0;i<Natom_in_path;++i){
        mol->atoms[path[i]].mark = 1;
        if ((i>=1)&&(mol->conmap.map[path[i-1]][path[i]] != '0')) 
           used_bond.map[path[i-1]][path[i]] = used_bond.map[path[i]][path[i-1]] = 1; 
      }
    }
     /* printf("#Natom_nonused %d\n", Natom_nonused);  */

  } while (Natom_nonused>0);

  setup_neighbors_for_UNITs(mol, HeadUnit);

  Free_UNSIGNED_CHAR2DMAP(&used_bond,mol->Natom);
  free(term_atom);
  free(path);

} /* end of Decompose_Atoms_into_BLOCKs() */






int set_XY_to_ring_UNIT(mol,unit)
  struct MOLECULE *mol;
  struct UNIT *unit;
{
  int i,j,ii,jj,Nround;
  struct UNIT HeadRing; 
  struct UNIT *rn,*sn; 
  int *path,*pathsmlt, Natom_pathsmlt, Nring,*ring2unit;
  char *mark, improper_ring;
  int Nshare,*ashare,ishare_r[12],ishare_s[12],Nring_add;
  struct CHAR2DMAP used_bond;

  ashare   = (int *)malloc(sizeof(int)*mol->Natom);
  path     = (int *)malloc(sizeof(int)*mol->Natom);
  pathsmlt = (int *)malloc(sizeof(int)*mol->Natom);
  ring2unit = (int *)malloc(sizeof(int)*mol->Natom);
  mark = (char *)malloc(sizeof(char)*mol->Natom);
  Malloc_CHAR2DMAP(&used_bond,mol->Natom);
  /*
  printf("#set_XY_to_ring_UNIT(mol,un: %d '%c' Natom %d)\n",unit->num,unit->type,unit->Natom);
   */
 
  if (unit->type != 'R') return(0);
 
 /*** [1] Find all the RING structures (like SSSR) in the unit ***/ 
  HeadRing.next = HeadRing.prev = NULL; HeadRing.num  = -1;  
  for (i=0;i<unit->Natom;++i) mol->atoms[unit->num_atoms[i]].mark = 0;

  for (i=0;i<unit->Natom;++i){
/*
    printf(">start_atom %d mark %d\n",mol->atoms[unit->num_atoms[i]].num_in_file,mol->atoms[unit->num_atoms[i]].mark);
*/
    if (mol->atoms[unit->num_atoms[i]].mark==0){
      for (j=0;j<unit->Natom;++j) mark[unit->num_atoms[j]] = 0;
      Natom_pathsmlt = -1;
      mark_neighbors_for_smallest_cycle_in_UNIT(&HeadRing,mol,unit,1,path,mark,&Natom_pathsmlt,pathsmlt,
                                            unit->num_atoms[i],unit->num_atoms[i]);
      if (Natom_pathsmlt>=3){
        /* printf("##find ring##(%s)\n",mol->name); */
        add_new_UNIT_to_tail(&HeadRing, Natom_pathsmlt,pathsmlt, 'R');
        for (j=0;j<Natom_pathsmlt;++j){
          mol->atoms[pathsmlt[j]].mark = 1;
          used_bond.map[pathsmlt[j]][pathsmlt[(j+1)%Natom_pathsmlt]] = '1';
          used_bond.map[pathsmlt[(j+1)%Natom_pathsmlt]][pathsmlt[j]] = '1';
        }
      }
    }
  }

  Nring = Number_of_UNITs(&HeadRing);
  /*
  printf("#Nring %d\n",Nring);
  show_UNITs(mol, &HeadRing, "list of RINGs");
   */

  unit->set_localXY = 'F';

  improper_ring = 0;
 
  /** If all the bonds are NOT covered by the (SSSR)-rings, this UNIT is "improper_ring" **/  
  /** It means this set of rings is NOT really SSSR !!  */
  for (i=0;i<unit->Natom;++i){
    ii = unit->num_atoms[i];
    for (j=i+1;j<unit->Natom;++j){
      jj = unit->num_atoms[j];
      if ((mol->conmap.map[ii][jj] != '0') && (used_bond.map[ii][jj] == '0')) improper_ring = 1;   
    }
  }


  /*** [2] Set local XY to the ring ***/
  /** (2-1) For single ring --> regular polyhedron **/
  if ((Nring==1)&&(HeadRing.next->Natom == unit->Natom)){
    rn = HeadRing.next;
    set_localXY_in_regular_polygon(rn);
    make_ring_to_unit_index(ring2unit,rn,unit);
    for (i=0;i<unit->Natom;++i){
      unit->pos[ring2unit[i]][0] = rn->pos[i][0];
      unit->pos[ring2unit[i]][1] = rn->pos[i][1];
      /*
      printf("ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f rn_i %d uint_i %d %f %f\n",
        mol->atoms[rn->num_atoms[i]].num_in_file,"CA","PEN",'-',"1",rn->pos[i][0], rn->pos[i][1],0.0,i,ring2unit[i],
      unit->pos[ring2unit[i]][0], unit->pos[ring2unit[i]][1]);
      */
    }
    unit->set_localXY = 'T';
    /*
    for (i=0;i<unit->Natom;++i){
      printf("ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f i %d anum %d\n",
        mol->atoms[unit->num_atoms[i]].num_in_file,"CA","PEN",'-',"1",
         unit->pos[i][0], unit->pos[i][1],0.0,i,unit->num_atoms[i]);
    }
    */
  }

  /*** (2-2) Checking proper (linear or tree-shape) multi rings, or not ***/
  if ((unit->set_localXY=='F')&&(Nring>1)){
    for (j=0;j<unit->Natom;++j) mark[unit->num_atoms[j]] = 0;
    rn = &HeadRing; sn = NULL;
    while (rn->next !=NULL){
      rn = rn->next;
      /* if (rn->Natom>8) {improper_ring = 1;}  if one of SSSR is more than 8-member, it is improper. */
      sn = rn;
      while (sn->next !=NULL){
        sn = sn->next;
        find_shared_atoms_bwn_two_rings(rn,sn,&Nshare,ashare,ishare_r,ishare_s);
        /* printf("#rn %d sn %d Nshare %d\n",rn->num,sn->num,Nshare); */
         /* If Nshare !=2, this is improper ring */
        if ((Nshare!=0)&&(Nshare!=2)){ 
          improper_ring = 1; 
          printf("#unit [%d %c Natom %d] is improper. because Nshare = %d.\n",unit->num,unit->type,unit->Natom,Nshare);
        }
        if (rn->num < sn->num) {for (j=0;j<Nshare;++j){ mark[ashare[j]] += 1; }}
      }
    }
    /* If one atom belongs to >2 rings, this unit is "cyclic"-ring structure */ 
    for (j=0;j<unit->Natom;++j){
      if (mark[unit->num_atoms[j]]>1) improper_ring = 1;
      /* printf("%d %d mark %d\n",unit->num_atoms[j],mol->atoms[unit->num_atoms[j]].num_in_file,mark[unit->num_atoms[j]]); */
    }
  }
  
  /*** (2-3) If proper multi ring, assemble rings one by one ***/
  if ((unit->set_localXY=='F')&&(Nring>1)&&(improper_ring==0)){

    printf("##(case:tree-ring) Add seed-connecting ring into the seed, one by one.\n");
    /** Set localXY into each ring **/
    rn = &HeadRing;
    while (rn->next !=NULL){
      rn = rn->next; 
      set_localXY_in_regular_polygon(rn);
      rn->set_localXY = 'T';
      rn->set_globalXY = 'F';
    }
    /** Take the first ring as seed **/
    sn = HeadRing.next;
    sn->set_globalXY = 'T';
    make_ring_to_unit_index(ring2unit,sn,unit);
    for (i=0;i<sn->Natom;++i){
      unit->pos[ring2unit[i]][0] = sn->pos[i][0];
      unit->pos[ring2unit[i]][1] = sn->pos[i][1];
    }
    /** Add seed-connecting ring into the seed, one by one. **/
    Nround = 0; 
    do{
      Nround += 1;
      Nring_add = 0;
      rn = &HeadRing;
      while ((rn->next != NULL)&&(Nring_add==0)){
        rn = rn->next;
        if (rn->set_globalXY != 'T'){
          sn = &HeadRing;
          while ((sn->next != NULL)&&(Nring_add==0)){
            sn = sn->next;
            if ((sn->set_globalXY=='T')&&(find_shared_atoms_bwn_two_rings(rn,sn,&Nshare,ashare,ishare_r,ishare_s)==2)){
               combine_new_ring_into_seed_ring(mol,unit,rn,sn,ishare_r,ishare_s);
               Nring_add += 1;
               rn->set_globalXY = 'T';
             } 
            }
         }
       }
      printf("#>ROUND[%d] Nring_add %d Nring %d\n",Nround,Nring_add,Nring);
    } while ((Nring_add>0)&&((Nround+1)<Nring));
  
    unit->set_localXY = 'T';
  }

  /** [5] Using distance geometry  **/
  if (unit->set_localXY=='F'){ set_XY_to_UNIT_by_distance_geometry(mol,unit,&HeadRing); }
 
  /** [6] Just place in the circle (better than nothing...) **/
  if (unit->set_localXY=='F'){ set_localXY_in_regular_polygon(unit); }
/*
  printf("#restraint_energy for the unit[%d] Natom %d ::  %f\n",unit->num,unit->Natom,restraint_energy(mol,unit,&HeadRing));
*/ 
  Free_CHAR2DMAP(&used_bond);
  free(ashare);
  free(path);
  free(pathsmlt);
  free(ring2unit);
  free(mark);
  return(1);

} /* end of set_XY_to_ring_UNIT() */






int set_hydrogen_XY(mol)
  struct MOLECULE *mol;
{
  int i,j,min_j,a,b,Nhydro_nei,Nheavy_nei, hydro_nei[64],heavy_nei[64];
  float bvec[2],rbvec[2],Lhbond,th,c_t,s_t,midp[3][2],score[3],min_score;
  /* printf("#set_hydrogen_XY(mol)\n"); */

 /*
 for (i=0;i<mol->Natom;++i){
   if (mol->atoms[i].one_char_ele=='H'){ 
    mol->atoms[i].Pos[0] = mol->atoms[i].Pos[1] = mol->atoms[i].Pos[2] = 0.0; 
   }
 }
 return(0);
  */

  Lhbond = 0.4*PAR.LengthOfBond2D;
  for (i=0;i<mol->Natom;++i){
    if ((mol->atoms[i].Pos[0] +  mol->atoms[i].Pos[1] + mol->atoms[i].Pos[2])==0.0){
       mol->atoms[i].mark = 0;
    }
    else mol->atoms[i].mark = 1;
  } 



 /*** [i-loop] : Check each heavy atom bound to hydrogens ***/

  for (i=0;i<mol->Natom;++i){

    /** (1) Check hydrogen neighbors, and heavy neighbors for atom[i] */
    Nhydro_nei = Nheavy_nei = 0;
    for (j=0;j<mol->Natom;++j){
      if (mol->conmap.map[i][j]!='0'){
        if (mol->atoms[j].one_char_ele=='H'){
          hydro_nei[Nhydro_nei] = j;
          Nhydro_nei += 1;


        }
        else {
          heavy_nei[Nheavy_nei] = j;
          Nheavy_nei += 1; 
        }
      }
    } /* j */

    /* printf("#i %d Nheavy_nei %d Nhydro_nei %d\n",i,Nheavy_nei,Nhydro_nei); */
    /** (2) Set XY of hydrogen atoms **/
    /** Nheavy_nei == 2 */
    if ((Nheavy_nei==2)&&(Nhydro_nei==1)){
       bvec[0] = mol->atoms[i].Pos[0] - (mol->atoms[heavy_nei[0]].Pos[0] + mol->atoms[heavy_nei[1]].Pos[0])/2.0;
       bvec[1] = mol->atoms[i].Pos[1] - (mol->atoms[heavy_nei[0]].Pos[1] + mol->atoms[heavy_nei[1]].Pos[1])/2.0;
       normalize_2Dvec(bvec);
       if (mol->atoms[hydro_nei[0]].mark == 0){
         mol->atoms[hydro_nei[0]].Pos[0] = mol->atoms[i].Pos[0] + Lhbond*bvec[0];
         mol->atoms[hydro_nei[0]].Pos[1] = mol->atoms[i].Pos[1] + Lhbond*bvec[1];
      }
    }
    else if ((Nheavy_nei==2)&&(Nhydro_nei==2)){
      bvec[0] = mol->atoms[i].Pos[0] - (mol->atoms[heavy_nei[0]].Pos[0] + mol->atoms[heavy_nei[1]].Pos[0])/2.0;
      bvec[1] = mol->atoms[i].Pos[1] - (mol->atoms[heavy_nei[0]].Pos[1] + mol->atoms[heavy_nei[1]].Pos[1])/2.0;
      normalize_2Dvec(bvec);
      for (j=0;j<2;++j){
        if (j==0) th =   M_PI * 40/180.0;
            else  th = - M_PI * 40/180.0;
        c_t = cos(th);
        s_t = sin(th);
        rbvec[0] = c_t * bvec[0] - s_t * bvec[1];
        rbvec[1] = s_t * bvec[0] + c_t * bvec[1];
        if (mol->atoms[hydro_nei[j]].mark == 0){
          mol->atoms[hydro_nei[j]].Pos[0] = mol->atoms[i].Pos[0] + Lhbond*rbvec[0];
          mol->atoms[hydro_nei[j]].Pos[1] = mol->atoms[i].Pos[1] + Lhbond*rbvec[1];
        }
      }
    }
    /** Nheavy_nei == 3 */
    else if ((Nheavy_nei==3)&&(Nhydro_nei==1)){

      a = 0; b = 0; 
      for (j=0;j<3;++j){
             if (j==0) {a = 0; b = 1;}
        else if (j==1) {a = 0; b = 2;}
        else if (j==2) {a = 1; b = 2;}
        bvec[0] = (mol->atoms[heavy_nei[a]].Pos[0] + mol->atoms[heavy_nei[b]].Pos[0])/2.0 - mol->atoms[i].Pos[0];
        bvec[1] = (mol->atoms[heavy_nei[a]].Pos[1] + mol->atoms[heavy_nei[b]].Pos[1])/2.0 - mol->atoms[i].Pos[1];
        normalize_2Dvec(bvec);
        midp[j][0] = mol->atoms[i].Pos[0] + Lhbond*bvec[0];
        midp[j][1] = mol->atoms[i].Pos[1] + Lhbond*bvec[1];
        score[j] = point_congestion_score(midp[j][0],midp[j][1],mol);
      } 
      min_j = 0; min_score = score[0];
      for (j=1;j<3;++j){ 
         if (score[j]<min_score){min_score = score[j]; min_j = j;} 
      }
      if (mol->atoms[hydro_nei[0]].mark == 0){
        mol->atoms[hydro_nei[0]].Pos[0] = midp[min_j][0];
        mol->atoms[hydro_nei[0]].Pos[1] = midp[min_j][1];
      }
    }

    /** Nheavy_nei == 1 */
    else if ((Nheavy_nei==1)&&(Nhydro_nei>0)){
      bvec[0] = mol->atoms[i].Pos[0] - mol->atoms[heavy_nei[0]].Pos[0];
      bvec[1] = mol->atoms[i].Pos[1] - mol->atoms[heavy_nei[0]].Pos[1];
      normalize_2Dvec(bvec);
      for (j=0;j<Nhydro_nei;++j){
        th =  M_PI*60.0/180.0;
        if (Nhydro_nei == 1) th =  M_PI*60.0/180.0;
        if (Nhydro_nei == 1) th =  -M_PI*60.0/180.0;
        else if (Nhydro_nei==2){
               if (j==0) th =  M_PI*60.0/180.0;
          else if (j==1) th = -M_PI*60.0/180.0;
        }
        else if (Nhydro_nei==3){
               if (j==0) th =  M_PI*80.0/180.0;
          else if (j==1) th =  0.0;
          else if (j==2) th = -M_PI*80.0/180.0;
        }

     
        c_t = cos(th);
        s_t = sin(th);
        rbvec[0] = c_t * bvec[0] - s_t * bvec[1];
        rbvec[1] = s_t * bvec[0] + c_t * bvec[1];

        if (mol->atoms[hydro_nei[j]].mark == 0){
          mol->atoms[hydro_nei[j]].Pos[0] = mol->atoms[i].Pos[0] + Lhbond*rbvec[0];
          mol->atoms[hydro_nei[j]].Pos[1] = mol->atoms[i].Pos[1] + Lhbond*rbvec[1];
        }
      }
    }
    /** Nheavy_nei == 0  (like H2_O, CH4, NH3)*/
    else if ((Nheavy_nei==0)&&(Nhydro_nei>0)){
      bvec[0] = 0.0; bvec[1] = 1.0;
      for (j=0;j<Nhydro_nei;++j){
        th =  M_PI*60.0/180.0;
        if (Nhydro_nei == 1) th =  0.0;
        else if (Nhydro_nei==2){
               if (j==0) th =   M_PI*60/180.0;
          else if (j==1) th =  -M_PI*60.0/180.0;
        }
        else if (Nhydro_nei==3){
               if (j==0) th =  M_PI*0/180.0;
          else if (j==1) th =  M_PI*120/180.0;
          else if (j==2) th =  M_PI*240/180.0;
        }
        else if (Nhydro_nei==4){
               if (j==0) th =  M_PI*0/180.0;
          else if (j==1) th =  M_PI*90/180.0;
          else if (j==2) th =  M_PI*180.0/180.0;
          else if (j==3) th =  M_PI*270.0/180.0;
        }
        c_t = cos(th);
        s_t = sin(th);
        rbvec[0] = c_t * bvec[0] - s_t * bvec[1];
        rbvec[1] = s_t * bvec[0] + c_t * bvec[1];
        if (mol->atoms[hydro_nei[j]].mark == 0){
          mol->atoms[hydro_nei[j]].Pos[0] = mol->atoms[i].Pos[0] + Lhbond*rbvec[0];
          mol->atoms[hydro_nei[j]].Pos[1] = mol->atoms[i].Pos[1] + Lhbond*rbvec[1];
        }
      }
    }
    /** other case: overlap onto the heavy atom (better than nothing ...) **/
    else{
      for (j=0;j<Nhydro_nei;++j){
       if (mol->atoms[hydro_nei[j]].mark == 0){
         mol->atoms[hydro_nei[j]].Pos[0] = mol->atoms[i].Pos[0];
         mol->atoms[hydro_nei[j]].Pos[1] = mol->atoms[i].Pos[1];
       }
      }
    }
 
  } /* i */

  return(1);
} /* end of set_hydrogen_XY() */




void set_localXY_in_regular_polygon(unit)
  struct UNIT *unit;
{
  float theta,R;
  int i;

  theta =  (2.0*M_PI)/unit->Natom;
  R = PAR.LengthOfBond2D*sqrt(1.0/2.0/(1.0-cos(theta)));
  for (i=0;i<unit->Natom;++i){
    unit->pos[i][0] = R * cos(theta*i);
    unit->pos[i][1] = R * sin(theta*i);
  }
  unit->set_localXY = 'T';
}


int combine_new_ring_into_seed_ring(mol,unit,addR,seedR,ishare_a,ishare_s)
  struct MOLECULE *mol;
  struct UNIT *unit;  /* parent UNIT including addR and seedR rings */
  struct UNIT *addR;  /* ring to be added */
  struct UNIT *seedR; /* seed ring */
  int    *ishare_a; /* index of the two shared atoms for addU */
  int    *ishare_s; /* index of the two shared atoms for seedU */
{
  float Os[2],Oa[2]; /* shared atom XY coordinates */
  float Bs[2],Ba[2]; /* shared bond vector coordinates */
  float Gs[2],Ga[2]; /* center of the ring */
  float dpos[2],th,c_t,s_t;
  float Pa[2],Qa[2],Ps[2],Qs[2];
  int i,*ring2unit,k;

  printf("#combine_new_ring_into_seed_ring(unit[%d] addR[%d] seedR[%d] add %d %d seed %d %d)\n",
    unit->num,addR->num,seedR->num,
mol->atoms[addR->num_atoms[ishare_a[0]]].num_in_file,
mol->atoms[addR->num_atoms[ishare_a[1]]].num_in_file,
mol->atoms[seedR->num_atoms[ishare_s[0]]].num_in_file,
mol->atoms[seedR->num_atoms[ishare_s[1]]].num_in_file);

  ring2unit = (int *)malloc(sizeof(int)*addR->Natom);

  /** (1) Checking chirarity of ring structure. */
  Ga[0] = Ga[1] = Gs[0] = Gs[1] = 0.0;
  for (i=0;i<addR->Natom;++i){ Ga[0] += addR->pos[i][0]; Ga[1] += addR->pos[i][1]; }
  Ga[0] /= addR->Natom; Ga[1] /= addR->Natom;
  
  for (i=0;i<seedR->Natom;++i){ Gs[0] += seedR->pos[i][0]; Gs[1] += seedR->pos[i][1]; }
  Gs[0] /= seedR->Natom; Gs[1] /= seedR->Natom;

  for (k=0;k<2;++k){  
    Pa[k] = addR->pos[ishare_a[0]][k] - Ga[k]; 
    Qa[k] = addR->pos[ishare_a[1]][k] - Ga[k]; 
    Ps[k] = seedR->pos[ishare_s[0]][k] - Gs[k]; 
    Qs[k] = seedR->pos[ishare_s[1]][k] - Gs[k]; 
  }

  /** If the chirarities of two rings are same, ring addR should be inverted. **/
  if ((Pa[0]*Qa[1]-Pa[1]*Qa[0])*(Ps[0]*Qs[1]-Ps[1]*Qs[0])>0.0){
    for (i=0;i<addR->Natom;++i){ 
      addR->pos[i][0] = -addR->pos[i][0]; 
    }
  }

  /** (2) Setup Oa and Ba for addR, Os and Bs for seedR */
  Oa[0] = addR->pos[ishare_a[0]][0];
  Oa[1] = addR->pos[ishare_a[0]][1];
  
  Ba[0] = addR->pos[ishare_a[1]][0] - Oa[0];
  Ba[1] = addR->pos[ishare_a[1]][1] - Oa[1];
  normalize_2Dvec(Ba);
 
  Os[0] = seedR->pos[ishare_s[0]][0];
  Os[1] = seedR->pos[ishare_s[0]][1];
 
  Bs[0] = seedR->pos[ishare_s[1]][0] - Os[0];
  Bs[1] = seedR->pos[ishare_s[1]][1] - Os[1];
  normalize_2Dvec(Bs);

 
  printf("Oa %f %f Ba %f %f Os %f %f Bs %f %f\n",Oa[0],Oa[1],Ba[0],Ba[1],Os[0],Os[1],Bs[0],Bs[1]); 
  /*** (3) Calculated rotated addR XY, and copy into mol->atoms[] **/
  /***  add_pos_new = R(theta) * (add_pos_old - Oa) + Os **/

  th = myacos(Bs[0]*Ba[0] + Bs[1]*Ba[1]);
  if ((Ba[0]*Bs[1]-Ba[1]*Bs[0])<0.0) th = -th;

  c_t = cos(th);
  s_t = sin(th);
  printf("th %f c_t %f s_t %f\n",th,c_t,s_t); 
  make_ring_to_unit_index(ring2unit,addR,unit);
  for (i=0;i<addR->Natom;++i){
    dpos[0] = addR->pos[i][0] - Oa[0];
    dpos[1] = addR->pos[i][1] - Oa[1];
    addR->pos[i][0] = (c_t * dpos[0] - s_t * dpos[1]) + Os[0];
    addR->pos[i][1] = (s_t * dpos[0] + c_t * dpos[1]) + Os[1];

    if ((i!=ishare_a[0])&&(i!=ishare_a[1])){
      unit->pos[ring2unit[i]][0] = addR->pos[i][0];
      unit->pos[ring2unit[i]][1] = addR->pos[i][1];
    }
  }

  free(ring2unit);
  return(1);

} /* end of combine_new_ring_into_seed_ring() */




void make_ring_to_unit_index(ring2unit,ring,unit)
  int *ring2unit;
  struct UNIT *ring;
  struct UNIT *unit;
{
  int i,j;

  for (i=0;i<ring->Natom;++i){ 
    for (j=0;j<unit->Natom;++j){ 
     if (ring->num_atoms[i] == unit->num_atoms[j]) ring2unit[i] = j;
    }
  } 
}



int find_shared_atoms_bwn_two_rings(rn,sn,Nshare,ashare,ishare_r,ishare_s)
  struct UNIT *rn,*sn;      /* rings */
  int *Nshare;              /* Number of common atoms */
  int *ashare;              /* array of atom number   for shared atoms */  
  int *ishare_r, *ishare_s; /* array of index in ring for shared atoms */
{
 int i,j;

  *Nshare = 0;
  for (i=0;i<rn->Natom;++i){
    for (j=0;j<sn->Natom;++j){
      if (rn->num_atoms[i]==sn->num_atoms[j]){
        ashare[*Nshare]   = rn->num_atoms[i];
        ishare_r[*Nshare] = i;
        ishare_s[*Nshare] = j;
        *Nshare += 1;
       }
    }
  }

 return(*Nshare);

} /* end of find_shared_atoms_bwn_two_rings() */





int set_XY_to_chain_UNIT(mol,unit)
  struct MOLECULE *mol;
  struct UNIT *unit;
{
  int i;
  float x[2];
  float th,c_t,s_t,Lbond;
/*
  printf("#set_XY_to_chain_UNIT(mol,unit[%d] type %c Natom %d)\n",unit->num,unit->type,unit->Natom);
 */
  Lbond = PAR.LengthOfBond2D;

  if (unit->type != 'C') return(0);

  unit->pos[0][0] = 0.0;
  unit->pos[0][1] = 0.0;

  if (unit->Natom ==1) return(1);
  
  if ((mol->atoms[unit->num_atoms[0]].Nnei_heavy==1)||(mol->atoms[unit->num_atoms[1]].Nnei_heavy==1))
     /* unit->pos[1][0] = Lbond * 0.7; */
     unit->pos[1][0] = Lbond; 
  else
     unit->pos[1][0] = Lbond;


  unit->pos[1][1] = 0.0;


  /** Zig-zag conformation **/
  for (i=2;i<unit->Natom;++i){
    x[0] = unit->pos[i-2][0] - unit->pos[i-1][0];
    x[1] = unit->pos[i-2][1] - unit->pos[i-1][1];
    normalize_2Dvec(x);
    /* if (mol->atoms[unit->num_atoms[i-1]].Nnei_heavy==4){ th = M_PI*160.0/180.0; }   */
     if (mol->atoms[unit->num_atoms[i-1]].Nnei_heavy==4){ th = M_PI*140.0/180.0; }   
    else th = M_PI*120.0/180.0;
    if ((i%2)==1){ th = -th;}

    c_t = cos(th);  s_t = sin(th);  
    unit->pos[i][0] =  Lbond*(c_t*x[0] - s_t*x[1]) + unit->pos[i-1][0];
    unit->pos[i][1] =  Lbond*(s_t*x[0] + c_t*x[1]) + unit->pos[i-1][1];
    /* printf("[%d] %f %f\n",i,unit->pos[i][0], unit->pos[i][1]); */
  }
  
  unit->set_localXY = 'T';
/*
  for (i=0;i<unit->Natom;++i){
    printf("#unit[%d %c %d] %d %s %f %f\n",unit->num,unit->type,unit->Natom,
    mol->atoms[unit->num_atoms[i]].num_in_file,
    mol->atoms[unit->num_atoms[i]].atomname,
    unit->pos[i][0],unit->pos[i][1]);
  }
*/
  return(1);

} /* end of set_XY_to_chain_UNIT() */




int mark_neighbors_for_smallest_cycle_in_UNIT(HeadRing,mol,unit,Natom_path,path,mark,Natom_pathsmlt, pathsmlt,ns,na)
  struct UNIT *HeadRing; 
  struct MOLECULE *mol;
  struct UNIT *unit;
  int   Natom_path;  /* Natom of the path */
  int  *path;        /* path of the atom number */
  char *mark;        /* mark of the atom number */
  int  *Natom_pathsmlt;  /* Number of atoms for the smallest ring */ 
  int  *pathsmlt;   /* path of the smallest ring */ 
  int  ns;        /* atom number of the start atom   (0..unit->Natom-1) */
  int  na;        /* atom number of the focused atom (0..unit->Natom-1) */
{
  int b,nb,i;
/*
  printf("#mark_neighbors_for_smallest_cycle_in_UNIT(s:%d a:%d Natom_path:%d Natom_pathsmlt:%d)\n", mol->atoms[ns].num_in_file,mol->atoms[na].num_in_file,Natom_path, *Natom_pathsmlt);
 */
  mark[na]       = 1;
  path[Natom_path-1] = na;

  if ((*Natom_pathsmlt<0)||(Natom_path <= *Natom_pathsmlt)){
    for (b=0;b<unit->Natom;++b){
      nb = unit->num_atoms[b];
      if ((na!=nb)&&(mol->conmap.map[na][nb] != '0')){
        if ((nb==ns)&&(Natom_path>=3)){
/*
          printf("#Find Smallest Ring!!\n");
          for (i=0;i<Natom_path;++i) { printf(" %d",mol->atoms[path[i]].num_in_file); } printf("\n");
*/
          *Natom_pathsmlt = Natom_path;
          for (i=0;i<Natom_path;++i) pathsmlt[i] = path[i]; 
        }
        else if (mark[nb]==0){
           mark_neighbors_for_smallest_cycle_in_UNIT(HeadRing,mol,unit,Natom_path+1,path,mark,Natom_pathsmlt, pathsmlt,ns,nb);
        }
      }
    }
  } 
  mark[na] = 0;
  return(0);

} /* end of mark_neighbors_for_smallest_cycle_in_UNIT() */


void add_new_UNIT_to_tail(HeadUnit, Natom_new, num_atoms_new, type_new)
  struct UNIT *HeadUnit;
  int Natom_new;
  int *num_atoms_new;
  char type_new;
{
  struct UNIT *rn;
  int i;

/*
  printf("#add_new_UNIT_to_tail(HeadUnit, Natom_new %d num_atoms_new, type_new %c)\n",Natom_new,type_new);
  for (i=0;i<Natom_new;++i){
    printf("#num_atoms_new %d %d\n",i,num_atoms_new[i]);
  }
*/
  rn = HeadUnit;
  while (rn->next != NULL){ 
    rn = rn->next;
  }

  rn->next = (struct UNIT*)malloc(sizeof(struct UNIT));
  rn->next->prev = rn;
  rn->next->next = NULL;
  rn->next->num = rn->num + 1;
  rn = rn->next;
  rn->Natom = Natom_new;
  rn->num_atoms = (int *)malloc(sizeof(int)*Natom_new);
  rn->Nneighbor = (int *)malloc(sizeof(int)*Natom_new);
  rn->neighbors = (int **)malloc(sizeof(int*)*Natom_new);
  rn->pos       = (float **)malloc(sizeof(float *)*Natom_new);

  for (i=0;i<Natom_new;++i){
    rn->num_atoms[i] = num_atoms_new[i];
    rn->neighbors[i] = (int *)malloc(sizeof(int)*Natom_new);
    rn->pos[i] = (float *)malloc(sizeof(float)*2);
    rn->pos[i][0] = 0.0;
    rn->pos[i][1] = 0.0;
  }
  rn->type = type_new;
  rn->set_localXY  = 'F';
  rn->set_globalXY = 'F';

} /* end of add_new_UNIT_to_tail() */


void setup_neighbors_for_UNITs(mol,HeadUnit)
  struct MOLECULE *mol;
  struct UNIT *HeadUnit;
{
  struct UNIT *unit;
  int i,j;
  
  unit = HeadUnit;

  while (unit->next != NULL){
    unit = unit->next;
    for (i=0;i<unit->Natom;++i){
      unit->Nneighbor[i] = 0;
      for (j=0;j<unit->Natom;++j){
        if (mol->conmap.map[unit->num_atoms[i]][unit->num_atoms[j]]!='0'){
          unit->neighbors[i][unit->Nneighbor[i]] = j;
          unit->Nneighbor[i] += 1;
        }
      }
   /* 
    printf("[%d] Nneighbor %d nei ",mol->atoms[unit->num_atoms[i]].num_in_file, unit->Nneighbor[i]);
    for (j=0;j<unit->Nneighbor[i];++j){ printf(" %d",mol->atoms[unit->num_atoms[unit->neighbors[i][j]]].num_in_file); }
    printf("\n");
   */
    }

  }



} /* end of setup_neighbors_for_UNIT() */



void show_UNITs(mol, HeadUnit,comment)
  struct MOLECULE *mol;
  struct UNIT *HeadUnit;
  char *comment;
{
  struct UNIT *rn;

  if (comment[0]!='\0'){ printf("#%s\n",comment); }
  printf("#");
  rn = HeadUnit;
  while (rn->next != NULL){
    rn = rn->next;
    show_UNIT(mol,rn,"");
  }
  printf("\n");

} /* end of show_UNITs() */


void show_UNIT(mol,un,comment)
  struct MOLECULE *mol;
  struct UNIT *un;
  char *comment;
{
  struct ATOM *an;
  int a,aa;
  printf("#>%s[%d %c] loc %c glob %c Natm %d (",comment,un->num,un->type,un->set_localXY,un->set_globalXY,un->Natom);
  for (a=0;a<un->Natom;++a){
    aa = un->num_atoms[a];
    if ((aa>=0)&&(aa<mol->Natom)){
      an = &(mol->atoms[aa]);
      printf(" %c%d",an->one_char_ele,an->num_in_file);
    }
    else {printf("#IMPROPOER a %d aa %d\n",a,aa);}
 }
 printf(")\n");
} /* end of show_UNIT() */



int Number_of_UNITs(HeadUnit)
  struct UNIT *HeadUnit;
{
  struct UNIT *rn;
  int n;

  n = 0;
  rn = HeadUnit;
  while (rn->next != NULL){
    rn = rn->next;
    n += 1;
  }
  return(n);
}


void normalize_2Dvec(x)
  float *x;
{
  float len;
  len = x[0]*x[0] + x[1]*x[1];
  if (len>0.0){ 
    len = sqrt(len);
    x[0] = x[0]/len; x[1] = x[1]/len;
   }
}


void congestion_score(mol,score,Natmpair_clash)
  struct MOLECULE *mol;
  float  *score;
  int    *Natmpair_clash;
{
  int i,j;
  float dx,dy,dd,Dclash_thre,DDclash_thre;

  *score = 0.0;
  *Natmpair_clash = 0;
 /*
  Dclash_thre = PAR.LengthOfBond2D * 0.5;
 */ 
  Dclash_thre = PAR.LengthOfBond2D * 0.6;
  DDclash_thre = Dclash_thre*Dclash_thre;

  for (i=0;i<mol->Natom;++i){
    if ((mol->atoms[i].one_char_ele != 'H')&&(mol->atoms[i].mark==1)){
      for (j=i+1;j<mol->Natom;++j){
        if ((mol->atoms[j].one_char_ele != 'H')&&(mol->atoms[j].mark==1)){
          dx =  mol->atoms[i].Pos[0] - mol->atoms[j].Pos[0];
          dy =  mol->atoms[i].Pos[1] - mol->atoms[j].Pos[1];
          dd = dx*dx + dy*dy;
          if (dd==0.0){ dd = 0.001;}
          *score += 1.0/dd; 
          if (dd <= DDclash_thre) *Natmpair_clash += 1;
       }
     } 
   }
 }

} /* end of congestion_score() */




float point_congestion_score(x,y,mol)
  float  x,y;   /* given point */
  struct MOLECULE *mol;
/*
 [CAUTION] Only calcuated for "mark==1" atoms !! 
 This is mainly for deciding position of hydrogens. 
*/
{
  int i;
  float dx,dy,dd,score;

  score = 0.0;

  for (i=0;i<mol->Natom;++i){
    if ((mol->atoms[i].one_char_ele != 'H')&&(mol->atoms[i].mark==1)){
      dx =  mol->atoms[i].Pos[0] - x;
      dy =  mol->atoms[i].Pos[1] - y;
      dd = dx*dx + dy*dy;
      if (dd==0.0){ dd = 0.001;}
      score += 1.0/dd; 
   }
 }
 return(score);
} /* end of point_congestion_score() */




float congestion_in_ring_for_dir(unit,inumO,dir_x,dir_y)
  struct UNIT *unit;
  int  inumO; /* atom index number for the origin atom */ 
  float dir_x, dir_y;
{
  int i;
  float op[2],dis,score;
  float prob_pos[2]; /* position of probe atom in (dir_x,dir_y) */

  score = 0.0;
  prob_pos[0] = unit->pos[inumO][0] + PAR.LengthOfBond2D * dir_x;
  prob_pos[1] = unit->pos[inumO][1] + PAR.LengthOfBond2D * dir_y;
  for (i=0;i<unit->Natom;++i){
    if (i!=inumO){
      op[0] = unit->pos[i][0] - prob_pos[0];
      op[1] = unit->pos[i][1] - prob_pos[1];
      dis = sqrt(op[0]*op[0] + op[1]*op[1]);
      normalize_2Dvec(op);
      score += 1.0/dis/dis; 
    } 
 }

 return(score);
} /* end of congestion_in_ring_for_dir() */





float congestion_in_marked_atoms_for_dir(mol,anumO,dir_x,dir_y,OnlyRingCrash)
  struct MOLECULE *mol;
  int    anumO; /* atom number for the origin atom */ 
  float dir_x, dir_y;
  char  OnlyRingCrash; /* if T, care only clashes with "small" ring (ring=='r') */
{
  int i;
  float op[2],dprod,dis,score,w;
  float prob_pos[2]; /* position of probe atom in (dir_x,dir_y) */
  float Wbonded_neighbor3;  /* weight for bonded neighbor atoms */
  struct RING *rn;
  char hit;


  if (OnlyRingCrash=='T'){
    for (i=0;i<mol->Natom;++i){ mol->atoms[i].subnum = 0;}
    rn = &(mol->HeadRing);
    while (rn->next != NULL){
      rn = rn->next;
      hit = 0;
      for (i=0;i<rn->Natom;++i){
        if (rn->num_atoms[i]==anumO){ hit = 1;}
      }
      if (hit==1){
        for (i=0;i<rn->Natom;++i){
          mol->atoms[rn->num_atoms[i]].subnum = 1;
        }
     }
    }
    /*
    for (i=0;i<mol->Natom;++i) {
      if (mol->atoms[i].subnum>0){
        printf("%s%d subnum %d\n",mol->atoms[i].atomname,mol->atoms[i].num_in_file,mol->atoms[i].subnum);
      }
    }
    */
  }
  

  Wbonded_neighbor3 = 2.0;
  score = 0.0;
  prob_pos[0] = mol->atoms[anumO].Pos[0] + PAR.LengthOfBond2D * dir_x;
  prob_pos[1] = mol->atoms[anumO].Pos[1] + PAR.LengthOfBond2D * dir_y;
  for (i=0;i<mol->Natom;++i){
    if ((i!=anumO) && (mol->atoms[i].mark==1) && ((OnlyRingCrash!='T')||(mol->atoms[i].subnum==1))){
      op[0] = mol->atoms[i].Pos[0] - prob_pos[0];
      op[1] = mol->atoms[i].Pos[1] - prob_pos[1];
      dis = sqrt(op[0]*op[0] + op[1]*op[1]);
      normalize_2Dvec(op);
      dprod = op[0]*dir_x + op[1]*dir_y;
      w = 1.0;
      if (mol->topodismap.map[anumO][i]<=3) w = Wbonded_neighbor3;
      score += w/dis/dis; 
    } 
  }

 return(score);
 
} /* end of congestion_in_marked_atoms_for_dir() */




void cal_mid_point(bdir,o,a,b,c)
  float *bdir;      /* [2] dir vector for mid point na-o-nb */
  float *o;            /* [2] xy of origin (given) */
  float *a,*b;      /* [2] xy of neighbors a,b,c (given)*/
{
  /* 
       m  
       | 
  b----o---a

  */

  /* printf("#cal_mid_point(theta, bdir,o,na,nb)\n"); */
  bdir[0] = (a[0]+b[0])/2.0 - o[0];
  bdir[1] = (a[1]+b[1])/2.0 - o[1];
  normalize_2Dvec(bdir);
} /* end of cal_mid_point() */




void cal_angular_mid_point(theta, bdir,o,a,b,c)
  float *theta;        /* andle between na-o-nb (to be calculated) */
  float *bdir;      /* [2] dir vector for mid point na-o-nb */
  float *o;            /* [2] xy of origin (given) */
  float *a,*b,*c;   /* [2] xy of neighbors a,b,c (given)*/
{
  float da[2],db[2],dc[2],th;
  /* 
  case :((a+b)/2-o) *(c-o) < 0
        a 
        /
  b----o
       |
       c 

  case :((a+b)/2-o) *(c-o) > 0
      a 
     /
    / 
   o---c
   |
   |
   b 
  */

  da[0] = a[0] - o[0];
  da[1] = a[1] - o[1];
  normalize_2Dvec(da);

  db[0] = b[0] - o[0];
  db[1] = b[1] - o[1];
  normalize_2Dvec(db);

  dc[0] = c[0] - o[0];
  dc[1] = c[1] - o[1];
  normalize_2Dvec(dc);

 /*
  bdir[0] = (a[0]+b[0])/2.0 - o[0];
  bdir[1] = (a[1]+b[1])/2.0 - o[1];
  */ 
  bdir[0] = (da[0]+db[0])/2.0;
  bdir[1] = (da[1]+db[1])/2.0;

  normalize_2Dvec(bdir);
  
  th = myacos(da[0]*db[0] + da[1]*db[1]);

  printf("#cal_angular_mid_point(theta, bdir,o,na,nb) bdir*dc %f\n",bdir[0]*dc[0]+bdir[1]*dc[1]);
  if ((bdir[0]*dc[0]+bdir[1]*dc[1])>0.0){
    th = 2.0*M_PI - th;
    bdir[0] = -bdir[0];
    bdir[1] = -bdir[1];
    printf("#reverse_direction\n");
  } 

  printf("#th %f o %f %f a %f %f b %f %f bdir %f %f\n",th,o[0],o[1],a[0],a[1],b[0],b[1],bdir[0],bdir[1]); 

  *theta = th;

} /* end of cal_angular_mid_point() */




void cal_angular_three_division_point(bdir,o,a,b)
  float bdir[4][2];   /* [4][2] dir vector for 1st three_division point na-o-nb */
  float *o;      /* [2] xy of origin (given) */
  float *a,*b;   /* [2] xy of neighbors a,b,c (given)*/
{
  float oa[2],ob[2],th,dth,c_t,s_t; 
  int i;
  /* 
   b2  b1  a 
       | /
  b----o
  */

  oa[0] = a[0] - o[0];
  oa[1] = a[1] - o[1];
  normalize_2Dvec(oa);
  ob[0] = b[0] - o[0];
  ob[1] = b[1] - o[1];
  normalize_2Dvec(ob);

  th = myacos(oa[0]*ob[0] + oa[1]*ob[1]);
  if ((oa[0]*ob[1]-oa[1]*ob[0])>0.0) th = -th; 
 
  for (i=0;i<4;++i){ 
         if (i==0) dth = th/3.0;
    else if (i==1) dth = 2.0*th/3.0;
    else if (i==2) dth = (-2*M_PI+th)/3.0;
    else if (i==3) dth = 2.0*(-2*M_PI+th)/3.0;
    else dth = th;
    c_t = cos(dth); 
    s_t = sin(dth); 
    bdir[i][0] = c_t * oa[0] - s_t * oa[1];
    bdir[i][1] = s_t * oa[0] + c_t * oa[1];
  }

} /* end of cal_angular_three_division_point() */









int Collision_Resolving(mol,HeadUnit)
  struct MOLECULE *mol; 
  struct UNIT *HeadUnit;
{
 int i,j,k,r,Nnonring_nei;
 struct UNSIGNED_CHAR2DMAP  ringbond; /* mark only bond in the "large"(>6-membered) ring */
 int NROUND,NROUND_MAX;
 int Nrepeat,Nrepeat_max;
 int *inv_bond_candO; /* Origin atom number fo inv_bond candidate */
 int *inv_bond_candE; /* Endatom number fo inv_bond candidate */
 int Ninv_bond_cand;
 float score0,score1;
 int Nclash0,Nclash1;
 float *initX,*initY; /* initial  XYs before Collision_Resolving */ 
 float *origX,*origY; /* original XYs befere inversion for each step */  
 float *newX, *newY;  /* new      XYs after  inversion for each step */
 char ofname[128];

 printf("#Collision_Resolving(mol)\n");
 /* [1] Finding clashing pairs, and assign candidate bonds for inversion **/
 Malloc_UNSIGNED_CHAR2DMAP(&ringbond,mol->Natom);
 inv_bond_candO = (int *)malloc(sizeof(int)*mol->Natom*mol->Natom);
 inv_bond_candE = (int *)malloc(sizeof(int)*mol->Natom*mol->Natom);
 
 initX = (float*)malloc(sizeof(float)*(mol->Natom+1));
 initY = (float*)malloc(sizeof(float)*(mol->Natom+1));
 newX  = (float*)malloc(sizeof(float)*(mol->Natom+1));
 newY  = (float*)malloc(sizeof(float)*(mol->Natom+1));
 origX = (float*)malloc(sizeof(float)*(mol->Natom+1));
 origY = (float*)malloc(sizeof(float)*(mol->Natom+1));

 /** [1] Preparing inv_bond map and  **/
 /* ringbond.map[i][j]=='r' only if bond(i,j) is in the "small"(<=6-membered) ring */
 /* ringbond.map[i][j]=='R' only if bond(i,j) is in the "large"(>6-membered) ring */
 
 for (i=0;i<mol->Natom;++i){
   for (j=0;j<mol->Natom;++j){
     ringbond.map[i][j] = '0';
     if (mol->conmap.map[i][j]!='0'){
            if ((mol->atoms[i].ring=='r')&&(mol->atoms[j].ring=='r')) ringbond.map[i][j] = 'r';
       else if ((mol->atoms[i].ring!=' ')&&(mol->atoms[j].ring!=' ')) ringbond.map[i][j] = 'R';
     }
   }
 }


 /** [2] Preparing inv_bond_candO[] and inv_bond_candE[]  **/
 Ninv_bond_cand = 0;
 for (i=0;i<mol->Natom;++i){
   for (j=0;j<mol->Natom;++j){
     if ((i!=j)&&(mol->conmap.map[i][j]!='0')){
       if ((ringbond.map[i][j]=='0')&&(mol->atoms[i].Nnei_heavy>=2)&&(mol->atoms[i].ring ==' ')){ 
         inv_bond_candO[Ninv_bond_cand] = i;
         inv_bond_candE[Ninv_bond_cand] = j;
         Ninv_bond_cand += 1;
       }
      if ((ringbond.map[i][j]=='R')&&(mol->atoms[i].Nnei_heavy>=3)&&(mol->atoms[i].ring=='R')&&(mol->atoms[j].ring=='R')){ 
          Nnonring_nei = 0;
          for (k=0;k<mol->Natom;++k){ 
             if ((k!=i)&&(k!=j)&&(mol->conmap.map[i][k]!='0')&&(ringbond.map[i][k]=='0')) Nnonring_nei += 1; 
           }
          if (Nnonring_nei>0){
            inv_bond_candO[Ninv_bond_cand] = i;
            inv_bond_candE[Ninv_bond_cand] = j;
            Ninv_bond_cand += 1;
   /*
            printf("inv O [%d %s] E [%d %s] conmap %c ringbond %c\n",
              mol->atoms[i].num_in_file, mol->atoms[i].atomname,
              mol->atoms[j].num_in_file, mol->atoms[j].atomname,mol->conmap.map[i][j],ringbond.map[i][j]);
     */
          }
       }
     }
   }
 }

 printf("#Ninv_bond_cand %d\n",Ninv_bond_cand);

 congestion_score(mol,&score0,&Nclash0);

 if (Ninv_bond_cand == 0) return(Nclash0);  /* no candidate for inversing bond !! */
 if (Nclash0==0) return(Nclash0);
 
 /*
 for (k=0;k<Ninv_bond_cand;++k){
   i = inv_bond_candO[k];
   j = inv_bond_candE[k];
   printf("#inv_bond_candidate O(%d %s)--- E(%d %s) '%c'\n",
     mol->atoms[i].num_in_file, mol->atoms[i].atomname,
     mol->atoms[j].num_in_file, mol->atoms[j].atomname,mol->conmap.map[i][j]);
 }
 */

 /* [3] Repeating inversion operations until collision resolving **/
 NROUND_MAX = 10;  
 Nrepeat_max = 10 * mol->Nheavyatom;
 /*
 Write_SDF_Molecule("before.sdf",mol,"");
 */

 copy_posXY_from_MOLECULE(initX,initY,mol);

 NROUND = Nrepeat = 0;
 while ((Nclash0>0)&&(NROUND<NROUND_MAX)){ 
   copy_posXY_to_MOLECULE(initX,initY,mol);
   Nrepeat = 0;
   while ((Nclash0>0)&&(Nrepeat<Nrepeat_max)){
     r = rand()%Ninv_bond_cand;
     congestion_score(mol,&score0,&Nclash0);
     inversion(newX,newY,mol,inv_bond_candO[r],inv_bond_candE[r],&ringbond);
     copy_posXY_from_MOLECULE(origX,origY,mol);
     copy_posXY_to_MOLECULE(newX,newY,mol);
     congestion_score(mol,&score1,&Nclash1);
     /* 
     printf("Nrepeat %d r %d bond %d %s --- %d %s ",Nrepeat,r,
       mol->atoms[inv_bond_candO[r]].num_in_file, mol->atoms[inv_bond_candO[r]].atomname, 
       mol->atoms[inv_bond_candE[r]].num_in_file, mol->atoms[inv_bond_candE[r]].atomname);
     printf("score0 %f score1 %f Nclash0 %d Nclash1 %d",score0,score1,Nclash0,Nclash1);
     */
    if ((Nrepeat%500)==0){ 
       printf("#ROUND %d rep %d score0 %f score1 %f Nclash0 %d Nclash1 %d\n",NROUND,Nrepeat,score0,score1,Nclash0,Nclash1);
       sprintf(ofname,"%d-%d.sdf",NROUND,Nrepeat);
       /* Write_SDF_Molecule(ofname,mol,"");    */
    }
     /*
     printf("inv O [%d %s] E [%d %s] Nclash0 %d Nclash1 %d score0 %f score1 %f\n",
       mol->atoms[inv_bond_candO[r]].num_in_file, mol->atoms[inv_bond_candO[r]].atomname,
       mol->atoms[inv_bond_candE[r]].num_in_file, mol->atoms[inv_bond_candE[r]].atomname,Nclash0,Nclash1,score0,score1);
     */

     /** ACCEPT **/ 
     if ((Nclash1<Nclash0)||(score1<score0)){
      score0 = score1;
      Nclash0 = Nclash1;
      /*
      if (check_bond_distance(mol)==0){
        printf("WOOPS!! inv O [%d %s] E [%d %s] Nclash0 %d Nclash1 %d score0 %f score1 %f\n",
         mol->atoms[inv_bond_candO[r]].num_in_file, mol->atoms[inv_bond_candO[r]].atomname,
         mol->atoms[inv_bond_candE[r]].num_in_file, mol->atoms[inv_bond_candE[r]].atomname,Nclash0,Nclash1,score0,score1);
      } 
       */
     }
     /** REJECT **/ 
     else{
       copy_posXY_to_MOLECULE(origX,origY,mol);
       congestion_score(mol,&score1,&Nclash1);
     }



     Nrepeat += 1;
   }

   NROUND += 1;
 }

 printf("#FINAL_ROUND %d rep %d score0 %f Nclash0 %d\n",NROUND,Nrepeat,score0,Nclash0);

 free(origX); 
 free(origY); 
 free(newX); 
 free(newY); 
 free(initX); 
 free(initY); 
 free(inv_bond_candO);
 free(inv_bond_candE);
 Free_UNSIGNED_CHAR2DMAP(&ringbond);

 return(Nclash0);

} /* end of Collision_Resolving() */




void copy_posXY_from_MOLECULE(posX,posY,mol)
  float *posX,*posY;
  struct MOLECULE *mol;
{
 int i;
 for (i=0;i<mol->Natom;++i){
   posX[i] = mol->atoms[i].Pos[0];
   posY[i] = mol->atoms[i].Pos[1];
 }
} /* end of copy_posXY_from_MOLECULE() */


void copy_posXY_to_MOLECULE(posX,posY,mol)
  float *posX,*posY;
  struct MOLECULE *mol;
{
 int i;
 for (i=0;i<mol->Natom;++i){
   mol->atoms[i].Pos[0] = posX[i];
   mol->atoms[i].Pos[1] = posY[i];
 }
 } /* end of copy_posXY_to_MOLECULE() */


int inversion(posnewX, posnewY, mol, Onum, Enum, ringbond)
  float *posnewX; /* [mol->Natom] */ 
  float *posnewY; /* [mol->Natom] */ 
  struct MOLECULE *mol;
  int Onum;   /* atom number for the Origin atom */
  int Enum;   /* atom number for the other terminal atom of the focused bond */
  struct UNSIGNED_CHAR2DMAP *ringbond;  /* [mol->Natom] if atoms i-j are connected in a ring, ringbond.map[i][j] = 'r' */
{
  int Nnei;
  int i,k,ii,nei0,nei1;
  float vecE[2],dpos[2],th,c_t,s_t,dotAE;
  float vecNei[10][2];
  int   anumNei[10];  /* atom number of neighboring heavy atoms */
  char ChangeType,OEbond; 

 /*
  printf("#inversion(O:%s%d E%s%d)\n",
    mol->atoms[Onum].atomname,mol->atoms[Onum].num_in_file,
    mol->atoms[Enum].atomname,mol->atoms[Enum].num_in_file);
  fflush(stdout);
 */

 /*
 
  Atoms O(ooh)  and E are fixed.  Atom O-E may be part of the ring.
 
  Atoms 0(zero) and 1 and 2 are candidates of rotations.
  Bond O-0,O-1,O-2 should be non-ring bond.

  2  1 
   \ |
     O-- 0
     |
     E
  */


  /** [1] Prepare anumNei[]  **/
  for (i=0;i<mol->Natom;++i){
    posnewX[i] = mol->atoms[i].Pos[0];
    posnewY[i] = mol->atoms[i].Pos[1];
  }
  if (ringbond->map[Onum][Enum]=='0') OEbond = '-'; else OEbond = 'R';

  Nnei = 0;
  for (k=0;k<10;++k) anumNei[k] = -1;
  
  for (i=0;i<mol->atoms[Onum].Nneighbor;++i){
    ii = mol->atoms[Onum].neighbors[i];
    if ((ii!=Enum)&&(mol->atoms[ii].one_char_ele != 'H')&&(ringbond->map[ii][Onum]=='0')){
      if (Nnei<3) anumNei[Nnei] = ii;
      Nnei += 1;
      if (Nnei>=10) return(0);
    }
  }

  for (k=0;k<2;++k){ vecE[k] = mol->atoms[Enum].Pos[k] - mol->atoms[Onum].Pos[k]; }
  normalize_2Dvec(vecE);
   
  for (i=0;i<Nnei;++i){
    if (anumNei[i]>=0){
      for (k=0;k<2;++k){
        vecNei[i][k] = mol->atoms[anumNei[i]].Pos[k] - mol->atoms[Onum].Pos[k];
      }
    }
    normalize_2Dvec(vecNei[i]);
  }

  /*
  printf("#inversion(mol,Onum:%d %s Enum:%d %s) Nnei %d\n",
  mol->atoms[Onum].num_in_file, mol->atoms[Onum].atomname,
  mol->atoms[Enum].num_in_file, mol->atoms[Enum].atomname,Nnei); fflush(stdout);
  */

  /** [2-1] Case for Nnei = 2 (Nnei_heavy == 3) (0 --> 1, 1 --> 0)**/
  /*
     1 
     |
     O-- 0
     |
     E
  */
  if (Nnei==2){  

    th = myacos(vecNei[0][0]*vecNei[1][0] + vecNei[0][1]*vecNei[1][1]);
    if ((vecNei[0][0]*vecNei[1][1]-vecNei[0][1]*vecNei[1][0])<0.0) th = -th;
    /* printf("#two_neighbor_rotation th %f\n",th); */
    c_t = cos(th);
    s_t = sin(th);
    for (i=0;i<mol->Natom;++i){
      if (mol->atoms[i].one_char_ele != 'H'){
        if (mol->topodismap.map[anumNei[0]][i]<mol->topodismap.map[Onum][i]){
/*        printf("#inverse2 %d %s atomO %d %s\n", mol->atoms[i].num_in_file, mol->atoms[i].atomname,mol->atoms[Onum].num_in_file, mol->atoms[Onum].atomname);  */
          dpos[0] = mol->atoms[i].Pos[0] - mol->atoms[Onum].Pos[0];
          dpos[1] = mol->atoms[i].Pos[1] - mol->atoms[Onum].Pos[1];
          posnewX[i] = ( c_t * dpos[0] - s_t * dpos[1]) + mol->atoms[Onum].Pos[0];
          posnewY[i] = ( s_t * dpos[0] + c_t * dpos[1]) + mol->atoms[Onum].Pos[1];
        }
        else if (mol->topodismap.map[anumNei[1]][i]<mol->topodismap.map[Onum][i]){
/*        printf("#inverse2 %d %s atomO %d %s\n", mol->atoms[i].num_in_file, mol->atoms[i].atomname,mol->atoms[Onum].num_in_file, mol->atoms[Onum].atomname);  */
          dpos[0] = mol->atoms[i].Pos[0] - mol->atoms[Onum].Pos[0];
          dpos[1] = mol->atoms[i].Pos[1] - mol->atoms[Onum].Pos[1];
          posnewX[i] = ( c_t * dpos[0] + s_t * dpos[1]) + mol->atoms[Onum].Pos[0];
          posnewY[i] = (-s_t * dpos[0] + c_t * dpos[1]) + mol->atoms[Onum].Pos[1];
        }
      }
    }
  }

  /** [2-2] Case for Nnei = 1 (Nnei_heavy == 2) (A --> A') **/
  /*
    (0') 
     |
     O-- 0
     |
     E
  */
 if (Nnei==1){  
    dotAE = vecNei[0][0]*vecE[0] + vecNei[0][1]*vecE[1];
    if ((fabs(dotAE)>0.01)&&(fabs(dotAE)<0.99)){
      vecNei[1][0] = 2*dotAE*vecE[0] - vecNei[0][0];
      vecNei[1][1] = 2*dotAE*vecE[1] - vecNei[0][1];

      th  = myacos(vecNei[0][0]*vecNei[1][0] + vecNei[0][1]*vecNei[1][1]);
      if ((vecNei[0][0]*vecNei[1][1]-vecNei[0][1]*vecNei[1][0])<0.0) th = -th;
      if (OEbond=='R') th = M_PI;

      /* printf("#one_neighbor_rotation th %f\n",th); */
      c_t = cos(th);
      s_t = sin(th);

      for (i=0;i<mol->Natom;++i){
        if (mol->atoms[i].one_char_ele != 'H'){
          if (mol->topodismap.map[anumNei[0]][i]<mol->topodismap.map[Onum][i]){
            dpos[0] = mol->atoms[i].Pos[0] - mol->atoms[Onum].Pos[0];
            dpos[1] = mol->atoms[i].Pos[1] - mol->atoms[Onum].Pos[1];
            posnewX[i] = ( c_t * dpos[0] - s_t * dpos[1]) + mol->atoms[Onum].Pos[0];
            posnewY[i] = ( s_t * dpos[0] + c_t * dpos[1]) + mol->atoms[Onum].Pos[1];
          }
        }
      }
    }
  }

  if (Nnei==3){ 
  /** [2-3] Case for Nnei = 3 (Nnei_heavy == 4) **/
  /*
     1 
     |
  2--O-- 0
     |
     E

  Three types of change are possible:
    ChangeType 0 : 0 <--> 1
    ChangeType 1 : 0 <--> 2
    ChangeType 2 : 1 <--> 2
  In this function, one of the three is randomly chosen.

  */
    ChangeType = rand()%3;
         if (ChangeType==0) {nei0 = 0; nei1 = 1;}
    else if (ChangeType==1) {nei0 = 0; nei1 = 2;}
    else                    {nei0 = 1; nei1 = 2;}


    th = myacos(vecNei[nei0][0]*vecNei[nei1][0] + vecNei[nei0][1]*vecNei[nei1][1]);
    if ((vecNei[nei0][0]*vecNei[nei1][1]-vecNei[nei0][1]*vecNei[nei1][0])<0.0) th = -th;
    /* printf("#three_neighbor_rotation th %f\n",th); */
    c_t = cos(th);
    s_t = sin(th);

    for (i=0;i<mol->Natom;++i){
      if (mol->atoms[i].one_char_ele != 'H'){
        if (mol->topodismap.map[anumNei[nei0]][i]<mol->topodismap.map[Onum][i]){
  /*      printf("#inverse2 %d %s atomO %d %s\n", mol->atoms[i].num_in_file, mol->atoms[i].atomname,mol->atoms[Onum].num_in_file, mol->atoms[Onum].atomname);  */
          dpos[0] = mol->atoms[i].Pos[0] - mol->atoms[Onum].Pos[0];
          dpos[1] = mol->atoms[i].Pos[1] - mol->atoms[Onum].Pos[1];
          posnewX[i] = ( c_t * dpos[0] - s_t * dpos[1]) + mol->atoms[Onum].Pos[0];
          posnewY[i] = ( s_t * dpos[0] + c_t * dpos[1]) + mol->atoms[Onum].Pos[1];
        }
        else if (mol->topodismap.map[anumNei[nei1]][i]<mol->topodismap.map[Onum][i]){
 /*       printf("#inverse2 %d %s atomO %d %s\n", mol->atoms[i].num_in_file, mol->atoms[i].atomname,mol->atoms[Onum].num_in_file, mol->atoms[Onum].atomname);  */
          dpos[0] = mol->atoms[i].Pos[0] - mol->atoms[Onum].Pos[0];
          dpos[1] = mol->atoms[i].Pos[1] - mol->atoms[Onum].Pos[1];
          posnewX[i] = ( c_t * dpos[0] + s_t * dpos[1]) + mol->atoms[Onum].Pos[0];
          posnewY[i] = (-s_t * dpos[0] + c_t * dpos[1]) + mol->atoms[Onum].Pos[1];
        }
      }
    } 
  }


  return(1);
} /* end of inversion() */







void Write_SDF_UNIT(ofname,unit,mol)
 char *ofname;
 struct UNIT    *unit;
 struct MOLECULE *mol;
{
 FILE *fpo;
 int i,j,Nbond;

 printf("#Write_SDF_UNIT()-->'%s'\n",ofname);
 fpo = fopen(ofname,"w");
 if (fpo==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }

 /** count Nbond **/
 Nbond = 0;
 for (i=0;i<unit->Natom;++i){
   for (j=i+1;j<unit->Natom;++j){
     if (mol->conmap.map[unit->num_atoms[i]][unit->num_atoms[j]]!='0'){
        Nbond += 1;
     }
   }
 }


 fprintf(fpo,"UNIT%d%c\n",unit->num,unit->type);
 fprintf(fpo,"  %s\n","UNIT");
 fprintf(fpo,"\n");
 /*  13 12  0     1  0  0  0  0  0999 V2000 */
 /*  39 41  0  0  1  0  0  0  0  0999 V2000 */
 /*  10  9  0     0  0  0  0  0  0999 V2000 */
 /* 107110  0     0  0  0  0  0  0999 V2000 */
 /*  30 31  0  0  1  0  0  0  0  0999 V2000 */
 fprintf(fpo,"%3d%3d  0  0  %c  0  0  0  0  0999 V2000\n",unit->Natom,Nbond,'0');
 /*    5.6720    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0 */

 for (i=0;i<unit->Natom;++i){
   fprintf(fpo,"%10.4f%10.4f%10.4f %-2s %2d  %d  %d  0  0  0  0  0  0  0  0  0\n",
     unit->pos[i][0],unit->pos[i][1],0.0,mol->atoms[unit->num_atoms[i]].element,
     mol->atoms[unit->num_atoms[i]].mass_diff,mol->atoms[unit->num_atoms[i]].char_charge,mol->atoms[unit->num_atoms[i]].stereo_parity);
 }
/*  1  6  1  0  0  0  0 */
 for (i=0;i<unit->Natom;++i){
   for (j=i+1;j<unit->Natom;++j){
     if (mol->conmap.map[unit->num_atoms[i]][unit->num_atoms[j]]!='0'){
       fprintf(fpo,"%3d%3d  %c  %c  0  0  0\n",
          i+1,j+1,mol->conmap.map[unit->num_atoms[i]][unit->num_atoms[j]], '0');
     }
   }
 }
 fprintf(fpo,"M  END\n");
 fclose(fpo);

} /* end of Write_SDF_UNIT() */



void Free_UNITs(HeadUnit)
 struct UNIT *HeadUnit;
{
 struct UNIT *un;
 int i;
 un = HeadUnit;
 if (un->next != NULL){
   while (un->next != NULL) un = un->next;
   while (un->prev!=NULL){
     un = un->prev;
     if (un->next != NULL){

       if (un->next->pos != NULL){
         for (i=0;i<un->next->Natom;++i){ free(un->next->pos[i]);}
         free(un->next->pos);
       } 

       if (un->next->neighbors != NULL){
         for (i=0;i<un->next->Natom;++i){ free(un->next->neighbors[i]);}
         free(un->next->neighbors);
       } 

       if (un->next->Nneighbor!=NULL) free(un->next->Nneighbor);
       if (un->next->num_atoms!=NULL) free(un->next->num_atoms);
       free(un->next);
     }
   }
   HeadUnit->next = NULL;
 }

} /* end of Free_UNITs() */


int check_bond_distance(mol)
  struct MOLECULE *mol;
{
 int i,j,k;
 float dpos[3],dis;
 for (i=0;i<mol->Natom;++i){
   if (mol->atoms[i].one_char_ele != 'H'){
     for (j=i+1;j<mol->Natom;++j){
       if ((mol->atoms[j].one_char_ele != 'H')&&(mol->conmap.map[i][j]!='0')){
          for (k=0;k<3;++k) dpos[k] = mol->atoms[i].Pos[k] - mol->atoms[j].Pos[k];
          dis = sqrt(dpos[0]*dpos[0] + dpos[1]*dpos[1] + dpos[2]*dpos[2]);
          if (dis>PAR.LengthOfBond2D*2){ return(0); } 
       }
    }
  }
 }
 return(1);
}



float myacos(dot_prod)
  float dot_prod;
{
       if (dot_prod >  1.0) return(acos( 1.0));
  else if (dot_prod < -1.0) return(acos(-1.0));
  else return(acos(dot_prod));
}



void Rescaling_Corrdinates_of_Molecule(mol)
  struct MOLECULE *mol;
{
  int i,j,k,Nhvy,Npair;
  float dif[3],gpos[3],dis,sumdis,aveLbond;

 /* [1] Calculate average distance bwn bonded heavy atoms  */
  Nhvy = Npair = 0;
  sumdis = aveLbond = 0.0;
  gpos[0] = gpos[1] = gpos[2] = 0.0;
  for (i=0;i<mol->Natom;++i){
    if (mol->atoms[i].one_char_ele != 'H'){
      for (k=0;k<3;++k) gpos[k] += mol->atoms[i].Pos[k];
      Nhvy += 1;
      for (j=i+1;j<mol->Natom;++j){
        if ((mol->atoms[j].one_char_ele != 'H')&&(mol->conmap.map[i][j] != '0')){
          sub_vec3(dif,mol->atoms[i].Pos,mol->atoms[j].Pos); 
          dis = sqrt(dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2]);
          sumdis += dis;
          Npair += 1;
        }
      }
    }
  }

  if (Nhvy >0){
    for (k=0;k<3;++k) gpos[k] = gpos[k]/Nhvy;
  }
  if (Npair>0) aveLbond = sumdis/Npair;

  printf("#gpos %f %f %f aveLbond %f\n",gpos[0],gpos[1],gpos[2],aveLbond);

 /** [2] Rescalig : xnew = LengthOfBond2D/L * (x-g)  */
 for (i=0;i<mol->Natom;++i){
   sub_vec3(dif,mol->atoms[i].Pos,gpos);
   for (k=0;k<3;++k){
     mol->atoms[i].Pos[k] = PAR.LengthOfBond2D * dif[k] /aveLbond + gpos[k];
   }
 }
} /* end of Rescaling_Corrdinates_of_Molecule() */
