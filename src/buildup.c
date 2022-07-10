/*

 <buildup.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 functions for "build-up" algorithm for MCS ditection

 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "buildup.h"
#include "MrgSrtMATCH.h"
#include "molpermu.h"


/** FUNCTIONS (GLOBAL) **/
void BuildUp();
void Incremental_Disconnected_BuildUp();
int  Number_of_Connected_Components_for_MATCH();
void Set_Number_of_Connected_Components_for_MATCHlist();
void Bubble_Sort_anumA_anumB_array_by_anumA();
void Bubble_Sort_anumA_anumB_array_by_anumB();
int  Check_NonRedundancy_of_new_MATCH_with_MATCHlist();
int  Check_Equivalence_bwn_Two_MATCHes();

/** FUNCTIONS (LOCAL) **/
static void Replace_the_last_node_with_new_MATCH();
static void Initialize_Malloced_MATCHlist();
static int Check_Atom_Overlap_bwn_Pair_and_MATCH();
static int Bond_Equivalence_for_Pair_and_MATCH();
static int Bond_Connection_for_Pair_and_MATCH();
static int Number_of_Connected_Components_for_Pair_and_MATCH();
static void Mark_Connected_Component();
static float Worst_select_dis_in_MATCH_list();
static int Check_Diff_ShortPath_bwn_Pair_and_MATCH();


void BuildUp(optMlist,pairMlist,molA,molB,Nkeep,ConnectGraph,maxDIFtopodis,maxNcomponent)
 struct MATCH *optMlist;     /* MATCHlist to be optmized (already Nkeep-MATCHes are malloced : Nmatch = smllerInt(molA.Natom,molB.Natom)) */
 struct MATCH *pairMlist; /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 struct MOLECULE *molA,*molB;
 int    Nkeep;
 char   ConnectGraph;     /* 'C'onnected,'D'isconnected,'T'opodis-constrained,'S'ubstructure,'I'somorphic */
 int    maxDIFtopodis;    /* accepted maximum difference of topological distance (number of bonds in the shortest path) */ 
 int    maxNcomponent;    /* Maximum Number of Cluster for Disconnected graph. 0:Don't care */
{
 int Npair,Npair_max;
 struct MATCH newMlist,*om,*pm,newm;
 int n,Nmatch_new,end,Ncomponent_new;
 float worst_select_dis;
 char overlap,equivalence,ring_consist,max_cluster_ok,shortest_path_ok;
 
 if (PAR.OutCalProcess=='T')
   printf("#BuildUp(Nkeep %d Len_pairMlist %d Len_optMlist %d ConnectGraph %c maxNcomponent %d)\n",
       Nkeep,Length_of_MATCHlist(pairMlist),Length_of_MATCHlist(optMlist),ConnectGraph,maxNcomponent); fflush(stdout);

 if (molA->Nheavyatom < molB->Nheavyatom) Npair_max = molA->Nheavyatom;
                                     else Npair_max = molB->Nheavyatom;
 /* printf("#Npair_max %d\n",Npair_max); */
 Malloc_MATCHlist(&newMlist, Npair_max, Nkeep);
 Malloc_MATCH(&newm,Npair_max);

 /** (1) Buildup procedures ***/
 if (optMlist->next->nodetype != 'E') Npair = optMlist->next->Npair;
 else Npair = 1;

 end = 0;
 while ((Npair<=Npair_max)&&(end==0)){
   /* printf(">round Npair %d/%d\n",Npair,Npair_max); fflush(stdout); */
   Initialize_Malloced_MATCHlist(&newMlist,Npair_max);
   worst_select_dis = -1.0; 
   Nmatch_new = 0;

   om = optMlist; /* om:old n-atom match */

   while ((om->next!=NULL)&&(om->next->nodetype!='E')){

     om = om->next;

     if ((Nmatch_new<Nkeep) || (om->select_dis <= worst_select_dis)) {

       pm = pairMlist; /* pm: pair (2-atom) match*/

       while (pm->next!=NULL){
         pm = pm->next;
         overlap = Check_Atom_Overlap_bwn_Pair_and_MATCH(pm,om);

         if (overlap==0){ 
           Ncomponent_new = equivalence =  ring_consist = max_cluster_ok = shortest_path_ok = 1; 

                if ((ConnectGraph=='D')||(ConnectGraph=='T')) 
                equivalence = Bond_Equivalence_for_Pair_and_MATCH(pm,om,molA,molB);
           else if ((ConnectGraph=='C')||(ConnectGraph=='S')||(ConnectGraph=='I')||(ConnectGraph=='t')) 
                equivalence = Bond_Equivalence_for_Pair_and_MATCH(pm,om,molA,molB)* 
                              Bond_Connection_for_Pair_and_MATCH(pm,om,molA,molB);

           if ((equivalence==1)&&((ConnectGraph=='D')||(ConnectGraph=='T')||(ConnectGraph=='t'))){ 
              if (maxNcomponent>0){
                Ncomponent_new = Number_of_Connected_Components_for_Pair_and_MATCH(pm,om,molA,molB);
                if (Ncomponent_new > maxNcomponent) max_cluster_ok = 0;
              }
              if (((ConnectGraph=='T')||(ConnectGraph=='t'))&&(maxDIFtopodis>=0)&&(maxNcomponent!=1)){
                shortest_path_ok = Check_Diff_ShortPath_bwn_Pair_and_MATCH(pm,om,molA,molB,maxDIFtopodis);
              }
          }

        }

         if ( (overlap==0) && (equivalence==1) && (ring_consist==1) && (max_cluster_ok==1) && (shortest_path_ok==1)){ 
            newm.Npair = om->Npair + 1;
            newm.Ncomponent = Ncomponent_new;
            for (n=0;n<om->Npair;++n){ 
              newm.anumA[n] = om->anumA[n]; 
              newm.anumB[n] = om->anumB[n];
            } 
            newm.anumA[om->Npair] = pm->anumA[0]; 
            newm.anumB[om->Npair] = pm->anumB[0]; 
            Bubble_Sort_anumA_anumB_array_by_anumA(&newm);
            newm.select_dis =  cal_select_dis_for_MATCH(&newm,molA,molB);
/*
            printf("#Npair %d Ncomponent_new %d Ncomponent_calcA %d B %d\n",newm.Npair,Ncomponent_new,
Number_of_Connected_Components_for_MATCH(newm.anumB,newm.Npair,molB), Number_of_Connected_Components_for_MATCH(newm.anumA,newm.Npair,molA));
*/
            if ( (Check_NonRedundancy_of_new_MATCH_with_MATCHlist(&newm,&newMlist,molA,molB)==1)  &&
                 ((ConnectGraph!='I') || (newm.select_dis==0.0)) &&
                 ((Nmatch_new<Nkeep) || (newm.select_dis <= worst_select_dis))   ){ 
                Replace_the_last_node_with_new_MATCH(&newMlist,&newm);
                Merge_Sort_Double_Linked_List_MATCH(&newMlist,'N');
                /* Print_MATCHlist(&newMlist,"newMlist(after_sort)");  */
                Nmatch_new += 1;
                worst_select_dis = Worst_select_dis_in_MATCH_list(&newMlist);
            }
         } 
       } /* pm */
      }
     } /* om */
   
   if (Nmatch_new>0){ 
     Copy_MATCHlist(optMlist,&newMlist);
   }
   else end = 1;
   Npair += 1;

 } /* while */


 Free_MATCHlist(&newMlist);
 Free_MATCH(&newm);

} /* end of BuildUp() */









void Incremental_Disconnected_BuildUp(Mlist,pairMlist,molA,molB,Nkeep,maxNcomponent)
 struct MATCH *Mlist;       /* MATCHlist to be optmized (already malloced) */
 struct MATCH *pairMlist;   /* MATCHlist of all the possible pairs (already malloced) */
 struct MOLECULE *molA,*molB;
 int    Nkeep;
 int    maxNcomponent;
{
 int r;

 if (PAR.OutCalProcess=='T')
   printf("#Incremental_Disconnected_Buildup(Nkeep %d maxNcomponent %d)\n",Nkeep,maxNcomponent);
 BuildUp(Mlist,pairMlist,molA,molB,Nkeep,'C',-1);
 /* Print_MATCHlist(Mlist,"optMlist(Connected)"); */

 if (maxNcomponent==-1) BuildUp(Mlist,pairMlist,molA,molB,Nkeep,'D',-1);
 else if (maxNcomponent>=2){
   for (r=2;r<=maxNcomponent;++r){
     BuildUp(Mlist,pairMlist,molA,molB,Nkeep,'D',r);
     /* Print_MATCHlist(Mlist,"optMlist(round:Disconnected)"); */
   }
 }

} /* end of Incremental_Disconnected_BuildUp(Mlist,pairMlist,molA,molB,Nkeep,maxNcomponent) */




void Initialize_Malloced_MATCHlist(M,Nmax_pair)
 struct MATCH *M;
 int Nmax_pair;
{
 struct MATCH *mn;
 int n;

 M->nodetype = 'I';
 mn = M;
 while (mn->next!=NULL){
   mn = mn->next;
   mn->nodetype = 'E';
   mn->select_dis = Nmax_pair * 1000.0;
   for (n=0;n<mn->Npair_malloc;++n) mn->anumA[n] = mn->anumB[n] = -1;
   mn->Npair = 0;
 }

} /* end of Initialize_Malloced_MATCHlist() */



void Replace_the_last_node_with_new_MATCH(Mlist,newm)
 struct MATCH *Mlist;
 struct MATCH *newm;
{
 struct MATCH *mn,*last;
 int n;

 /*
 printf("#Replace_the_last_node_with_new_MATCH(Mlist,newm)\n");
 */
 /** [1] Find the last node ("last") **/
 mn = Mlist;
 while ((mn->next!=NULL) && (mn->next->nodetype != 'E')){
   mn = mn->next;
 }

 if (mn->next==NULL) {
    last = mn; 
  }
 else { 
  last = mn->next; 
  last->nodetype = 'N';
  if (last->next !=NULL){last->next->nodetype = 'E';}
 }

 /** [2] Replace the last node with new MATCH() */
 last->select_dis = newm->select_dis;
 last->Dneiatm  = newm->Dneiatm;
 last->Dextcon  = newm->Dextcon;
 last->Dtopodis = newm->Dtopodis;

 for (n=0;n<PAR.max_atomtype;++n){
   last->Natomtype[n] = newm->Natomtype[n];
 }

 last->Npair  = newm->Npair;
 last->Ncomponent  = newm->Ncomponent;
 last->nodetype  = newm->nodetype;
 for (n=0;n<newm->Npair;++n){
   last->anumA[n] = newm->anumA[n];
   last->anumB[n] = newm->anumB[n];
 } 

} /* end of Replace_the_last_node_with_new_MATCH() */



void Bubble_Sort_anumA_anumB_array_by_anumA(m)
 struct MATCH *m;
{
 /*
   Sort anumA[] and anumB[] in anumA order.
                               ^^^^^
  ex) [(1,0), (2,2), (3,1), (0,5)] -> [(0,5), (1,0), (2,2), (3,1)]
 */ 
 int Nchange,n,a0,b0;
 
 /*
 printf("#Bubble_Sort_anumA_anumB_array(m)\n");
 printf("BEFORE SORT");
 for (n=0;n<m->Npair;++n){ printf("(%d %d)",m->anumA[n],m->anumB[n]); }
 printf("\n");
 */

 do{ 
   Nchange = 0;
   for (n=0;n<(m->Npair-1);++n){
     if (m->anumA[n]>m->anumA[n+1]){
       a0 = m->anumA[n];
       b0 = m->anumB[n];
       m->anumA[n] = m->anumA[n+1];
       m->anumB[n] = m->anumB[n+1];
       m->anumA[n+1] = a0;
       m->anumB[n+1] = b0;
       ++Nchange;
      }
   }
 } while (Nchange >0);

 /*
 printf("AFTER SORT");
 for (n=0;n<m->Npair;++n){ printf("(%d %d)",m->anumA[n],m->anumB[n]); }
 printf("\n");
 */
} /* end of Bubble_Sort_anumA_anumB_array_by_anumA() */



void Bubble_Sort_anumA_anumB_array_by_anumB(m)
 struct MATCH *m;
{
 /*
   Sort anumA[] and anumB[] in anumB order.
                               ^^^^^
  ex) [(1,0), (2,2), (3,1), (0,5)] -> [(1,0), (3,1), (2,2), (0,5)]
 */ 
 int Nchange,n,a0,b0;
 
 /*
 printf("#Bubble_Sort_anumA_anumB_array(m)\n");
 printf("BEFORE SORT");
 for (n=0;n<m->Npair;++n){ printf("(%d %d)",m->anumA[n],m->anumB[n]); }
 printf("\n");
 */

 do{ 
   Nchange = 0;
   for (n=0;n<(m->Npair-1);++n){
     if (m->anumB[n]>m->anumB[n+1]){
       a0 = m->anumA[n];
       b0 = m->anumB[n];
       m->anumA[n] = m->anumA[n+1];
       m->anumB[n] = m->anumB[n+1];
       m->anumA[n+1] = a0;
       m->anumB[n+1] = b0;
       ++Nchange;
      }
   }
 } while (Nchange >0);

/*
 printf("AFTER SORT");
 for (n=0;n<m->Npair;++n){ printf("(%d %d)",m->anumA[n],m->anumB[n]); }
 printf("\n");
 */

} /* end of Bubble_Sort_anumA_anumB_array_by_anumB() */



int Check_Atom_Overlap_bwn_Pair_and_MATCH(pm,om)
  struct MATCH *pm;  /* one pair match (to be added to om)*/
  struct MATCH *om;  /* match in previous stage */ 
{ int n;
  for (n=0;n<om->Npair;++n){
    if ((pm->anumA[0]==om->anumA[n])||(pm->anumB[0]==om->anumB[n]))  return(1);
  }
  return(0);
} /* end of Check_Atom_Overlap_bwn_Pair_and_MATCH() */


int Bond_Equivalence_for_Pair_and_MATCH(pm,om,molA,molB)
  struct MATCH *pm;  /* one pair match (to be added to om)*/
  struct MATCH *om;  /* match in previous stage */ 
  struct MOLECULE *molA,*molB;
{ int n,ap,bp,ao,bo;
  
  ap = pm->anumA[0];
  bp = pm->anumB[0];
  
  if (PAR.bondtype_class == 'X'){
    for (n=0;n<om->Npair;++n){
      ao = om->anumA[n]; 
      bo = om->anumB[n]; 
      if ( ((molA->conmap.map[ap][ao] == '0') && (molB->conmap.map[bp][bo] != '0')) 
            ||  
           ((molA->conmap.map[ap][ao] != '0') && (molB->conmap.map[bp][bo] == '0'))) return(0); 
    }
  }

  else {
    for (n=0;n<om->Npair;++n){
      ao = om->anumA[n]; 
      bo = om->anumB[n]; 
      if (molA->conmap.map[ap][ao] != molB->conmap.map[bp][bo]) return(0); 
    }
  }


  return(1);

} /* end of Bond_Equivalence_for_Pair_and_MATCH() */



int Bond_Connection_for_Pair_and_MATCH(pm,om,molA,molB)
  struct MATCH *pm;  /* one pair match (to be added to om)*/
  struct MATCH *om;  /* match in previous stage */ 
  struct MOLECULE *molA,*molB;
{ int n,ap,bp,ao,bo;
  int connect;

  ap = pm->anumA[0];
  bp = pm->anumB[0];
  if ((ap<0)||(ap>=molA->Natom)||(bp<0)||(bp>=molB->Natom)) { printf("#ERROR:ap %d bp %d are out of range.\n",ap,bp); exit(1);}
  connect = 0; 
  for (n=0;n<om->Npair;++n){
     ao = om->anumA[n]; 
     bo = om->anumB[n];  
     /*
     if ((molA->conmap.map[ap][ao] != '0') && (molB->conmap.map[bp][bo] != '0')) connect = 1;
     */
     if ((ao>=0)&&(ao<molA->Natom)&&(bo>=0)&&(bo<molB->Natom) &&
        (molA->conmap.map[ap][ao] != '0') && (molB->conmap.map[bp][bo] != '0') ) { connect = 1; return(connect);}
  }

  return(connect);

} /* end of Bond_Connection_for_Pair_and_MATCH() */




int Number_of_Connected_Components_for_Pair_and_MATCH(pm,om,molA,molB)
  struct MATCH *pm;  /* one pair match (to be added to om)*/
  struct MATCH *om;  /* match in previous stage */ 
  struct MOLECULE *molA,*molB;
{ int a,n,ap,bp,ao,bo;
  int NedgeA, NedgeB,Ncluster;
 
  if (om->Ncomponent <=0){
     printf("#ERROR:Ncomponent %d Npair %d\n",om->Ncomponent,om->Npair); exit(1);
  }
 
 /** [1] Counring Nedge bwn pm and om **/ 
  ap = pm->anumA[0]; 
  bp = pm->anumB[0];
  NedgeA = NedgeB = 0;
  for (n=0;n<om->Npair;++n){
     ao = om->anumA[n]; 
     bo = om->anumB[n]; 
     if (molA->conmap.map[ao][ap]!='0')  NedgeA += 1;
     if (molB->conmap.map[bo][bp]!='0')  NedgeB += 1;
  }

 /** [2] (1) If Nedge==0, return (om->Ncomponent + 1) **/ 
 if ((NedgeA==0) && (NedgeB==0)){ ++PAR.M; return(om->Ncomponent + 1);}
 
 /**     (2) If Nedge==1, return (om->Ncomponent) **/ 
 else if ((NedgeA==1) && (NedgeB==1)) { ++PAR.N; return(om->Ncomponent);}

 /**     (3) If Nedege>=2, do reccursive count of connected component **/
 else if (NedgeA==NedgeB) { 
    ++PAR.O;
    /** initialize **/
  
    for (a=0;a<molA->Natom;++a) molA->atoms[a].mark = 0;

    for (n=0;n<om->Npair;++n){
       ao = om->anumA[n]; 
       molA->atoms[ao].mark = 1; 
    } 

   ap = pm->anumA[0];
   molA->atoms[ap].mark = 1; 

   /* Recurrsive Mark_Connected_Component */
   Ncluster = 1;
    for (a=0;a<molA->Natom;++a){
      if (molA->atoms[a].mark==1){
        Ncluster += 1;
        Mark_Connected_Component(molA,a,Ncluster);
      }
    } 
    return(Ncluster-1);
 }
 else{
  printf("#ERROR:NedgeA:%d NedgeB:%d\n",NedgeA,NedgeB); exit(1);
 }

} /* end of Number_of_Connected_Components_for_Pair_and_MATCH() */



int Number_of_Connected_Components_for_MATCH(anum,Npair,mol)
  unsigned char *anum;
  int Npair; 
  struct MOLECULE *mol; 
{ int a,n;
  int Ncluster;

  /* printf("#Number_of_Connected_Components_for_MATCH(anum,Npair %d,mol)\n",Npair); fflush(stdout); 
   */
  if (Npair==0) return(0);
  
  /** initialize **/
  for (a=0;a<mol->Natom;++a){mol->atoms[a].mark = 0;}
  
  for (n=0;n<Npair;++n){mol->atoms[anum[n]].mark = 1; }


  mol->atoms[anum[0]].mark = 1; 

  /* Recurrsive Mark_Connected_Component */
  Ncluster = 1;
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].mark==1){
      Ncluster += 1;
      Mark_Connected_Component(mol,a,Ncluster);
    }
  } 
  return(Ncluster-1);

} /* end of Number_of_Connected_Components_for_MATCH() */



void Set_Number_of_Connected_Components_for_MATCHlist(Mlist,molA)
  struct MATCH *Mlist;
  struct MOLECULE *molA; 
{ 
  struct MATCH *m;
  m = Mlist;
  while ((m->next !=NULL)&&(m->next->nodetype!='E')){
    m = m->next;
    m->Ncomponent = Number_of_Connected_Components_for_MATCH(m->anumA,m->Npair,molA);
 }

} /* end of Set_Number_of_Connected_Components_for_MATCHlist() */




void Mark_Connected_Component(mol,focus_a,Ncluster)
 struct MOLECULE *mol;
 int    focus_a;
 int    Ncluster;
{
  int b;
  mol->atoms[focus_a].mark = Ncluster;
  for (b=0;b<mol->Natom;++b){  
    if ((b!=focus_a) && (mol->atoms[b].mark==1) &&(mol->conmap.map[focus_a][b] != '0')){
      Mark_Connected_Component(mol,b,Ncluster); 
    }
  }
} /* end of Mark_Connected_Component() */





int Check_NonRedundancy_of_new_MATCH_with_MATCHlist(newm,Mlist,molA,molB)
  struct MATCH *newm;  /* new MATCH (to be added )*/
  struct MATCH *Mlist;  /* old MATCHlist */ 
  struct MOLECULE *molA,*molB;
{ 
 struct MATCH *om;
 int n,a;
 char same,nsame; 
 om = Mlist;

 while ((om->next!=NULL) &&(om->next->nodetype!='E')){
   om = om->next;
   same = 1;
   for (n=0;n<om->Npair;++n){
     if ((newm->anumA[n]!=om->anumA[n]) || (newm->anumB[n]!=om->anumB[n])) same = 0;
   }
   
   if (same==1) return(0);  /* Redundant !! */
  
   /** Check of Symmetric Redundancy **/ 
   if ((PAR.GenEquivAtomPermu=='B')||(PAR.NonRedSelectDis=='T')){
     if ((newm->Dneiatm==om->Dneiatm)&&(newm->Dextcon==om->Dextcon)
         &&(newm->Dtopodis==om->Dtopodis)){
       nsame = 1;
       for (a=0;a<PAR.max_atomtype;++a){
         if (newm->Natomtype[a]!=om->Natomtype[a]) nsame = 0;
       }
       if ((PAR.NonRedSelectDis=='T')&&(nsame==1)) return(0);
       if ((PAR.GenEquivAtomPermu=='B')&&(nsame==1)
          &&(Check_Equivalence_bwn_Two_MATCHes_by_Permutation(newm,om,molA,molB)==1)) return(0); 
     }
  }
 }
 
 return(1); /* Nonredundant !! */

} /* end of Check_NonRedundancy_of_new_MATCH_with_MATCHlist() */


int Check_Equivalence_bwn_Two_MATCHes(mP,mQ)
  struct MATCH *mP,*mQ;
{
 char same;
 int n;

 same = 1;
 for (n=0;n<mP->Npair;++n){
   if ((mP->anumA[n]!=mQ->anumA[n])||(mP->anumB[n]!=mQ->anumB[n]))
     { same = 0; return(0);}  /* Not equvalent !! */
 }
 if (same==1) return(1);  /* Equivalnet !! */
 return(0);
} /* end of Check_Equivalence_bwn_Two_MATCHes() */




float Worst_select_dis_in_MATCH_list(Mlist)
 struct MATCH *Mlist;
{
 float worst_select_dis;
 struct MATCH *mn;
 char init;

 worst_select_dis = 0.0;
 mn = Mlist;
 init = 1;
 while ((mn->next!=NULL) && (mn->next->nodetype != 'E')){
   mn = mn->next;
   /* printf("#select_dis %f\n",mn->select_dis); */
   if (init==1) {worst_select_dis = mn->select_dis; init = 0;}
   else{
    if (mn->select_dis > worst_select_dis) worst_select_dis = mn->select_dis;
   } 
 }

 /* printf("#-->worst_select_dis %f\n",worst_select_dis); */
 return(worst_select_dis);

} /* end of Worst_select_dis_in_MATCH_list() */



int Check_Diff_ShortPath_bwn_Pair_and_MATCH(pm,om,molA,molB,maxDIFtopodis)
  struct MATCH *pm;  /* one pair match (to be added to om)*/
  struct MATCH *om;  /* match in previous stage */ 
  struct MOLECULE *molA,*molB;
  int    maxDIFtopodis;
{ 
 int i;
  for (i=0;i<om->Npair;++i){
    if (abs(molA->topodismap.map[pm->anumA[0]][om->anumA[i]] - molB->topodismap.map[pm->anumB[0]][om->anumB[i]])>maxDIFtopodis) return(0);
  }
  return(1);
} /* end of Check_Diff_ShortPath_bwn_Pair_and_MATCH() */

