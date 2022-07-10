/*
  <c_clique.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  Bron-Kerbosch algorithm for connected MCS detection.
  based on:

  I. Koch.  Enumerating all connected maximal common subgraphs in two graphs.
  Theoretical Computer Science, 250, (2001), 1-30.

  F.Cazals, C.Karande.  An algorithm for reporting maximal c-cliques.
  Theoretical Computer Science, 349, (2005), 484-490.

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
#include "MrgSrtMATCH.h"
#include "clique.h"
#include "options.h"

/*** FUNCTIONS (GLOBAL) ***/
int  C_clique_BronKerbosch_Version1();
int  C_clique_Cazals_Karande();
void Evaluate_Connected_Degree_of_Atom_Pair();

/*** FUNCTIONS (LOCAL) ***/
static int recurrent_func_C_clique_BronKerbosch_Version1();
static int recurrent_func_C_clique_Cazals_Karande();
static int evaluate_c_edge_connection();

static int Nrecurrent;


int C_clique_BronKerbosch_Version1(optMlist,pairlist,Npairlist,molA,molB)
 struct MATCH *optMlist;    /* MATCHlist to be optmized (already malloced: Nmatch = smllerInt(NatmA,NatmB)) */
 struct MATCH *pairlist;    /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 int    Npairlist;           /* Number of pair in *pairlist */ 
 struct MOLECULE *molA,*molB;
{
 /*
  Based on Koch's algotithm (2001)
 */ 
 unsigned char  uanumA,uanumB,vanumA,vanumB;
 struct MATCH *u,*v;
 struct MATCH CS;
 char   *CAC_CAD_NOTstr;    /* {'c','d','n',' '}-array 
                                for Connected Candidate(CAC),Disconnected Candidate (CAD) and Not(NOT) set */ 
 char   *Tstr;              /* {'t',' '}-array. 't':already searched node */ 


 printf("#C_clique_BronKerbosch_Version1(optMlist,pairlist,Npairlist:%d molA,molB)\n",Npairlist);
 Malloc_MATCH(&CS,1);
 CS.Npair = 1;
 CAC_CAD_NOTstr = (char *)malloc(sizeof(char)*Npairlist); 
 
 Tstr = (char *)malloc(sizeof(char)*Npairlist); 
 Init_Bit_String(Tstr,' ',Npairlist);

 u = pairlist;
 while (u->next != NULL){
   u = u->next;
   uanumA = u->anumA[0];
   uanumB = u->anumB[0];
   CS.anumA[0] = uanumA;
   CS.anumB[0] = uanumB;

   Init_Bit_String(CAC_CAD_NOTstr,' ',Npairlist);
   
   v = pairlist;
   while (v->next != NULL){
     v = v->next;
     vanumA = v->anumA[0];
     vanumB = v->anumB[0];
     if ((uanumA != vanumA) && (uanumB != vanumB)){
       /* 'c'-edge connection */
       if (evaluate_c_edge_connection(molA,uanumA,vanumA,molB,uanumB,vanumB)==1){
         if (Tstr[v->num]=='t') CAC_CAD_NOTstr[v->num] = 'n';
         else                   CAC_CAD_NOTstr[v->num] = 'c';
 
       }
       /* 'd'-edge connection */
       if ((molA->conmap.map[uanumA][vanumA]=='0') && (molB->conmap.map[uanumB][vanumB]=='0')){
          CAC_CAD_NOTstr[v->num] = 'd'; 
       }
     }
     /* printf("v [%3d] (%3d %3d) -->'%c'\n",v->num,vanumA,vanumB,CAC_CAD_NOTstr[v->num]); */
   }

   /* printf("FIRST_TRY %d (%d %d)\n",u->num,uanumA+1,uanumB+1);  */
   recurrent_func_C_clique_BronKerbosch_Version1(optMlist,&CS,pairlist,Npairlist,CAC_CAD_NOTstr,molA,molB);
   Tstr[u->num] = 't'; 
   /* printf("FIRST_done %d(%d %d)\n",u->num,uanumA+1,uanumB+1); */
 } 

 free(CAC_CAD_NOTstr);
 Free_MATCH(&CS);
 return(1);
} /* end of C_clique_BronKerbosch_Version1() */








int recurrent_func_C_clique_BronKerbosch_Version1(optMlist,CS,pairlist,Npairlist,CAC_CAD_NOTstr,molA,molB)
 struct MATCH *optMlist;    /* MATCHlist to be optmized (already malloced: Nmatch = smllerInt(NatmA,NatmB)) */
 struct MATCH *CS;          /* MATCH of complete subgraph (matching atom pairs, Nmax_match=smallerInt(NatmA,NatomB) */ 
 struct MATCH *pairlist;    /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 int    Npairlist;          /* Number of pair list */ 
 char   *CAC_CAD_NOTstr;    /* {'c','d','n',' '}-array 
                                for Connected Candidate(CAC),Disconnected Candidate (CAD) and Not(NOT) set */ 
 struct MOLECULE *molA,*molB;
{
 /*
  Based on Koch's algotithm (2001)
 */ 
 int    N_CAC, N_CAD,N_NOT;  /* Nmember of CAC, CAD, and NOT **/ 
 int    i,unum;
 unsigned char  uanumA,uanumB,vanumA,vanumB;
 struct MATCH *u,*v;
 struct MATCH CSnew;
 char   *CAC_CAD_NOTstr_new,connect;

 /* 
 printf("int recurrent_func_C_clique_BronKerbosch_Version1() CS.Npair %d\n",CS->Npair); 
 Show_Bit_String(CAC_CAD_NOTstr,Npairlist,"");
 fflush(stdout); 
*/
 
 N_CAC  = Sum_Bit_String(CAC_CAD_NOTstr,'c',Npairlist);
 N_CAD  = Sum_Bit_String(CAC_CAD_NOTstr,'d',Npairlist);
 N_NOT  = Sum_Bit_String(CAC_CAD_NOTstr,'n',Npairlist);
 /*
 printf("N_CAC %d N_CAD %d N_NOT %d\n",N_CAC,N_CAD,N_NOT); 
 Print_MATCH(CS,"CS");
 */
  
 /** 'E':ENDING PROCEDURE **/
 if ((N_CAC==0) && (N_NOT>0)) return(0);
 
 if ((N_CAC==0) && (N_NOT==0)){
   CS->select_dis = CS->Npair; 
   Add_MATCH_to_MATCHlist_Keep_Only_Best(CS,optMlist); 
   /* printf("return 1\n"); */
   return(1);
 } 

 /** 'R':RECURSIVE CALLING of 'C_clique()' function **/
 while (N_CAC > 0){

    /** [1] delete 'u' from 'CAC' **/
    unum  = Get_First_Target_Symbol_In_Bit_String(CAC_CAD_NOTstr,'c',Npairlist);
    CAC_CAD_NOTstr[unum] = ' '; 
    u = pairlist;
    while ((u->next != NULL) && (u->num != unum)){ u = u->next;}
    if (u->num != unum) {printf("#ERROR:unum %d is not found.\n",unum); exit(1);} 
   
    uanumA = u->anumA[0]; 
    uanumB = u->anumB[0]; 
    /* printf("add [%d](%d %d)\n",cnum,canumA,canumB); fflush(stdout);  */
    /** making CSnew (CS + c)**/
    Malloc_MATCH(&CSnew,CS->Npair+1); 
    Copy_MATCH(&CSnew,CS);
    CSnew.Npair = CS->Npair + 1;
    CSnew.anumA[CS->Npair] = uanumA;
    CSnew.anumB[CS->Npair] = uanumB;

    /** making CAC_CAD_NOTstr_new[] **/
    CAC_CAD_NOTstr_new = (char *)malloc(sizeof(char)*Npairlist); 
    for (i=0;i<Npairlist;++i) CAC_CAD_NOTstr_new[i] = CAC_CAD_NOTstr[i];
    
    /** for v in (CADnew), if (c and v) are c_edged, then add v into CADnew. **/
    v = pairlist;
    while (v->next != NULL){
      v = v->next;
      vanumA = v->anumA[0];
      vanumB = v->anumB[0];
      connect = ' ';
      if ((uanumA != vanumA) && (uanumB !=vanumB)){
        if (evaluate_c_edge_connection(molA,uanumA,vanumA,molB,uanumB,vanumB)==1) connect = 'c';
        if ((molA->conmap.map[uanumA][vanumA]=='0') && (molB->conmap.map[uanumB][vanumB]=='0')) connect = 'd';
      }

      if ((CAC_CAD_NOTstr_new[v->num]=='d')&&(connect=='c')) CAC_CAD_NOTstr_new[v->num] = 'c'; 
      if (connect==' ') CAC_CAD_NOTstr_new[v->num] = ' ';  
      /* printf("c[%3d](%3d %3d) v [%3d](%3d %3d) -->'%c'\n",cnum,canumA,canumB,v->num,vanumA,vanumB,CAC_CAD_NOTstr[v->num]); fflush(stdout);  */
    }

    /** [4] Do BronKerbosch(), recursively. **/
    /* printf("add %d (%d %d)\n",unum,uanumA,uanumB); */
    recurrent_func_C_clique_BronKerbosch_Version1(optMlist,&CSnew,pairlist,Npairlist,CAC_CAD_NOTstr_new, molA,molB);
  
    free(CAC_CAD_NOTstr_new); 
    Free_MATCH(&CSnew);

    /** [5] Add c to 'NOT' **/
    CAC_CAD_NOTstr[unum] = 'n';
    /* printf("done %d (%d %d)\n",unum,uanumA,uanumB); */
    N_CAC  = Sum_Bit_String(CAC_CAD_NOTstr,'c',Npairlist);

  } /* while */

  return(0);

} /* end of recurrent_func_C_clique_BronKerbosch_Version1() */




int C_clique_Cazals_Karande(optMlist,pairlist,Npairlist,molA,molB,ConnectType,maxDIFtopodis)
 struct MATCH *optMlist;    /* MATCHlist to be optmized (already malloced: Nmatch = smllerInt(NatmA,NatmB)) */
 struct MATCH *pairlist;    /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 int    Npairlist;           /* Number of pair in *pairlist */ 
 struct MOLECULE *molA,*molB;
 char   ConnectType;        /* 'C'onnected, 't'opo_dis-constrained MCS */
 int    maxDIFtopodis;
{
 /*
  Based on Cazals-Karande algotithm (2005)
 */ 
 unsigned char  uanumA,uanumB,vanumA,vanumB;
 struct MATCH *u,*v;
 struct MATCH CS;
 char   *CAC_CAD_NOC_NODstr;    /* {'c','d','C','D'}-array 
         for Connected Candidate(CAC),Disconnected Candidate (CAD) and Not_connected(NOC), Not_disconnected(NOD) set */ 
 char   *Tstr;              /* {'t',' '}-array. 't':already searched node */ 
 char  topodis_ok;

 printf("#C_clique_Cazals_Karande(optMlist,pairlist,Npairlist:%d molA,molB)\n",Npairlist);
 Malloc_MATCH(&CS,1);
 CS.Npair = 1;
 CAC_CAD_NOC_NODstr = (char *)malloc(sizeof(char)*Npairlist); 
 
 Tstr = (char *)malloc(sizeof(char)*Npairlist); 
 Init_Bit_String(Tstr,' ',Npairlist);
 Nrecurrent = 0;

 u = pairlist;
 while (u->next != NULL){
   u = u->next;
   uanumA = u->anumA[0];
   uanumB = u->anumB[0];
   CS.anumA[0] = uanumA;
   CS.anumB[0] = uanumB;

   Init_Bit_String(CAC_CAD_NOC_NODstr,' ',Npairlist);
   
   v = pairlist;
   while (v->next != NULL){
     v = v->next;
     vanumA = v->anumA[0];
     vanumB = v->anumB[0];
     if ((uanumA != vanumA) && (uanumB != vanumB)){
       /* 'c'-edge connection */
       if (evaluate_c_edge_connection(molA,uanumA,vanumA,molB,uanumB,vanumB)==1) { 
         if (Tstr[v->num]=='t') CAC_CAD_NOC_NODstr[v->num] = 'C';
         else                   CAC_CAD_NOC_NODstr[v->num] = 'c';
 
       }
       /* 'd'-edge connection */
       if ((molA->conmap.map[uanumA][vanumA]=='0') && (molB->conmap.map[uanumB][vanumB]=='0')){

            topodis_ok = 1;
            if ((ConnectType=='t')
               &&(abs(molA->topodismap.map[uanumA][vanumA] - molB->topodismap.map[uanumB][vanumB])>maxDIFtopodis))
             {topodis_ok = 0;}

         if (topodis_ok==1){
           if (Tstr[v->num]=='t') CAC_CAD_NOC_NODstr[v->num] = 'D';
           else                   CAC_CAD_NOC_NODstr[v->num] = 'd';
         }
       }
     }
     /* printf("v [%3d] (%3d %3d) -->'%c'\n",v->num,vanumA,vanumB,CAC_CAD_NOTstr[v->num]); */
   }

   /* printf("FIRST_TRY %d (%d %d)\n",u->num,uanumA+1,uanumB+1);  */
   recurrent_func_C_clique_Cazals_Karande(optMlist,&CS,pairlist,Npairlist,CAC_CAD_NOC_NODstr,molA,molB,ConnectType);
   Tstr[u->num] = 't'; 
   /* printf("FIRST_done %d(%d %d)\n",u->num,uanumA+1,uanumB+1); */

 } 

 free(CAC_CAD_NOC_NODstr);
 Free_MATCH(&CS);
 return(1);

} /* end of C_clique_Cazals_Karande() */







int recurrent_func_C_clique_Cazals_Karande(optMlist,CS,pairlist,Npairlist,CAC_CAD_NOC_NODstr,molA,molB,ConnectType,maxDIFtopodis)
 struct MATCH *optMlist;    /* MATCHlist to be optmized (already malloced: Nmatch = smllerInt(NatmA,NatmB)) */
 struct MATCH *CS;          /* MATCH of complete subgraph (matching atom pairs, Nmax_match=smallerInt(NatmA,NatomB) */ 
 struct MATCH *pairlist;    /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 int    Npairlist;          /* Number of pair list */ 
 char   *CAC_CAD_NOC_NODstr;    /* {'c','d','C','D'}-array 
                                for Connected Candidate(CAC),Disconnected Candidate (CAD) and Not(NOT) set */ 
 struct MOLECULE *molA,*molB;
 char   ConnectType;        /* 'C'onnected, 't'opo_dis-constrained MCS */
 int    maxDIFtopodis;
{
 /*
  Based on Kazals_Karande algotithm (2005)
 */ 
 int    N_CAC, N_CAD,N_NOC, N_NOD;  /* Nmember of CAC, CAD, NOC, NOD **/ 
 int    i,unum,ret;
 unsigned char  uanumA,uanumB,vanumA,vanumB;
 struct MATCH *u,*v;
 struct MATCH CSnew;
 char   *CAC_CAD_NOC_NODstr_new,connect;

 /* 
 printf("int recurrent_func_C_clique_Cazals_Karande() CS.Npair %d\n",CS->Npair);
 Show_Bit_String(CAC_CAD_NOC_NODstr,Npairlist,"");
 fflush(stdout); 
 */
 
 N_CAC  = Sum_Bit_String(CAC_CAD_NOC_NODstr,'c',Npairlist);
 N_CAD  = Sum_Bit_String(CAC_CAD_NOC_NODstr,'d',Npairlist);
 N_NOC  = Sum_Bit_String(CAC_CAD_NOC_NODstr,'C',Npairlist);
 N_NOD  = Sum_Bit_String(CAC_CAD_NOC_NODstr,'D',Npairlist);
 /*
 printf("N_CAC %d N_CAD %d N_NOT %d\n",N_CAC,N_CAD,N_NOT); 
 Print_MATCH(CS,"CS");
 */

 
 /** EXIT() for TIME-UP !!**/
 if (PAR.MAX_COMP_TIME_SEC>0.0){
   Nrecurrent += 1;
   if (Nrecurrent==10000){
     Set_END_DATE(); 
     if (PAR.COMP_TIME_SEC>PAR.MAX_COMP_TIME_SEC){
       printf("#WARNING:TIME-UP. COMP_TIME_SEC(%lf) is over MAX (%lf)\n",PAR.COMP_TIME_SEC, PAR.MAX_COMP_TIME_SEC); 
       /* exit(1); */
       return(-1);
     }
    Nrecurrent = 0;
   } 
 }

 /** 'E':ENDING PROCEDURE **/
 if ((N_CAC==0) && (N_NOC>0)){ return(0); }
 
 if ((N_CAC==0) && (N_NOC==0)){
   CS->select_dis = CS->Npair; 
   Add_MATCH_to_MATCHlist_Keep_Only_Best(CS,optMlist); 
   return(1);
 } 


 /** 'R':RECURSIVE CALLING of 'C_clique()' function **/
 while (N_CAC > 0){

    /** [1] delete 'u' from 'CAC' **/
    unum  = Get_First_Target_Symbol_In_Bit_String(CAC_CAD_NOC_NODstr,'c',Npairlist);
    CAC_CAD_NOC_NODstr[unum] = ' '; 
    u = pairlist;
    while ((u->next != NULL) && (u->num != unum)){ u = u->next;}
    if (u->num != unum) {printf("#ERROR:unum %d is not found.\n",unum); exit(1);} 
   
    uanumA = u->anumA[0]; 
    uanumB = u->anumB[0]; 
    /* printf("add [%d](%d %d)\n",cnum,canumA,canumB); fflush(stdout);  */
    /** making CSnew (CS + c)**/
    Malloc_MATCH(&CSnew,CS->Npair+1); 
    Copy_MATCH(&CSnew,CS);
    CSnew.Npair = CS->Npair + 1;
    CSnew.anumA[CS->Npair] = uanumA;
    CSnew.anumB[CS->Npair] = uanumB;

    /** making CAC_CAD_NOTstr_new[] **/
    CAC_CAD_NOC_NODstr_new = (char *)malloc(sizeof(char)*Npairlist); 
    for (i=0;i<Npairlist;++i) CAC_CAD_NOC_NODstr_new[i] = ' ';
    
    /** for v in (CADnew), if (c and v) are c_edged, then add v into CADnew. **/
    v = pairlist;
    while (v->next != NULL){
      v = v->next;
      vanumA = v->anumA[0];
      vanumB = v->anumB[0];
      connect = ' ';
      if ((uanumA != vanumA) && (uanumB !=vanumB)){
        if (evaluate_c_edge_connection(molA,uanumA,vanumA,molB,uanumB,vanumB)==1) connect = 'c';
        if ((molA->conmap.map[uanumA][vanumA]=='0') && (molB->conmap.map[uanumB][vanumB]=='0')){
            if (ConnectType=='t'){
              if ((abs(molA->topodismap.map[uanumA][vanumA] - molB->topodismap.map[uanumB][vanumB])<=maxDIFtopodis)) connect = 'd'; 
            } 
            else connect = 'd'; 
        }
      }

      if ((CAC_CAD_NOC_NODstr[v->num]=='d')&&(connect=='d')) CAC_CAD_NOC_NODstr_new[v->num] = 'd';
 
      if ((CAC_CAD_NOC_NODstr[v->num]=='c')&&(connect!=' ')) CAC_CAD_NOC_NODstr_new[v->num] = 'c'; 
      if ((CAC_CAD_NOC_NODstr[v->num]=='d')&&(connect=='c')) CAC_CAD_NOC_NODstr_new[v->num] = 'c'; 
      
      if ((CAC_CAD_NOC_NODstr[v->num]=='D')&&(connect=='d')) CAC_CAD_NOC_NODstr_new[v->num] = 'D';
      
      if ((CAC_CAD_NOC_NODstr[v->num]=='C')&&(connect!=' ')) CAC_CAD_NOC_NODstr_new[v->num] = 'C';
      if ((CAC_CAD_NOC_NODstr[v->num]=='D')&&(connect=='c')) CAC_CAD_NOC_NODstr_new[v->num] = 'C';

      if (connect==' ') CAC_CAD_NOC_NODstr_new[v->num] = ' ';  
      /* printf("c[%3d](%3d %3d) v [%3d](%3d %3d) -->'%c'\n",cnum,canumA,canumB,v->num,vanumA,vanumB,CAC_CAD_NOTstr[v->num]); fflush(stdout);  */
    }

    /** [4] Do BronKerbosch(), recursively. **/
    /* printf("add %d (%d %d)\n",unum,uanumA,uanumB); */
    ret  = recurrent_func_C_clique_Cazals_Karande(optMlist,&CSnew,pairlist,Npairlist,CAC_CAD_NOC_NODstr_new, molA,molB,ConnectType);
  
    free(CAC_CAD_NOC_NODstr_new); 
    Free_MATCH(&CSnew);
    if (ret==-1) return(-1);
    /** [5] Add c to 'NOT' **/
    CAC_CAD_NOC_NODstr[unum] = 'C';
    /* printf("done %d (%d %d)\n",unum,uanumA,uanumB); */
    N_CAC  = Sum_Bit_String(CAC_CAD_NOC_NODstr,'c',Npairlist);

  } /* while */

  return(0);

} /* end of recurrent_func_C_clique_Cazals_Karande() */




void Evaluate_Connected_Degree_of_Atom_Pair(Mlist,molA,molB)
 struct MATCH    *Mlist;
 struct MOLECULE *molA,*molB;
{
 int pA,pB,qA,qB;
 struct MATCH *pn,*qn;

 /*
   Evaluated connected degree is stored in 'pn->select_dis'.
 */

 pn = Mlist;
 while (pn->next != NULL){
   pn = pn->next;
   pA = pn->anumA[0];
   pB = pn->anumB[0];
   pn->select_dis = 0.0;
   qn = Mlist;
   while (qn->next != NULL){
     qn = qn->next;
     qA = qn->anumA[0];
     qB = qn->anumB[0];
     if ((pA != qA) &&  ( pB != qB)){ 
       if (evaluate_c_edge_connection(molA,pA,qA,molB,pB,qB)==1)  pn->select_dis += 1.0;
    }
   }
 }

} /* end of Evaluate_Connected_Degree_of_Atom_Pair() */





int evaluate_c_edge_connection(molA,uA,vA,molB,uB,vB)
  struct MOLECULE *molA;
  int    uA,vA; /* atom number for molA */
  struct MOLECULE *molB;
  int    uB,vB; /* atom number for molB */
{
  if (PAR.bondtype_class == 'X'){
    if ((molA->conmap.map[uA][vA]!='0') && (molB->conmap.map[uB][vB]!='0')) return(1);
    else return(0);
  }
  else{
    if ((molA->conmap.map[uA][vA]!='0') && (molB->conmap.map[uB][vB]!='0') &&
        (molA->conmap.map[uA][vA] == molB->conmap.map[uB][vB]) ) return(1);
    else return(0);
  }
}

