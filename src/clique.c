/*
  <clique.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


   functions for Bron-Kerbosch algorithm for clique detection 


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "MrgSrtMATCH.h"
#include "options.h"

/*** FUNCTIONS (GLOBAL) ***/
int BronKerbosch_Version1_BitString();
int BronKerbosch_Pivot_BitString();
void Evaluate_Degree_of_Atom_Pair();
int  Delete_Small_Degree_of_Atom_Pair();
int  Keep_Only_Maximum_Npair_MATCH();
int  Keep_Only_First_K_MATCH();
int  Renumbering_MATCHlist();
void Init_Bit_String();
int  Sum_Bit_String();
void Show_Bit_String();
int Get_First_Target_Symbol_In_Bit_String();
int Get_i_th_Target_Symbol_In_Bit_String();

/*** FUNCTIONS (LOCAL) ***/

static int evaluate_edge_connection();
static int Nrecurrent;



int BronKerbosch_Version1_BitString(optMlist,CS,pairlist,Npair,CA_NOTstring,molA,molB,ConnectType,maxDIFtopodis)
 struct MATCH *optMlist;    /* MATCHlist to be optmized (already malloced: Nmatch = smllerInt(NatmA,NatmB)) */
 struct MATCH *CS;          /* MATCH of complete subgraph (matching atom pairs, Nmax_match=smallerInt(NatmA,NatomB) */ 
 struct MATCH *pairlist;    /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 int    Npair;              /* Number of pair list */ 
 char   *CA_NOTstring;     /* {'c','n',' '}-array for Candidate(CA) and Not(NOT) set */ 
 struct MOLECULE *molA,*molB;
 char   ConnectType;        /* 'D'isconnected, 'T'opo_dis-constrained MCS */
 int   maxDIFtopodis;      /* only for ConnectType=='T' */
{
 int    N_CA, N_NOT;  /* Nmember of CA and NOT **/ 
 int    cnum, ret;
 unsigned char  canumA, canumB, topodis_ok;
 struct MATCH *v;
 struct MATCH CSnew;
 char   *CA_NOTstring_new;
 /*
 printf("#BronKerbosch_Version1_BitString(optMlist,CS,pairlist,Npair %d CAstring,NOTstring,molA,molB)\n",Npair);
 fflush(stdout); 
 */
 
 N_CA   = Sum_Bit_String(CA_NOTstring,'c',Npair);
 N_NOT  = Sum_Bit_String(CA_NOTstring,'n',Npair);
 /* 
 printf("N_CA %d N_NOT %d\n",N_CA,N_NOT); 
 Print_MATCH(CS,"CS");
 */


 /** EXIT() for TIME-UP !!**/
 if (N_CA==Npair) Nrecurrent = 0;
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
 if ((N_CA==0) && (N_NOT>0)) { return(0); }
 
 if ((N_CA==0) && (N_NOT==0)){
   CS->select_dis = CS->Npair; 
   Add_MATCH_to_MATCHlist_Keep_Only_Best(CS,optMlist);
   /* Print_MATCHlist(optMlist,"temp");  */
   return(1);
 } 

 /** 'R':RECURSIVE CALLING of 'BronKerbosch()' function **/
 while (N_CA > 0){

    /** [1] delete 'c' from 'CAlist' **/
    N_CA   = Sum_Bit_String(CA_NOTstring,'c',Npair);
    cnum  = Get_First_Target_Symbol_In_Bit_String(CA_NOTstring,'c',Npair);
    CA_NOTstring[cnum] = ' '; 
    v = pairlist;
    while ((v->next != NULL) && (v->num != cnum)){ v = v->next;}
    if (v->num != cnum) {printf("#ERROR:cnum %d is not found.\n",cnum); exit(1);} 
   
    canumA = v->anumA[0]; 
    canumB = v->anumB[0]; 

    /** making CSnew (CS + c)**/
    Malloc_MATCH(&CSnew,CS->Npair+1); 
    Copy_MATCH(&CSnew,CS);
    CSnew.Npair = CS->Npair + 1;
    CSnew.anumA[CS->Npair] = canumA;
    CSnew.anumB[CS->Npair] = canumB;

    /** [2] making CA_NOTlist_new (atom pairs connecting to the atom pair 'c' in CAlist) **/
    /** [3] making NOTlist_new (atom pairs connecting to the atom pair 'c' in NOTlist) **/
    CA_NOTstring_new = (char *)malloc(sizeof(char)*Npair); 
    Init_Bit_String(CA_NOTstring_new,' ',Npair); 
    v = pairlist;
    while (v->next != NULL){
      v = v->next;
      if ( ((CA_NOTstring[v->num]=='c')||(CA_NOTstring[v->num]=='n')) &&
            (canumA != v->anumA[0]) && (canumB != v->anumB[0]) &&
            (evaluate_edge_connection(molA,canumA,v->anumA[0],molB,canumB,v->anumB[0])==1)) {
          topodis_ok = 1;
          if ((ConnectType=='T')
             &&(abs(molA->topodismap.map[canumA][v->anumA[0]]-molB->topodismap.map[canumB][v->anumB[0]])>maxDIFtopodis)) 
           {topodis_ok = 0;}

          if (topodis_ok==1) CA_NOTstring_new[v->num] = CA_NOTstring[v->num];
      } 
    }

    /** [4] Do BronKerbosch(), recursively. **/
    ret = BronKerbosch_Version1_BitString(optMlist,&CSnew,pairlist,Npair,CA_NOTstring_new, molA,molB, ConnectType,maxDIFtopodis);
    Free_MATCH(&CSnew);
    free(CA_NOTstring_new); 
    if (ret==-1) return(-1);
    /** [5] Add c to 'NOT' **/
    CA_NOTstring[cnum] = 'n';
    N_CA   = Sum_Bit_String(CA_NOTstring,'c',Npair);

  } /* while */

  return(0);

} /* end of BronKerbosch_Version1_BitString() */















int BronKerbosch_Pivot_BitString(optMlist,CS,pairlist,Npair,CA_NOTstring,molA,molB,ConnectType,maxDIFtopodis)
 struct MATCH *optMlist;    /* MATCHlist to be optmized (already malloced: Nmatch = smllerInt(NatmA,NatmB)) */
 struct MATCH *CS;          /* MATCH of complete subgraph (matching atom pairs, Nmax_match=smallerInt(NatmA,NatomB) */ 
 struct MATCH *pairlist;    /* MATCHlist of all the possible pairs (already malloced: Nmatch = 1 ) */
 int    Npair;              /* Number of pair list */ 
 char   *CA_NOTstring;     /* {'c','n',' '}-array for Candidate(CA) and Not(NOT) set */ 
 struct MOLECULE *molA,*molB;
 char  ConnectType;        /* 'D'isconnected, 'T'opo_dis-constrained MCS */
 int   maxDIFtopodis;      /* only for ConnectType=='T' */
{
 int    N_CA, N_NOT;  /* Nmember of CA and NOT **/ 
 int    cnum,pnum,ret;
 unsigned char canumA,canumB,panumA,panumB,topodis_ok;
 int    i;
 unsigned char   connect_p_and_c;
 struct MATCH *v;
 struct MATCH CSnew;
 char   *CA_NOTstring_new;

 /*
 printf("#BronKerbosch_Version1_BitString(optMlist,CS,pairlist,Npair %d CAstring,NOTstring,molA,molB)\n",Npair);
 fflush(stdout); 
 printf("CS.Npair %d\n",CS->Npair); fflush(stdout); 
 */ 
 N_CA   = Sum_Bit_String(CA_NOTstring,'c',Npair);
 N_NOT  = Sum_Bit_String(CA_NOTstring,'n',Npair);
 /*
 if ((PAR.N%10000)==0){
  printf("#BronKerbosch_Pivot(CS.Npair %2d Ntry %d Nfree %d Nkept %d Nstop %d Nout %d)\n",
  CS->Npair,PAR.N,PAR.M,PAR.N-PAR.M,PAR.O,PAR.P); fflush(stdout);
 }
 */
 
 /* 
 printf("N_CA %d N_NOT %d\n",N_CA,N_NOT); 
 Print_MATCH(CS,"CS");
 */

 /** EXIT() for TIME-UP !!**/
 if (N_CA==Npair) Nrecurrent = 0;
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
 if ((N_CA==0) && (N_NOT>0)) { 
  PAR.O += 1; 
  /* 
  printf("END\n"); fflush(stdout);
  Print_MATCH(CS,"end");
  */
  return(0); }
 
  if ((N_CA==0) && (N_NOT==0)){
    CS->select_dis = CS->Npair; 
    PAR.P += 1;
    Add_MATCH_to_MATCHlist_Keep_Only_Best(CS,optMlist); 
    /* PAR.P += Add_MATCH_to_MATCHlist_Keep_Only_Best(CS,optMlist); */
    /* Print_MATCHlist(optMlist,"temp");  */
    /*
    printf("OUT\n"); fflush(stdout);
    Print_MATCH(CS,"out");
    */
    return(1);
  } 

  /*** choosing pivot head of CAstring ***/
  pnum  = Get_First_Target_Symbol_In_Bit_String(CA_NOTstring,'c',Npair);
  if (pnum<0) {printf("ERROR: no more '1' in CAstring for p"); exit(1);}
  v = pairlist;
  while ((v->next != NULL) && (v->num != pnum)){ v = v->next;}
  if (v->num != pnum) {printf("#ERROR:pnum %d is not found.\n",pnum); exit(1);} 
  panumA = v->anumA[0]; 
  panumB = v->anumB[0]; 


  /* printf("#pivot_num %d A %d B %d\n",pnum,panumA,panumB); fflush(stdout); */
  /** 'R':RECURSIVE CALLING of 'BronKerbosch()' function **/
  i = 0;
  while (Get_i_th_Target_Symbol_In_Bit_String(CA_NOTstring,i,'c',Npair)>=0){ 
    /** [1] delete 'c' from 'CAlist' **/
    N_CA   = Sum_Bit_String(CA_NOTstring,'c',Npair);
    cnum  = Get_i_th_Target_Symbol_In_Bit_String(CA_NOTstring,i,'c',Npair); 
    if (cnum<0) {printf("ERROR: no more '1' in CAstring in %d for c.(num1 for CAstring:%d)",i,Sum_Bit_String(CA_NOTstring,'c',Npair)); exit(1);}
    v = pairlist;
    while ((v->next != NULL) && (v->num != cnum)){ v = v->next;}
    if (v->num != cnum) {printf("#ERROR:pnum %d is not found.\n",cnum); exit(1);} 
    canumA = v->anumA[0]; 
    canumB = v->anumB[0]; 

    /* printf("#c_num %d A %d B %d\n",cnum,canumA,canumB); fflush(stdout); */
    /** checking connection 'c' and 'p' **/
    if ((canumA == panumA) || (canumB == panumB)) connect_p_and_c = 0;
    else if (evaluate_edge_connection(molA,canumA,v->anumA[0],molB,canumB,v->anumB[0])==1) {
      connect_p_and_c = 1;
      if ((ConnectType=='T')
         &&(abs(molA->topodismap.map[canumA][panumA]-molB->topodismap.map[canumB][panumB])>maxDIFtopodis))connect_p_and_c = 0;
    } 
    else connect_p_and_c = 0;

    
    if (connect_p_and_c == 0){
      PAR.N += 1;
      CA_NOTstring[cnum] = ' ';
      /** making CSnew (CS + c)**/
      Malloc_MATCH(&CSnew,CS->Npair+1); 
      Copy_MATCH(&CSnew,CS);
      CSnew.Npair = CS->Npair + 1;
      CSnew.anumA[CS->Npair] = canumA;
      CSnew.anumB[CS->Npair] = canumB;

      /** [2] making CA_NOTlist_new (atom pairs connecting to the atom pair 'c' in CAlist) **/
      /** [3] making NOTlist_new (atom pairs connecting to the atom pair 'c' in NOTlist) **/
      CA_NOTstring_new = (char *)malloc(sizeof(char)*Npair); 
      Init_Bit_String(CA_NOTstring_new,' ',Npair); 
      v = pairlist;
      while (v->next != NULL){
        v = v->next;
        if ( ((CA_NOTstring[v->num]=='c') || (CA_NOTstring[v->num]=='n')) &&
             (canumA != v->anumA[0]) && (canumB != v->anumB[0]) &&
             (evaluate_edge_connection(molA,canumA,v->anumA[0],molB,canumB,v->anumB[0])==1) ) {
            topodis_ok = 1;
            if ((ConnectType=='T')
               &&(abs(molA->topodismap.map[canumA][v->anumA[0]]-molB->topodismap.map[canumB][v->anumB[0]])>maxDIFtopodis)) 
             {topodis_ok = 0;}

          if (topodis_ok==1) CA_NOTstring_new[v->num] = CA_NOTstring[v->num];
        } 
      }
 
      /** [4] Do BronKerbosch(), recursively. **/
      ret = BronKerbosch_Pivot_BitString(optMlist,&CSnew,pairlist,Npair,CA_NOTstring_new,molA,molB,ConnectType,maxDIFtopodis);
      Free_MATCH(&CSnew);
      free(CA_NOTstring_new); 
      PAR.M += 1;
      if (ret==-1) return(-1);
      /** [5] Add c to 'NOT' **/
      CA_NOTstring[cnum] = 'n';
      N_CA   = Sum_Bit_String(CA_NOTstring,'c',Npair);
    }

    if (connect_p_and_c == 1)  i=i+1;

  } /* while */

  return(0);

} /* end of BronKerbosch_Pivot_BitString() */




void Evaluate_Degree_of_Atom_Pair(Mlist,molA,molB,ConnectType,maxDIFtopodis)
 struct MATCH    *Mlist;
 struct MOLECULE *molA,*molB;
 char   ConnectType;        /* 'D'isconnected, 'T'opo_dis-constrained MCS */
 int   maxDIFtopodis;      /* only for ConnectType=='T' */
{
 int pA,pB,qA,qB;
 struct MATCH *pn,*qn;
 char topodis_ok;

 /*
   Evaluated degree is stored in "pn->selct_dis".
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
     if ( (pA != qA) &&  ( pB != qB) &&
           (evaluate_edge_connection(molA,pA,qA,molB,pB,qB)==1) ) {
        topodis_ok = 1;
        if ((ConnectType=='T')
           &&(abs(molA->topodismap.map[pA][qA]-molB->topodismap.map[pB][qB])>maxDIFtopodis)) 
         {topodis_ok = 0;}
        
        if (topodis_ok==1) pn->select_dis += 1.0;
    }
   } 
 }

} /* end of  Evaluate_Degree_of_Atom_Pair() */



int Delete_Small_Degree_of_Atom_Pair(Mlist,thre_degree)
 struct MATCH    *Mlist;
 int    thre_degree; /* if pn->degree < thre_degree, delete the node "pn". */
{
 struct MATCH *pn;
 int Npair_rest;
 float max_degree,thre_degreeF;

 /*
   * Evaluated degree should be stored in "pn->selct_dis".
     by the function "Evaluate_Degree_of_Atom_Pair()".

   * If thre_degree < max_delete, then
     delete a node pn with pn->degree < thre_degree.

 */ 
 /*
 printf("#void Delete_Small_Degree_of_Atom_Pair(Mlist,molA,molB,thre_degree:%d)\n",thre_degree);
 */ 
 thre_degreeF = (float)thre_degree;
 
 Npair_rest = 0;
 max_degree = 0.0;
 pn = Mlist;
 while (pn->next != NULL){
   pn = pn->next;
   Npair_rest += 1;
   if (pn->select_dis > max_degree) max_degree = pn->select_dis;
 }
 

 if (thre_degree < max_degree){
   Npair_rest = 0;
   pn = Mlist;
   while (pn->next != NULL){
     pn = pn->next;
     if (pn->select_dis < thre_degreeF){
       if (pn->prev != NULL) pn->prev->next = pn->next;
       if (pn->next != NULL) pn->next->prev = pn->prev;
     }
     else { Npair_rest += 1;} 
   }
 }
 return(Npair_rest);
} /* end of Delete_Small_Degree_of_Atom_Pair() */




int Keep_Only_Maximum_Npair_MATCH(Mlist)
 struct MATCH    *Mlist;
{
 struct MATCH *pn,*qn;
 int Npair_max,N;

 printf("#Keep_Only_Maximum_Npair_MATCH(Mlist,molA,molB,min_degree)\n");
 /** (1) Find Npair_max **/
 pn = Mlist;
 Npair_max = 0;
 while (pn->next != NULL){
   pn = pn->next;
   if (pn->Npair>Npair_max) Npair_max = pn->Npair;
 }
 printf("#Npair_max %d\n",Npair_max);
 /** (2) Delete pn with pn->Npair< Npair_max **/
 pn = Mlist;
 N = 0;
 while (pn->next != NULL){
   pn = pn->next;
   if (pn->Npair<Npair_max) {
     qn = pn;
     if (pn->prev != NULL) pn->prev->next = pn->next;
     if (pn->next != NULL) pn->next->prev = pn->prev;
     pn = pn->prev;
     Free_MATCH(qn);
   }
   else N = N + 1;
 }
 return(N);
} /* end of Keep_Only_Maximum_Npair_MATCH() */



int Keep_Only_First_K_MATCH(Mlist,K)
 struct MATCH    *Mlist;
 int K;
{
 struct MATCH *pn,*qn;
 int N,Nkeep;

 printf("#Keep_Only_Maximum_Npair_MATCH(Mlist,molA,molB,min_degree)\n");
 pn = Mlist;
 N = Nkeep = 0; 
 while (pn->next != NULL){
   pn = pn->next;
   N = N + 1; 
   if (N>K) {
     qn = pn;
     if (pn->prev != NULL) pn->prev->next = pn->next;
     if (pn->next != NULL) pn->next->prev = pn->prev;
     pn = pn->prev;
     Free_MATCH(qn);
   }
  else {Nkeep += 1;}
 }

 return(Nkeep);
} /* end of Keep_Only_Maximum_Npair_MATCH() */


int Renumbering_MATCHlist(Mlist)
 struct MATCH    *Mlist;
{
 struct MATCH *pn;
 int N;

 pn = Mlist;
 Mlist->num = -1;
 N = 0; 
 while (pn->next != NULL){
   pn = pn->next;
   pn->num = N;
   N = N + 1; 
 }
 return(N);
} /* end of Renumbering_MATCHlist() */







void Init_Bit_String(bitstr,init_symbol,len)
  char *bitstr;
  char init_symbol;
  int  len;
{
  int i;
  for (i=0;i<len;++i) bitstr[i] = init_symbol;
} /* end of Init_Bit_String() */


int Sum_Bit_String(bitstr,target,len)
  char *bitstr;
  char target;  /* target_symbol */
  int  len;
{
  int i,N;
  N = 0;
  for (i=0;i<len;++i){
    if (bitstr[i]==target) N += 1;
  }
  return(N);
} /* end of Sum_Bit_String() */



void Show_Bit_String(bitstr,len,comment)
  char *bitstr;
  int  len;
  char *comment;
{
  int i;
  printf("#%s:",comment); 
  for (i=0;i<len;++i){ 
    if (bitstr[i]==' ') printf("%d_ ",i);
    else printf("%d%c ",i,bitstr[i]);
 }
  printf("\n");
} /* end of Show_Bit_String() */



int Get_First_Target_Symbol_In_Bit_String(bitstr,target,len)
  char *bitstr;
  char target; 
  int  len;
/*
 return the first appearing non zero item
 ex) bitstr = '001101', return(2)  
*/
{
  int i;
  for (i=0;i<len;++i){ 
    if (bitstr[i]==target) return(i);
  }
  return(-1);
} /* end of Get_First_Target_Symbol_In_Bit_String() */


int Get_i_th_Target_Symbol_In_Bit_String(bitstr,i,target,len)
  char *bitstr;
  int  i;
  char target;
  int  len;
/*
 return the i-th appearing non zero item
 ex) bitstr = '00110101' and k = 2, return(3)  
*/
{
  int j,k;
  k = 0;
  for (j=0;j<len;++j){ 
    if (bitstr[j]==target) { 
      if (k==i) return(j);
      k += 1;
      }
  }
  return(-1);

} /* end of Get_i_th_Target_Symbol_In_Bit_String() */




int evaluate_edge_connection(molA,uA,vA,molB,uB,vB)
  struct MOLECULE *molA;
  int    uA,vA; /* atom number for molA */
  struct MOLECULE *molB;
  int    uB,vB; /* atom number for molB */
{
   if ((molA->conmap.map[uA][vA]=='0') && (molB->conmap.map[uB][vB]=='0')) return(1);

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

