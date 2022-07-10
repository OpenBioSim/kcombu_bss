/*

 <MrgSrtLIBMOL.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
 functions for MergeSort for linked linear list. 

 coded by Takeshi Kawabata.

 Reference :
 Yoshiyuki Kondo. "Algorithms and Data Structures for C programmers".
  Soft Bank Publishing (1998). written in Japanese.  

 <Standard Definition of node 'LIBMOL'> 

 struct LIBMOL{
  struct LIBMOL *next;
  struct LIBMOL *prev;
  float  score_for_sort;
  float  alph_score;
 };


 <Standard Usage> 
 
 struct LIBMOL HEAD,*bn;
 
 bn = Merge_Sort_List_LIBMOL(Head.next,'-'); 
 Head.next = bn;
   or 
 Merge_Sort_Double_Linked_List_LIBMOL(head,type)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2DMAP.h" 
#include "molecule.h" 
#include "ioLib.h" 

/*** FUNCTIONS (GLOBAL) ***/
struct LIBMOL *Merge_Sort_List_LIBMOL();
void Merge_Sort_Double_Linked_List_LIBMOL();
void Insertion_Sort_Double_Linked_List_LIBMOL();

/*** FUNCTIONS (LOCAL) ***/
static void Set_Prev_Pointer();
static struct LIBMOL *Merge_List();
static void Show_NODE();


void Show_NODE(head)
 struct LIBMOL *head;
{ 
 struct LIBMOL *bn;      

 /* bn = head->next;   */
 bn = head;
 while (bn != NULL){ 
   printf("%f\n",bn->score_for_sort); 
   bn = bn->next;
 }

} /* end of Show_NODE() */



void Merge_Sort_Double_Linked_List_LIBMOL(head,type)
 struct LIBMOL *head;
 char type; /* if 'D'ecreasing, large to small, otherwise small to large */ 
{
 struct LIBMOL *bn;
 bn = Merge_Sort_List_LIBMOL(head->next,type); 
 head->next = bn;
 Set_Prev_Pointer(head);

} /* end of Merge_Sort_Double_Linked_List_LIBMOL() */



struct LIBMOL *Merge_List(a,b,type)
 struct LIBMOL *a,*b;
 char type; /* if 'D', decreasing order large to small, otherwise small to large */ 
{
 struct LIBMOL head, *p;

 /* Initialize */ 
 p = &head;

 /* Repeat until lists 'a' and 'b' become empty */

 /* The smaller node is attached to the list 'p' */
 if (type != 'D') 
 while ((a!=NULL)&&(b!=NULL)){
       if ( (a->score_for_sort < b->score_for_sort)||
            ((a->score_for_sort == b->score_for_sort)&&(a->alph_score < b->alph_score)))
        {p->next = a; p = a; a = a->next; }
  else  {p->next = b; p = b; b = b->next; }
 }  

 else
 
 /* The larger node is attached to the list 'p' */
 while ((a!=NULL)&&(b!=NULL)){
    if ((a->score_for_sort >  b->score_for_sort) || 
        ((a->score_for_sort == b->score_for_sort)&&(a->alph_score > b->alph_score)))
         {p->next = a; p = a; a = a->next; }
   else  {p->next = b; p = b; b = b->next; }
 }  


 /* The rest nodes are attached to the list 'p'. */
       if ((a == NULL) && (b !=NULL)) p->next = b; 
 else if  ((b == NULL) && (a !=NULL)) p->next = a; 
   
 return(head.next);

} /* end of  Merge_List() */






struct LIBMOL *Merge_Sort_List_LIBMOL(x,type)
 struct LIBMOL *x;
 char type;
{
 struct LIBMOL *a, *b, *p;


 /* if x has no nodes, or one nodes */
 if ((x==NULL) || (x->next ==NULL)) return(x); 

 /* a is the 1st node, and 'b' is the 2nd or 3rd node */
 a = x;
 b = x->next;
 if (b->next !=NULL) b = x->next;

 /* 
   Progress a by 1 node, and 'b' by two nodes. 
   if 'b' reaches to the end, 'a' reaches to the midpoint.
 */

 while (b->next !=NULL){
  a = a->next;
  b = b->next;
  if (b->next !=NULL) b = b->next;
  }

 /* Truncate the list on the midpoint */
 p = a->next; a->next = NULL;

 return(Merge_List(Merge_Sort_List_LIBMOL(x,type),Merge_Sort_List_LIBMOL(p,type),type));


} /* end of *Merge_Sort_List_LIBMOL() */




void Set_Prev_Pointer(head)
 struct LIBMOL *head;
{
 struct LIBMOL *bn,*pn;
 int N;

 N = 0;
 bn = head;
 while (bn->next !=NULL){
   pn = bn;
   bn = bn->next;
   /* bn->num = N; */
   N += 1;
   bn->prev = pn; 
  }

} /* end of Set_Prev_Pointer() */


void Insertion_Sort_Double_Linked_List_LIBMOL(head,Ntop,type)
 struct LIBMOL *head;
 int    Ntop;    /* number of item for sorting */ 
 char   type;    /* if 'D'ecreasing, large to small, otherwise small to large */ 
{
 struct LIBMOL *x; /* focus */
 struct LIBMOL *h; /* head */
 struct LIBMOL *t; /* top */
 struct LIBMOL *b; /* buffer */
 int n,Nnode,ntop;
 float opt_value;

 printf("#Insertion_Sort_Double_Linked_List_LIBMOL(head,Ntop:%d, type:%c)\n",Ntop,type);

 /* [0] Check number of all the node < Ntop */
 Nnode = 0;
 x = head;
 while (x->next != NULL){
   x = x->next;
   Nnode += 1;
 }
 ntop = Ntop; 
 if (Ntop>Nnode){ ntop = Nnode;}

 /* [1] Do insertion_sort */
 h = head;  
 for (n=0;n<ntop;++n){ 
   x = h;
   t = NULL;
   opt_value = 0.0;

   while (x->next != NULL){
     x = x->next;
     if ((t==NULL)||
         ((type=='D')&&(x->score_for_sort > opt_value)) ||
         ((type!='D')&&(x->score_for_sort < opt_value)) ){
       t = x;
       opt_value = x->score_for_sort;
     }
   }
   /* printf("[%d] opt_value %f\n",n,opt_value);  */
   b = h->next;
   t->prev->next = t->next;
   if (t->next != NULL){
     t->next->prev = t->prev;
   } 
 
   h->next = t;
   t->next = b;
   t->prev = h;
   if (b != NULL){ b->prev = t;}
   h = t;    
 }

 /* 
 x = head;
 n = 0;
 while ((x->next != NULL) && (n < ntop*10)){
   x = x->next;
   printf("[%d] %f\n",n,x->score_for_sort);
   ++n;
 }
 */

} /* end of Insertion_Sort_Double_Linked_List_LIBMOL() */
