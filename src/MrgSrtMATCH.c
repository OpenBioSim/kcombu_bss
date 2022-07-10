/*

 <MrgSrtMATCH.c>
 
 functions for MergeSort for linked linear list. 

 coded by Takeshi Kawabata.

 Reference :
 Yoshiyuki Kondo. "Algorithms and Data Structures for C programmers".
  Soft Bank Publishing (1998). written in Japanese.  

 <Standard Definition of node 'MATCH'> 

 struct MATCH{
  struct MATCH *next;
  struct MATCH *prev;
  float  select_dis;
 };


 <Standard Usage> 
 
 struct MATCH HEAD,*bn;
 
 bn = Merge_Sort_List_MATCH(Head.next,'-'); 
 Head.next = bn;
   or 
 Merge_Sort_Double_Linked_List_MATCH(head,type)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "molecule.h" 

/*** FUNCTIONS (GLOBAL) ***/
struct MATCH *Merge_Sort_List_MATCH();
void Merge_Sort_Double_Linked_List_MATCH();

/*** FUNCTIONS (LOCAL) ***/
static void Set_Prev_Pointer();
static struct MATCH *Merge_List();
static void Show_NODE();


void Show_NODE(head)
 struct MATCH *head;
{ 
 struct MATCH *bn;      

 /* bn = head->next;   */
 bn = head;
 while (bn != NULL){ 
   printf("%f\n",bn->select_dis); 
   bn = bn->next; 
 }

} /* end of Show_NODE() */



void Merge_Sort_Double_Linked_List_MATCH(head,type)
 struct MATCH *head;
 char type; /* if 'D', large to small, otherwise small to large */ 
{
 struct MATCH *bn;
 bn = Merge_Sort_List_MATCH(head->next,type); 
 head->next = bn;
 Set_Prev_Pointer(head);

} /* end of Merge_Sort_Double_Linked_List_MATCH() */



struct MATCH *Merge_List(a,b,type)
 struct MATCH *a,*b;
 char type; /* if 'D', decreasing order large to small, otherwise small to large */ 
{
 struct MATCH head, *p;

 /* Initialize */ 
 p = &head;

 /* Repeat until lists 'a' and 'b' become empty */

 /* The smaller node is attached to the list 'p' */
 if (type != 'D') 
 while ((a!=NULL)&&(b!=NULL)){
  if (a->select_dis <= b->select_dis){p->next = a; p = a; a = a->next; }
                       else  {p->next = b; p = b; b = b->next; }
 }  

 else
 
 /* The larger node is attached to the list 'p' */
 while ((a!=NULL)&&(b!=NULL)){
  if (a->select_dis >= b->select_dis){p->next = a; p = a; a = a->next; }
                       else  {p->next = b; p = b; b = b->next; }
 }  


 /* The rest nodes are attached to the list 'p'. */
       if ((a == NULL) && (b !=NULL)) p->next = b; 
 else if  ((b == NULL) && (a !=NULL)) p->next = a; 
   
 return(head.next);

} /* end of  Merge_List() */






struct MATCH *Merge_Sort_List_MATCH(x,type)
 struct MATCH *x;
 char type;
{
 struct MATCH *a, *b, *p;


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

 return(Merge_List(Merge_Sort_List_MATCH(x,type),Merge_Sort_List_MATCH(p,type),type));


} /* end of *Merge_Sort_List_MATCH() */




void Set_Prev_Pointer(head)
 struct MATCH *head;
{
 struct MATCH *bn,*pn;

 bn = head;
 while (bn->next !=NULL){
   pn = bn;
   bn = bn->next;
   bn->prev = pn; 
  }

} /* end of Set_Prev_Pointer() */

