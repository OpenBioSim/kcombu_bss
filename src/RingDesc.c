/*

 <RingDesc.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 for calcularing ring descriptor


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "molprop.h"
#include "molring.h"


/** FUNCTIONS (GLOBAL) **/
void Set_Ring_Descriptor();
void show_ring_descriptor();
float similarity_bwn_ring_descriptors();
float similarity_bwn_oneatom_ring_descriptors();
int Number_of_ring_block_type();

/** FUNCTIONS (LOCAL) **/


/*** VARIABLES (LOCAL) **/

void Set_Ring_Descriptor(mol)
  struct MOLECULE *mol;
{
  struct RING *rn;
  int a,i,size,Natom_no_block;

  mol->ring_descriptor = (unsigned char *)malloc(sizeof(unsigned char) * PAR.max_ring_descriptor);  

  for (i=0;i<PAR.max_ring_descriptor;++i) mol->ring_descriptor[i] = 0;

  Natom_no_block = 0; 
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].ringblock == ' ')  Natom_no_block += 1;
  }

  /** (1) [0]:acyclic atoms **/
  if (Natom_no_block<256) mol->ring_descriptor[0] = Natom_no_block;
  else mol->ring_descriptor[0] = 255;
    
  /** (2) [1..max_ring-2]: rings **/
  rn = &(mol->HeadRing); 
  while (rn->next != NULL){
    rn = rn->next;
    size = rn->Natom;
    if (size > PAR.max_ring_size) size = PAR.max_ring_size;
    if (size>=3){
      mol->ring_descriptor[size-2] += 1;
    }
  }

  /** (3) [max_ring-1.. max_ring + max_block-4]: blocks **/
  rn = &(mol->HeadBlock); 
  while (rn->next != NULL){
    rn = rn->next;
    size = rn->Natom;
    if (size > PAR.max_block_size) size = PAR.max_block_size;
    if (size>=3){
      mol->ring_descriptor[size + PAR.max_ring_size - 4] += 1;
    }
  }

} /* end of Set_Ring_Descriptor() */







int Number_of_ring_block_type(s)
 char *s;
 /*
  "a"->0 "r3"->1 "r4"->2,... "r[R]"-> R-2,
  "b3"->R-1,"b4"->R,... "b[B]: B+R-4,
   otherwise return(-1).  
  */ 
{
  int L,size;
  char tail[32];

  L = strlen(s);
  if ((L==1)&&(s[0]=='a')) return(0);
  if (L>=2){
    Get_Part_Of_Line(tail,s,1,L-1);     
    size = atoi(tail);
    if ((s[0]=='r') && (size>=3) && (size<=PAR.max_ring_size))  return(size-2);
    if ((s[0]=='b') && (size>=3) && (size<=PAR.max_block_size)) return(size + PAR.max_ring_size-4);
  }
  return(-1);
} /* end of Number_of_ring_block_type() */


void show_ring_descriptor(desc)
  unsigned char *desc;
{
 int i;
 printf(" ring_desc:maxde %d maxring %d maxblk %d:",PAR.max_ring_descriptor,PAR.max_ring_size, PAR.max_block_size);
 for (i=0;i<PAR.max_ring_descriptor;++i){
    if (desc[i]>0){
            if (i==0) printf(" a ");
       else if (i<(PAR.max_ring_size-1)) printf(" r%d ",i+2);
       else printf(" b%d ",i-PAR.max_ring_size+4);
       printf("%d",desc[i]);
     }
  }
 printf("\n");
} /* end of show_ring_descriptor() */



float similarity_bwn_ring_descriptors(desA, desB)
   unsigned char *desA; /* ring descriptor for molecule A [PAR.max_ring_descriptor] */
   unsigned char *desB; /* ring descriptor for molecule B [PAR.max_ring_descriptor] */
{
  int i,common,A,B,a,b,size;
  float sim;

  A = B = common = 0;
  for (i=0;i<PAR.max_ring_descriptor;++i){
    if (i==0) size = 1;
    else if (i<(PAR.max_ring_size-1))     size = i + 2; 
    else if (i<(PAR.max_ring_descriptor)) size = i - PAR.max_ring_size + 4; 
    /* printf("i %d size %d\n",i,size); */
    a = desA[i] * size;
    b = desB[i] * size;
    A += a;
    B += b;
    if (a<b) common += a; else common += b;
  }
  if ((A + B - common)>0)
    sim = (float)common/((float)A + (float)B - (float)common);
  else
      sim = 0.0;
  return(sim);
} /* end of similarity_bwn_ring_descriptors() */



float similarity_bwn_oneatom_ring_descriptors(oneA,ringA,oneB,ringB,unitype)
   unsigned char *oneA;  /* one_atom descriptor for molecule A */
   unsigned char *ringA; /* ring     descriptor for molecule A */
   unsigned char *oneB;  /* one_atom descriptor for molecule B */
   unsigned char *ringB; /* ring     descriptor for molecule B */
   char unitype; /* 'A'verage, 'I':minimum */
{
  /*
   < number of one atom descriptor > 
   C:0 O:1 N:2 S:3 P:4
   C@:5 O@:6 N@:7 S@:8 
   O1:9 N1:10 other:11
  */
  int i,common,A,B,a,b,size;
  float sim;
  int Ncommon_acycle,Ncommon_cycle;
  int NatomA, NatomB, NcycleA, NcycleB;
  int Ncommon_cycle_one;
  int Ncycle_ringA,Ncycle_ringB, Ncommon_cycle_ring;
  int Ncycle_blockA,Ncycle_blockB,Ncommon_cycle_block;
  char rbtype;
/*
  printf("#float similarity_bwn_oneatom_ring_descriptors(oneA,ringA,oneB,ringB)\n");
*/
  /** (1) count common atom for one atom descriptor **/
  NatomA = NatomB = Ncommon_cycle = Ncommon_acycle = Ncommon_cycle_one = 0;
  NcycleA = NcycleB = 0; 
  for (a=0;a<PAR.max_atomtype;++a){
    NatomA += oneA[a];
    NatomB += oneB[a];
    if ((a<=4)||(a>=9)){
      if (oneA[a]<oneB[a]) Ncommon_acycle += oneA[a];
                   else    Ncommon_acycle += oneB[a];
    } 
    else{ 
      if (oneA[a]<oneB[a]) Ncommon_cycle_one += oneA[a];
                   else    Ncommon_cycle_one += oneB[a];
      NcycleA += oneA[a];
      NcycleB += oneB[a];
    } 
 }


  /** (2) count common atom for ring descriptor **/
  A = B = common = 0;
  Ncycle_ringA = Ncycle_ringB = Ncycle_blockA = Ncycle_blockB = 0;
  Ncommon_cycle_ring = Ncommon_cycle_block = 0; 
  for (i=1;i<PAR.max_ring_descriptor;++i){
         if (i<(PAR.max_ring_size-1))     {size = i + 2;  rbtype = 'r';}
    else if (i<(PAR.max_ring_descriptor)) {size = i - PAR.max_ring_size + 4;  rbtype = 'b';}
    /* printf("i %d size %d\n",i,size); */
    a = ringA[i] * size;
    b = ringB[i] * size;
    if (rbtype=='r'){
      Ncycle_ringA += a;
      Ncycle_ringB += b;
      if (a<b) Ncommon_cycle_ring += a; else Ncommon_cycle_ring += b;
    }
    else if (rbtype=='b'){
      Ncycle_blockA += a;
      Ncycle_blockB += b;
      if (a<b) Ncommon_cycle_block += a; else Ncommon_cycle_block += b;
    }
  }
/*
  if ((NcycleA>0) && (NcycleB>0)){
       scale = (float)(NcycleA+NcycleB)/(float)(Ncycle_ringA+Ncycle_ringB);
  }
  printf("NcycleA %d NcycleB %d Ncycle_ringA %d Ncycle_ringB %d Ncommon_ring_orig %d\n",NcycleA,NcycleB,Ncycle_ringA, Ncycle_ringB,Ncommon_cycle_ring);
  printf("Ncommon_acycle %d Ncommon_cycle %d %d %d scale %f\n",Ncommon_acycle,Ncommon_cycle_one,Ncommon_cycle_ring,Ncommon_cycle_block,scale); 
*/
  if (Ncommon_cycle_ring>Ncommon_cycle_one) Ncommon_cycle_ring = Ncommon_cycle_one;
 
  if (unitype=='A'){
    Ncommon_cycle = (Ncommon_cycle_one + Ncommon_cycle_ring + Ncommon_cycle_block)/3;
  }

  if (unitype=='I'){
    if ((Ncommon_cycle_one <= Ncommon_cycle_ring) && (Ncommon_cycle_one <= Ncommon_cycle_block))
      Ncommon_cycle = Ncommon_cycle_one;
    if ((Ncommon_cycle_ring <= Ncommon_cycle_one) && (Ncommon_cycle_ring <= Ncommon_cycle_block))
      Ncommon_cycle = Ncommon_cycle_ring;
    if ((Ncommon_cycle_block <= Ncommon_cycle_one) && (Ncommon_cycle_block <= Ncommon_cycle_ring))
      Ncommon_cycle = Ncommon_cycle_block;
  }


  common = Ncommon_acycle + Ncommon_cycle;
 
  if ((NatomA + NatomB - common)>0)
    sim = (float)common/((float)NatomA + (float)NatomB - (float)common);
  else
    sim = 0.0;
  
  /* printf("#similatity %f\n",sim); */
   
  return(sim);
} /* end of similarity_bwn_oneatom_ring_descriptors() */

