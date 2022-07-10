/*

 <DualAtmPair.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  for calculating dual atom pair descriptors

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
#include "moltopodis.h"
#include "ringblock.h"
#include "RingDesc.h"

/**** FUNCTIONS (GLOBAL) ***/
void Set_DualAtomPair_Descriptor();
extern void get_param_from_index_of_dual_atompair_descriptor();


/**** FUNCTIONS (LOCAL) ***/
static void count_cyclic_and_acyclic_bond_on_the_path();


void Set_DualAtomPair_Descriptor(mol)
  struct MOLECULE *mol;
{
  /*
 
   Desctiptor : [atom_type(a)] - [Ncyclic_bond(c)] - [Nacyclic-bond(d)] - [atom_type(b)]
         
    (Note that a >= b.)
   a = [0..A-1]
   b = [0..A-1]
   c = [0..C]
   d = [0..D]
  
   A = PAR.max_atomtype
   C = PAR.max_separation_cyclic 
   D = PAR.max_separation_acyclic 

   Length_of_descriptor = 1/2 * A(A+1)(C+1)(D+1)


   
   Apair = A*(A+1)/2

                    
   index = Apair*(C+1*)*d + Apair*c + a*(a+1)/2 + b


*/
  int a, b, c, d, aa, bb, i, j, k, ind,Apair;
  int *path, Natom_in_path,Ncyc,Nacyc;

 /*
  printf("#void Set_DualAtomPair_Descriptor(%s)\n",mol->name);
  */ 
  path = (int *)malloc(sizeof(int)*mol->Natom); 


  Apair = (PAR.max_atomtype  * (PAR.max_atomtype + 1))/2;
  PAR.max_atompair_descriptor   = Apair * (PAR.max_separation_cyclic + 1) * (PAR.max_separation_acyclic + 1);

 /*
  printf("#max_atomtype %d max_separation %d max_atomtype_descriptor %d\n",
     PAR.max_atomtype, PAR.max_separation, PAR.max_atompair_descriptor);

 *   printf("#max_atomtype %d Apair %d PAR.max_separation %d  max_atompair_descriptor %d\n",PAR.max_atomtype,Apair,PAR.max_separation, PAR.max_atompair_descriptor);  
 *    */

  mol->atompair_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor);

  for (k=0;k<PAR.max_atompair_descriptor;++k) mol->atompair_descriptor[k] = 0;

  if (mol->topodismap.malloced!=1) {printf("#ERROR:topodismap is not malloced.PROGRAM_TYPE %c\n",PAR.PROGRAM_TYPE); exit(1);}

  for (i=0;i<mol->Natom;++i){
    if (mol->atoms[i].one_char_ele != 'H'){
      aa =  Number_of_atomtype(mol->atoms[i].atomtype);
      for (j=i+1;j<mol->Natom;++j){
        if (mol->atoms[j].one_char_ele != 'H'){
          bb =  Number_of_atomtype(mol->atoms[j].atomtype);
          if (aa > bb){a = aa; b = bb;}else {a = bb; b = aa;}

          set_shortest_path_between_two_atoms(mol,i,j,path,&Natom_in_path);
          count_cyclic_and_acyclic_bond_on_the_path(mol,path,Natom_in_path,&Ncyc,&Nacyc);
          if ((Ncyc <= PAR.max_separation_cyclic) && (Nacyc <= PAR.max_separation_acyclic)){
            c = Ncyc;
            d = Nacyc;
            ind = Apair * (PAR.max_separation_cyclic +1) * d + Apair *c + a*(a+1)/2 + b;
            get_param_from_index_of_dual_atompair_descriptor(ind,&a,&b,&c,&d);
            if ((ind<PAR.max_atompair_descriptor)&&(mol->atompair_descriptor[ind]<255)){
               mol->atompair_descriptor[ind] += 1;
            }
          }
        }
      } /* j */
    }
  } /* i */


  free(path);

} /* end of Set_DualAtomPair_Descriptor() */



void count_cyclic_and_acyclic_bond_on_the_path(mol,path,Natom_in_path,Ncyc,Nacyc)
  struct MOLECULE *mol;
  int *path;
  int Natom_in_path;
  int *Ncyc;   /* Number of cyclic  bond (to be calculated) */
  int *Nacyc;  /* Number of acyclic bond (to be calculated) */
{
  int k,a0,a1;

  *Ncyc = 0;
  *Nacyc = 0;

  for (k=1;k<Natom_in_path;++k){
    a0 = path[k-1]; 
    a1 = path[k]; 
    if ((mol->atoms[a0].ringblock!=' ')&&(mol->atoms[a0].ringblock == mol->atoms[a1].ringblock)) *Ncyc = *Ncyc + 1;
    else *Nacyc = *Nacyc + 1;
  }

} /* end count_cyclic_and_acyclic_bond_on_the_path() */





void get_param_from_index_of_dual_atompair_descriptor(ind,a,b,c,d)
  int ind;    /* input index */
  int *a,*b;  /* atom type [0..PAR.max_atomtype-1]  */
  int *c;     /* Nbond_cyclic  [0 .. PAR.max_separation_cyclic]  */
  int *d;     /* Nbond_acyclic [0 .. PAR.max_separation_acyclic]  */
{
  int Apair,ApairC,P,Q;
  float a0;

/*
  ind = Apair*(C+1*)*d + Apair*c + a*(a+1)/2 + b
*/
  Apair = (PAR.max_atomtype  * (PAR.max_atomtype + 1))/2;
  ApairC = Apair * (PAR.max_separation_cyclic+1);
  
  *d = ind / ApairC;
  P  = ind - ApairC * (*d);
  *c = P/Apair;
  Q = P - Apair * (*c);  
  a0 = (-1.0 + sqrt(8.0*Q+1.0))/2.0;
  *a = (int)floor(a0);
  *b = Q - (*a)*((*a)+1)/2;
  /*
  printf("Apair %d ApairC %d index %d abcd %d %d %d %d\n",Apair,ApairC, ind,*a,*b,*c,*d);
  */

} /* end of get_param_from_index_of_dual_atompair_descriptor() */

