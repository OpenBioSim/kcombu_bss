/*

 <ringblock.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 for assingning "ring" and "block" of ring structure for the moluecule.


 The definition of the ring is based on the concept SSSR (Smallest Set of the Smallest Rings).

 The algorithm is based on the following article:

 Bo Tao Fan, Annick Panaye, Jean-Pierre Doucet, Alain Barbu.
 Ring percetion. A new algorithm for directly finding the smallest set of smallest rings
 from a connection table.  J.Chem.Inf.Comput.Sci., 1993, 33, 657-662.

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
void Set_Ring_Block_and_SSSR_Ring();

/** FUNCTIONS (LOCAL) **/
static int set_mark0_for_open_acyclic_atoms();
static void set_ring_blocks();
static void set_SSSR_rings_for_each_block();
static void Set_Property_using_BLOCK_and_RING();
static void set_NC_in_the_block();
static void delete_atom_mark_for_N2_atoms_in_smallest_ring();
static void add_new_RING_to_tail();
static int mark_neighbors_for_all_the_rings();
static int mark_neighbors_for_smallest_ring();
static void restore_hyding_atoms();
static void show_ringatms_with_NC();
static int is_metal_element();


/*** VARIABLES (LOCAL) **/
static int  *NATOM_IN_PATH;  /* [mol->Natom] (malloc later) : for finding circle */
static int  *NC;             /* [mol->Natom] (malloc later) : Number of connectivity */
static int  *BLOCK_NUM;      /* [mol->Natom] (malloc later) : for finding SSSR in the block*/
static char *ATOM_MARK;      /* [mol->Natom] (malloc later) : 1:in the block, -1: out of the block 0:deleted */
static int *EQUIV_STA_ATM;   /* [mol->Natom] (malloc later) (for merging overlapping rings for finding blocks */

static int N_RECURRENT;
static int N_RECURRENT_MAX;


void Set_Ring_Block_and_SSSR_Ring(mol)
 struct MOLECULE *mol;
{
  /* 
  printf("#Set_Ring_Block_and_SSSR_Ring('%s')\n",mol->filename); fflush(stdout);
  printf("#Nheavyatom %d Nheavybond %d Vmax %d\n",mol->Nheavyatom, mol->Nheavybond, 2*mol->Nheavyatom - mol->Nheavybond);
  */
  N_RECURRENT_MAX = 1000000;
 
  N_RECURRENT = 0;
  set_mark0_for_open_acyclic_atoms(mol);

 /*
  printf("#N_RECURRENT %d\n",N_RECURRENT);
  for (i=0;i<mol->Natom;++i){
    printf("[%d] mark %d\n",i,mol->atoms[i].mark);
  }
 */

  N_RECURRENT = 0;
  set_ring_blocks(mol);
  
  N_RECURRENT = 0;
  set_SSSR_rings_for_each_block(mol);
/*
  printf("#N_RECURRENT %d N_RECURRENT_MAX %d\n",N_RECURRENT,N_RECURRENT_MAX);
*/
 
  Set_Property_using_BLOCK_and_RING(mol);

/*
  Show_Ring_Structure(mol,&(mol->HeadBlock),"BLOCK");
  Show_Ring_Structure(mol,&(mol->HeadRing),"RING");
*/

} /* end Set_Ring_Block_and_SSSR_Ring() */




void Set_Property_using_BLOCK_and_RING(mol)
 struct MOLECULE *mol;
{
 int i,j;
 struct RING *rn, *blk;
 char abcstr[100],metal;
 sprintf(abcstr,"abcdefghijklmnopqrstuvwxyz1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZ");

  for (i=0;i<mol->Natom;++i){
    mol->atoms[i].ring     = '-';
    mol->atoms[i].aromatic = '-';
    mol->atoms[i].ringblock = ' ';
  }

  /** [1] Using the Blocks, Set mol->atoms[].ring and mol->atoms[].block **/
 
  blk = &(mol->HeadBlock);
  while (blk->next != NULL){
    blk = blk->next;
/*
    printf("blk %d Natom %d\n",blk->num,blk->Natom);
*/
    for (i=0;i<blk->Natom;++i){
      mol->atoms[blk->num_atoms[i]].ring = 'r';
      mol->atoms[blk->num_atoms[i]].ringblock = abcstr[blk->num%strlen(abcstr)];
    } 

  }


 /** [2] Using the SSSR Rings, Set 'r'ing type bond  and mol->atoms[].aromatic**/
 /*
   Condition of aromatic ring: 
    (i) Number of atom <= PAR.maxNatom_in_aromatic_ring (about eight)
    (ii) Not contain any metal atoms.
         ( this is mainly for comparison between HEM and TBV (open HEM without iron) ).
 */ 
  rn = &(mol->HeadRing);
   while (rn->next != NULL){
     rn = rn->next;
     metal = 0;
     for (i=0;i<rn->Natom;++i){ 
        j = rn->num_atoms[i];
        /*
        printf(" (%s %d %d)", mol->atoms[j].element, mol->atoms[j].num_in_file,
         is_metal_element(mol->atoms[j].element));
        */
        if (is_metal_element(mol->atoms[j].element)==1) metal = 1;
     }

     if ((metal==0)&&(rn->Natom <= PAR.maxNatom_in_aromatic_ring)){
       for (i=0;i<rn->Natom;++i){ 
         if (mol->atoms[rn->num_atoms[i]].aromatic == '-') {mol->atoms[rn->num_atoms[i]].aromatic = 'A';}
       }
     }
  }

} /* end of Set_Property_using_BLOCK_and_RING() */







int set_mark0_for_open_acyclic_atoms(mol)
  struct MOLECULE *mol;
{
 /*
  > GOAL of this function 
    atoms[].mark = 0 :for open acyclic atoms and non heavy atoms.
    atoms[].mark = 1 :for other atoms
 
  >> CAUTION <<
   "mol->atoms[].Nnei_heavy"  and "mol->atoms[].neighbors[]" should be assigned before
   perfoming this function.

 */
 int a,i,n,nround,Nacyclic,Nmark0;
 int *Nnei_trim; /* Nneighbor for trimed molecule */

 /* printf("#set_mark0_for_open_acyclic_atoms()\n"); */


 Nnei_trim = (int *)malloc(sizeof(int)*mol->Natom);

 /* [1] Initialization (Nnei_trim = Nnei_heavy )*/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].one_char_ele != 'H'){
     mol->atoms[a].mark = 1;
     Nnei_trim[a] =  mol->atoms[a].Nnei_heavy;
   }
 }


 /* [2] Repeat triming atoms with Nneighbor <=1  */
 Nacyclic = 0;
 nround = 0; 
 do{
   Nmark0  = 0;
   
   /* (2-1) Assign mark:=0  for Nnei_trim <= 1*/
   for (a=0;a<mol->Natom;++a){
     if ((mol->atoms[a].one_char_ele!='H')&&(mol->atoms[a].mark==1) && (Nnei_trim[a]<=1)){
       mol->atoms[a].mark = 0; 
       Nmark0 += 1;
       Nacyclic += 1;
     }
  }

   /* (2-2) Recalculate Nneighbor and stored in Nnei_trim */
   for (a=0;a<mol->Natom;++a){
     Nnei_trim[a] = 0;
     for (i=0;i<mol->atoms[a].Nneighbor;++i){
       n = mol->atoms[a].neighbors[i];
       if ((mol->atoms[n].one_char_ele!='H') && (mol->atoms[n].mark==1)){Nnei_trim[a] += 1;} 
     }
   }
   /* printf("#round [%d] Nmark0 %d\n",nround,Nmark0); */
   nround += 1;
 } while (Nmark0>0);

 free(Nnei_trim);
 return(Nacyclic);

} /* end of set_mark0_for_open_acyclic_atoms() */





void set_ring_blocks(mol)
  struct MOLECULE *mol;
{
 /*
  ** Initial state
    atoms[].mark = 0 : for open acyclic atoms and non heavy atoms.
    atoms[].mark = 1 : for other atoms

  ** NATOM_IN_PATH[]      is for storing number of atoms in the path.  
 
  >> CAUTION <<
   The function "set_mark0_for_open_acyclic_atoms()" 
            and "set_mark0_for_closed_acyclic_atoms()" should be perfomed 
   before executing this function.
 */
 int a,b;
 int Nringatms,Nringatms_smlt;
 int *ringatms;   /* [mol->Natom] (malloc later)  */
 int *ringatms_smlt;   /* [mol->Natom] (malloc later)  */
 int *STA_ATM_BLOCK;   /* [mol->Natom] (malloc later)  */
 int overlap_min_sta; 


 /* printf("#set_ring_blocks()\n"); */

 /** [1] initialize various arrays **/
 NATOM_IN_PATH      = (int *)malloc(sizeof(int)*mol->Natom);
 EQUIV_STA_ATM       = (int *)malloc(sizeof(int)*mol->Natom);
 STA_ATM_BLOCK = (int *)malloc(sizeof(int)*mol->Natom);
 ringatms           = (int *)malloc(sizeof(int)*mol->Natom);
 ringatms_smlt      = (int *)malloc(sizeof(int)*mol->Natom);
 ATOM_MARK          = (char *)malloc(sizeof(char)*mol->Natom);
 for (a=0;a<mol->Natom;++a){ 
   STA_ATM_BLOCK[a] = -1;
   ATOM_MARK[a] = 1;
 }
 mol->HeadBlock.next = NULL;
 mol->HeadBlock.prev = NULL;
 mol->HeadBlock.num  = -1;
 mol->HeadBlock.Natom  = 0;
 mol->HeadBlock.num_atoms  = NULL;
 mol->Nringblock = 0;

 
 /** [2] Calculate all the rings throuth atom a, and set EQUIV_STA_ATM[x]=a  **/
 for (a=0;a<mol->Natom;++a){
   /*
   printf("[%d] mark %d start_atom_in_ring %d\n",mol->atoms[a].num_in_file,mol->atoms[a].mark, EQUIV_STA_ATM[a]);
   */
   if ((mol->atoms[a].mark==1) && (STA_ATM_BLOCK[a]==-1)){
 
     /* (2-A) Finding the smallest ring through atom 'a', using the recurrent function **/

     Nringatms_smlt = 0;
     for (b=0;b<mol->Natom;++b){ NATOM_IN_PATH[b] = 0; EQUIV_STA_ATM[b] = -1;}
     N_RECURRENT = 0;


     /* (2-B) Set STA_ATM_BLOCK[] and EQUIV_STA_ATM[].
 
        For the ring atom r starting from the atom a, STA_ATM_BLOCK[r] = a
        If the ring atoms starting from the atom a are overlaped with the ring staring from the atom b,
        then EQUIV_STA_ATM[b] := a.
 
    */
     mark_neighbors_for_all_the_rings(mol,a,a,1,ringatms);
     overlap_min_sta  = a;
     for (b=0;b<mol->Natom;++b){ 
       if ((EQUIV_STA_ATM[b]>=0)&&(STA_ATM_BLOCK[b]>=0)){
         if (STA_ATM_BLOCK[b]<overlap_min_sta) overlap_min_sta = STA_ATM_BLOCK[b];
         EQUIV_STA_ATM[STA_ATM_BLOCK[b]] = a;
       }
       if (EQUIV_STA_ATM[b]>=0) STA_ATM_BLOCK[b] = a;
     }

/*
     for (c=0;c<mol->Natom;++c){ 
       printf("%d sta_atm_num %d final_sta_atm %d\n",c,EQUIV_STA_ATM[c],STA_ATM_BLOCK[c]);
      }
*/

/*
     mark_neighbors_for_smallest_ring(mol,-1,a,a,1,ringatms,&Nringatms_smlt,ringatms_smlt);

     EQUIV_STA_ATM[a] = a;
     overlap_min_sta  = a;
     for (c=0;c<Nringatms_smlt;++c){ 
       r = ringatms_smlt[c];
       if (STA_ATM_BLOCK[r]>=0){
          if (STA_ATM_BLOCK[r]<overlap_min_sta) overlap_min_sta = STA_ATM_BLOCK[r];
          EQUIV_STA_ATM[STA_ATM_BLOCK[r]] = a;
       }
       STA_ATM_BLOCK[r] = a;
     }
*/



     /*
       (2-C) Making unified STA_ATM_BLOCK[]  

        If EQUIV_STA_ATM[b] ==a,  the ring starting from the atom b should be unified into 
        the ring staring from the atom a. 

     */
     for (b=0;b<mol->Natom;++b){
       if ((STA_ATM_BLOCK[b]>0) && (EQUIV_STA_ATM[STA_ATM_BLOCK[b]]==a)){
       STA_ATM_BLOCK[b] = overlap_min_sta;
       }
     }

   }
 }

 /** [3] Making "blocks" as the list of structure RING **/
 mol->HeadBlock.next = NULL;
 mol->HeadBlock.num = -1;
 for (a=0;a<mol->Natom;++a){
   if (STA_ATM_BLOCK[a] == a){
     Nringatms = 0;
     for (b=0;b<mol->Natom;++b){
       if (STA_ATM_BLOCK[b] == a){
         ringatms[Nringatms] = b;
         Nringatms += 1;
       }
     }
 
    add_new_RING_to_tail(&(mol->HeadBlock), Nringatms, ringatms);
   }
 }

 mol->Nringblock = Number_of_RING(&(mol->HeadBlock));

 free(ATOM_MARK);
 free(ringatms_smlt);
 free(ringatms);
 free(STA_ATM_BLOCK);
 free(EQUIV_STA_ATM);
 free(NATOM_IN_PATH);

} /* end of set_ring_blocks() */






void set_SSSR_rings_for_each_block(mol)
  struct MOLECULE *mol;
{
  struct RING *blk;
  int a,i,min_a,minNC;
  int *ringatms, *ringatms_smlt;
  int Nringatms_smlt; 
  int Natom_blk, Nedge_blk,Nround,Nround_max;
  char enddowhile; 
  /*
  * ATOM_MARK[] for {deleted:0,non-deleted:1, out_of_the_block:-1}  
  */

/*
  printf("#set_SSSR_rings_for_each_block(%s Nringblock %d)\n",mol->name,mol->Nringblock);
*/

  /** [1] Initialize variables **/
  ringatms       = (int *)malloc(sizeof(int)*mol->Natom);
  ringatms_smlt  = (int *)malloc(sizeof(int)*mol->Natom);
  NC             = (int *)malloc(sizeof(int)*mol->Natom);
  NATOM_IN_PATH  = (int *)malloc(sizeof(int)*mol->Natom);
  BLOCK_NUM      = (int *)malloc(sizeof(int)*mol->Natom);
  ATOM_MARK      = (char *)malloc(sizeof(char)*mol->Natom);

  /* ATOM_MARK[a] : 1:inside of the blocks, 0: outside of the blocks, 
                    2:hyding in the even round, restoring in the next even round.
                    3:hyding in the odds round, restoring in the next odds round.
  */

  /** [1] Initializing **/
  for (a=0;a<mol->Natom;++a) {
    BLOCK_NUM[a] = -1;
    ATOM_MARK[a] = 0;
  }
  blk = &(mol->HeadBlock);
  while (blk->next != NULL){
    blk = blk->next;
    for (i=0;i<blk->Natom;++i){
      a = blk->num_atoms[i];
      BLOCK_NUM[a] = blk->num;
      ATOM_MARK[a] = 1;
    }
  }

  mol->HeadRing.next = NULL;
  mol->HeadRing.prev = NULL;
  mol->HeadRing.num_atoms = NULL; 
  mol->HeadRing.Natom  = 0;
  mol->HeadRing.num    = -1;

  /** [2] Making the smallest rings for each block blk ***/
  Nround = 0;
  blk = &(mol->HeadBlock);
  while (blk->next != NULL){
    blk = blk->next;
    Nround_max = blk->Natom;
    /*
    printf(">>>BLOCK [%d]",blk->num);
    for (i=0;i<blk->Natom;++i){printf(" %d",mol->atoms[blk->num_atoms[i]].num_in_file); }
    printf("\n");
    */

    enddowhile = 0;
    Nround = 0;
    do{

      Nround += 1;
      restore_hyding_atoms(mol,Nround,ATOM_MARK);
      set_NC_in_the_block(mol,blk,ATOM_MARK,NC,&Natom_blk,&Nedge_blk,&min_a,&minNC);
  
      /* 
      printf("[round %d] blk %d Natom %d Nedge %d Nring %d minNC %d min_a %d\n",
         Nround,blk->num,Natom_blk, Nedge_blk, Nedge_blk - Natom_blk+1,minNC,mol->atoms[min_a].num_in_file);
    */

      /** (2-A) : If no atom is left --> enddowhile **/
      if ((Natom_blk==0)||(minNC==0)){
        /* printf("#no non-hyding candidate atom is left.\n"); */
        enddowhile = 1;
      } 
      /** (2-B) : If minimum NC atoms has NC=1 --> delete it. **/
      else if (minNC<2){
        /* printf("#Atom[%d] with NC %d has no ring.\n",mol->atoms[min_a].num_in_file,minNC);  */
        ATOM_MARK[min_a] = 0;
      }
      /** (2-C) : Otherwise, finding the smallest ring by starting min_a atom. **/
      else{
        Nringatms_smlt = 0;
        for (a=0;a<mol->Natom;++a) {NATOM_IN_PATH[a] = 0;}
  
        mark_neighbors_for_smallest_ring(mol,blk->num,min_a,min_a,1,ringatms,&Nringatms_smlt, ringatms_smlt);

        /* show_ringatms_with_NC(mol, Nringatms_smlt, ringatms_smlt,NC,"smallest_ring"); */
        if (Nringatms_smlt==0){
          /* printf("#No ring is found !!\n");  */
          ATOM_MARK[min_a] = 0;
        }
        else{
          add_new_RING_to_tail(&(mol->HeadRing), Nringatms_smlt, ringatms_smlt);
          delete_atom_mark_for_N2_atoms_in_smallest_ring(mol,blk,ATOM_MARK, Nringatms_smlt, ringatms_smlt,NC,Nround);
        }
      }
  /*
      set_NC_in_the_block(mol,blk,ATOM_MARK,NC,&Natom_blk,&Nedge_blk,&min_a,&minNC);
      printf("[after round %d] blk %d Natom %d Nedge %d Nring %d\n",Nround,blk->num,Natom_blk, Nedge_blk, Nedge_blk - Natom_blk +  1);
  */
    } while ((Natom_blk>0)&&(enddowhile==0)&&(Nround<Nround_max)); 

   if (Nround==Nround_max){
     printf("#WARNING(%s block %d/%d):set_SSSR_rings stops because of too many round over %d\n",
       mol->filename,blk->num,mol->Nringblock,Nround);
   }
 
 } /* blk */


 mol->Nring = Number_of_RING(&(mol->HeadRing));
 free(ATOM_MARK);
 free(BLOCK_NUM);
 free(NATOM_IN_PATH);
 free(NC);
 free(ringatms_smlt);
 free(ringatms);

} /* end of set_SSSR_rings_for_each_block() */





void set_NC_in_the_block(mol,blk,atom_mark,nc, Natom, Nedge, min_a, minNC)
  struct MOLECULE *mol; 
  struct RING *blk; 
  char *atom_mark;  /* [mol->Natom] 1:inside, 0:outside, 2:hyding */
  int *nc;          /* number of connectivity [mol->Natom] */
  int *Natom;       /* (to be calculated) */
  int *Nedge;       /* (to be calculated) */
  int *min_a;       /* atom number with minimum NC (to be calculated) */
  int *minNC;       /* minimum NC (to be calculated) */
{
  int i,j,a,b;
  char init;
 
  /** (1) Calulate NC (number of connectivity)  **/
  for (a=0;a<mol->Natom;++a) {nc[a] = 0;}
  (*Natom) = (*Nedge) = 0;
  for (i=0;i<blk->Natom;++i){
    a = blk->num_atoms[i];
    if (atom_mark[a]==1){
      (*Natom) += 1;
      for (j=i+1;j<blk->Natom;++j){
        b = blk->num_atoms[j];
        if ((atom_mark[b]==1)&&(mol->conmap.map[a][b]!='0')){
          NC[a] += 1;
          NC[b] += 1;
          (*Nedge) += 1;
        }
      }
    }
  }

 /** (2) Calulate minNC and min_a **/
 init = 1; *minNC = 0; *min_a = 0;
 for (i=0;i<blk->Natom;++i){
   a = blk->num_atoms[i];
   if (atom_mark[a]==1){
     if ((init==1)||(nc[a] < *minNC)){ 
       *minNC = nc[a]; *min_a = a; init = 0;
     }
   }
 }

} /* end of set_NC_in_the_block() */





void delete_atom_mark_for_N2_atoms_in_smallest_ring(mol,blk,atom_mark, Nringatms, ringatms,nc,Nround)
  struct MOLECULE *mol;
  struct RING *blk;
  char *atom_mark;   /* [mol->Natom]  (to be changed ) */
  int Nringatms;     /* number of atoms in the smallest ring */ 
  int *ringatms;     /* array  of atoms in the smallest ring (circular) */ 
  int *nc;           /* number of connectivity [mol->Natom] */
  int Nround;   
{
  int i,as,ae,N2,ip,i0,in; 
  char state;
  int *sta_sg2, *end_sg2,*len_sg2, Nsg,s,max_len,max_s;
  /*

  Basically, the atom a with NC=2 are deleted in this procedure, in other words,
  atom_mark[a] with a with NC=2 is changed from "1" into "0".
  Its detailed procedure is as follow:
  [1] "2-segments" data (sta_sg2[],end_sg2[],len_sg2[]) is generated from the
      circular atom number sequence "ringatms[]". Because of circularity, sta_sg2[] is
      sometimes larger than end_sg2[].
  
  >> example of "2"-sequenctial segment << 
   "222333" -> "xxx333"
   sta_sg2[0] = 0, end_sg2[0] = 2;
 
   "322323" -> "3xx323" 
   sta_sg2[0] = 1, end_sg2[0] = 2;
   sta_sg2[1] = 4, end_sg2[1] = 4;
   
  "222223" -> "xxxxx3" 
   sta_sg2[0] = 0, end_sg2[0] = 4;

  "233222" - > "x33xxx" 
   sta_sg2[0] = 3, end_sg2[0] = 0;

  [2] Delete (atom_mark[]-->0) for the atoms in the largest 2-segment. 

  [3] If there is no atom with NC=2, hyde (atom_mark[] --> 1+Nround%2)the starting atoms.  

  */


/*
  printf("#void delete_atom_mark_for_N2_atoms_in_smallest_ring(mol,blk,atom_mark, Nringatms, ringatms,nc)\n");
*/
  /** [1] make '2'-sequential segments (sta_sg2[], end_sg2[], len_sg2[]) **/
  sta_sg2 = (int *)malloc(sizeof(int)*Nringatms);
  end_sg2 = (int *)malloc(sizeof(int)*Nringatms);
  len_sg2 = (int *)malloc(sizeof(int)*Nringatms);

  /*
  for (i=0;i<Nringatms;++i){ printf("%d",nc[ringatms[i]]); }
  printf("\n");
  */

  Nsg = 0;
  state = '-';
  N2 = 0;
  for (i=0;i<Nringatms;++i){
    if (nc[ringatms[i]]==2) N2 += 1;
    ip = (i-1 + Nringatms) % Nringatms;
    i0 = i % Nringatms;
    in = (i+1) % Nringatms;
    if ((nc[ringatms[ip]]!=2) && (nc[ringatms[i0]]==2)){sta_sg2[Nsg] = i0; state = '2';}
    if ((state=='2')&&(nc[ringatms[i0]]==2) && (nc[ringatms[in]]!=2)){end_sg2[Nsg] = i0; Nsg += 1; state = '-';}
  }
  if (state=='2'){
    while ((state=='2')&&(i<2*Nringatms)){
      i0 = i % Nringatms;
      in = (i+1) % Nringatms;
      if ((state=='2')&&(nc[ringatms[i0]]==2) && (nc[ringatms[in]]!=2)){end_sg2[Nsg] = i0; Nsg += 1; state = '-';}
      i += 1;
    }
  }

  /* If all the atoms have NC=2 */
  if (N2==Nringatms){
    sta_sg2[0] = 0; end_sg2[0] = Nringatms-1; len_sg2[0] = Nringatms;
    Nsg = 1;
  }  

  /* Find the maximum length 2-segment */ 
  max_len = 0; max_s = 0;
  for (s=0;s<Nsg;++s){
    len_sg2[s] = end_sg2[s] - sta_sg2[s] + 1;
    if (len_sg2[s]<0) { len_sg2[s] = end_sg2[s] + Nringatms - sta_sg2[s] + 1;}
    if (len_sg2[s]>max_len) {max_s = s; max_len = len_sg2[s];}
    /* printf("sta %d end %d len %d\n",sta_sg2[s],end_sg2[s],len_sg2[s]); */
  }

  /* 
  printf("#Nsg %d\n",Nsg);
  printf("#max_len %d max_s %d\n",max_len,max_s);
  */

  /** [2] delete (atom_mark[]-->0) for the max length segment **/
  if (N2>0){
     as = sta_sg2[max_s];
     ae = end_sg2[max_s];
     if (ae<as) ae = ae + Nringatms;
     for (i=as;i<=ae;++i){
       i0 = i % Nringatms;
       atom_mark[ringatms[i0]] = 0; 
        /* printf("delete %s%d\n",mol->atoms[ringatms[i0]].element,mol->atoms[ringatms[i0]].num_in_file); */
     }
  }

  /** [3] If no N2 atom, hyde (atom_mark[]-->2+Nround%2) for the starting atom **/
  if (N2==0){
 /*
    printf("#(%s) (round %d) No N2 atoms in this ring. --> hyde atom %d ",mol->filename,Nround,mol->atoms[ringatms[0]].num_in_file);
    printf(" atom_mark %d ",atom_mark[ringatms[0]]);
    printf(" --> %d\n", atom_mark[ringatms[0]]);
  */ 
    atom_mark[ringatms[0]] =  2 + Nround%2;
  }

  free(len_sg2);
  free(end_sg2);
  free(sta_sg2);

} /* end of delete_atom_mark_for_N2_atoms_in_smallest_ring() */




void add_new_RING_to_tail(HeadRing, Nringatms, ringatms)
  struct RING *HeadRing;
  int Nringatms;
  int *ringatms; 
{
  struct RING *rn;
  int i;

  rn = HeadRing; 
  while (rn->next != NULL) rn = rn->next;
  rn->next = (struct RING*)malloc(sizeof(struct RING));
  rn->next->prev = rn;
  rn->next->next = NULL;
  rn->next->num = rn->num + 1;
  rn = rn->next;
  rn->Natom = Nringatms;
  rn->num_atoms = (int *)malloc(sizeof(int)*Nringatms);
  for (i=0;i<Nringatms;++i) rn->num_atoms[i] = ringatms[i];

} /* end of add_new_RING_to_tail() */



int mark_neighbors_for_all_the_rings(mol,s,a,Nringatms,ringatms)
  struct MOLECULE *mol;
  int s;             /* number of start   atoms (0..,mol->Natom-1) */
  int a;             /* number of focused atoms (0..,mol->Natom-1) */
  int Nringatms;     /* number of member atoms in the path */
  int *ringatms;     /* [0,1,..,Nringatms-1] array of member atoms in the ring */
{
  int i,b,c;
 /*
 *  printf("#mark_neighbors_for_all_the_rings(s %d a %d Nringatms:%d\n",s,a,Nringatms); 
 *   */
  N_RECURRENT += 1;
  NATOM_IN_PATH[a] = Nringatms;
  ringatms[Nringatms-1] = a;

  for (i=0;i<mol->atoms[a].Nneighbor;++i){
    b = mol->atoms[a].neighbors[i];
    if ((a!=b) && (mol->atoms[b].mark==1) && (NATOM_IN_PATH[b]==0)&&
        (mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H')){
         if (N_RECURRENT<N_RECURRENT_MAX){
           mark_neighbors_for_all_the_rings(mol,s,b,Nringatms+1,ringatms);
         }
    }
    else if ((b==s)&&(Nringatms>2)){
/*
      printf("#Find Circle for blocks !! b=s:%d Nringatms %d\n",mol->atoms[b].num_in_file,Nringatms);
*/
      for (c=0;c<Nringatms;++c){ EQUIV_STA_ATM[ringatms[c]] = s; }
      NATOM_IN_PATH[ringatms[Nringatms-1]] = 0;
      return(1);
    }
  }

 NATOM_IN_PATH[a] = 0;
 return(0);

} /* end of mark_neighbors_for_all_the_rings() */









int mark_neighbors_for_smallest_ring(mol,blk_num,s,a,Nringatms,ringatms,Nringatms_smlt, ringatms_smlt)
  struct MOLECULE *mol;
  int blk_num;         /* block number for focusing blocks. if blk_num == -1, focus all of the atoms. */ 
  int s;               /* number of start   atoms (0..,mol->Natom-1) */
  int a;               /* number of focused atoms (0..,mol->Natom-1) */
  int Nringatms;       /* number of member atoms in the path */
  int *ringatms;       /* [0,1,..,Nringatms-1] array of member atoms in the ring */
  int *Nringatms_smlt; /* Nringatms for the smallest ring through atom blk->num_atoms[s] */
  int *ringatms_smlt;  /* [0,1,..,Nringatms_smlt] array of member atoms in the smallest ring */
{
  int i,b,c;
  
  N_RECURRENT += 1; 
  NATOM_IN_PATH[a]   = Nringatms;
  ringatms[Nringatms-1] = a;

 /*
  printf("#mark_neighbors_for_smallest_ring(s %d a %d Nringatms %d Nringatms_smlt %d)",
   mol->atoms[s].num_in_file, mol->atoms[a].num_in_file,Nringatms, *Nringatms_smlt);
  for (i=0;i<Nringatms;++i){
    printf(" %d",mol->atoms[ringatms[i]].num_in_file);
  }
  printf("\n");
 */


  for (i=0;i<mol->atoms[a].Nneighbor;++i){
    b = mol->atoms[a].neighbors[i];
    if ((a!=b) && (ATOM_MARK[b]==1) && ((blk_num==-1)||(BLOCK_NUM[b]==blk_num)) && (NATOM_IN_PATH[b]==0) && 
          ((*Nringatms_smlt==0)||(Nringatms <= *Nringatms_smlt))  ){
       if (N_RECURRENT<N_RECURRENT_MAX){
         mark_neighbors_for_smallest_ring(mol,blk_num,s,b,Nringatms+1,ringatms, Nringatms_smlt, ringatms_smlt); 
        }
       else{
/*
         printf("#N_RECURRENT IS OVER N_RECCURENT_MAX %d. Nringatms %d Nringatms_smlt %d\n",N_RECURRENT_MAX,Nringatms,*Nringatms_smlt);
         for (c=0;c<Nringatms;++c){ printf(" %d ",mol->atoms[ringatms[c]].num_in_file); }
         printf("\n");
  */
       }
    }
    else if ((b==s)&&(Nringatms>2)){
      /*
      printf("#Find Circle for smallest !! b=s:%d Nringatms %d Nringatms_smlt %d\n",mol->atoms[b].num_in_file,Nringatms,*Nringatms_smlt);
      */
      if ((*Nringatms_smlt==0)||(Nringatms < *Nringatms_smlt)){
        *Nringatms_smlt = Nringatms;
        for (c=0;c<Nringatms;++c){ ringatms_smlt[c] = ringatms[c]; }
      }
      NATOM_IN_PATH[ringatms[Nringatms-1]] = 0;
      return(1);
    }   
  }

 NATOM_IN_PATH[a] = 0;
 return(0);

} /* end of mark_neighbors_for_smallest_ring() */




void restore_hyding_atoms(mol,Nround,atom_mark)
  struct MOLECULE *mol;
  int Nround;
  char *atom_mark; /* [Natom] */
{
  int a;
  /*
   If (Nround is even), atom_mark[a]:: 2 --> 1.
   If (Nround is odds), atom_mark[a]:: 3 --> 1.
 */

  if ((Nround%2)==0){ 
      for (a=0;a<mol->Natom;++a){
        if (atom_mark[a]==2) {
           atom_mark[a] = 1; 
           /* printf("#(round %d)restore %d 2->1\n",Nround,mol->atoms[a].num_in_file); */
        }
      }
    }
      
   if ((Nround%2)==1){ 
     for (a=0;a<mol->Natom;++a){
       if (atom_mark[a]==3) {
          atom_mark[a] = 1; 
          /* printf("#(round %d)restore %d 3->1\n",Nround,mol->atoms[a].num_in_file); */
       }
     }
   }
} /* end of restore_hyding_atoms() */




void show_ringatms_with_NC(mol, Nringatms, ringatms,NC,comment)
  struct MOLECULE *mol;
  int Nringatms;
  int *ringatms;
  int *NC;
  char *comment;
{
  int i,a;
  printf("#%s(%d)",comment,Nringatms);
  for (i=0;i<Nringatms;++i){
    a = ringatms[i];
    printf(" [%d %s %d]",mol->atoms[a].num_in_file,mol->atoms[a].atomname,NC[a]);
  }
  printf("\n");
}





int is_metal_element(ele)
  char *ele;
{
  if ((strcmp(ele,"Li")==0)||(strcmp(ele,"LI")==0)) return(1);
  if ((strcmp(ele,"Mg")==0)||(strcmp(ele,"MG")==0)) return(1);
  if ((strcmp(ele,"Na")==0)||(strcmp(ele,"NA")==0)) return(1);
  if ((strcmp(ele,"Al")==0)||(strcmp(ele,"AL")==0)) return(1);
  if (strcmp(ele,"K")==0) return(1);
  if ((strcmp(ele,"Ca")==0)||(strcmp(ele,"CA")==0)) return(1);
  if ((strcmp(ele,"Ti")==0)||(strcmp(ele,"TI")==0)) return(1);
  if (strcmp(ele,"V")==0) return(1);
  if ((strcmp(ele,"Cr")==0)||(strcmp(ele,"CR")==0)) return(1);
  if ((strcmp(ele,"Mn")==0)||(strcmp(ele,"MN")==0)) return(1);
  if ((strcmp(ele,"Fe")==0)||(strcmp(ele,"FE")==0)) return(1);
  if ((strcmp(ele,"Co")==0)||(strcmp(ele,"CO")==0)) return(1);
  if ((strcmp(ele,"Ni")==0)||(strcmp(ele,"NI")==0)) return(1);
  if ((strcmp(ele,"Cu")==0)||(strcmp(ele,"CU")==0)) return(1);
  if ((strcmp(ele,"Zn")==0)||(strcmp(ele,"ZN")==0)) return(1);
  if ((strcmp(ele,"Ga")==0)||(strcmp(ele,"GA")==0)) return(1);
  if ((strcmp(ele,"Mo")==0)||(strcmp(ele,"MO")==0)) return(1);
  if ((strcmp(ele,"Ru")==0)||(strcmp(ele,"RU")==0)) return(1);
  if ((strcmp(ele,"Rh")==0)||(strcmp(ele,"RH")==0)) return(1);
  if ((strcmp(ele,"Pd")==0)||(strcmp(ele,"PD")==0)) return(1);
  if ((strcmp(ele,"Ag")==0)||(strcmp(ele,"AG")==0)) return(1);
  if ((strcmp(ele,"Cd")==0)||(strcmp(ele,"CD")==0)) return(1);
  if ((strcmp(ele,"Pt")==0)||(strcmp(ele,"PT")==0)) return(1);
  if ((strcmp(ele,"Au")==0)||(strcmp(ele,"AU")==0)) return(1);
  if ((strcmp(ele,"Hg")==0)||(strcmp(ele,"HG")==0)) return(1);
  return(0);
}
