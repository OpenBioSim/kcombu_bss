/*

 <ioSMILES.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 * String of SMILES is repesented by the linear list of LUNIT.
 
   LUNIT has five types : 'A'tom, 'N'umber, 'B'ond, 'P'arenthesis, '.':period

   CC(C)C1CNCC1C

   [C] [-] [(] [-] [C] [)] [-] [C] [-] [1] [-] [C] [-] [N] [=] [C] [-] [C] [-] [1] [-] [C]
    A   B   P   B   A   P   B   A   B   N   B   A   B   A   B   A   B   A   B   N   B   A

 * Production rules 
    (1) A[N]...B[N]   --> A...B    [bond A-B] 
        A=[N]...B=[N] --> A...B    [bond A=B] 
        A#[N]...B#[N] --> A...B    [bond A#B] 
        A=[N]...B[N]  --> A...B    [bond A=B] 
        A#[N]...B[N]  --> A...B    [bond A#B] 
        A[N]...B=[N]  --> A...B    [bond A=B] 
        A[N]...B#[N]  --> A...B    [bond A#B] 
    (2) A - B [end] --> A [end]  [bond A-B] 
        A = B [end] --> A [end]  [bond A=B] 
        A # B [end] --> A [end]  [bond A#B] 
        A - B )     --> A )      [bond A-B] 
        A = B )     --> A )      [bond A=B] 
        A # B )     --> A )      [bond A#B] 
    (3) A (- B)     --> A        [bond A-B] 
        A (= B)     --> A        [bond A=B] 
        A (# B)     --> A        [bond A#B] 
    (4) .A[end]     --> [end]    [do nothing] 


** Stereo-chemical descriptions are ignored.

  example)
    d-alanine:  N[C@H](C)C(=O)O
    l-alanine:  N[C@@H](C)C(=O)O
 
   trans-2-butene: C/C=C/C
     cis-2-butene: C/C=C\C

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
#include "ioSMILES.h"
#include "moltopodis.h"
#include "molring.h"
#include "ringblock.h"


/** FUNCTIONS (GLOBAL) **/
int Translate_Lines_into_SMILES_Molecule();
int make_LUNITs_from_Molecule_by_DFS();
void make_string_from_LUNITs();
int length_string_from_LUNITs();
void Write_SMILES_Molecule();
void randomize_atom_order_of_MOLECULE();
void Set_Unique_Extended_Connectivity();
int  Set_Unique_Rank_of_Atoms();
void Set_Topological_Distance_FW_with_Order_Priority();
int  Number_of_LUNITs();
void Free_LUNITs();
void Set_Connected_Structure();
void Keep_Only_Largest_Connected_Structure();
void Set_Aromatic_Bond_Atoms();

/** FUNCTIONS (LOCAL) **/
static int  make_linear_list_of_LUNIT_from_string();
static void make_atoms_by_linear_list_of_LUNIT();
static void make_bonds_by_production_rules();
static void SMILEbond_str_from_CONMAP();
static char CONMAPbond_char_from_SMILEbond();
static void renumber_LUNITs_and_lunit_num();
static struct LUNIT* add_new_LUNIT_tail();
static struct LUNIT* add_new_LUNIT();
static struct LUNIT* insert_new_LUNIT();
static void remove_LUNIT();
static void show_LUNITs();
static int number_of_bond();
static int priority_atom_pairs();
static int priority_bond_bwn_o_and_a0_a1();
static int recursive_generate_tree_by_DFS();

static int Mark_Neighbors_for_Connected();
static void Bubble_Sort_Index();

static void output_rank_atom_rasmol_script();

static int PRIME[600] = {
   2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
 947, 953, 967, 971, 977, 983, 991, 997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,
1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,
1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,
1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,
1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,1619,1621,1627,1637,1657,
1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,
1823,1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,1951,1973,1979,1987,
1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,2069,2081,2083,2087,2089,2099,2111,2113,2129,
2131,2137,2141,2143,2153,2161,2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,
2293,2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,2399,2411,2417,2423,
2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,2543,2549,2551,2557,2579,2591,2593,2609,2617,
2621,2633,2647,2657,2659,2663,2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,
2749,2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,2879,2887,2897,2903,
2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,3011,3019,3023,3037,3041,3049,3061,3067,3079,
3083,3089,3109,3119,3121,3137,3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,
3259,3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,3389,3391,3407,3413,
3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,3527,3529,3533,3539,3541,3547,3557,3559,3571,
3581,3583,3593,3607,3613,3617,3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,
3733,3739,3761,3767,3769,3779,3793,3797,3803,3821,3823,3833,3847,3851,3853,3863,3877,3881,3889,3907,
3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,4001,4003,4007,4013,4019,4021,4027,4049,4051,4057,
4073,4079,4091,4093,4099,4111,4127,4129,4133,4139,4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,
4241,4243,4253,4259,4261,4271,4273,4283,4289,4297,4327,4337,4339,4349,4357,4363,4373,4391,4397,4409
};


int Translate_Lines_into_SMILES_Molecule(line,mol)
  char *line; 
  struct MOLECULE *mol;
{
  struct LUNIT HeadLunit;
  char *line_smiles;
  int  Nword, Wsta[100],Wend[100];

/*
>> SMILES exmples generated by OpenBabel <<

NC(C(=O)O)Cc1ccc(cc1)O	TYR
c1nc(c2c(n1)cc(c(c2)OCCCN1CCOCC1)OC)Nc1cc(c(cc1)F)Cl	D01977
*/


  line_smiles = (char *)malloc(sizeof(char)*(strlen(line)+1));
  Split_to_Words(line,' ',&Nword,Wsta,Wend,100);
  Get_Part_Of_Line(line_smiles,line,Wsta[0],Wend[0]);
  if (Nword>=2){
    Get_Part_Of_Line(mol->name,line,Wsta[1],Wend[1]);
  }
  else{
    Get_Part_Of_Line(mol->name,line_smiles,0,127);
  }
  printf("#Translate_Lines_into_SMILES_Molecule('%s' '%s')\n",line_smiles,mol->name);

  mol->Natom = make_linear_list_of_LUNIT_from_string(&HeadLunit,line_smiles);
  /* 
  show_LUNITs(&HeadLunit); 
  exit(1);
  */
  Malloc_MOLECULE(mol,mol->Natom,0); 
  make_atoms_by_linear_list_of_LUNIT(mol,&HeadLunit);
  make_bonds_by_production_rules(mol,&HeadLunit);
  /* show_LUNITs(&HeadLunit); */
  if (Number_of_LUNITs(&HeadLunit)>1){
    printf("#ERROR(translate_SMILES): %d LUNITs are left for '%s'.\n",Number_of_LUNITs(&HeadLunit),line);
    /* show_LUNITs(&HeadLunit); */
    exit(1); 
    return(0); 
  }
  Free_LUNITs(&HeadLunit);
  Set_one_char_ele(mol);
  free(line_smiles);
 /*
  free(HeadLunit.next);
  Write_SDF_Molecule("temp.sdf",mol,"FROM SMILES"); 
  */


  return(1); 
} /* end of Translate_Lines_into_SMILES_Molecule(mol) */



 
int make_linear_list_of_LUNIT_from_string(HeadLunit,line0)
  struct LUNIT *HeadLunit;
  char *line0; 
{
  int Natom;
  int i,j,Lline,atom_num;
  char str[32],s,s_next,s_prev,type,type_prev;
  char *line;
/* 
  printf("#make_linear_list_of_LUNIT_from_string(HeadLunit,'%s')\n",line);
*/

  /************************************************/
  /** [1] Making 'clean' and non-Hydrogen line[] **/
  /************************************************/

  Lline = strlen(line0);
  if (line0[Lline-1]=='\n') {line0[Lline-1] = '\0'; Lline = strlen(line0);}
  line = (char *)malloc(sizeof(char)*(Lline+1));

/* 
    N[C@H](C)C(=O)O  --> N[C](C)C(=O)O
    N[C@@H](C)C(=O)O --> N[C](C)C(=O)O
 
   trans-butene: C/C=C/C --> CC=CC
     cis-butene: C/C=C\C --> CC=CC
ZINC65365034:
 C[NH+]1CCC(CC1)Oc2ccc(cc2)NC(=O)[C@@H]3CCCN(C3)S(=O)(=O)C --> C[N]1CCC(CC1)Oc2ccc(cc2)NC(=O)[C]3CCCN(C3)S(=O)(=O)C

2022/06/10:

[H][C@@]12OC3=C(O)C=CC4=C3[C@@]11CCN(C)[C@]([H])(C4)[C@]1(O)CCC2=O 
     -> [][C]12OC3=C(O)C=CC4=C3[C]11CCN(C)[C]([])(C4)[C]1(O)CCC2=O
     -> [C]12OC3=C(O)C=CC4=C3[C]11CCN(C)[C]()(C4)[C]1(O)CCC2=O
     -> [C]12OC3=C(O)C=CC4=C3[C]11CCN(C)[C](C4)[C]1(O)CCC2=O

*/
  /* sprintf(line,"%s",line0); */
   /* (1-1) remove improper characters and 'H' from the line */
  printf("#make_clean '%s'",line0); 
  j = 0;
  for (i=0;i<Lline;++i){
    if ((isalnum(line0[i]))||(line0[i]=='-')||(line0[i]=='=')||(line0[i]=='#') || (line0[i]=='.') /* modified (2017/11/06) */
         ||(line0[i]=='(')||(line0[i]==')')||(line0[i]=='[')||(line0[i]==']')) {
      if (line0[i] != 'H') { 
        line[j] = line0[i]; 
	++j;
      }
    }
  }
  line[j] = '\0'; Lline = j;

  /* (1-2) remove empty '[]' and empty '()' [2022/06/10] */
  for (i=0;i<Lline;++i){ line0[i] = line[i]; }
  j = 0;
  for (i=0;i<Lline;++i){
    if ( ((line0[i]=='(') && (line0[i+1]==')')) || ((line0[i]=='[') && (line0[i+1]==']')) ){
      i += 1;
    }
    else{
      line[j] = line0[i];
      ++j;
    }
  }
  line[j] = '\0'; Lline = j;

  /* (1-3) remove again empty '[]' and empty '()' again [2022/06/10] */
  for (i=0;i<Lline;++i){ line0[i] = line[i]; }
  j = 0;
  for (i=0;i<Lline;++i){
    if ( ((line0[i]=='(') && (line0[i+1]==')')) || ((line0[i]=='[') && (line0[i+1]==']')) ){
      i += 1;
    }
    else{
      line[j] = line0[i];
      ++j;
    }
  }
  line[j] = '\0'; Lline = j;


  printf("-->'%s'\n",line); 
  /*************************************/
  /** [2] Making linear list of LUNIT **/
  /*************************************/

  Natom = 0; 

  HeadLunit->next = HeadLunit->prev = NULL;
  HeadLunit->num = -1;
  HeadLunit->type = ' '; HeadLunit->str[0] = '\0';
 
  Natom = 0;
  i = 0;
  type = type_prev = s_next = s_prev = ' ';
  while (i<Lline){
    s = line[i];
    if (i<(Lline-1)) s_next = line[i+1];
               else  s_next = ' ';
 
     atom_num = -1;

    /** 'A' tom **/
    if (isalpha(s)!=0){
      type = 'A';
      atom_num = Natom;
      Natom += 1;
      if ((isalpha(s_next)!=0)&&(islower(s_next)!=0)&&(s_next!='c')&&(s_next!='n')&&(s_next!='o')&&(s_next!='s')){
        sprintf(str,"%c%c",s,s_next);
        i += 2;
      }
      else{
        sprintf(str,"%c",s);
        i += 1;
      }
    }
    /** 'N'umber **/
    else if (isdigit(s)!=0){
      sprintf(str,"%c",s);
      type = 'N';
      i += 1;
    }
    /** 'P'arenthesis **/
    else if ((line[i]=='(')||(line[i]==')')){
      sprintf(str,"%c",s);
      type = 'P';
      i += 1;
    }
    /** 'B'ond **/
    else if ((line[i]=='-')||(line[i]=='=')||(line[i]=='#')){ 
      sprintf(str,"%c",s);
      type = 'B';
      i += 1;
    }
    /** '[':blacket, such as in 'C[N+](C)(C)C' **/
    else if (line[i] == '['){
      j = i;
      while ((line[j] != ']' )&&(j<Lline)){ 
        j += 1;
      }
      if (line[j]!=']'){
         printf("#ERROR(translate SMILES):%d-th '%c' does not have corresponding ']' in '%s'.\n",i,line[i],line);
         exit(1);
      }
      /** the case for '[]'  (modified in 2022/06/10) */
      if ((line[i]=='[') && (line[i+1] == ']')){
        type = ' ';
        str[0] = '\0';
      }
      else{
        type = 'A';
        atom_num = Natom;
        Natom += 1;
        Get_Part_Of_Line(str,line,i+1,j-1);
      }
      /* printf("i %d line[i] '%c' j %d line[j] '%c' type '%c' str '%s'\n",i,line[i],j,line[j],type,str);  */
      i = j + 1;
    }
    /*** '%': numbers >=10, such as in '1CC2CC3CCC4CC5CCC6CCC7CCC8CC%10CC(C1)C2C%11C3C4C9C5C6C7C8C9C%10%11' */
    else if ((line[i]=='%')&&(Lline>(i+1))&&(isdigit(line[i+1])!=0)){ 
      j = i+1;
      while ((isdigit(line[j])!=0)&&(j<Lline)){ j += 1;}
      if ((j-i-1)<1){
         printf("#ERROR(translate SMILES):%d-th '%c' is improper in '%s'.\n",i,line[i],line);
         exit(1);
      }
      Get_Part_Of_Line(str,line,i+1,j-1);
      type = 'N';
      i = j;
    }
   /**** '.':period **/
    else if (line[i]=='.'){
      sprintf(str,".");
      type = '.';
      i += 1;
    }
    else{
      i += 1;
      type = ' ';
      str[0] = '\0';
    }

    /*** add default hidden single bond ***/
    if ((type_prev=='A')&&(type=='A')){ add_new_LUNIT_tail(HeadLunit,'B',"-",-1);}
    if ((type_prev=='N')&&(type=='A')){ add_new_LUNIT_tail(HeadLunit,'B',"-",-1);}
    if ((type_prev=='P')&&(type=='A')){ add_new_LUNIT_tail(HeadLunit,'B',"-",-1);}

    if (type != ' '){
      /* printf("#str '%s' type '%c'\n",str,type); */
      add_new_LUNIT_tail(HeadLunit,type,str,atom_num);
    }
    type_prev = type;
  }

  /*
  show_LUNITs(HeadLunit); 
  exit(1);  
  */

  free(line);
  return(Natom);

} /* end of make_linear_list_of_LUNIT_from_string() */


void make_atoms_by_linear_list_of_LUNIT(mol,HeadLunit)
  struct MOLECULE *mol;
  struct LUNIT    *HeadLunit;
{
  struct LUNIT *xn;
  int n,L;
  char str[8];
  xn = HeadLunit;
  while (xn->next != NULL){
    xn = xn->next;
    if (xn->type=='A'){
      n = xn->atom_num;
      /* 
   examples of str of 'A' Lunit:
   'C','c','N','n','Cl','Mg','Fe'.
 
   str may be included charge, such as 'N+','C-','Fe+2','Fe++' */
      sprintf(str,"%s",xn->str);
      L = strlen(str);

      if ((L>=3)&&(str[L-2]=='+')&&(str[L-1]=='2')) {str[L-2] = '\0'; mol->atoms[n].char_charge = '2';}
      if ((L>=3)&&(str[L-2]=='+')&&(str[L-1]=='+')) {str[L-2] = '\0'; mol->atoms[n].char_charge = '2';}
      if ((L>=2)&&(str[L-2]!='+')&&(str[L-1]=='+')) {str[L-1] = '\0'; mol->atoms[n].char_charge = '3';}
      if ((L>=2)&&(str[L-2]!='-')&&(str[L-1]=='-')) {str[L-1] = '\0'; mol->atoms[n].char_charge = '5';}
      if ((L>=3)&&(str[L-2]=='-')&&(str[L-1]=='2')) {str[L-2] = '\0'; mol->atoms[n].char_charge = '6';}
      if ((L>=3)&&(str[L-2]=='-')&&(str[L-1]=='-')) {str[L-2] = '\0'; mol->atoms[n].char_charge = '6';}
      sprintf(mol->atoms[n].element,"%s",str);
      Change_String_ToUpper(mol->atoms[n].element);
      sprintf(mol->atoms[n].atomname,"%s",mol->atoms[n].element);
      mol->atoms[n].num_in_file = n + 1;  
    }
  }

} /* end of make_atoms_by_linear_list_of_LUNIT() */


void make_bonds_by_production_rules(mol,HeadLunit)
  struct MOLECULE *mol;
  struct LUNIT    *HeadLunit;
{
  struct LUNIT *a,*b,*c,*d,*e,*f;
  int Nproduct,i,j,nbond; 
  char bondtype,end; 

  /** Production rules **/
  /*  (1) A[N]...B[N]   --> A...B    [bond A-B] */
  /*      A=[N]...B=[N] --> A...B    [bond A=B] */
  /*      A#[N]...B#[N] --> A...B    [bond A#B] */
  /*      A=[N]...B[N]  --> A...B    [bond A=B] */
  /*      A#[N]...B[N]  --> A...B    [bond A#B] */
  /*      A[N]...B=[N]  --> A...B    [bond A=B] */
  /*      A[N]...B#[N]  --> A...B    [bond A#B] */
  /*  (2) A - B [end] --> A [end]  [bond A-B] */
  /*      A = B [end] --> A [end]  [bond A=B] */
  /*      A # B [end] --> A [end]  [bond A#B] */
  /*      A - B )     --> A )      [bond A-B] */
  /*      A = B )     --> A )      [bond A=B] */
  /*      A # B )     --> A )      [bond A#B] */
  /*  (3) A (- B)     --> A        [bond A-B] */
  /*      A (= B)     --> A        [bond A=B] */
  /*      A (# B)     --> A        [bond A#B] */
  /*  (4) .A[end]     --> [end]    [do nothing] */

  /* printf("#make_bonds_by_production_rules()\n"); */
  mol->Nbond = 0;
  mol->chiral_flag = '0';

  do {
    Nproduct = 0; 

  /*** [I] PRODUCTION FOR RULE (1) ***/
  /*  (1) A[N]...B[N]   --> A...B    [bond A-B] */
  /*      A=[N]...B=[N] --> A...B    [bond A=B] */
  /*      A#[N]...B#[N] --> A...B    [bond A#B] */
  /*      A=[N]...B[N]  --> A...B    [bond A=B] */
  /*      A#[N]...B[N]  --> A...B    [bond A#B] */
  /*      A[N]...B=[N]  --> A...B    [bond A=B] */
  /*      A[N]...B#[N]  --> A...B    [bond A#B] */
    a = HeadLunit;
    while (a->next!=NULL){
      a = a->next;
      b = a->next;
      if (b != NULL) c = b->next; else c = NULL;
      if ((((a->type=='A')&&(b!=NULL)&&(b->type=='N'))) ||
          (((a->type=='A')&&(b!=NULL)&&(b->type=='B')&&(c!=NULL)&&(c->type=='N'))) ){
         d = b;
         end = 0;
         while ((d->next != NULL)&&(end==0)){
           d = d->next;
           e = d->next;
           if (e != NULL) f = e->next; else f = NULL;
  /*      a b    c d  */
  /*  (1) A[N]...B[N]   --> A...B    [bond A-B] */
           if ((a->type=='A')&&(b!=NULL)&&(b->type=='N')&&
               (d->type=='A')&&(e!=NULL)&&(e->type=='N')&&(strcmp(b->str,e->str)==0)){
              mol->conmap.map[a->atom_num][d->atom_num] = mol->conmap.map[d->atom_num][a->atom_num] = '1';
              remove_LUNIT(e); remove_LUNIT(b);
              Nproduct += 1;
              mol->Nbond += 1;
              end = 1;
           }
  /*      ab c    d e  */
  /*      A=[N]...B[N]  --> A...B    [bond A=B] */
  /*      A#[N]...B[N]  --> A...B    [bond A#B] */
           if ((a->type=='A')&&(b!=NULL)&&(b->type=='B')&&(c!=NULL)&&(c->type=='N')&&
               (d->type=='A')&&(e!=NULL)&&(e->type=='N')&&(strcmp(c->str,e->str)==0)){
              mol->conmap.map[a->atom_num][d->atom_num] = mol->conmap.map[d->atom_num][a->atom_num] 
                = CONMAPbond_char_from_SMILEbond(b->str[0]);
              remove_LUNIT(e); remove_LUNIT(c); remove_LUNIT(b);
              Nproduct += 1;
              mol->Nbond += 1;
              end = 1;
           }

  /*      a b    de f  */
  /*      A[N]...B=[N]  --> A...B    [bond A=B] */
  /*      A[N]...B#[N]  --> A...B    [bond A#B] */
           if ((a->type=='A')&&(b!=NULL)&&(b->type=='N')&&(d->type=='A')&&(e!=NULL)&&(e->type=='B')
               &&(f!=NULL)&&(f->type=='N')&&(strcmp(b->str,f->str)==0)){
              mol->conmap.map[a->atom_num][d->atom_num] = mol->conmap.map[d->atom_num][a->atom_num] 
                = CONMAPbond_char_from_SMILEbond(e->str[0]);
              remove_LUNIT(f); remove_LUNIT(e); remove_LUNIT(b);
              Nproduct += 1;
              mol->Nbond += 1;
              end = 1;
           }
  /*      ab c    de f   */
  /*      A=[N]...B=[N] --> A...B    [bond A=B] */
  /*      A#[N]...B#[N] --> A...B    [bond A#B] */
           if ((a->type=='A')&&(b!=NULL)&&(b->type=='B')&&(c!=NULL)&&(c->type=='N')&&
               (d->type=='A')&&(e!=NULL)&&(e->type=='B')&&(f!=NULL)&&(f->type=='N')&&
               (strcmp(b->str,e->str)==0)&&(strcmp(c->str,f->str)==0)){
              mol->conmap.map[a->atom_num][d->atom_num] = mol->conmap.map[d->atom_num][a->atom_num] 
                = CONMAPbond_char_from_SMILEbond(b->str[0]);
              remove_LUNIT(f); remove_LUNIT(e); remove_LUNIT(c); remove_LUNIT(b);
              Nproduct += 1;
              mol->Nbond += 1;
              end = 1;
           }

         } /* d */
      }  
    } /* a */

  /*** [II] PRODUCTION FOR RULE (2) and (3) ***/
  /*  (2) A - B [end] --> A [end]  [bond A-B]  */
  /*      A = B [end] --> A [end]  [bond A=B]  */
  /*      A # B [end] --> A [end]  [bond A#B]  */
  /*      A - B )     --> A )      [bond A-B]  */
  /*      A = B )     --> A )      [bond A=B]  */
  /*      A # B )     --> A )      [bond A#B]  */
  /*  (3) A (- B)     --> A        [bond A-B]  */
  /*      A (= B)     --> A        [bond A=B]  */
  /*      A (# B)     --> A        [bond A#B]  */

    a = HeadLunit;
    while (a->next != NULL){
      a = a->next;
      b = c = d = e = NULL;
      if (a->next != NULL){
        b = a->next;
        if (b->next != NULL){
          c = b->next;
          if (c->next != NULL){
             d = c->next;
            if (d->next != NULL) e = d->next;
          }
        }
      }
      /*  (2) A - B [end] --> A [end]  [bond A-B]  */
      /*      A = B [end] --> A [end]  [bond A=B]  */
      /*      A # B [end] --> A [end]  [bond A#B]  */
      /*      A - B )     --> A )      [bond A-B]  */
      /*      A = B )     --> A )      [bond A=B]  */
      /*      A # B )     --> A )      [bond A#B]  */
      if ((a != NULL) && (b!=NULL) && (c!=NULL) && 
          (a->type=='A')&&(b->type=='B')&&(c->type=='A')&&
          ((d==NULL) || ((d!=NULL)&&(d->type=='P')&&(d->str[0]==')')) )){
/*
          printf("#rule (2) A bond B ) or [end] --> A [A bond B] for '%s' '%s' '%s'\n",a->str,b->str,c->str);  
*/
              if (b->str[0]=='=') bondtype = '2';
         else if (b->str[0]=='#') bondtype = '3';
         else bondtype = '1';
         mol->conmap.map[a->atom_num][c->atom_num] = mol->conmap.map[c->atom_num][a->atom_num] = bondtype;
         remove_LUNIT(c);
         remove_LUNIT(b);
        /* show_LUNITs(HeadLunit);   */
         Nproduct += 1;
         mol->Nbond += 1;
      }

      /*  (3) A (- B)     --> A [bond A-B]  */
      /*      A (= B)     --> A [bond A=B]  */
      /*      A (# B)     --> A [bond A#B]  */
      if ((a != NULL) && (b!=NULL) && (c!=NULL) && (d!=NULL) && (e!=NULL) &&
          (a->type=='A')&&(b->type=='P')&&(c->type=='B')&&(d->type=='A')&&(e->type=='P')&&
          (b->str[0]=='(')&&(e->str[0]==')')){
/*
          printf("#rule (3) A (bond B) for '%s:%d' '%s' '%s' '%s:%d' '%s'\n",
    a->str,mol->atoms[a->atom_num].num_in_file,b->str,c->str,d->str,mol->atoms[d->atom_num].num_in_file,e->str);  
*/
              if (c->str[0]=='=') bondtype = '2';
         else if (c->str[0]=='#') bondtype = '3';
         else bondtype = '1';
         mol->conmap.map[a->atom_num][d->atom_num] = mol->conmap.map[d->atom_num][a->atom_num] = bondtype;
         remove_LUNIT(e);
         remove_LUNIT(d);
         remove_LUNIT(c);
         remove_LUNIT(b);
         /* show_LUNITs(HeadLunit);  */
         Nproduct += 1;
         mol->Nbond += 1;
   

      }
       /*  (4) .A[end]        --> [end]       [do nothing] */
      if ((a != NULL) && (b!=NULL) &&  (a->type=='.')&&(b->type=='A') && (c==NULL)){
         remove_LUNIT(a);
         remove_LUNIT(b);
         /* show_LUNITs(HeadLunit);  */
         Nproduct += 1;
      }

    }

  } while (Nproduct > 0);


  /*** [III] make BONDs ***/

 if (mol->Nbond>0){
   mol->bonds = (struct BOND*)malloc(sizeof(struct BOND)*mol->Nbond);
   nbond = 0;
   for (i=0;i<mol->Natom;++i){
     for (j=i+1;j<mol->Natom;++j){
       if (mol->conmap.map[i][j]!='0'){
         mol->bonds[nbond].num   =  nbond;
         mol->bonds[nbond].anum1 = i;
         mol->bonds[nbond].anum2 = j;
         mol->bonds[nbond].bondtype = mol->conmap.map[i][j];
         mol->bonds[nbond].bondstereo = '0';
         nbond += 1;
       }
     }
   } 
 }


} /* end of make_bonds_by_production_rules() */




void Write_SMILES_Molecule(ofname,mol,KeepOnlyLargestConGraph)
  char *ofname;
  struct MOLECULE *mol;
  char KeepOnlyLargestConGraph; /* If 'T', 'XXXX.YYY' --> 'XXXX' */
{
  struct LUNIT HeadLunit;
  FILE *fpo;
  char *string;

  /*** (1) make LUNITs from Molecule ***/
  /*
  make_LUNITs_from_Molecule(&HeadLunit,mol);
  */

  make_LUNITs_from_Molecule_by_DFS(&HeadLunit,mol);
  /* show_LUNITs(&HeadLunit);  */

  string = (char *)malloc(sizeof(char)*(length_string_from_LUNITs(mol,&HeadLunit)));
  make_string_from_LUNITs(mol,string,&HeadLunit,KeepOnlyLargestConGraph); 
  /* printf("#%s\n",string); */
  if (ofname[0]=='-'){ fpo = stdout; }
  else{
     printf("#Write_SMILES_Molecule()-->'%s'\n",ofname);
     fpo = fopen(ofname,"w");
  }

  /* fprintf(fpo,"%s\n",string); */
  fprintf(fpo,"%s %s\n",string,mol->name); 
  if (ofname[0]!='-'){ fclose(fpo); }
  
  Free_LUNITs(&HeadLunit);
  free(string);
  
} /* end of Write_SMILES_Molecule() */














void SMILEbond_str_from_CONMAP(smilestr,conmap_bond)
  char *smilestr;
  char conmap_bond;
{
       if (conmap_bond=='3') smilestr[0] = '#';
  else if (conmap_bond=='2') smilestr[0] = '=';
  else smilestr[0] = '-';
  smilestr[1] = '\0';
}


char CONMAPbond_char_from_SMILEbond(smilechar)
  char smilechar;
{
       if (smilechar=='#') return('3');
  else if (smilechar=='=') return('2');
  else return('1');
}





void renumber_LUNITs_and_lunit_num(mol,HeadLUNIT,lunit_num)
  struct MOLECULE *mol;
  struct LUNIT *HeadLUNIT;
  int *lunit_num;  /* array[mol->Natom]. 
                      If atom[a] is included in LUNITs,lunit_num[a] = lunit_number,otherwise lunit_num[a] = -1; */
{
  struct LUNIT *xn;
  int i, n;

  for (i=0;i<mol->Natom;++i) lunit_num[i] = -1;
  HeadLUNIT->num = -1;
  xn = HeadLUNIT;
  n = 0;
  while (xn->next != NULL){
    xn = xn->next;
    xn->num = n;
    if ((xn->type=='A') && (xn->atom_num>=0)) { 
      lunit_num[xn->atom_num] = n; n += 1; 
    }
  }

} /* end of renumber_LUNITs_and_lunit_num() */




struct LUNIT* add_new_LUNIT_tail(HeadLUNIT,new_type,new_str,atom_num)
  struct LUNIT *HeadLUNIT;
  char  new_type; 
  char  new_str[8]; 
  int   atom_num;
{
  struct LUNIT *new,*last;
  last = HeadLUNIT;
  while (last->next != NULL) last = last->next;
  last->next = (struct LUNIT *)malloc(sizeof(struct LUNIT));
  new = last->next;
  new->prev = last;
  new->next = NULL;
  new->type = new_type;
  sprintf(new->str,"%s",new_str);
  new->atom_num = atom_num;
  return(new);
} /* end of add_new_LUNIT_tail() */





struct LUNIT* add_new_LUNIT(prev_unit,new_type,new_str,atom_num)
  struct LUNIT *prev_unit; /* unit just before the new adding unit */
  char  new_type; 
  char  new_str[8]; 
  int   atom_num;
{
  /* this function is also used for the 'insertion' */
  struct LUNIT *head_next,*new;

  head_next = prev_unit->next;

  prev_unit->next = (struct LUNIT *)malloc(sizeof(struct LUNIT));
  new = prev_unit->next; 
  new->prev = prev_unit;
  new->next = head_next;
  if (head_next != NULL) head_next->prev = new;

  new->type = new_type;
  sprintf(new->str,"%s",new_str);
  new->atom_num = atom_num;
  return(new);
} /* end of add_new_LUNIT() */



struct LUNIT* insert_new_LUNIT(HeadLUNIT,atom_num_inserted,new_type,new_str,atom_num)
  struct LUNIT *HeadLUNIT;
  int   atom_num_inserted;
  char  new_type; 
  char  new_str[8]; 
  int   atom_num;
{
  struct LUNIT *xn,*ins,*new;

  /*
  printf("# insert_new_LUNIT(HeadLUNIT,atom_num_inserted,new_type %c ,new_str %s ,atom_num %d)\n",new_type,new_str,atom_num);
  */

  xn = HeadLUNIT;
  ins = NULL;

  while ((xn->next != NULL)&&(ins==NULL)){
    xn = xn->next;
    if ((xn->type=='A') && (xn->atom_num == atom_num_inserted)) ins = xn;
  }
  if (ins==NULL){
    printf("#ERROR:(insert_new_LUNIT):atom_num_inserted %d is not found.\n",atom_num_inserted);
    exit(1);
  }

  new = add_new_LUNIT(ins,new_type,new_str,atom_num);
  return(new);

} /* end of insert_new_LUNIT() */










void remove_LUNIT(rmv_lunit)
  struct LUNIT *rmv_lunit;
{
/*
  printf("#remove_LUNIT(%c '%s' %d)\n",rmv_lunit->type,rmv_lunit->str,rmv_lunit->atom_num);
*/
  if (rmv_lunit->prev != NULL){
    rmv_lunit->prev->next = rmv_lunit->next;
  }

  if (rmv_lunit->next != NULL){
    rmv_lunit->next->prev = rmv_lunit->prev;
  }
  free(rmv_lunit);

} /* end of remove_LUNIT() */



void show_LUNITs(HeadLUNIT)
  struct LUNIT *HeadLUNIT;
{
  struct LUNIT *xn;
/* 
  show LUNIT for debugging 
 
  LUNIT has five types : 'A'tom, 'N'umber, 'B'ond, 'P'arenthesis, '.':period

  'A'tom is shown by ['element''atom_number'], such as '[C1]' or [N2]' .

*/
  xn = HeadLUNIT;
  while (xn->next!=NULL){
    xn = xn->next;
    /*
    printf("%c '%s' %d\n",xn->type,xn->str,xn->atom_num);
    printf(" %s",xn->str);
    */ 
    if (xn->type=='A'){ printf("[%s%d]",xn->str,xn->atom_num+1);}
    else printf("%s",xn->str);
  } 
  printf("\n");
} /* end of show_LUNITs() */



void make_string_from_LUNITs(mol, string,HeadLUNIT,KeepOnlyLargestConGraph)
  struct MOLECULE *mol;
  char *string;
  struct LUNIT *HeadLUNIT;
  char KeepOnlyLargestConGraph; /* If 'T', 'XXXX.YYY' --> 'XXXX' */
{
  struct LUNIT *xn;
  int n,i,L;
  
  n = 0; 
  string[n] = '\0';
  xn = HeadLUNIT;
  while (xn->next!=NULL){
    xn = xn->next;
    L = strlen(xn->str);
    /* printf("n %d '%s' string '%s'\n",n,xn->str,string); */
    if (xn->type=='A'){ 
      if (L>1) { string[n] = '['; n += 1;}
      if (mol->atoms[xn->atom_num].aromatic == 'A'){
        for (i=0;i<L;++i){ 
          string[n] = tolower(xn->str[i]); 
          n += 1;
        }
      }
      else{
        for (i=0;i<L;++i){ 
          if (i==0) string[n] = xn->str[i];  
               else string[n] = tolower(xn->str[i]); 
          n += 1;
        }
      }
      if (L>1) { string[n] = ']'; n += 1;}
    }
    else if (xn->type=='B'){ 
      if (xn->str[0]!='-'){
        string[n] =  xn->str[0];
        n += 1;
      }
    }
    else if (xn->type=='P'){ 
      string[n] = xn->str[0];
      n += 1;
    }
    else if (xn->type=='N'){
      if (L==1){  
        string[n] = xn->str[0]; n += 1;
      }
      else{
        string[n] = '%'; n += 1;
        for (i=0;i<L;++i){ string[n] = xn->str[i]; n += 1;}
      }
      string[n] = '\0';
    }
    else if (xn->type=='.'){ 
      string[n] = '.';
      n += 1;
    }

  } /* while */

  string[n] = '\0'; 

  if (KeepOnlyLargestConGraph=='T'){
     for (i=0;i<n;++i){
       if (string[i]=='.') string[i] = '\0';
     }
  }

} /* end of make_string_from_LUNITs() */






int length_string_from_LUNITs(mol, HeadLUNIT)
  struct MOLECULE *mol;
  struct LUNIT *HeadLUNIT;
{
  struct LUNIT *xn;
  int n,i,L;
  
  n = 0; 
  xn = HeadLUNIT;
  while (xn->next!=NULL){
    xn = xn->next;
    L = strlen(xn->str);
    if (xn->type=='A'){ 
      if (L>1) { n += 1;}
      if (mol->atoms[xn->atom_num].aromatic == 'A'){
        for (i=0;i<L;++i){ n += 1; }
      }
      else{
        for (i=0;i<L;++i){ n += 1; }
      }
      if (L>1) {  n += 1;}
    }
    else if (xn->type=='B'){ 
      if (xn->str[0]!='-'){ n += 1; }
    }
    else if (xn->type=='P'){ n += 1; }
    else if (xn->type=='.'){ n += 1; }
    else if (xn->type=='N'){
      if (L==1){  n += 1; }
      else{
        n += 1;
        for (i=0;i<L;++i){  n += 1;}
      }
    }
  }
  return(n+1);
} /* end of length_string_from_LUNITs() */






int Number_of_LUNITs(HeadLUNIT)
  struct LUNIT *HeadLUNIT;
{
  struct LUNIT *xn;
  int n, nPs,nPe;
  n = 0;
  nPs = nPe = 0;
  xn = HeadLUNIT;
  while (xn->next!=NULL){ 
    xn = xn->next;
    if ((xn->type=='P') && (xn->str[0]=='('))  nPs += 1;
    if ((xn->type=='P') && (xn->str[0]==')'))  nPe += 1;
    n += 1;
  } 

  /* printf("nPs %d nPe %d\n",nPs,nPe); */
  return(n);
} /* end of Number_of_LUNITs() */


void randomize_atom_order_of_MOLECULE(mol)
  struct MOLECULE *mol;
{
  int *new_frm_old;
  int b,r,i0,i1,j0,j1,k;
  struct MOLECULE molO;
  struct RING *or,*nr;
  printf("#randomize_atom_order_of_MOLECULE()\n");
  Initialize_MOLECULE(&molO);
  Malloc_MOLECULE(&molO,mol->Natom,mol->Nbond);
  Copy_MOLECULE(&molO,mol); 
  /* [1] make new_frm_old and randomly replace new_frm_old[i] and new_frm_old[j] */
  new_frm_old = (int *)malloc(sizeof(int)*(mol->Natom));
  for (i0=0;i0<mol->Natom;++i0){ new_frm_old[i0] = i0; }

  for (r=0;r<mol->Natom;++r){
    i0 = rand()%mol->Natom;
    j0 = rand()%mol->Natom;
    b = new_frm_old[i0];
    new_frm_old[i0] = new_frm_old[j0];
    new_frm_old[j0] = b;
  }
  /*
  for (i0=0;i0<mol->Natom;++i0){
    i1 = new_frm_old[i0];
    printf("i %d %s new_frm_old %d %s\n",
    mol->atoms[i0].num_in_file, mol->atoms[i0].atomname,
    mol->atoms[i1].num_in_file, mol->atoms[i1].atomname);
  }
  */

  /* [2] Copy ATOMs by new_frm_old */
  for (i0=0;i0<mol->Natom;++i0){
    i1 = new_frm_old[i0];
    Copy_ATOM(&(mol->atoms[i1]),&(molO.atoms[i0]));
    for (k=0;k<molO.atoms[i0].Nneighbor;++k){
      mol->atoms[i1].neighbors[k] = new_frm_old[molO.atoms[i0].neighbors[k]];
    }
  } 

  /* [3] Copy BONDs by new_frm_old */
  for (b=0;b<mol->Nbond;++b){
    i0 = molO.bonds[b].anum1; 
    i1 = new_frm_old[i0];
    j0 = molO.bonds[b].anum2; 
    j1 = new_frm_old[j0];

    mol->bonds[b].anum1    = i1;
    mol->bonds[b].anum2    = j1;

  /*
    printf("BOND[%d] old i %d %s -- %d %s  new %d %s -- %d %s\n",b,
      molO.atoms[i0].num_in_file, molO.atoms[i0].atomname,
      molO.atoms[j0].num_in_file, molO.atoms[j0].atomname,
      mol->atoms[i1].num_in_file,mol->atoms[i1].atomname,
      mol->atoms[j1].num_in_file,mol->atoms[j1].atomname);
    */
  }
  

  /* [4] Copy 2DMAP by new_frm_old */
  for (i0=0;i0<mol->Natom;++i0){
    i1 = new_frm_old[i0];
    for (j0=0;j0<mol->Natom;++j0){
      j1 = new_frm_old[j0];
      if (mol->conmap.malloced==1)
        mol->conmap.map[i1][j1]     = molO.conmap.map[i0][j0];
      if (mol->topodismap.malloced==1)
        mol->topodismap.map[i1][j1] = molO.topodismap.map[i0][j0];
      if (mol->pathmap.malloced==1)
        mol->pathmap.map[i1][j1]    = new_frm_old[molO.pathmap.map[i0][j0]];
      if (mol->dismap.malloced==1)
        mol->dismap.map[i1][j1]     = molO.dismap.map[i0][j0];
    }
  }

  /** [5] Copy Ring/Block structure **/
  nr = &(mol->HeadRing); 
  or = &(molO.HeadRing);
  while (or->next != NULL){
    nr = nr->next;
    or = or->next;
    for (k=0;k<or->Natom;++k){ 
        nr->num_atoms[k] = new_frm_old[or->num_atoms[k]]; 
   /*
      printf("ring[%d] new %d %s old %d %s\n",k,
       mol->atoms[nr->num_atoms[k]].num_in_file, mol->atoms[nr->num_atoms[k]].atomname,
       molO.atoms[or->num_atoms[k]].num_in_file, molO.atoms[or->num_atoms[k]].atomname);
    */
    }
  }
  
  nr = &(mol->HeadBlock); 
  or = &(molO.HeadBlock);
  while (or->next != NULL){
    nr = nr->next;
    or = or->next;
    for (k=0;k<or->Natom;++k){ 
      nr->num_atoms[k] = new_frm_old[or->num_atoms[k]]; 
      /*
      printf("block[%d] new %d %s old %d %s\n",k,
       mol->atoms[nr->num_atoms[k]].num_in_file, mol->atoms[nr->num_atoms[k]].atomname,
       molO.atoms[or->num_atoms[k]].num_in_file, molO.atoms[or->num_atoms[k]].atomname);
      */

    }
  }

 /* if SDF, rename num_in_file */
 if (mol->filetype=='S'){
   for (i1=0;i1<mol->Natom;++i1){
     mol->atoms[i1].num_in_file = i1 + 1;
   }
 }


 
  Show_Ring_Structure(mol,&(mol->HeadRing),"RING");
  Show_Ring_Structure(mol,&(mol->HeadRing),"BLOCK");
  printf("HeadRing %d Natom %d\n",molO.HeadRing.num,molO.HeadRing.Natom);
  Free_MOLECULE(&molO);
  free(new_frm_old);
} /* end of randomize_atom_order_of_MOLECULE() */





void Set_Unique_Extended_Connectivity(mol,start_anum,constr_num)
  struct MOLECULE *mol;
  int   start_anum;  /* atom number of the starting atom. if (start_anum <0), not specified */
  int   constr_num;   /* connected structure number */
{
 /*
  based on 
    D.Weininger, A.Weininger, J.L.Weininger. 
    SMILES. 2. Algorithm for Generation of Unique SMILES Notation.
    J.Chem.Inf.Comput.Sci, 1989, 29, 97-101.
 
  * ECuniq[a] for atom a not in the connected structure (constr_num) is set to 0. 

  * Finaly set up mol->atoms[a].rank.
 */
 int a,b,r, maxNround,Nround_same_uniq,maxNround_same_uniq;
 int Natom_uniq1,Natom_uniq0;
 char end,string[32];
 long *ECuniq,*ECuniq1;

 ECuniq  = (long *)malloc(sizeof(long)*mol->Natom); 
 ECuniq1 = (long *)malloc(sizeof(long)*mol->Natom); 
 /*
 printf("#Set_Unique_Extended_Connectivity(start_anum:%d,constr_num:%d OutRankForSMILES:%c)\n",start_anum,constr_num, PAR.OutRankForSMILES);
 */ 
 /** [1] Initialize (uniqEC[] := degree(Nneighbor_heavyatom) )**/
 for (a=0;a<mol->Natom;++a){
   ECuniq[a]  =  0;
   if ((mol->atoms[a].one_char_ele!='H')&&(mol->atoms[a].constr_num == constr_num)){
     for (b=0;b<mol->Natom;++b){
       if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H') && (mol->atoms[b].constr_num==constr_num)){
         ECuniq[a] += 1;
       }
     }
   }
 }

 Natom_uniq0 = Set_Unique_Rank_of_Atoms(mol,ECuniq,start_anum,constr_num);

 if ((PAR.OutRankForSMILES == 'T')&&(start_anum==-1)){
   sprintf(string,"rank_init.ras");
   output_rank_atom_rasmol_script(string,mol,ECuniq);
 }

 /** [2] Repeat multiplying neighbors of PRIME[rank] while (uniqEC1>uniqEC0)**/
  maxNround = mol->Nheavyatom;
  maxNround_same_uniq = 1; 
  r = 1;
  Nround_same_uniq = 0;
  end = 0;
  while (end==0){ 

   for (a=0;a<mol->Natom;++a){
     ECuniq1[a] = 0;
     if ((mol->atoms[a].one_char_ele!='H')&&(mol->atoms[a].constr_num == constr_num)){
      ECuniq1[a] = 1;
       for (b=0;b<mol->Natom;++b){
         if ((mol->conmap.map[a][b] != '0') && (mol->atoms[b].one_char_ele != 'H') &&(mol->atoms[b].constr_num == constr_num)){
           ECuniq1[a] = ECuniq1[a] * PRIME[mol->atoms[b].rank];
          }
        }
      }
    }
    
    Natom_uniq1 = Set_Unique_Rank_of_Atoms(mol,ECuniq1,start_anum,constr_num);

    if (Natom_uniq1 == Natom_uniq0) { Nround_same_uniq += 1;} else{ Nround_same_uniq = 0;}
   /*
    printf("#round %d uniqEC1 %d (%d) / %d Nround_same_uniq %d\n",r,uniqEC1,uniqEC0,mol->Nheavyatom,Nround_same_uniq);
   */
    Natom_uniq0 = Natom_uniq1;

    for (a=0;a<mol->Natom;++a){ ECuniq[a] = ECuniq1[a];}

    if ((r>maxNround)||(Nround_same_uniq>=maxNround_same_uniq)) { end = 1;}

    r += 1;

    if (PAR.OutRankForSMILES == 'T'){
      sprintf(string,"%s_s%d_round%d.ras",mol->name,start_anum,r);
      output_rank_atom_rasmol_script(string,mol,ECuniq);
    }

  } /* while */

 /*
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].one_char_ele != 'H'){
      printf("%s%d EC %d rank %d\n",mol->atoms[a].atomname,mol->atoms[a].num_in_file,(int)ECuniq[a],mol->atoms[a].rank);
    }
  }
  */

  if (PAR.OutRankForSMILES == 'T'){
    sprintf(string,"%s_rank%d.ras",mol->name,start_anum);
    output_rank_atom_rasmol_script(string,mol,ECuniq);
  }


  free(ECuniq1);
  free(ECuniq);
} /* end of Set_Unique_Extended_Connectivity() */


void output_rank_atom_rasmol_script(ofname,mol,ECuniq)
  char *ofname;
  struct MOLECULE *mol;
  long *ECuniq;
{
  FILE *fpo;
  int a;

  fpo = fopen(ofname,"w"); 
  printf("#-->'%s'\n",ofname);
  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].one_char_ele != 'H'){
      fprintf(fpo,"select atomno==%d\nlabel \"%d:%ld\"\n",mol->atoms[a].num_in_file,mol->atoms[a].rank,ECuniq[a]);
    }
  }
  fclose(fpo);
}



int Set_Unique_Rank_of_Atoms(mol,ECuniq,start_anum,constr_num)
  struct MOLECULE *mol;
  long  *ECuniq;  /* ECuniq[mol->Natom] : */ 
  int   start_anum;  /* atom number of the starting atom. if (start_anum <0), not specified */
  int   constr_num;  /* connected structure number */
{
 int a,b,h;
 int *index, *elenum, *bondnum;
 int Nchange,rank,Nheavyatom;
/*
  printf("#Set_Unique_Rank_of_Atoms(start_anum:%d constr_num:%d)\n",start_anum,constr_num);  
*/
 /*
  
   Atom with  smaller rank has a higher priority  
  
  (1) smaller [ECuniq] : unique extended connectivity  
  (2) smaller [elenum] : 
         2 * element  number (C:0,H:1, HE:2,LI:3,BE:4,B:5,N:6,O:7) + 1 * (aromatics=='A')
      Non-aromatic has more priority than Aromatic !!
  (3) smaller [bondnum]: maximum bond valence number. (-:1,=:2,#:3) 
 
  * rank[a] for atom a not in the connected structure (constr_num) is set to -1. 
 
 */
 index   = (int *)malloc(sizeof(int)*mol->Natom);
 elenum  = (int *)malloc(sizeof(int)*mol->Natom);
 bondnum = (int *)malloc(sizeof(int)*mol->Natom);
 /** [1] Setup index[], elenum[] and bondnum[] **/
 /*
 for (a=0;a<mol->Natom;++a){
   printf("%s%d EC %d\n",mol->atoms[a].element,mol->atoms[a].num_in_file,ECuniq[a]);
 }
 */

 Nheavyatom = 0;
 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].rank = -1;
   if (mol->atoms[a].constr_num == constr_num){ 
     elenum[a]  = 2*Number_of_Element(mol->atoms[a].element);
     if (mol->atoms[a].aromatic == 'A') elenum[a] += 1;  /* Non-aromatic has more priority than Aromatic. */
     if ((start_anum>=0)&&(a==start_anum)) ECuniq[a] = 0; /* ECuniq[start_anum] should be zero !!! */
     bondnum[a] = 1;
     if (mol->atoms[a].one_char_ele != 'H'){ 
       index[Nheavyatom]   = a;
       for (b=0;b<mol->Natom;++b){
         if (mol->conmap.map[a][b]=='2') bondnum[a] = 2;
         if (mol->conmap.map[a][b]=='3') bondnum[a] = 3;
       }
       Nheavyatom += 1;
     }
  }
 }
 /*
  printf("#Nheavyatom %d\n",Nheavyatom);
 */

 /** [2] Bubble Sort **/
 do{
   Nchange = 0;
   for (h=1;h<Nheavyatom;++h){
       if ( (ECuniq[index[h-1]] >  ECuniq[index[h]])|| 
          ( (ECuniq[index[h-1]] == ECuniq[index[h]])&&(elenum[index[h-1]]>elenum[index[h]])) ||
          ( (ECuniq[index[h-1]] == ECuniq[index[h]]) &&
            (elenum[index[h-1]]==elenum[index[h]]) && (bondnum[index[h-1]]>bondnum[index[h]]))
        ){ 
       b = index[h];
       index[h] = index[h-1];
       index[h-1] = b; 
       Nchange += 1;
     }
   }
 } while (Nchange > 0);

 rank = 0;
 /** [3] Set up rank **/
 mol->atoms[index[0]].rank = rank;
 for (h=1;h<Nheavyatom;++h){
     if ( (ECuniq[index[h-1]] <  ECuniq[index[h]])|| 
          ((ECuniq[index[h-1]] == ECuniq[index[h]])&&(elenum[index[h-1]]<elenum[index[h]])) ||
          ((ECuniq[index[h-1]] == ECuniq[index[h]])
          &&(elenum[index[h-1]]==elenum[index[h]]) && (bondnum[index[h-1]]<bondnum[index[h]]))
      ){  rank += 1;}
     mol->atoms[index[h]].rank = rank;
 }

 /*
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].one_char_ele != 'H'){
     printf("%-2s%3d ECuniq %10d ele %d bond %d rank %3d\n",
     mol->atoms[a].atomname, mol->atoms[a].num_in_file, (int)ECuniq[a],elenum[a],bondnum[a],mol->atoms[a].rank);
   }
 }
 printf("#Nrank %4d\n",rank+1);
 */

 free(bondnum);
 free(elenum);
 free(index);
 return(rank+1);

} /* end of Set_Unique_Rank_of_Atoms() */




void Set_Topological_Distance_FW_with_Order_Priority(mol, includeH)
  struct MOLECULE *mol;
  char   includeH; /* 'H': including hydrogen atoms */
{
 /*
  This function is based on "Floyd-Warshall" algorithm. 

  mol->topodismap :  topological distance map (number of bonds on the shortest path) 
  mol->pathmap    :  path for topological distance p[i][j]:the last node just before j on the path i->j


  */

 int a,b,c,newdis;

 printf("#Set_Topological_Distance_FW_with_Order_Priority(mol)\n"); 
 /** (1) Initialize topodismap[][] **/

 for (a=0;a<mol->Natom;++a){
   for (b=0;b<mol->Natom;++b){
     if (mol->conmap.map[a][b]!='0')  mol->topodismap.map[a][b] = 1;
                                else  mol->topodismap.map[a][b] = 255;
     mol->pathmap.map[a][b] = a;
   }
   mol->topodismap.map[a][a] = 0;
 }


 /** (2) Calculate length of shortest path **/
 /*     a    */
 /*   /  \   */
 /*  b....c  */
 /*  if (Dbc > Dba + Dac ): Dbc := Dba + Dac  */
 /*                         Pbc := Pac        */

 for (a=0;a<mol->Natom;++a){
   if ((mol->atoms[a].one_char_ele != 'H') || (includeH=='H')){
     for (b=0;b<mol->Natom;++b){
       if ((b!=a)&&((mol->atoms[b].one_char_ele != 'H')||(includeH=='H'))){
         for (c=0;c<mol->Natom;++c){
           if ((c!=a)&&((mol->atoms[c].one_char_ele != 'H')||(includeH=='H'))){
              newdis = mol->topodismap.map[b][a] +  mol->topodismap.map[a][c];
              if (mol->topodismap.map[b][c]>newdis){
                mol->topodismap.map[b][c] =newdis;
                mol->pathmap.map[b][c]    = mol->pathmap.map[a][c];
              }
              else if ((mol->topodismap.map[b][c]==newdis)&&(mol->pathmap.map[a][c] != mol->pathmap.map[b][c])&&
                  (priority_bond_bwn_o_and_a0_a1(mol,c,mol->pathmap.map[a][c],mol->pathmap.map[b][c])==0)){
/*
                printf("a %s%d b %s%d c %s%d disbc %d newdis %d bc_path_prec %s%d bc_path_prec_new_cand %s%d\n",
                 mol->atoms[a].atomname,mol->atoms[a].num_in_file,
                 mol->atoms[b].atomname,mol->atoms[b].num_in_file,
                 mol->atoms[c].atomname,mol->atoms[c].num_in_file,mol->topodismap.map[b][c],newdis,
                 mol->atoms[mol->pathmap.map[b][c]].atomname,mol->atoms[mol->pathmap.map[b][c]].num_in_file,
                 mol->atoms[mol->pathmap.map[a][c]].atomname,mol->atoms[mol->pathmap.map[a][c]].num_in_file
                );
*/
                mol->topodismap.map[b][c] =newdis;
                mol->pathmap.map[b][c]    = mol->pathmap.map[a][c];
              }
           }
         }
      }
    }
  }
 }

} /* end of Set_Topological_Distance_FW_with_Order_Priority() */



int number_of_bond(bondchar)
  char bondchar;
{
       if (bondchar=='0') return(0);
  else if (bondchar=='r') return(2);
  else if (bondchar=='2') return(3);
  else if (bondchar=='3') return(4);
  else return(1); 
}




int priority_atom_pairs(mol,lunit_num,An,Bn,Ap,Bp,symmetric)
  struct MOLECULE *mol;
  int     *lunit_num;  /* array[mol->Natom]. 
                       If atom[a] is included in LUNITs,lunit_num[a] = lunit_number,otherwise lunit_num[a] = -1; */
  int    An,Bn; /* atom numbers for candidate(-1).  */
  int    Ap,Bp; /* atom numbers for candidate(+1).  */
  char  symmetric;  /* If 'S', assign higher priority atom as 'A', another as 'B'. otherwise: atom 'A' has higher priority than 'B' */
 /*
   return(-1): if (an,bn) has a more priority than (ap,bp). 
   return(+1): if (ap,bp) has a more priority than (an,bn). 
  */
{
  int an,bn,ap,bp;

  if (symmetric=='S'){
    if (mol->atoms[An].rank < mol->atoms[Bn].rank) {an = An; bn = Bn;} else {an = Bn; bn = An;}
    if (mol->atoms[Ap].rank < mol->atoms[Bp].rank) {ap = Ap; bp = Bp;} else {ap = Bp; bp = Ap;}
  }
  else {an = An; bn = Bn; ap = Ap; bp = Bp; }


         if (mol->atoms[an].rank <   mol->atoms[ap].rank){ return(-1);}
   else  if (mol->atoms[an].rank >   mol->atoms[ap].rank){ return(+1);}
   else  if (mol->atoms[an].rank ==  mol->atoms[ap].rank){
          if (mol->atoms[bn].rank <   mol->atoms[bp].rank){ return(-1);}
     else if (mol->atoms[bn].rank >   mol->atoms[bp].rank){ return(+1);}
     else if (mol->atoms[bn].rank == mol->atoms[bp].rank){
          if ((lunit_num[an]>=0)&&(lunit_num[ap]>=0)){
            if (lunit_num[an] <  lunit_num[ap]){ return(-1);}
            if (lunit_num[an] >  lunit_num[ap]){ return(+1);}
          }
          if ((lunit_num[bn]>=0)&&(lunit_num[bp]>=0)){
            if (lunit_num[bn] <  lunit_num[bp]){ return(-1);}
            if (lunit_num[bn] >  lunit_num[bp]){ return(+1);}
          }
     } 
   }
  return(0); /* If two bonds are really equivalent, return(0) */
} /* end of priority_atom_pairs() */




int priority_bond_bwn_o_and_a0_a1(mol,o,a0,a1)
  struct MOLECULE *mol;
  int o;  /* atom number of focused origin atom */ 
  int a0; /* atom number of connected atom candidate0 */ 
  int a1; /* atom number of connected atom candidate1 */ 
 /*
   return(0): if a0 has a more priority than a1. 
   return(1): if a1 has a more priority than a0. 
  */
{
  int b0,b1;

       if (mol->atoms[a0].rank < mol->atoms[a1].rank){ return(0); }
  else if (mol->atoms[a1].rank < mol->atoms[a0].rank){ return(1);}
  else if (mol->atoms[a1].rank == mol->atoms[a0].rank){
/*
      printf("#bond_decision o %s%d a0 %s%d a1 %s%d b0 %d b1 %d\n",
     mol->atoms[o].atomname, mol->atoms[o].num_in_file,
     mol->atoms[a0].atomname, mol->atoms[a0].num_in_file,
     mol->atoms[a1].atomname, mol->atoms[a1].num_in_file, b0,b1);
 */
      b0 = number_of_bond(mol->conmap.map[o][a0]);
      b1 = number_of_bond(mol->conmap.map[o][a1]);
           if (b0<b1) {return(0);}
      else if (b1<b0) {return(1);}
  }
  return(0); /* If two bonds are really equivalent, return(0), aribitrarily!! */
} /* end of priority_bond_bwn_o_and_a0_a1() */



void Set_Aromatic_Bond_Atoms(mol)
  struct MOLECULE *mol;
{
  struct RING *rn;
  int i,j,ai,aj;
  int Nbond,Nbond1,Nbond2,NbondR,Naddaroma;
  int C,N,O,S;  
  int Npi,sumNpi; /* Number of pi electron */

 /*  printf("#Set_Aromatic_Bond_Atoms()\n"); */
  /*  
  maxNatom_in_ring = 8;
  Set_Ring_Structure_with_maxNatom_in_ring(mol,&HeadRing,maxNatom_in_ring);
  Show_Ring_Structure(mol,&HeadRing,"SSSR-Ring");
  */

  for (i=0;i<mol->Natom;++i){
    mol->atoms[i].aromatic = ' ';
  }

  /* rn = &HeadRing;  */
  rn = &(mol->HeadRing);  
  while (rn->next != NULL){
    rn = rn->next;
    rn->aromatic = ' ';
  }


  /*** [1] Decide aromaticity for each RING ***/
  do{
    Naddaroma = 0;
    rn = &(mol->HeadRing); 
    /* rn = &HeadRing;  */
    
    while (rn->next != NULL){
      rn = rn->next;
      if (rn->aromatic == ' '){
        sumNpi = 0;
        C = N = O = S = 0;
        for (i=0;i<rn->Natom;++i){
          
          Npi = Nbond = Nbond1 = Nbond2 = NbondR = 0;
          ai = rn->num_atoms[i];

          if (mol->atoms[ai].one_char_ele == 'C') C += 1;
          if (mol->atoms[ai].one_char_ele == 'S') S += 1;
          if (mol->atoms[ai].one_char_ele == 'N') N += 1;
          if (mol->atoms[ai].one_char_ele == 'O') O += 1;

          for (aj=0;aj<mol->Natom;++aj){
            if ((ai!=aj)&&(mol->atoms[aj].one_char_ele != 'H')){
              if (mol->conmap.map[ai][aj]!='0') Nbond  += 1;
                   if (mol->conmap.map[ai][aj]=='1') Nbond1 += 1;
              else if (mol->conmap.map[ai][aj]=='2') Nbond2 += 1;
              else if (mol->conmap.map[ai][aj]=='r') NbondR += 1;
            }
          } /* aj */

          if (mol->atoms[ai].one_char_ele=='C'){
             if ((Nbond==2)||(Nbond==3)){
                    if ((Nbond1>=1)&&(Nbond2==1)) Npi = 1;
               else if  (NbondR>=1)  Npi = 1;
             }
          } 
          else if (mol->atoms[ai].one_char_ele=='N'){
                  if ((Nbond==2)&&(Nbond2==1)) Npi = 1;
             else if ((Nbond==2)&&(NbondR==2)) Npi = 1;
             else if ((Nbond==2)&&(Nbond1==2)) Npi = 2;
             else if ((Nbond==2)&&(Nbond1==2)) Npi = 2;
             else if ((Nbond==3)&&(NbondR==2)&&(Nbond1==1)) Npi = 1;
             else if ((Nbond==3)&&(Nbond2==1)&&(Nbond1==2)) Npi = 2;
             else if ((Nbond==3)&&(Nbond1==3)) Npi = 2;
          } 
          else if ((mol->atoms[ai].one_char_ele=='O')||(mol->atoms[ai].one_char_ele=='S')){
                  if ((Nbond==2) && (Nbond1==2)) Npi = 2;
             else if ((Nbond==2) && (Nbond==2)) Npi = 2;
          } 

         /*
         printf("%s%d Nbond %d [1]:%d [2]:%d [r]:%d Npi %d\n",
           mol->atoms[ai].atomname,mol->atoms[ai].num_in_file,Nbond,Nbond1,Nbond2,NbondR,Npi);
          */
         sumNpi += Npi;

        } /* i */

        /*
        printf("#RING[%d] Natom %d sumNpi %d\n",rn->num,rn->Natom,sumNpi);     
        */
        if ((rn->Natom==5)||(rn->Natom==6)){ 
          if ( ((C+S+N+O)==rn->Natom) &&
               ((sumNpi==6)||(sumNpi==10)||(sumNpi==14)||(sumNpi==18)) ) rn->aromatic = 'T';

          /** for pyrrole **/
          if ( (rn->Natom==5) && (C==4) && (N==1) && ((sumNpi==5)||(sumNpi==6)) )rn->aromatic = 'T';
        
          /** for pyrrole **/
          if ( (rn->Natom==6) && (C==5) && (N==1) && ((sumNpi==6)||(sumNpi==7)) )rn->aromatic = 'T';
        }

        if (rn->aromatic == 'T'){
          Naddaroma += 1;
          for (i=0;i<rn->Natom;++i){
            ai = rn->num_atoms[i];
            mol->atoms[ai].aromatic = 'A';
              for (j=0;j<rn->Natom;++j){
                aj = rn->num_atoms[j];
                if (mol->conmap.map[ai][aj] != '0') mol->conmap.map[ai][aj] = 'r';
              }
          }
        }
     }
   } /* rn */

 }while (Naddaroma > 0);


/*
  for (i=0;i<mol->Natom;++i){
    printf("#%s%d aromatic '%c'\n",mol->atoms[i].atomname,mol->atoms[i].num_in_file,mol->atoms[i].aromatic);
  }
 Free_RING(&HeadRing);
 */
} /* end of Set_Aromatic_Bond_Atoms() */


void Free_LUNITs(HeadUnit)
 struct LUNIT *HeadUnit;
{
 struct LUNIT *un;
 /* printf("#Free_LUNITs(HeadUnit)\n"); fflush(stdout); */
 un = HeadUnit;
 if (un->next != NULL){
   while (un->next != NULL) un = un->next;
   while (un->prev!=NULL){
     un = un->prev;
     if (un->next != NULL){
       free(un->next);
     }
   }
 }
 HeadUnit->next = NULL;
 HeadUnit->prev = NULL;
} /* end of Free_LUNITs() */




int make_LUNITs_from_Molecule_by_DFS(HeadLunit,mol)
  struct LUNIT *HeadLunit;
  struct MOLECULE *mol;
{
  int i,j,sa,min_i,min_j,c; 
  int *lunit_num;  /* array[mol->Natom]. 
                      If atom[a] is included in LUNITs,lunit_num[a] = lunit_number,otherwise lunit_num[a] = -1; */
  struct CHAR2DMAP bond_mark;
  int Nadd,Ncycle,Nunmark_bond,Nbond;
  char buff1[8],buff2[8];
  struct LUNIT *tail;
 
  /* printf("#make_LUNITs_from_Molecule_by_DFS(Nheavyatom %d)\n",mol->Nheavyatom); */
  
  Set_Aromatic_Bond_Atoms(mol);
  Set_Connected_Structure(mol);
 /*
  printf("#Nconstr %d\n",mol->Nconstr);
  for (i=0;i<mol->Natom;++i){
    printf("%s%d constr %d\n",mol->atoms[i].element,mol->atoms[i].num_in_file,mol->atoms[i].constr_num);
  } 
 */

  HeadLunit->next = NULL; HeadLunit->prev = NULL; HeadLunit->num  = -1;
  lunit_num = (int *)malloc(sizeof(int)*(mol->Natom));
  Ncycle = 0;

  /** Repeat making LUNITs for each connected structure **/

  for (c=0;c<mol->Nconstr;++c){
    /* printf(">>>>cluster %d\n",c); */
    if (c>0) add_new_LUNIT_tail(HeadLunit,'.',".",-1);
    /* show_LUNITs(HeadLunit); */
    /*** [0] Initial preparation ***/  

    Set_Unique_Extended_Connectivity(mol,-1,c);

    Malloc_UNSIGNED_CHAR2DMAP(&bond_mark,mol->Natom);

    for (i=0;i<mol->Natom;++i){
      for (j=0;j<mol->Natom;++j){ bond_mark.map[i][j] = 0; }
    }

    if (mol->Nheavyatom==1){
      for (i=0;i<mol->Natom;++i){
        if (mol->atoms[i].one_char_ele != 'H'){
          add_new_LUNIT_tail(HeadLunit,'A',mol->atoms[i].element,i);
        }
      }
     return(Number_of_LUNITs(HeadLunit));
    }
    renumber_LUNITs_and_lunit_num(mol,HeadLunit,lunit_num);

    /** [1] find the lowest (0)-ranked atom */
    sa = -1;
    for (i=0;i<mol->Natom;++i){
      mol->atoms[i].mark = 0;
      if (mol->atoms[i].rank == 0){ sa = i;}
    }
    if (sa==-1) {printf("#ERROR:no rank==0 atom in the molecule.\n"); exit(1);}
/*  
    printf("#starting_atom:%s%d rank %d\n",mol->atoms[sa].atomname,mol->atoms[sa].num_in_file,mol->atoms[sa].rank);  
 */  
    Set_Unique_Extended_Connectivity(mol,sa,c); 
/* 
    printf("#starting_atom:%s%d rank %d\n",mol->atoms[sa].atomname,mol->atoms[sa].num_in_file,mol->atoms[sa].rank);  
 */
    /** [2] Making LUNITs tree-structure recursively */
    tail = HeadLunit;
    while (tail->next != NULL) {tail = tail->next;}
    recursive_generate_tree_by_DFS(sa,tail,mol,&bond_mark,1,c);

    Number_of_LUNITs(HeadLunit);
    /** [3] Count Nunmark_bond */

    renumber_LUNITs_and_lunit_num(mol,HeadLunit,lunit_num);
    Nbond = 0;
    Nunmark_bond = 0;
    for (i=0;i<mol->Natom;++i){
      if ((mol->atoms[i].one_char_ele != 'H')&&(mol->atoms[i].constr_num==c)){
        for (j=i+1;j<mol->Natom;++j)
          if ((mol->atoms[j].one_char_ele != 'H')&&(mol->atoms[j].constr_num==c)){
            if ((mol->conmap.map[i][j] != '0') && (bond_mark.map[i][j] ==0))  {Nunmark_bond += 1;}
            if (mol->conmap.map[i][j] != '0') {Nbond += 1;}
          }
      }
    }

    /* printf("#Nbond %d Nunmark_bond %d\n",Nbond,Nunmark_bond); */
   /*** [3] Add numbers (1,2,3,..) for  cycle-edges ***/
    /* printf("#Nunmark_bond %d\n",Nunmark_bond); */
    if (Nunmark_bond>0){
      do{
        Nadd = 0;
        min_i = -1; min_j = -1;
        for (i=0;i<mol->Natom;++i){
          if ((mol->atoms[i].one_char_ele != 'H') && (mol->atoms[i].constr_num==c) && (lunit_num[i]>=0)){
            for (j=i+1;j<mol->Natom;++j){
              /*** The connecting atoms min_i and min_j are symmetric.  ***/
              if ((mol->atoms[j].one_char_ele != 'H') && (mol->atoms[j].constr_num==c) && (lunit_num[j]>=0) 
                   && (mol->conmap.map[i][j]!='0') && (bond_mark.map[i][j]==0)){ 
                   if (((min_i<0)&&(min_j<0))|| 
                      ((min_i>=0)&&(min_j>=0)&&(priority_atom_pairs(mol,lunit_num,min_i,min_j,i,j,'S')==1))){
                     min_i = i; min_j = j; 
                   }
              }
           }
          }
        }

       if ((min_i>=0)&&(min_j>=0)){
          Ncycle += 1;
          sprintf(buff1,"%d",Ncycle);
          SMILEbond_str_from_CONMAP(buff2,mol->conmap.map[min_i][min_j]);
          insert_new_LUNIT(HeadLunit,min_i,'N',buff1,-1);
          insert_new_LUNIT(HeadLunit,min_i,'B',buff2,-1);
          insert_new_LUNIT(HeadLunit,min_j,'N',buff1,-1);
          insert_new_LUNIT(HeadLunit,min_j,'B',buff2,-1);
          /* show_LUNITs(HeadLunit); */
          bond_mark.map[min_i][min_j] = bond_mark.map[min_j][min_i] = 1;
          Nunmark_bond -= 1;
          Nadd += 1;
       }
      }while (Nadd>0); 
    }
 
  } /* c (loop for connected structure) */


  Number_of_LUNITs(HeadLunit);
  free(lunit_num);
  Free_UNSIGNED_CHAR2DMAP(&bond_mark);
  return(1);
} /* end of make_LUNITs_from_Molecule_by_DFS() */








int recursive_generate_tree_by_DFS(focus_anum,preLunit,mol,bond_mark,init_path,constr_num)
  int    focus_anum;     /* focusing atom number (0,1,2,...) */
  struct MOLECULE *mol;
  struct LUNIT *preLunit;
  struct CHAR2DMAP *bond_mark;
  char   init_path;   /* if (init_path==1), this is in the initial path starting from the starting atom */
  int    constr_num;   /* connected structure number */
{
 struct LUNIT *focusLunit,*newLunit;
 int i,index[10],Nchange,buff,Nneighbor_heavy,Nneighbor_add;
 char init_path_next,bstr[8];

 /*
 printf("#DFS(focus:%s%d,rank %d init_path:%d)\n",
   mol->atoms[focus_anum].atomname, mol->atoms[focus_anum].num_in_file,mol->atoms[focus_anum].rank, init_path);
 */

 /** [1] add focus_atom into LUNITs  **/
 focusLunit =  add_new_LUNIT(preLunit,'A',mol->atoms[focus_anum].element,focus_anum);
 mol->atoms[focus_anum].mark = 1;
  
 /** [2] set index[] of neighbors **/
 Nneighbor_heavy = 0;
 for (i=0;i<mol->Natom;++i){
   if ((focus_anum!=i)&&(mol->atoms[i].one_char_ele !='H')&&(mol->atoms[i].constr_num==constr_num)
      &&(mol->atoms[i].mark==0)&&(mol->conmap.map[focus_anum][i]!='0')){ 
     index[Nneighbor_heavy] = i;
     Nneighbor_heavy += 1;  
   }
 }

 /** [3] bubble sort for index[] by rank **/
 do{
   Nchange = 0;
   for (i=1;i<Nneighbor_heavy;++i){
     if (mol->atoms[index[i-1]].rank > mol->atoms[index[i]].rank){
       buff = index[i]; 
       index[i] = index[i-1]; 
       index[i-1] = buff; 
       Nchange += 1;
     }
   }
 } while (Nchange>0); 
 
 /** [4] repeat executing "recursive_generate_tree_by_DFS()" for each neighbor  **/
 init_path_next = init_path;
 Nneighbor_add = 0;
 for (i=0;i<Nneighbor_heavy;++i){
   /*
   printf("#(%s%d --> %s%d)\n",
     mol->atoms[focus_anum].atomname, mol->atoms[focus_anum].num_in_file,
     mol->atoms[index[i]].atomname, mol->atoms[index[i]].num_in_file);
   */
   if ((mol->atoms[index[i]].mark==0)&&(mol->atoms[index[i]].constr_num==constr_num)){
      bond_mark->map[focus_anum][index[i]] = bond_mark->map[index[i]][focus_anum] = 1; 
      SMILEbond_str_from_CONMAP(bstr,mol->conmap.map[focus_anum][index[i]]);
     if (i>0){ 
       newLunit =  add_new_LUNIT(focusLunit,'P',"(",-1); 
       init_path_next = 0;
       newLunit = add_new_LUNIT(newLunit,'B',bstr,-1);
       recursive_generate_tree_by_DFS(index[i],newLunit,mol,bond_mark,init_path_next,constr_num);
     }
     else{
       newLunit = add_new_LUNIT(focusLunit,'B',bstr,-1);
       recursive_generate_tree_by_DFS(index[i],newLunit,mol,bond_mark,init_path_next,constr_num);
     }
     Nneighbor_add += 1;
   }
 }

 if ((Nneighbor_add==0)&&(init_path==0)){
     add_new_LUNIT(focusLunit,'P',")",-1);
 }

 return(1);

} /* end of recursive_generate_tree_by_DFS() */





void Set_Connected_Structure(mol)
 struct MOLECULE *mol;
{
 int a,i,noheavy_cluster_num;
 int *Natom_constr;  /* [0..mol->Nconstr-1] (malloc_later): number of heavyatoms for connected str. */
 int *new_from_orig; /* [0..mol->Nconstr-1] (malloc_later): new constr number from original constr number */
 int *Selem_constr;  /* [0..mol->Nconstr-1] (malloc_later): Sum of atomic element number (H:1, He:2, .., C:6, N:7, O:8,)  */


 /** [1] initialize (constr_num, mark) **/
 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].constr_num = -1;
 }

 mol->Nconstr = 0;

 /** [2] Mark_Neighbors_for_Connected for each atom **/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].constr_num == -1){
     Mark_Neighbors_for_Connected(mol,a,mol->Nconstr);
     mol->Nconstr += 1;
    }
 }

 /** output **/
 /*
  for (a=0;a<mol->Natom;++a){
    printf("#constr %d (%d %s %s)\n", mol->atoms[a].constr_num,mol->atoms[a].num_in_file,mol->atoms[a].element,mol->atoms[a].atomname);
   }
 */

 /** [3] Malloc Natom_constr[], Selem_constr[] and new_from_orig[]   **/
 Natom_constr   = (int *)malloc(sizeof(int)*mol->Nconstr);
 Selem_constr   = (int *)malloc(sizeof(int)*mol->Nconstr);
 new_from_orig  = (int *)malloc(sizeof(int)*mol->Nconstr);

 for (i=0;i<mol->Nconstr;++i){ 
   Natom_constr[i] = Selem_constr[i] = 0;
 }

 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].constr_num>=0){
     /* printf("%d %s constr_num %d\n",mol->atoms[a].num_in_file,mol->atoms[a].element,mol->atoms[a].constr_num); */
     if (mol->atoms[a].one_char_ele != 'H'){
       Natom_constr[mol->atoms[a].constr_num] += 1;
       Selem_constr[mol->atoms[a].constr_num] += Number_of_Element(mol->atoms[a].element);
     }
   }
 }

 /** [3] Sort constr_num by Natom_constr[] **/
 Bubble_Sort_Index(Natom_constr,Selem_constr,new_from_orig,mol->Nconstr);

 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].constr_num = new_from_orig[mol->atoms[a].constr_num];
 }

 printf("#N_CONNECTED_STRUCTURE %d\n",mol->Nconstr);
 for (i=0;i<mol->Nconstr;++i){
   printf("#CONNECTED_STR %d N_ATOM %d S_ELEM %d\n",i,Natom_constr[i],Selem_constr[i]);
 }

 /** [4] Omit non-heavy-atom (only H) cluster, by decreasing mol->Nconstr **/ 
 noheavy_cluster_num = -1;
 for (i=mol->Nconstr-1;i>=0;--i){
   if (Natom_constr[i]==0){ 
      noheavy_cluster_num = i;
   }
 } 

 if (noheavy_cluster_num >=0){ 
   mol->Nconstr = noheavy_cluster_num;
 }

 free(Natom_constr);
 free(Selem_constr);
 free(new_from_orig);

} /* end of Set_Connected_Structure() */





int Mark_Neighbors_for_Connected(mol,a,Nconstr)
  struct MOLECULE *mol;
  int a;        /* index of the focused atom */
  int Nconstr;  /* number of connected graph */
{
 int b;

 mol->atoms[a].constr_num = Nconstr;

 for (b=0;b<mol->Natom;++b){
  if ((a!=b) && (mol->conmap.map[a][b] != '0')){
    if (mol->atoms[b].constr_num==-1) Mark_Neighbors_for_Connected(mol,b,Nconstr);
    }

 } 



  return(1);

} /* end of Mark_Neighbors_for_Connected() */




void Bubble_Sort_Index(X,Y, new_from_orig,N)
  int *X;           /* array of values         to be sorted [0..N-1] (larger has a priority )*/
  int *Y;           /* array of another values to be sorted [0..N-1] (smaller has a priority )*/
  int *new_from_orig; /* 0,1,2,...,N-1  */
  int  N;     /* Number of items */
{
 int i,Nchange,buff;
 int *ind; /* index */

 ind = (int *)malloc(sizeof(int)*N);
 for (i=0;i<N;++i) ind[i] = i;

 do{
   Nchange = 0;
   for (i=1;i<N;++i){
     if ( (X[i-1]<X[i]) || ((X[i-1]==X[i])&&(Y[i-1]>Y[i])) ){
       buff = X[i-1];  X[i-1] = X[i];  X[i] = buff;
       buff = Y[i-1];  Y[i-1] = Y[i];  Y[i] = buff;
       buff = ind[i-1]; ind[i-1] = ind[i]; ind[i] = buff;
       Nchange += 1;
     }
   }
 } while (Nchange>0);

 for (i=0;i<N;++i){ new_from_orig[ind[i]] = i;}

 free(ind);
} /* end of Bubble_Sort_Index() */


void Keep_Only_Largest_Connected_Structure(mol)
  struct MOLECULE *mol;
/**
 added by T.Kawabata (2018/09/07).
 
 The function  'Set_Connected_Structure(mol)' has to be done before this function.

 Following variable has to be already assigend.
   mol->Nconstr
   mol->atoms[a].constr_num

 Keep only ATOMs with mol->atoms[a].constr_num == 0.

**/
{
  struct MOLECULE newmol;
  int a,b,aa,bb,B,BB;
  int *atomnum_new_frm_org; 
  int *bondnum_new_frm_org; 


  Initialize_MOLECULE(&newmol);
  newmol.Natom = 0;
  newmol.Nbond = 0;

  atomnum_new_frm_org = (int *)malloc(sizeof(int)*mol->Natom);
  bondnum_new_frm_org = (int *)malloc(sizeof(int)*mol->Nbond);

  /* [1] count newmol.Natom */
  aa = 0;
  for (a=0;a<mol->Natom;++a){
    atomnum_new_frm_org[a] = -1;
    if (mol->atoms[a].constr_num==0){
       atomnum_new_frm_org[a] = aa;
       aa += 1;
    }

  }
  newmol.Natom = aa;

  /* [2] count newmol.Nbond */
  BB = 0;
  for (B=0;B<mol->Nbond;++B){
    bondnum_new_frm_org[B] = -1;
    if ((mol->atoms[mol->bonds[B].anum1].constr_num==0) && (mol->atoms[mol->bonds[B].anum2].constr_num==0)){
      bondnum_new_frm_org[B] = BB;
      BB += 1;
    }
  }
  newmol.Nbond = BB;

  printf("#Keep_Only_Largest_Connected_Structure() Natom %d --> %d  Nbond %d --> %d\n", mol->Natom, newmol.Natom, mol->Nbond, newmol.Nbond);

  if (newmol.Natom < mol->Natom){
    printf("#KPLCN_NATOM_ORIG %d\n", mol->Natom);
    printf("#KPLCN_NATOM_NEW  %d\n", newmol.Natom);

    Malloc_MOLECULE(&newmol,newmol.Natom,newmol.Nbond);

 
    /* [3] Copy_ATOM from a to atomnum_new_frm_org[a] */
    for (a=0;a<mol->Natom;++a){
      aa = atomnum_new_frm_org[a];
      if ((aa>=0)&&(aa<newmol.Natom)){
        Copy_ATOM(&(newmol.atoms[aa]),&(mol->atoms[a]));
        newmol.atoms[aa].num = aa;  /* <= this is very important !! ATOM.num is defined in [0..Natom-1] */
    
      }
    }


    /* [4] Set newmol.bonds[] */
    for (B=0;B<mol->Nbond;++B){
      BB = bondnum_new_frm_org[B];
      if ((BB>=0)&&(BB<newmol.Nbond)){ 
        newmol.bonds[BB].num   = BB;
        newmol.bonds[BB].anum1 = atomnum_new_frm_org[mol->bonds[B].anum1];
        newmol.bonds[BB].anum2 = atomnum_new_frm_org[mol->bonds[B].anum2];
        newmol.bonds[BB].bondtype   = mol->bonds[B].bondtype;
        newmol.bonds[BB].bondstereo = mol->bonds[B].bondstereo;
      }
    }


    /* [5] Set newmol.conmap */
    for (a=0;a<mol->Natom;++a){
      aa = atomnum_new_frm_org[a];
      if ((aa>=0)&&(aa<newmol.Natom)){
        for (b=0;b<mol->Natom;++b){
          bb = atomnum_new_frm_org[b];
          if ((bb>=0)&&(bb<newmol.Natom)){
            newmol.conmap.map[aa][bb] = mol->conmap.map[a][b];
          }
        }
      }
    }


    Set_MOLECULE(&newmol,'B');
    Show_Ring_Structure(&newmol,&(newmol.HeadRing),"newmol.HeadRing");
    Show_Ring_Structure(&newmol,&(newmol.HeadBlock),"newmol.HeadBlock");
    Copy_MOLECULE(mol,&newmol);
/*
    Show_Ring_Structure(mol,&(mol->HeadRing),"mol.HeadRing");
    Show_Ring_Structure(mol,&(mol->HeadBlock),"mol.HeadBlock");
 */
    Free_MOLECULE(&newmol); 
  }

  free(bondnum_new_frm_org);
  free(atomnum_new_frm_org);

} /* end of Keep_Only_Largest_Connected_Structure() */
