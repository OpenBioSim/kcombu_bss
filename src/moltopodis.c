/*

 <moltopodis.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  for calculating topological distance (number of bonds on the shortest path between atoms)

  Three algorithms are prepared:
   Set_Topological_Distance_MBF();  "Moore-Bellman-Ford" algorithm
   Set_Topological_Distance_FW();   "Floyd-Warshall" algorithm
   Set_Topological_Distance_Original();  "Moore-Bellman-Ford"-like algorithm discovered by myself. 

  Actually, computational times for these 3 algorithm are very similar.

'kcombu -A SupLIG/C7M_2c8wB -B SupLIG/ESH_2cf8H -fA P -fB P -mtd 1 -wtd 1.0'
   Original:#COMP_TIME  0.006417 seconds
   MBF     :#COMP_TIME  0.006444 seconds
   FW      :#COMP_TIME  0.006486 seconds

'kcombu -A SupLIG/BAY_1vijA -B SupLIG/DMP_1bveA -fA P -fB P -mtd 1 -wtd 1.0'
   Original:#COMP_TIME  0.016864 seconds
   MBF     :#COMP_TIME  0.044874 seconds
   FW      :#COMP_TIME  0.017924 seconds


 ** Floyd-Warshall algorithm has a little advantage, it does not use
   mol->atoms[a].Nneighbor and mol->atoms[a].neighbors[].

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

/** FUNCTIONS (GLOBAL) **/
void Set_Topological_Distance();
int MaxDiff_Topological_Distance();
void set_shortest_path_between_two_atoms();
void show_shortest_pathes_between_atoms();
void show_pathmap();

/** FUNCTIONS (LOCAL) **/
extern void Set_Topological_Distance_MBF();
extern void Set_Topological_Distance_FW();
extern void Set_Topological_Distance_Original();
static void show_topodismap();

void Set_Topological_Distance(mol, includeH)
  struct MOLECULE *mol;
  char   includeH; /* 'H': including hydrogen atoms */
{
  /* printf("Set_Topological_Distance()\n"); */
  Set_Topological_Distance_FW(mol,includeH);
/*
  show_topodismap(mol,"topodismap",includeH);
  show_shortest_pathes_between_atoms(mol);
  show_shortest_pathes_between_atoms(mol);
  show_pathmap(mol,"pathmap",includeH);
*/

/*
  printf("#Set_Topological_Distance(mol)\n");
  Set_Topological_Distance_MBF(mol,includeH);
  Set_Topological_Distance_Original(mol,includeH);
 */
}




void Set_Topological_Distance_MBF(mol,includeH)
  struct MOLECULE *mol;
  char   includeH; /* 'H': including hydrogen atoms */
{
 /*
  This function is based on "Moore-Bellman-Ford" algorithm. 
 
   >> CAUTION << 
     The function Set_Nneighbor_CHNOSP(mol) should be executed before executing this funcion.
     Nneighbor, neighbors[] are neccesary.
 */

 int i,a,b,c,newdis,Nchange;

 printf("#Set_Topological_Distance_MBF(mol)\n");
 /** (1) Initialize topodismap[][] **/ 
 for (a=0;a<mol->Natom;++a){
   for (b=0;b<mol->Natom;++b){ mol->topodismap.map[a][b] = mol->Natom; }
   mol->topodismap.map[a][a] = 0;
 }

 /** (2) Calculate length of shortest path **/ 
 for (a=0;a<mol->Natom;++a){ 
   if ((mol->atoms[a].one_char_ele != 'H')||(includeH=='H')){
 
     do{
       Nchange = 0;
       for (b=0;b<mol->Natom;++b){ 
         if ((mol->atoms[b].one_char_ele != 'H')||(includeH=='H')){
           for (i=0;i<mol->atoms[b].Nneighbor;++i){
             c  = mol->atoms[b].neighbors[i];
             newdis = mol->topodismap.map[a][b] + 1;
             if (mol->topodismap.map[a][c]>newdis) { mol->topodismap.map[a][c] = newdis; Nchange += 1;}
           }
         } 
       }
     } while (Nchange>0);

   }
 } 
} /* end of Set_Topological_Distance_BMH() */



void Set_Topological_Distance_FW(mol, includeH)
  struct MOLECULE *mol;
  char   includeH; /* 'H': including hydrogen atoms */
{
 /*
  This function is based on "Floyd-Warshall" algorithm. 
 */

 int a,b,c,newdis;

 /* printf("#Set_Topological_Distance_FW(mol)\n"); */
 /** (1) Initialize topodismap[][] **/ 

 for (a=0;a<mol->Natom;++a){
   for (b=0;b<mol->Natom;++b){
     if (mol->conmap.map[a][b]!='0'){  mol->topodismap.map[a][b] = 1;}
                                else{  mol->topodismap.map[a][b] =  mol->Natom;}
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
           }
         }
      }
    }
  }
 }
 
} /* end of Set_Topological_Distance_FW() */




void Set_Topological_Distance_Original(mol,includeH)
  struct MOLECULE *mol;
  char   includeH; /* 'H': including hydrogen atoms */
{
 /*
  >> CAUTION << 
  The function Set_Nneighbor_CHNOSP(mol) should be executed before executing this funcion.
   Nneighbor, neighbors[] are neccesary.
 */

 int a,b,c,i,path_min,Nassign_path;
 char replace,first;

 printf("#Set_Topological_Distance_Original(mol)\n");
 /** (1) Initialize topodismap[][] **/
 for (a=0;a<mol->Natom;++a){
   for (b=0;b<mol->Natom;++b){
     mol->topodismap.map[a][b] = 0;
   }
 }

 /** (2) Calculate length of shortest path **/
 for (a=0;a<mol->Natom;++a){
   if (mol->atoms[a].one_char_ele != 'H'){
     for (b=0;b<mol->Natom;++b) { mol->atoms[b].mark = -1; }
     mol->atoms[a].mark = 0;

     do{
        Nassign_path = 0;
        for (b=0;b<mol->Natom;++b){
          replace = first = 0;
          if ((b!=a)&&((mol->atoms[b].one_char_ele != 'H')||(includeH=='H'))){
            path_min = mol->topodismap.map[a][b];
            for (i=0;i<mol->atoms[b].Nneighbor;++i){
              c  = mol->atoms[b].neighbors[i];
              if ((c!=b)&&((mol->atoms[c].one_char_ele != 'H')||(includeH=='H'))
                   &&(mol->atoms[c].mark>=0)){
                      if (path_min==0) { path_min = mol->atoms[c].mark + 1; first = 1;}
                 else if (path_min>(mol->atoms[c].mark+1)){
                   path_min = mol->atoms[c].mark + 1;
                   replace = 1;
                   /*
 *                    printf("##a %d b %d path_min %d\n",mol->atoms[a].num_in_file,mol->atoms[b].num_in_file,path_min);
 *                                       */
                }
              }
            }
            if ((first==1) || (replace==1)){
              mol->atoms[b].mark = path_min;
              /* mol->topodismap.map[a][b] = mol->topodismap.map[b][a] = path_min; */
              mol->topodismap.map[a][b] = path_min;
              /*
 *               printf(">>a %d b %d path_min %d\n",mol->atoms[a].num_in_file,mol->atoms[b].num_in_file,path_min);
 *                             */
              Nassign_path += 1;
            }
          }
        }
      } while (Nassign_path>0);

    }
  }
} /* end of Set_Topological_Distance_Original() */




void show_shortest_pathes_between_atoms(mol)
  struct MOLECULE *mol;
{
  int i,j,k;
  int *path,Natom_in_path;

  path = (int *)malloc(sizeof(int)*mol->Natom);
  for (i=0;i<mol->Natom;++i){
    for (j=0;j<mol->Natom;++j){
      if (i!=j){
         set_shortest_path_between_two_atoms(mol,i,j,path,&Natom_in_path);
        printf("[%d-%d](%d)",mol->atoms[i].num_in_file,mol->atoms[j].num_in_file,Natom_in_path); 
        for (k=0;k<Natom_in_path;++k){
          printf(" %s%d",mol->atoms[path[k]].atomname,mol->atoms[path[k]].num_in_file); 
        }
        printf("\n");
      }
    }
  }

 free(path);

} /* end of show_shortest_pathes_between_atoms() */


void set_shortest_path_between_two_atoms(mol,s,e,path,Natom_in_path)
  struct MOLECULE *mol;
  int    s,e;            /* atom number for start and end points */
  int    *path;          /* [mol->Natom]: array of atom numbers on the path. Length of path[] should be mol->Natom. */
  int    *Natom_in_path; /* length of path[].  It should be topodismap[s][e] + 1. */ 
{
  int k;
/*
  printf("#void set_shortest_path_between_two_atoms(mol,s:%d e:%d path[s][e]:%d)\n",s,e,mol->pathmap.map[s][e]);
*/ 
  k = s;
  path[0] = k;
  *Natom_in_path = 1;
  while ((k!=e)&& ((*Natom_in_path)<mol->Natom)){
    k = mol->pathmap.map[e][k];
    if (*Natom_in_path>=mol->Natom){printf("#ERROR:Natom_in_path is over mol->Natom %d.\n",mol->Natom); exit(1);}
    path[*Natom_in_path] = k;   
    (*Natom_in_path) += 1;
  }

} /* end of set_shortest_pathes_between_atoms() */






void show_topodismap(mol,comment,includeH)
  struct MOLECULE *mol;
  char *comment;
  char   includeH; /* 'H': including hydrogen atoms */
{
 int i,j;

 printf(">%s\n",comment);
 printf("     ");
 for (j=0;j<mol->Natom;++j) {
   if ((mol->atoms[j].one_char_ele != 'H')||(includeH=='H')){
    printf("%c%2d ",mol->atoms[j].one_char_ele,mol->atoms[j].num_in_file);
   }
 }
 printf("\n");

 for (i=0;i<mol->Natom;++i){
   if ((mol->atoms[i].one_char_ele != 'H')||(includeH=='H')){
    printf("%c%2d:",mol->atoms[i].one_char_ele,mol->atoms[i].num_in_file);
    for (j=0;j<mol->Natom;++j){
      if ((mol->atoms[j].one_char_ele != 'H')||(includeH=='H')){
        printf(" %3d",mol->topodismap.map[i][j]);
      }
    }
    printf("\n");
   }
 }

} /* end of show_topodismap() */


void show_pathmap(mol,comment,includeH)
  struct MOLECULE *mol;
  char *comment;
  char   includeH; /* 'H': including hydrogen atoms */
{
 int i,j;

 printf(">%s\n",comment);
 printf("     ");
 for (j=0;j<mol->Natom;++j) {
   if ((mol->atoms[j].one_char_ele != 'H')||(includeH=='H')){
    printf("%c%2d ",mol->atoms[j].one_char_ele,mol->atoms[j].num_in_file);
   }
 }
 printf("\n");

 for (i=0;i<mol->Natom;++i){
   if ((mol->atoms[i].one_char_ele != 'H')||(includeH=='H')){
    printf("%c%2d:",mol->atoms[i].one_char_ele,mol->atoms[i].num_in_file);
    for (j=0;j<mol->Natom;++j){
      if ((mol->atoms[j].one_char_ele != 'H')||(includeH=='H')){
        printf(" %3d",mol->atoms[mol->pathmap.map[i][j]].num_in_file);
      }
    }
    printf("\n");
   }
 }

} /* end of show_pathmap() */




int MaxDiff_Topological_Distance(M,molA,molB)
  struct MATCH    *M;
  struct MOLECULE *molA,*molB;
{
  int i,j,ai,bi,aj,bj,diff,maxdiff;

  maxdiff = 0;
  for (i=0;i<M->Npair;++i){
    ai = M->anumA[i];
    bi = M->anumB[i];
    for (j=i+1;j<M->Npair;++j){
      aj = M->anumA[j];
      bj = M->anumB[j];
      diff = abs(molA->topodismap.map[ai][aj]-molB->topodismap.map[bi][bj]);
      if (diff>maxdiff) maxdiff = diff;
    }
  }

  return(maxdiff);

} /* end of MaxDiff_Topological_Distance() */


