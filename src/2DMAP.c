/*
 
 <2DMAP.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
*/ 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include "globalvar.h"
#include "ioLib.h" 
#include "2DMAP.h" 
#include "molecule.h" 
#include "options.h" 

/** FUNCTIONS (GLOBAL) **/
int Malloc_CHAR2DMAP();
int Malloc_UNSIGNED_CHAR2DMAP();
int Malloc_FLOAT2DMAP();
int Malloc_INT2DMAP();
int Free_CHAR2DMAP();
int Free_UNSIGNED_CHAR2DMAP();
int Free_FLOAT2DMAP();
int Free_INT2DMAP();
void Write_Distance_FLOAT2DMAP_in_Phylip_format();
void show_conmap();

/** FUNCTIONS (LOCAL) **/



int Malloc_CHAR2DMAP(M,N)
 struct CHAR2DMAP *M;
 int N;
{
 int i,j;

 /* printf("#Malloc_CHAR2DMAP(M,N %d)\n",N); */
 M->N = N;

 M->map = (char **)malloc(sizeof(char *) * N);
 for (i=0;i<N;++i){
   M->map[i] = (char *)malloc(sizeof(char) * N);
 }

 for (i=0;i<N;++i)
  for (j=0;j<N;++j) M->map[i][j] = '0';

 M->malloced   = 1;
 M->calculated = 0;
 return(1);
} /* end of Malloc_CHAR2DMAP() */


int Malloc_UNSIGNED_CHAR2DMAP(M,N)
 struct UNSIGNED_CHAR2DMAP *M;
 int N;
{
 int i,j;

 /* printf("#Malloc_CHAR2DMAP(M,N %d)\n",N); */
 M->N = N;

 M->map = (unsigned char **)malloc(sizeof(unsigned char *) * N);
 for (i=0;i<N;++i){
   M->map[i] = (unsigned char *)malloc(sizeof(unsigned char) * N);
 }

 for (i=0;i<N;++i)
  for (j=0;j<N;++j) M->map[i][j] = 0;

 M->malloced   = 1;
 M->calculated = 0;
 return(1);
} /* end of Malloc_UNSIGNED_CHAR2DMAP() */




int Malloc_FLOAT2DMAP(M,N)
 struct FLOAT2DMAP *M;
 int N;
{
 int i,j;

 M->N = N;

 M->map = (float **)malloc(sizeof(float *) * N);
 for (i=0;i<N;++i){
   M->map[i] = (float *)malloc(sizeof(float ) * N);
 }

 for (i=0;i<N;++i){
   for (j=0;j<N;++j){ M->map[i][j] = 0.0;}
 }
 M->malloced = 1;
 M->calculated = 0;
 return(1);
} /* end of Malloc_FLOAT2DMAP() */



int Malloc_INT2DMAP(M,N)
 struct INT2DMAP *M;
 int N;
{
 int i,j;

 M->N = N;

 M->map = (int **)malloc(sizeof(int *) * N);
 for (i=0;i<N;++i){
   M->map[i] = (int *)malloc(sizeof(int ) * N);
 }

 for (i=0;i<N;++i){
   for (j=0;j<N;++j){ M->map[i][j] = 0.0;}
 }
 M->malloced = 1;
 M->calculated = 0;
 return(1);
} /* end of Malloc_INT2DMAP() */



int Free_CHAR2DMAP(M)
 struct CHAR2DMAP *M;
{
 int i;

 /* printf("int Free_CHAR2DMAP(M->malloced %d N %d)\n",M->malloced,M->N); */
 if (M->malloced==1){
  for (i=M->N-1;i>=0;--i){ 
    free(M->map[i]); 
  }
  free(M->map); 
 }
 M->N = 0;
 M->malloced = 0;
 M->calculated = 0;
 M->map = NULL;
 return(1);
} /* end of Free_CHAR2DMAP() */


int Free_UNSIGNED_CHAR2DMAP(M)
 struct UNSIGNED_CHAR2DMAP *M;
{
 int i;

 /* printf("int Free_UNSIGEND_CHAR2DMAP(M->malloced %d N %d)\n",M->malloced,M->N); */
 if (M->malloced==1){
  for (i=M->N-1;i>=0;--i){ 
    free(M->map[i]); 
  }
  free(M->map); 
 }
 M->N = 0;
 M->malloced = 0;
 M->calculated = 0;
 M->map = NULL;
 return(1);
} /* end of Free_UNSIGNED_CHAR2DMAP() */


int Free_FLOAT2DMAP(M)
 struct FLOAT2DMAP *M;
{
 int i;
 /*
 printf("#Free_FLOAT2DMAP(M) N %d malloced %d calculated %d\n",M->N,M->malloced, M->calculated); fflush(stdout);
 */
 if (M->malloced==1){
   for (i=0;i<M->N;++i){
      /* printf("i %d/%d\n",i,M->N); fflush(stdout);  */
      free(M->map[i]);
   } 
   free(M->map);
 }
 M->map = NULL;
 M->N = 0;
 M->malloced = 0;
 M->calculated = 0;
 return(1);
} /* end of Free_FLOAT2DMAP() */

int Free_INT2DMAP(M)
 struct INT2DMAP *M;
{
 int i;

 if (M->malloced==1){
  for (i=M->N-1;i>=0;--i){ 
    free(M->map[i]); 
  }
  free(M->map); 
 }
 M->N = 0;
 M->malloced = 0;
 M->calculated = 0;
 M->map = NULL;
 return(1);
} /* end of Free_INT2DMAP() */





void show_conmap(mol,comment,includeH)
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
        printf("  %c ",mol->conmap.map[i][j]);
      }
    }
    printf("\n");
   }
 }

} /* end of show_conmap() */







void Write_Distance_FLOAT2DMAP_in_Phylip_format(ofname,M,List,type,WithPropType)
  char *ofname;
  struct FLOAT2DMAP *M;  /* Score Matrix */
  struct LIBMOL *List;
  char   type;  /* 'S'core, 'D'istance*/
  char   WithPropType;  /* 'T':output with properties */
{
 FILE *fp;
 struct LIBMOL *nn;
 int Nlist,i,j;
 char name[128];
 
 printf("#Write_Distance_FLOAT2DMAP_in_Phylip_format(type '%c' WithPropType '%c')-->'%s'\n", type,WithPropType,ofname);
 fp = fopen(ofname,"w");
 if (fp==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }
 Nlist = 0;
 nn = List;
 while (nn->next!=NULL) {nn = nn->next; Nlist+=1;}

 fprintf(fp,"#COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"#START_DATE %s\n",PAR.START_DATE);
 fprintf(fp,"#END_DATE   %s\n",Get_Date_String());
 fprintf(fp,"%d\n",Nlist);
 nn = List;
 for (i=0;i<Nlist;++i){
   nn = nn->next;
   if (WithPropType=='T'){
     sprintf(name,"%s_%s",nn->mol->core_molname,nn->class);
     fprintf(fp,"%-20s ",name);
    } 
   else fprintf(fp,"%-20s ",nn->mol->core_molname); 

   for (j=0;j<Nlist;++j){
     if (type=='S') fprintf(fp," %8.3f",M->map[i][j]);
     if (type=='D') fprintf(fp," %8.3f",1.0 - M->map[i][j]);
   }
   fprintf(fp,"\n");
 }
 fclose(fp); 

} /* end of Write_Distance_FLOAT2DMAP_in_Phylip_format() */
