/*

 <AtmPairDesc.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

  for calculating atom pair descriptors

*/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "molprop.h"
#include "ringblock.h"
#include "RingDesc.h"
#include "DualAtmPair.h"
#include "options.h"
#include "ioLINE.h"

#define MAX_CHAR_LINE 1024


/** FUNCTIONS (GLOBAL) **/
void Set_AtomPair_Descriptor();
void Write_Headers_of_AtomPair_Descriptor();
int  Read_Headers_of_AtomPair_Descriptor();
void Write_Footers_of_AtomPair_Descriptor();
void Write_AtomPair_Descriptor();
void write_atompair_descriptor_for_one_molecule();
float similarity_bwn_atompair_descriptors();
float similarity_bwn_atompair_descriptors_oneatom_normalized();
int substructure_condition_atompair_descriptors();
int isomorphic_condition_atompair_descriptors();
int Read_AtomPair_Descriptor_into_LIBMOL();
struct LIBMOL* Read_AtomPair_Descriptor_into_LIBMOL_One_by_One();


/** FUNCTIONS (LOCAL) **/
static int Get_Two_Integers_at_Head_of_Line();
static void get_param_from_index_of_atompair_descriptor();
static void get_param_from_index_of_ringblock_atompair_descriptor();


void Set_AtomPair_Descriptor(mol)
  struct MOLECULE *mol;
{
  /*

   Desctiptor : [atom_type(a)] x [atom_type(b)] x [separation(s)]
    
   (Note that a >= b.)

   A = PAR.max_atomtype
   S = PAR.max_separation 
   R = PAR.ring_atomtype (R<=A) 

   Apair = A*(A+1)/2
   
   if (ringblock_atomtype != 'T')
     index = Apair * (s-1) + a*(a+1)/2 + b

   if (ringblock_atomtype == 'T')
     if (not ringblock)
       index = Apair * (s-1) + a*(a+1)/2 + b                
     else if (ringblock) && (a<R) && (b<R):
       index = Apair * S  + Rpair * (s-1) + a*(a+1)/2 + b

  */
  int a, b, s, aa, bb, i, j, k, ind, Apair,Rpair;

/*
  printf("#void Set_AtomPair_Descriptor(max_atomtype %d max_separation %d max_atompair_descriptor %d)\n",PAR.max_atomtype, PAR.max_separation, PAR.max_atompair_descriptor); fflush(stdout);
 */
  Apair = (PAR.max_atomtype  * (PAR.max_atomtype + 1))/2;
  Rpair = (PAR.ring_atomtype * (PAR.ring_atomtype + 1))/2;

  if (PAR.type_atompair=='R'){
    PAR.max_atompair_descriptor   = Apair * PAR.max_separation + Rpair * PAR.max_separation;
  }
  else{
    PAR.max_atompair_descriptor   = Apair * PAR.max_separation;
  }


 /*
  printf("#max_atomtype %d max_separation %d max_atomtype_descriptor %d\n",
     PAR.max_atomtype, PAR.max_separation, PAR.max_atompair_descriptor);

  printf("#max_atomtype %d Apair %d PAR.max_separation %d  max_atompair_descriptor %d\n",PAR.max_atomtype,Apair,PAR.max_separation, PAR.max_atompair_descriptor);  
 */

  mol->atompair_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor);
 
  for (k=0;k<PAR.max_atompair_descriptor;++k){ mol->atompair_descriptor[k] = 0;}
  
  if (mol->topodismap.malloced!=1) {
   printf("#ERROR(Set_AtomPair_Descriptor):topodismap is not malloced for '%s' '%s'\n",mol->name,mol->filename); 
   exit(1);
   }
  


  for (i=0;i<mol->Natom;++i){
    if (mol->atoms[i].one_char_ele != 'H'){
      aa =  Number_of_atomtype(mol->atoms[i].atomtype);
      for (j=i+1;j<mol->Natom;++j){
        if (mol->atoms[j].one_char_ele != 'H'){
          bb =  Number_of_atomtype(mol->atoms[j].atomtype);
          if (aa > bb){a = aa; b = bb;}else {a = bb; b = aa;}
           if (mol->topodismap.map[i][j] <= PAR.max_separation){ 
            s =  mol->topodismap.map[i][j];
            if (s>PAR.max_separation){ s = PAR.max_separation; }
            if (s<1){
               printf("#ERROR(Set_AtomPair_Descriptor()):topodismap[%d][%d] is %d. Something is wrong...\n",i,j,s);
               exit(1);
            }
            ind = Apair * (s-1) + a*(a+1)/2 + b;

            if (PAR.type_atompair=='R'){
               if ((mol->atoms[i].ringblock == mol->atoms[j].ringblock) && (a<PAR.ring_atomtype) && (b<PAR.ring_atomtype)){
                 ind = Apair * PAR.max_separation  + Rpair * (s-1) + a*(a+1)/2 + b;
               }
            }

           if ((ind<PAR.max_atompair_descriptor)&&(mol->atompair_descriptor[ind]<255)){
              /* printf("#ind %d s %d a %d b %d max_atompair_descriptor %d\n",ind, s,a,b,PAR.max_atompair_descriptor); */
              mol->atompair_descriptor[ind] += 1; 
           }   
         }
        }
      }
    }  
  }
} /* end of Set_AtomPair_Descriptor() */



void Write_Headers_of_AtomPair_Descriptor(fpo, Nlist,idirlib,des_filemode,HeadLibFile,lib_filetype)
  FILE *fpo; 
  int  Nlist;
  char *idirlib;
  char des_filemode;
  struct LINENODE *HeadLibFile;
  char  lib_filetype; /* 'S'df,'2':mol2, 'P'db */
{
  struct LINENODE *ln;
  int n;
  fprintf(fpo,"#>>File_for_AtomPair_Descriptor<<\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE    %s\n",PAR.START_DATE);
  fprintf(fpo,"#ATOMTYPE_CLASS  %c\n",PAR.atomtype_class);
  fprintf(fpo,"#MAX_ATOMTYPE    %d\n",PAR.max_atomtype);
  fprintf(fpo,"#TYPE_ATOMPAIR  %c\n",PAR.type_atompair);
  fprintf(fpo,"#MAX_SEPARATION  %d\n",PAR.max_separation);
  fprintf(fpo,"#MAX_SEPARATION_CYCLIC  %d\n",PAR.max_separation_cyclic);
  fprintf(fpo,"#MAX_SEPARATION_ACYCLIC %d\n",PAR.max_separation_acyclic);
  fprintf(fpo,"#MAX_ATOMPAIR_DESCRIPTOR   %d\n",PAR.max_atompair_descriptor);
  if (idirlib[0] !='\0') fprintf(fpo,"#LIBRARY_DIRECTORY %s\n",idirlib);
  fprintf(fpo,"#LIBRARY_FILETYPE %c\n",lib_filetype);
  if (des_filemode=='B')
    fprintf(fpo,"#FILE_MODE BINARY\n");
  else
    fprintf(fpo,"#FILE_MODE TEXT\n");

  n = 0;
  ln = HeadLibFile;
  while (ln->next != NULL){
    ln = ln->next;
    n += 1;
  }

  if (n>0){
    fprintf(fpo,"#N_LIBRARY_FILE %d\n",n);
    ln = HeadLibFile;
    while (ln->next != NULL){
      ln = ln->next;
      fprintf(fpo,"#LIBRARY_FILE_%d %s\n",ln->num,ln->line);
    }
  }
 
  fprintf(fpo,"#>[num_mol] [compound_filename] [num_file] [file_offset] [Nheavyatom] [molecular_formula] [oneatom_descriptor]\n");
} /* end of Write_Headers_of_AtomPair_Descriptor() */





int Read_Headers_of_AtomPair_Descriptor(ifname,idirlib,des_filemode,HeadLibFile,lib_filetype)
  char *ifname;
  char *idirlib;
  char *des_filemode;
  struct LINENODE *HeadLibFile;
  char *lib_filetype;
{
  FILE *fp;
  char line[MAX_LENGTH_LINE],end,buff[MAX_LENGTH_LINE];
  int Nword,Wsta[50],Wend[50];
  int Nmollib,Apair;
  
  printf("#Read_Headers_of_AtomPair_Descriptor('%s')\n",ifname);
  fp = fopen(ifname,"r");
  if (fp==NULL){printf("#ERROR:Can't open descriptor file '%s'\n",ifname); exit(1);}
  end = 0;
  Nmollib = 0;

/** [1] Read head of the file **/
/*
#>>File_for_AtomPair_Descriptor<<
#COMMAND dkcombu -M D -idml ligand_multi -ode ligand_multi.des -ltail -01 -nml L
#DATE    Nov 20,2011 10:2:52
#ATOMTYPE_CLASS  K
#MAX_ATOMTYPE    12
#TYPE_ATOMPAIR  F
#MAX_SEPARATION  10
#MAX_SEPARATION_CYCLIC  5
#MAX_SEPARATION_ACYCLIC 5
#MAX_ATOMPAIR_DESCRIPTOR   780
#LIBRARY_FILETYPE 2
#FILE_MODE TEXT
#N_LIBRARY_FILE 3
#LIBRARY_FILE_0 ligand_multi/namiki201010_03.mol2
#LIBRARY_FILE_1 ligand_multi/namiki201010_01.mol2
#LIBRARY_FILE_2 ligand_multi/namiki201010_02.mol2
#>[num_mol] [compound_filename] [num_file] [file_offset] [Nheavyatom] [molecular_formula] [oneatom_descriptor]
>0 00143744-01 0 0 28 C20_H22_N2_O4_?2 C@ 11 O@ 1 C 9 O 2 O1 1 N1 2 X 2
  0  11 C@ 1 C@
  1   2 O@ 1 C@
*/


  while ((feof(fp)==0) && (end==0)){  
    line[0] = '\0';    
    fgets(line,MAX_LENGTH_LINE,fp);
    Split_to_Words(line,' ',&Nword,Wsta,Wend,50);
    
    if (strncmp(line,"#LIBRARY_FILE_",14)==0){ 
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      Add_string_to_LINENODEs(buff,HeadLibFile);
    }

    if ((strncmp(line,"#LIBRARY_DIRECTORY",18)==0)&&(idirlib[0]=='\0')){ 
      Get_Part_Of_Line(idirlib,line,Wsta[1],Wend[1]);
    }
    if (strncmp(line,"#FILE_MODE BINARY",17)==0) *des_filemode = 'B'; 
    else if (strncmp(line,"#FILE_MODE TEXT",  15)==0) *des_filemode = 'T'; 

    else if (strncmp(line,"#NMOL_IN_LIBRARY", 16)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      Nmollib = atoi(buff); 
    }
  
    else if (strncmp(line,"#LIBRARY_FILETYPE", 17)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      *lib_filetype= buff[0]; 
    }

   else if (strncmp(line,"#ATOMTYPE_CLASS",15)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      PAR.atomtype_class = buff[0]; 
    }

    else if (strncmp(line,"#TYPE_ATOMPAIR",14)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      PAR.type_atompair = buff[0]; 
    }
    else if (strncmp(line,"#MAX_SEPARATION_CYCLIC",22)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      PAR.max_separation_cyclic = atoi(buff); 
    }

    else if (strncmp(line,"#MAX_SEPARATION_ACYCLIC",23)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      PAR.max_separation_acyclic = atoi(buff); 
    }

    else if (strncmp(line,"#MAX_SEPARATION",15)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      PAR.max_separation = atoi(buff); 
    }
 
    if (line[0]=='>') end = 1;

  } /* while */

/*
** [2] Read head of the file **
ex.)
713   2 C  10 O@
723   2 N  10 C@
747   2 O1 10 C@
770   3 X  10 N@
//
#NMOL_IN_LIBRARY 207875
#DATE_START Nov 20,2011 10:2:52
#DATE_END   Nov 20,2011 10:4:56
#COMP_TIME  123.647238 seconds
*/
  fseek(fp,-240,SEEK_END);
  while (feof(fp)==0){  
    line[0] = '\0';    
    fgets(line,MAX_LENGTH_LINE,fp);
    Split_to_Words(line,' ',&Nword,Wsta,Wend,50);
    
    if (strncmp(line,"#NMOL_IN_LIBRARY", 16)==0){
      Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
      Nmollib = atoi(buff); 
    }

  }

  fclose(fp); 
  printf("#Nmollib %d\n",Nmollib);
  Set_max_atomtype_by_atomtype_class();
  Apair = (PAR.max_atomtype * (PAR.max_atomtype + 1))/2;
  if (PAR.type_atompair=='R'){  
    PAR.max_atompair_descriptor   = Apair * PAR.max_separation + (PAR.ring_atomtype)*(PAR.ring_atomtype + 1)/2 * PAR.max_separation;
  }
  else if (PAR.type_atompair=='D'){  
    PAR.max_atompair_descriptor   = Apair * (PAR.max_separation_cyclic+1)*(PAR.max_separation_acyclic + 1);
  }
  else{
    PAR.max_atompair_descriptor   = Apair * PAR.max_separation;
  } 
  printf("#atomtype_class %c max_separation %d max_atompair_descriptor %d\n",PAR.atomtype_class, PAR.max_separation,PAR.max_atompair_descriptor);
  return(Nmollib);
} /* end of Read_Headers_of_AtomPair_Descriptor() */




void Write_Footers_of_AtomPair_Descriptor(fpo,Nmollib)
  FILE *fpo;
  int  Nmollib;
{
  Set_END_DATE();
  fprintf(fpo,"#NMOL_IN_LIBRARY %d\n",Nmollib); 
  fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
  fprintf(fpo,"#DATE_END   %s\n",PAR.END_DATE);
  fprintf(fpo,"#COMP_TIME_PROGRESS  %lf seconds\n",PAR.COMP_TIME_SEC_PROGRESS);
  fprintf(fpo,"#COMP_TIME           %lf seconds\n",PAR.COMP_TIME_SEC);
} /* end of Write_Footers_of_AtomPair_Descriptor() */





void Write_AtomPair_Descriptor(fpo,mol,num_mol,molname,class,output_mode,num_libfile,offset_libfile)
  FILE *fpo; 
  struct MOLECULE *mol;
  int  num_mol;    /* number of molecule */
  char *molname;
  char *class;
  char output_mode;     /* 'D':with description, 'R'educed_notation (only number) 'B'inary */
  int  num_libfile;     /* number of libfile  */
  long offset_libfile;  /* file offset for library file */
{
  int s,r,a,b,c,d;
  char sep_str[16]; 
  unsigned short k;
  unsigned char ch;

  /*
  fprintf(fp,">%s %d %s",molname,mol->Nheavyatom,mol->molform);
  */
  fprintf(fpo,">%d %s %d %ld %d %s",num_mol,molname,num_libfile,offset_libfile,mol->Nheavyatom,mol->molform);
  /* printf(">%d %s %d %ld %d %s",num_mol,molname,num_libfile,offset_libfile,mol->Nheavyatom,mol->molform); */


  for (a=0;a<PAR.max_atomtype;++a){
     if (mol->oneatom_descriptor[a]>0){
      fprintf(fpo," %s %d",atomtype_string_from_number(a),mol->oneatom_descriptor[a]);
      /* printf(" %s %d",atomtype_string_from_number(a),mol->oneatom_descriptor[a]); */
    }
  }
  fprintf(fpo,"\n");

  /* printf("\n"); */

  /*
  if (class != NULL){ fprintf(fp,"#class %s\n",class); }
  */

  for (k=0;k<PAR.max_atompair_descriptor;++k){
    if (mol->atompair_descriptor[k]>0){
      if (PAR.type_atompair=='R'){
        get_param_from_index_of_ringblock_atompair_descriptor(k,&r,&s,&a,&b);
        if (r==1) sprintf(sep_str,"%db",s); else sprintf(sep_str,"%d",s);
      }
      else if (PAR.type_atompair=='D'){
        get_param_from_index_of_dual_atompair_descriptor(k,&a,&b,&c,&d);
        sprintf(sep_str,"%d-%d",c,d);
      }
      else {
        get_param_from_index_of_atompair_descriptor(k,&s,&a,&b);
        sprintf(sep_str,"%d",s);
      }

      if (output_mode=='R'){ fprintf(fpo,"%3d %3d\n",k,mol->atompair_descriptor[k]); }

      else if (output_mode=='B'){ 
        fwrite(&k,1,2,fpo);
        ch = mol->atompair_descriptor[k];
        fwrite(&ch,1,1,fpo);
      }
      else{
        fprintf(fpo,"%3d %3d %-2s %s %-2s\n",
          k,mol->atompair_descriptor[k],atomtype_string_from_number(a),sep_str,atomtype_string_from_number(b));
      } 
     }
  }

  if (output_mode=='B'){ 
    k = 65535;
    c = '\n';
    fwrite(&k,1,2,fpo);
    fwrite(&c,1,1,fpo);
  }

  else{  fprintf(fpo,"//\n"); }

} /* end of Write_AtomPair_Descriptor() */




void write_atompair_descriptor_for_one_molecule(ofname,desc,comment)
  char *ofname;       /* if '-', stdout */
  unsigned char *desc;
  char *comment;
{
  int k,r,s,a,b,c,d;
  char sep_str[32];
  FILE *fpo;

  if (ofname[0]=='-')  fpo = stdout;
  else{
    fpo = fopen(ofname,"w");
    printf("#write_atompair_descriptor_for_one_molecule()-->'%s'\n",ofname);
    if (fpo==NULL){printf("#ERROR:Can't write to '%s'\n",ofname); exit(1);} 
  }
 
  fprintf(fpo,">%s\n",comment);
  for (k=0;k<PAR.max_atompair_descriptor;++k){
    if (desc[k]>0){
      sep_str[0] = '\0';
      if (PAR.type_atompair=='R'){
        get_param_from_index_of_ringblock_atompair_descriptor(k,&r,&s,&a,&b);
        if (r==1) sprintf(sep_str,"%db",s); else sprintf(sep_str,"%d",s);
      }
      else if (PAR.type_atompair=='D'){
        get_param_from_index_of_dual_atompair_descriptor(k,&a,&b,&c,&d);
        sprintf(sep_str,"%d-%d",c,d);
      }
      else {
        get_param_from_index_of_atompair_descriptor(k,&s,&a,&b);
        sprintf(sep_str,"%d",s);
      }
     fprintf(fpo,"%3d %3d %-2s %s %-2s\n",
         k,desc[k],atomtype_string_from_number(a),sep_str,atomtype_string_from_number(b));
     }
  }
  if (ofname[0]!='-'){ fclose(fpo); }

} /* end of write_atompair_descriptor_for_one_molecule() */






float similarity_bwn_atompair_descriptors(desA, desB)
   unsigned char *desA; /* descriptor for molecule A [PAR.max_atompair_descriptor] */
   unsigned char *desB; /* descriptor for molecule B [PAR.max_atompair_descriptor] */
{
  int k,common,A,B,a,b;
  float sim;

  A = B = common = 0; 
  sim = 0.0;
  for (k=0;k<PAR.max_atompair_descriptor;++k){
    a = desA[k];
    b = desB[k];
    A += a;
    B += b;
    if (a<b) common += a; else common += b;
  }
  if (PAR.simtype_descriptor=='T'){  /* Tanimoto */
    if ((A + B - common)>0)
      sim = (float)common/((float)A + (float)B - (float)common);
    else
      sim = 0.0;
  }
  else if (PAR.simtype_descriptor=='C'){   /* Cahart */
    sim = 2.0*(float)common/((float)A + (float)B);
  }
  return(sim);
}







float similarity_bwn_atompair_descriptors_oneatom_normalized(desA, desB)
   unsigned char *desA; /* descriptor for molecule A [PAR.max_atompair_descriptor] */
   unsigned char *desB; /* descriptor for molecule B [PAR.max_atompair_descriptor] */
{
  int k,common,A,B,a,b;
  float sim,Ncommon,NA,NB;
  
  A = B = common = 0; 
  sim = 0.0;
  for (k=0;k<PAR.max_atompair_descriptor;++k){
    a = desA[k];
    b = desB[k];
    A += a;
    B += b;
    if (a<b) common += a; else common += b;
  }

  Ncommon = (1.0 + sqrt(1+8.0*(float)common))/2.0;
  NA      = (1.0 + sqrt(1+8.0*(float)A))/2.0;
  NB      = (1.0 + sqrt(1+8.0*(float)B))/2.0;

  if (PAR.simtype_descriptor=='T'){
    if ((NA+NB-Ncommon)>0.0)
      sim = Ncommon/(NA + NB - Ncommon);
    else
      sim = 0.0;
  }
  else if (PAR.simtype_descriptor=='C'){ sim = 2.0*Ncommon/(NA + NB); }
  return(sim);
}







int substructure_condition_atompair_descriptors(desA, desB)
   unsigned char *desA; /* descriptor for molecule A [PAR.max_atompair_descriptor] */
   unsigned char *desB; /* descriptor for molecule B [PAR.max_atompair_descriptor] */
 /*
   Check the substructure required condition for molA is a part of molB.
 */ 
{
  int k;
  for (k=0;k<PAR.max_atompair_descriptor;++k){
    if (desA[k]>desB[k]) return(0);
  }
  return(1);

} /* substructure_condition_atompair_descriptors() */



int isomorphic_condition_atompair_descriptors(desA, desB)
   unsigned char *desA; /* descriptor for molecule A [PAR.max_atompair_descriptor] */
   unsigned char *desB; /* descriptor for molecule B [PAR.max_atompair_descriptor] */
 /*
   Check the substructure required condition for molA is a part of molB.
 */ 
{
  int k;
  for (k=0;k<PAR.max_atompair_descriptor;++k){
    if (desA[k]!=desB[k]) return(0);
  }
  return(1);
} /* isomorphic_condition_atompair_descriptors() */







int Read_AtomPair_Descriptor_into_LIBMOL(fname,HeadList,text_bin_mode,name_desc_mode)
 char   *fname;
 struct LIBMOL *HeadList;
 char   text_bin_mode;    /* 'T'ext, 'R'educe_text, 'B'inary */
 char   name_desc_mode;   /* '>' read only molecular name line */
{
/**
####  AtomPair Descriptor file in text format ###
#>[compound_filename] [file_offset] [Nheavyatom] [molecular_formula] [oneatom_descriptor] [ring_descriptor] 
>GOL 0 6 C3_O3 C 3 O1 3 a 6
#[descriptor_num] [Ncount_descriptor] [atomtype1] [separation] [atomtype2]
  0   2 C   1 C
 45   3 O1  1 C
 78   1 C   2 C
123   4 O1  2 C
201   2 O1  3 C
210   2 O1  3 O1
288   1 O1  4 O1
//
>XDL 0 8 C5_N_O2 C@ 5 N@ 1 O1 2 a 2 r6 1 b6 1
 20   4 C@  1 C@
 33   2 N@  1 C@
 50   2 O1  1 C@
 98   4 C@  2 C@
111   2 N@  2 C@
128   3 O1  2 C@
130   1 O1  2 N@
176   2 C@  3 C@
189   1 N@  3 C@
206   3 O1  3 C@
208   1 O1  3 N@
210   1 O1  3 O1
284   2 O1  4 C@
//
>TRP 0 15 C11_N2_O2 C 3 C@ 8 N@ 1 O1 2 N1 1 a 6 r5 1 r6 1 b9 1 6-0-5 6
  0   2 C   1 C
 15   1 C@  1 C
 20   8 C@  1 C@
:
518   4 O1  7 C@
528   1 N1  7 C@
596   2 O1  8 C@
//
>ATP 0 31 C10_N5_O13_P3 C 1 O 3 P 3 C@ 9 O@ 1 N@ 4 O1 9 N1 1 a 17 r5 2 r6 1 b5 1 b9 1 6-0-5 6 5-1-5 6
  1   1 O   1 C
 11   5 P   1 O
 15   1 C@  1 C
 20   5 C@  1 C@
 26   2 O@  1 C@
 33   9 N@  1 C@
 49   7 O1  1 P
 50   2 O1  1 C@
 60   1 N1  1 C@
:
756   3 O1 10 O1
761   1 N1 10 P
//

####  AtomPair Descriptor file in binary format ###
#>[compound_filename] [file_offset] [Nheavyatom] [molecular_formula] [oneatom_descriptor]
>GOL 0 6 C3_O3 C 3 O1 3 a 6
[pattern_num(2byte)][count(1byte)][pattern_num(2byte)][count(1byte)]..[65535]['\n']
>XDL 8 C5_N_O2 C@ 5 N@ 1 O1 2 a 2 r6 1 b6 1
[pattern_num(2byte)][count(1byte)][pattern_num(2byte)][count(1byte)]..[65535]['\n']
*/

 FILE *fp;
 char line[MAX_CHAR_LINE],buff[MAX_CHAR_LINE],buff2[MAX_CHAR_LINE];
 int  i,a,j,r,len,Nword,Wsta[100],Wend[100];
 struct LIBMOL *nn;
 int Nlist,  Ncount;
 char start_atompair_descriptor; 
 unsigned short k;
 unsigned char c;
 long file_offset;

 printf("#Read_AtomPair_Descriptor_into_LIBMOL('%s',text_bin_mode '%c' name_desc_mode '%c')\n",
  fname,text_bin_mode, name_desc_mode);
 /*
 printf("#max_atompair_descriptor %d\n",PAR.max_atompair_descriptor);
 */ 
 HeadList->next = NULL;
 HeadList->num  = -1;
 nn = HeadList;

 start_atompair_descriptor = 0;

 fp = fopen(fname,"r");
 if (fp==NULL){printf("#ERROR:Can't open descriptor file '%s'\n",fname); exit(1);}
 Nlist = 0;
 while (feof(fp)==0){
   if (start_atompair_descriptor==0){
     line[0] = '\0';
     file_offset = ftell(fp); 
     fgets(line,MAX_CHAR_LINE,fp);
     len = strlen(line);
     if (line[len-1]=='\n') { line[len-1]='\0'; len = len -1;}

      /*** [1] read ">" line ***/
     if ((len>0)&&(line[0]=='>')&&(isascii(line[1])!=0)){
       /**(1A) malloc LIBMOL varialble  ***/
       nn->next = (struct LIBMOL*)malloc(sizeof(struct LIBMOL));
       nn->next->next = NULL;
       nn->next->prev = nn;
       nn = nn->next;
       nn->offset_descfile = file_offset;
       nn->num = Nlist;
       nn->class = NULL;
       nn->mol = NULL;
 
       /**(1B) read name, Nheavyatom, molform   ***/
       Split_to_Words(line,' ',&Nword,Wsta,Wend,100);

       Get_Part_Of_Line(buff,line,Wsta[0]+1,Wend[0]);
       nn->num = atoi(buff);

       Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]); 
       nn->name = (char *)malloc(sizeof(char)*(Wend[1]-Wsta[1]+2));
       sprintf(nn->name,"%s",buff);

       Get_Part_Of_Line(buff,line,Wsta[2],Wend[2]); 
       nn->num_file = atoi(buff);

       Get_Part_Of_Line(buff,line,Wsta[3],Wend[3]);
       nn->offset_libfile = (long)atoi(buff);
       
       Get_Part_Of_Line(buff,line,Wsta[4],Wend[4]); 
       nn->Nheavyatom = atoi(buff);

       Get_Part_Of_Line(buff,line,Wsta[5],Wend[5]); 
       nn->molform = (char *)malloc(sizeof(char)*(Wend[5]-Wsta[5]+2));
       sprintf(nn->molform,"%s",buff);
      
       /*** (1C) interpret oneatom descriptor and ring descriptor */ 
       nn->oneatom_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atomtype);
       for (a=0;a<PAR.max_atomtype;++a) nn->oneatom_descriptor[a] = 0;

       nn->ring_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_ring_descriptor);
       for (a=0;a<PAR.max_ring_descriptor;++a) nn->ring_descriptor[a] = 0;

       i = 6; 
       while(i<Nword){
         Get_Part_Of_Line(buff, line,Wsta[i],Wend[i]); 
         Get_Part_Of_Line(buff2,line,Wsta[i+1],Wend[i+1]); 
         j = atoi(buff2); 
         if (j>256) j = 255; 
         if (isupper(buff[0])!=0){
           a = Number_of_atomtype(buff);
           if (a<PAR.max_atomtype) nn->oneatom_descriptor[a] = j;
         }
         else if (islower(buff[0])!=0){
           r = Number_of_ring_block_type(buff);
           if (r>=0) nn->ring_descriptor[r] = j;
         }
         i += 2; 
       }

       /*** (1D) malloc atom pair descriptor  */ 
       if (name_desc_mode != '>'){
         nn->atompair_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor);
         for (k=0;k<PAR.max_atompair_descriptor;++k) nn->atompair_descriptor[k] = 0;
       }

       Nlist += 1;
       start_atompair_descriptor = 1;
     }
   }
  
    /*** [2] read atom pair descriptor ***/
   else if (start_atompair_descriptor==1){
     if (text_bin_mode=='B'){
       fread(&k,1,2,fp);
       fread(&c,1,1,fp);
       if (k==65535) { start_atompair_descriptor = 0;}
       else if ((k<PAR.max_atompair_descriptor)&&(name_desc_mode != '>')){
         nn->atompair_descriptor[k] = c;
         /* printf("%d %d\n",k,c); */
       }
     }
     else{
       line[0] = '\0';
       fgets(line,MAX_CHAR_LINE,fp);
       len = strlen(line);
       if (line[len-1]=='\n') { line[len-1]='\0'; len = len -1;}
       if ((line[0]!='/') && (name_desc_mode !='>')){
         Get_Two_Integers_at_Head_of_Line(line,&k,&Ncount);
         if (Ncount>255) Ncount = 255;
         if (k<PAR.max_atompair_descriptor) nn->atompair_descriptor[k] = Ncount;
         else {printf("#WARNING k (%d) is over max_atompair\n",k);} 
       }
       if (line[0]=='/'){start_atompair_descriptor = 0;}
     }
   }
 }
 fclose(fp);
 printf("Nlist %d\n",Nlist);
 return(Nlist);

} /* end of Read_AtomPair_Descriptor_into_LIBMOL() */










struct LIBMOL* Read_AtomPair_Descriptor_into_LIBMOL_One_by_One(fp,atompair_descriptor,text_bin_mode,name_desc_mode)
 FILE  *fp; 
 unsigned char    *atompair_descriptor;
 char   text_bin_mode;    /* 'T'ext, 'R'educed_text,'B'inary */
 char   name_desc_mode;   /* 'D':read only descriptor */
/*
 Instead of "nn->atompair_descriptor[k]", use "*descriptor" given from the outside.
 This is for reducing memory usage.
*/
{
/**
####  AtomPair Descriptor file in text format ###
#>[compound_filename] [num_file] [file_offset] [Nheavyatom] [molecular_formula] [oneatom_descriptor] [ring_descriptor]
>GOL 0 0 6 C3_O3 C 3 O1 3
#[descriptor_num] [Ncount_descriptor] [atomtype1] [separation] [atomtype2]
 14   2 C  1 C
 49   3 O1 1 C
 92   1 C  2 C
127   4 O1 2 C
205   2 O1 3 C
210   2 O1 3 O1
288   1 O1 4 O1
//
>XDL 0 0 8 C5_N_O2 C@ 5 N@ 1 O1 2
  0   4 C@ 1 C@
  3   2 N@ 1 C@
 45   2 O1 1 C@
 78   4 C@ 2 C@
 81   2 N@ 2 C@
123   3 O1 2 C@
125   1 O1 2 N@
156   2 C@ 3 C@
159   1 N@ 3 C@
201   3 O1 3 C@
203   1 O1 3 N@
210   1 O1 3 O1
279   2 O1 4 C@
//
>00000001_00005000/NS-00002599.sdf 0 15 C12_N_O2 C@ 8 N@ 1 C 4 O1 2 %0 0:1:2:3:4:5
####  AtomPair Descriptor file in binary format ###
#>[compound_filename] [num_file] [file_offset] [Nheavyatom] [molecular_formula] [oneatom_descriptor]
>GOL 0 0 6 C3_O3 C 3 O1 3 a 6
[pattern_num(2byte)][count(1byte)][pattern_num(2byte)][count(1byte)]..[65535]['\n']
>XDL 0 0 8 C5_N_O2 C@ 5 N@ 1 O1 2 a 2 r6 1 b6 1
[pattern_num(2byte)][count(1byte)][pattern_num(2byte)][count(1byte)]..[65535]['\n']
*/

 char line[MAX_CHAR_LINE],buff[MAX_CHAR_LINE],buff2[MAX_CHAR_LINE];
 int  r,a,i,j,n,len,Nword,Wsta[100],Wend[100];
 struct LIBMOL *nn;
 int  Ncount;
 char start_atompair_descriptor; 
 unsigned short k;
 unsigned char c;

 /* printf("#Read_AtomPair_Descriptor_into_LIBMOL_One_by_One(fp,descriptor,text_bin_mode '%c')\n",text_bin_mode); */

 start_atompair_descriptor = 0;
 n = 0;
 nn = NULL;
 
 while (feof(fp)==0){
   if (start_atompair_descriptor==0){
     line[0] = '\0';
     fgets(line,MAX_CHAR_LINE,fp);
     len = strlen(line);
     /* printf("%d:%s\n",n,line); */
     n += 1;
     if (line[len-1]=='\n') { line[len-1]='\0'; len = len -1;}
      /*** [1] read ">" line ***/
     if ((len>0)&&(line[0]=='>')&&(isascii(line[1])!=0)){
       if (name_desc_mode != 'D'){

       /**(1A) malloc LIBMOL varialble  ***/
         nn = (struct LIBMOL*)malloc(sizeof(struct LIBMOL));
         nn->num = 0; 
         nn->score_for_sort =  nn->select_dis = 0.0; 
         nn->tanimoto_mcs = nn->tanimoto_oneatom =  nn->tanimoto_atompair =  0.0;
         nn->Npair = 0;
         nn->molform  = nn->class   = NULL;
         nn->mol = NULL;
         nn->property0 = NULL;
         nn->oneatom_descriptor = NULL; 
         nn->atompair_descriptor = NULL; 
         /**(1B) read name, Nheavyatom, molform   ***/
         Split_to_Words(line,' ',&Nword,Wsta,Wend,100);

         Get_Part_Of_Line(buff,line,Wsta[0]+1,Wend[0]); 
         nn->num = atoi(buff);
 
         Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]); 
         nn->name = (char *)malloc(sizeof(char)*(Wend[1]-Wsta[1]+2));
         sprintf(nn->name,"%s",buff);
         
         Get_Part_Of_Line(buff,line,Wsta[2],Wend[2]); 
         nn->num_file = atoi(buff);
         
         Get_Part_Of_Line(buff,line,Wsta[3],Wend[3]); 
         nn->offset_libfile = (long)atoi(buff);
 
         Get_Part_Of_Line(buff,line,Wsta[4],Wend[4]); 
         nn->Nheavyatom = atoi(buff);

         Get_Part_Of_Line(buff,line,Wsta[5],Wend[5]);
         nn->molform = (char *)malloc(sizeof(char)*(strlen(buff)+1));
         sprintf(nn->molform,"%s",buff);

         /*** (1C) interpret oneatom descriptor and ring descriptor */ 
         nn->oneatom_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atomtype);
         for (a=0;a<PAR.max_atomtype;++a) nn->oneatom_descriptor[a] = 0;
       
         nn->ring_descriptor = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_ring_descriptor);
         for (a=0;a<PAR.max_ring_descriptor;++a) nn->ring_descriptor[a] = 0;
         
         /* printf("'%s'\n",line); fflush(stdout); */
         
         i = 6; 
         while(i<(Nword-1)){
           Get_Part_Of_Line(buff, line,Wsta[i],Wend[i]); 
           Get_Part_Of_Line(buff2,line,Wsta[i+1],Wend[i+1]); 
           j = atoi(buff2); 
           /* printf("#'%s' %d/%d '%s' '%s'\n",line,i,Nword,buff,buff2);  */
           if (j>256) j = 255; 
           if (isupper(buff[0])!=0){
             a = Number_of_atomtype(buff);
             /* printf("#buff '%s' a %d\n",buff,a); */
             if (a<PAR.max_atomtype) nn->oneatom_descriptor[a] = j;
           }
           else if ((buff[0]=='%') && (buff[1]=='0')){
             nn->property0 = (char *)malloc(sizeof(char)*(strlen(buff2)+1));
             sprintf(nn->property0,"%s",buff2); 
             /* printf("property0 '%s'\n",nn->property0); */
           }
           else if (islower(buff[0])!=0){
             r = Number_of_ring_block_type(buff);
             if (r>=0) nn->ring_descriptor[r] = j;
           }
           i += 2; 
         }

       /*** (1D) malloc atom pair descriptor  */ 
        for (k=0;k<PAR.max_atompair_descriptor;++k) atompair_descriptor[k] = 0;
        start_atompair_descriptor = 1;
     }
    } 
   } /* start_atompair_descriptor == 0 */
    /*** [2] read atom pair descriptor ***/
   else if (start_atompair_descriptor==1){
     if (text_bin_mode=='B'){
       fread(&k,1,2,fp);
       fread(&c,1,1,fp);
       if (k==65535) { start_atompair_descriptor = 0; return(nn);}
       else if (k<PAR.max_atompair_descriptor) atompair_descriptor[k] = c;
       /* printf("k %d c %d\n",k,c);  */
      }
     else{
       line[0] = '\0';
       fgets(line,MAX_CHAR_LINE,fp);
       len = strlen(line);
       if (line[len-1]=='\n') { line[len-1]='\0'; len = len -1;}
       if (line[0]!='/'){
         Get_Two_Integers_at_Head_of_Line(line,&k,&Ncount);
         if (Ncount>255) Ncount = 255;
         if (k<PAR.max_atompair_descriptor) atompair_descriptor[k] = Ncount;
         else {printf("#WARNING k (%d) is over max_atompair(%d)\n",k,PAR.max_atompair_descriptor);} 
       }
       if (line[0]=='/'){start_atompair_descriptor = 0; return(nn);}
     }
   }
 } /* while */
 return(NULL);

} /* end of Read_AtomPair_Descriptor_into_LIBMOL_One_by_One() */








int Get_Two_Integers_at_Head_of_Line(line,int1,int2)
  char *line;
  int  *int1,*int2;
{
/*
   0   2 C   1 C
  45   3 O1  1 C
  78   1 C   2 C
 123   4 O1  2 C
*/
  /*
  status:ii11sss2eeeeeeee 
  string:  45   3 O1  1 C
  */

  int i,j;
  char end,state,buff[32],isnum,sym;

  i = j =  0;
  end = 0;
  state = 'i'; 
  *int1 = -1; *int2 = -1;

  while (end==0){
    isnum = 0;
    sym = line[i];
    if ((48<=sym)&&(sym<=57)) isnum = 1;

    if ((isnum==0) && (sym!=' ')&&(sym!='\n')&&(sym!='\0')) return(0);

    if (isnum==1){
      if (state=='i'){
        state = '1';
        buff[0] = sym;
        j = 1;
      }
      else if (state=='s'){
        state = '2';
        buff[0] = sym;
        j = 1;
      }
      else if ((state=='1')||(state=='2')){
        buff[j] = sym;
        j += 1;
      }
    }
    else if ((state=='1') && (line[i]==' ')) {
     buff[j] = '\0';
     *int1 = atoi(buff);
     state = 's';
    }
    else if ((state=='2') && ((line[i]==' ')||(line[i]=='\n')||(line[i]=='\0'))) {
     buff[j] = '\0';
     *int2 = atoi(buff);
     state = 's';
     return(1);
    }
    if (line[i]=='\0') end = 1;
    i += 1;
  }

  return(0);
} /* end of Get_Two_Integers_at_Head_of_Line() */
        


void get_param_from_index_of_atompair_descriptor(ind,s,a,b)
  int ind;    /* input index */
  int *s,*a,*b;  /* to be calculated */
{
  int Apair,M;
  float a0;

/*
  index = Apair * (s-1) + a*(a+1)/2 + b
*/
  Apair = (PAR.max_atomtype  * (PAR.max_atomtype + 1))/2;

  *s = ind/Apair + 1;
  M  = ind % Apair;
  a0 = (-1.0 + sqrt(8.0*M+1.0))/2.0; 
  *a = (int)floor(a0);
  *b = M - (*a)*((*a)+1)/2;
}


void get_param_from_index_of_ringblock_atompair_descriptor(ind,r,s,a,b)
  int ind;    /* input index */
  int *r,*s,*a,*b;  /* to be calculated */
{
  int Apair,Rpair,M;
  float a0;
  
  Apair = (PAR.max_atomtype  * (PAR.max_atomtype + 1))/2;
  Rpair = (PAR.ring_atomtype * (PAR.ring_atomtype + 1))/2;

  if (ind<Apair*PAR.max_separation){
    *r = 0;
    *s = ind/Apair + 1;
    M  = ind % Apair;
    a0 = (-1.0 + sqrt(8.0*M+1.0))/2.0; 
    *a = (int)floor(a0);
    *b = M - (*a)*((*a)+1)/2;
   }
  else{
    *r = 1;
    *s = (ind - Apair * PAR.max_separation) / Rpair + 1;
    M  = (ind - Apair * PAR.max_separation) % Rpair;
    a0 = (-1.0 + sqrt(8.0*M+1.0))/2.0; 
    *a = (int)floor(a0);
    *b = M - (*a)*((*a)+1)/2;
  }
}
