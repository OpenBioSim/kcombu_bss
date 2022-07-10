/*

 <io_match.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 functions for input/output with "structure MATCH"

 
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
#include "buildup.h"
#include "MrgSrtMATCH.h"
#include "options.h"
#include "moltopodis.h"
#include "gen2D.h"
#include "stereo_check.h"
#include "molprop.h"



/** FUNCTIONS (GLOBAL) **/
void Print_MATCH();
void Print_MATCHlist();
void Write_MATCH_in_rasmol_script();
void Write_MATCH_in_SDF();
void Read_MATCH();
void Write_MATCHlist();
void Paste_Coordinates_into_matched_molA_atoms();
void Write_MATCHlist_in_Horizontal_NumStr();
void Write_MATCHlist_in_Horizontal_Element();
void Write_MATCHlist_in_Horizontal_Pairwise();
void Write_Superimposed_MATCH_in_rasmol_script();
int Length_of_MATCHlist();
int Nmark_in_MATCHlist();
void Assign_Aligned_AtomNumber_to_tFactor();
int atom_num_from_num_in_file();

/** FUNCTIONS (LOCAL) **/
static void write_molA_molB_number_head_lines();
static void mark1_for_comparison_target_atoms();


void Print_MATCH(M,comment)
  struct MATCH *M;
  char *comment;
{
 int n;
 /** CAUTION: This functions output anumA + 1 and anumB + 1  !! **/
 /**          This is simply for that human favors the "natural" number **/

 printf("#MATCH '%s' Npair %d ",comment,M->Npair);
 for (n=0;n<M->Npair;++n){
    printf("(%d %d) ",M->anumA[n]+1,M->anumB[n]+1);
   }
  printf("\n");

} /* end of Print_MATCHlist() */


void Print_MATCHlist(Mlist,comment)
  struct MATCH *Mlist;
  char *comment;
{
 int n,N,Nmatch;
 struct MATCH *mn;
 
 /** CAUTION: This functions output anumA + 1 and anumB + 1  !! **/
 /**          This is simply for that human favors the "natural" number **/

 Nmatch = Length_of_MATCHlist(Mlist);
 printf("#>Nmatch %d:%s\n",Nmatch,comment);
 mn = Mlist;
 N = 0;
 while ((mn->next != NULL) && (mn->next->nodetype != 'E')){
  mn = mn->next;
  printf("#[%3d] %8.3f  %d: ",N,mn->select_dis,mn->Npair);
  for (n=0;n<mn->Npair;++n){
     /* printf("(%d %d)",mn->anumA[n]+1,mn->anumB[n]+1); */
     printf("|%d %d",mn->anumA[n]+1,mn->anumB[n]+1); 
   }
  printf("\n");
  N += 1;
 }
} /* end of Print_MATCHlist() */







void Write_MATCH_in_rasmol_script(ofname,m,molA,molB,moltypeAB,coltypeMAB,labeltype)
  char *ofname;
  struct MATCH *m;
  struct MOLECULE *molA, *molB;
  char   moltypeAB;  /* 'A' or 'B' */
  char   coltypeMAB; /* colored by 'M'atch_num, atom_num of mol'A', atom_num of mol 'B' */ 
  char   labeltype;  /* 'L':label on, otherwise:label false */ 
{
 FILE *fp;
 char color[10][16];
 int Ncolor,i,num;
 char abcstr[100];

 /* sprintf(abcstr,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"); */
 sprintf(abcstr,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789");
 
 Ncolor = 9;
 sprintf(color[0],"red");    sprintf(color[1],"redorange");
 sprintf(color[2],"orange"); sprintf(color[3],"pink");
 sprintf(color[4],"yellow"); sprintf(color[5],"green");
 sprintf(color[6],"cyan");   sprintf(color[7],"purple");
 sprintf(color[8],"blue");

 fp = fopen(ofname,"w");
 printf("#Write_MATCH_in_rasmol_script(moltypeAB %c)-->'%s'\n",moltypeAB,ofname);
  if (fp==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }
 fprintf(fp,"echo \"#matching atoms for molecule %c\"\n",moltypeAB);
 fprintf(fp,"echo \"#COMMAND %s\"\n",PAR.COMMAND);
 fprintf(fp,"background white\n");
 fprintf(fp,"set specular true\n");
 fprintf(fp,"select all\n");
 fprintf(fp,"color gray\n");
 fprintf(fp,"wireframe 20\n");
 fprintf(fp,"spacefill false\n");
 fprintf(fp,"label false\n");
 for (i=0;i<m->Npair;++i){
   num = 0;
   if (moltypeAB=='A') num = molA->atoms[m->anumA[i]].num_in_file;
   if (moltypeAB=='B') num = molB->atoms[m->anumB[i]].num_in_file;
   fprintf(fp,"select atomno==%d\n",num);
   fprintf(fp,"# PAIR molA %d %s (%s) <--> molB %d %s (%s)\n",
molA->atoms[m->anumA[i]].num_in_file, molA->atoms[m->anumA[i]].element,molA->atoms[m->anumA[i]].atomtype,
molB->atoms[m->anumB[i]].num_in_file, molB->atoms[m->anumB[i]].element,molB->atoms[m->anumB[i]].atomtype);
   /* fprintf(fp,"spacefill 100\n"); */
   /* fprintf(fp,"spacefill 75\n");  */
   fprintf(fp,"spacefill 100\n"); 
   if (labeltype=='L'){ 
    /* fprintf(fp,"label \"%c\"\n",abcstr[i%62]); */
    fprintf(fp,"label \"%d\"\n",i+1); 
   }
        if (coltypeMAB=='A')
     fprintf(fp,"color %s\n",color[(int)floor((float)Ncolor*molA->atoms[m->anumA[i]].num_heavy/(float)molA->Nheavyatom)]);
   else if (coltypeMAB=='B')
     fprintf(fp,"color %s\n",color[(int)floor((float)Ncolor*molB->atoms[m->anumB[i]].num_heavy/(float)molB->Nheavyatom)]);
   else
     fprintf(fp,"color %s\n",color[(int)floor((float)Ncolor*i/(float)m->Npair)]);
 }
 fprintf(fp,"select all\n");
 if (labeltype=='L'){
   fprintf(fp,"set fontstroke 1\n");
   fprintf(fp,"set fontsize 12\n");
 }
 fclose(fp); 

} /* end of Write_MATCH_in_rasmol_script() */




void Write_MATCH_in_SDF(ofname,m,molA,molB)
  char *ofname;
  struct MATCH *m;
  struct MOLECULE *molA, *molB;
{
 FILE *fpo;
 int i,j,a,b,Nbond;
 
 printf("#Write_MATCH_in_SDF()-->'%s'\n",ofname);

 /** [0] Count Nbond **/
 Nbond = 0;
 for (i=0;i<m->Npair;++i){
   a = m->anumA[i];
   for (j=i+1;j<m->Npair;++j){
     b = m->anumA[j];
     if (molA->conmap.map[a][b] != '0') Nbond += 1;
    }
 } 

 /** [1] Write MATCH in SDF format **/
 fpo = fopen(ofname,"w");
 if (fpo==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }
 fprintf(fpo,"%s\n",molA->name);
 fprintf(fpo,"  common_sub_structure\n");
 fprintf(fpo,"\n");
 /*  13 12  0     1  0  0  0  0  0999 V2000 */
 fprintf(fpo,"%3d%3d  0     1  0  0  0  0  0999 V2000\n",m->Npair,Nbond);
 /*    5.6720    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0 */
 for (i=0;i<m->Npair;++i){
   a = m->anumA[i];
   fprintf(fpo,"%10.4f%10.4f%10.4f %s   0  0  0  0  0  0  0  0  0  0  0  0\n",
    molA->atoms[a].Pos[0],molA->atoms[a].Pos[1],molA->atoms[a].Pos[2],molA->atoms[a].element);
 }

 for (i=0;i<m->Npair;++i){
   a = m->anumA[i];
   for (j=i+1;j<m->Npair;++j){
     b = m->anumA[j];
     if (molA->conmap.map[a][b] != '0')
       fprintf(fpo,"%3d%3d  %c  0  0  0  0\n", i+1,j+1,molA->conmap.map[a][b]);
    }
 } 
 fprintf(fpo,"M  END\n");
 fclose(fpo); 

} /* end of Write_MATCH_in_SDF() */








void Write_Superimposed_MATCH_in_rasmol_script(ofname,m,molA,molB)
  char *ofname;
  struct MATCH *m;
  struct MOLECULE *molA, *molB;
{
 FILE *fp;
 char color[10][16];
 int Ncolor,i,num;
 
 Ncolor = 9;
 sprintf(color[0],"red");    sprintf(color[1],"redorange");
 sprintf(color[2],"orange"); sprintf(color[3],"pink");
 sprintf(color[4],"yellow"); sprintf(color[5],"green");
 sprintf(color[6],"cyan");   sprintf(color[7],"purple");
 sprintf(color[8],"blue");

 fp = fopen(ofname,"w");
 if (fp==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }

 printf("#Write_MATCH_in_rasmol_script()-->'%s'\n",ofname);
 fprintf(fp,"echo \"#matching atoms for superimposed molecules\"\n");
 fprintf(fp,"echo \"#COMMAND %s\"\n",PAR.COMMAND);
 fprintf(fp,"background white\n");
 fprintf(fp,"set specular true\n");
 fprintf(fp,"select all\n");
 fprintf(fp,"color gray\n");
 fprintf(fp,"wireframe 20\n");
 fprintf(fp,"spacefill false\n");
 for (i=0;i<m->Npair;++i){
   num = molA->atoms[m->anumA[i]].num_in_file;
   fprintf(fp,"select atomno==%d && *:A\n",num);
   fprintf(fp,"spacefill 100\n");
   fprintf(fp,"color %s\n",color[(int)floor((float)Ncolor*i/(float)m->Npair)]);
 }

 for (i=0;i<m->Npair;++i){
   num = molB->atoms[m->anumB[i]].num_in_file;
   fprintf(fp,"select atomno==%d && *:B\n",num);
   fprintf(fp,"spacefill 100\n");
   fprintf(fp,"color %s\n",color[(int)floor((float)Ncolor*i/(float)m->Npair)]);
 }


 fprintf(fp,"select all\n");
 fclose(fp); 

} /* end of Write_Superimposed_MATCH_in_rasmol_script() */













void Write_MATCHlist(ofname,M,molA,molB,OutType,RankMatchOutput,HeadComment)
  char *ofname;
  struct MATCH *M;
  struct MOLECULE *molA, *molB;
  char   OutType; /* 'O'ne,'A'll */
  int    RankMatchOutput;   /* 1,2,... */
  struct LINENODE *HeadComment;
{
 FILE *fp;
 int i,j,a,b,numA,numB,N,Nnei_diff;
 struct MATCH *mn;
 int Nlen_match;
 char outmatch;
 struct LINENODE *ln;
 int    rankMatchOutput;   /* 1,2,... */

 if (RankMatchOutput>=1){
   rankMatchOutput = RankMatchOutput;
 }
 else{
   rankMatchOutput = 1;
 }

/***** FORMAT EXAMPLE ****
>1
#[num_in_fileA/B] [numA/B]   [atomnameA/B][atomtype] [EC_A/B] [ECdiff] [Nnei_diff]
2345  2495  5     24     PB   PN  P    24 22  2  0
2346  2496  6     25     O1B  O1N O1    6  6  0  0
2347  2497  7     26     O2B  O2N O1    6  6  0  0
2348  2498  8     27     O3B  O5D O    11 11  0  2
2349  2472  9     1      PA   PA  P    22 22  0  0
2350  2473  10    2      O1A  O1A O1    6  6  0  0
2351  2474  11    3      O2A  O2A O1    6  6  0  0
2352  2494  12    23     O3A  O3  O    12 12  0  0
2353  2475  13    4      O5'  O5B O    11 11  0  0
2354  2476  14    5      C5'  C5B C    13 13  0  0
2355  2477  15    6      C4'  C4B C@   18 18  0  0
2356  2478  16    7      O4'  O4B O@   15 15  0  0
:
*/

 
 printf("#Write_MATCHlist(OutType '%c' rank %d)-->'%s'\n",OutType, rankMatchOutput,ofname);
 fp = fopen(ofname,"w");
 if (fp==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }
 Nlen_match = Length_of_MATCHlist(M);
 fprintf(fp,"#>> Atom_Number MATCHing file <<\n");
 fprintf(fp,"#COMMAND '%s'\n",PAR.COMMAND);
 Set_END_DATE();
 fprintf(fp,"#DATE_START %s\n",PAR.START_DATE);
 fprintf(fp,"#DATE_END   %s\n",PAR.END_DATE);
 fprintf(fp,"#COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);
 fprintf(fp,"#AlgoType %c\n",PAR.AlgoType);
 fprintf(fp,"#ConnectGraphType %c\n",PAR.ConnectGraphType);
 fprintf(fp,"#Weight Wneiatm %.2f Wextcon %.2f Wtopodis %.2f\n",PAR.Wneiatm, PAR.Wextcon, PAR.Wtopodis);
 fprintf(fp,"#CalcFinished %c\n",PAR.calc_finish);
 fprintf(fp,"#MoleculeA %s\n",molA->filename);
 fprintf(fp,"#MoleculeB %s\n",molB->filename);
 fprintf(fp,"#FiletypeA %c\n",molA->filetype);
 fprintf(fp,"#FiletypeB %c\n",molB->filetype);
 fprintf(fp,"#NatomA %d\n",molA->Natom);
 fprintf(fp,"#NatomB %d\n",molB->Natom);
 fprintf(fp,"#NheavyatomA %d\n",molA->Nheavyatom);
 fprintf(fp,"#NheavyatomB %d\n",molB->Nheavyatom);
 fprintf(fp,"#TotalNatompair %d\n",PAR.TotalNatompair);
 fprintf(fp,"#NpermuA %d\n",molA->Npermu);
 fprintf(fp,"#NpermuB %d\n",molB->Npermu);
 fprintf(fp,"#Len_of_MATCHlist %d\n",Nlen_match);
 if (OutType=='O'){
   fprintf(fp,"#RankMatchOutput %d\n",rankMatchOutput);
 }

 ln = HeadComment;
 while (ln->next != NULL){
   ln = ln->next;
   fprintf(fp,"#%s\n",ln->line);
 }
/*
 fprintf(fp,"#[numA(1)] [num_in_fileA(2)] [atomnameA(3)] --- [numB(5)] [num_in_fileB(6)] [atomnameB(7)]\n");
 fprintf(fp,"#[numpair(8)] [atomtype(9)] [ECA(10)] [ECB(11)] [ECdiff(12)] [Nnei_diff(13)]\n");
 fprintf(fp,"#[num_in_fileA(1)] [num_in_fileB(2)] [numA(3)] [numB(4)] [atomnameA(5)] [atomnameB(6)]\n");
 fprintf(fp,"#[atomtype(7)] [EC_A(8)] [EC_B(9)] [ECdiff(10)] [Nnei_diff(11)]\n");
*/

 mn = M;
 N = 0;
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){ 
   mn = mn->next;
   N += 1;

   outmatch = 0;
   if (OutType=='A'){outmatch = 1;}   

   if ((OutType=='O')&&(N==rankMatchOutput)){ outmatch = 1; }
 

   if (outmatch==1){
     fprintf(fp,">%d\n",N);
     fprintf(fp,"#Npair_atom  %d\n",mn->Npair);
     fprintf(fp,"#tanimoto    %f\n",(float)mn->Npair/(float)(molA->Nheavyatom + molB->Nheavyatom - mn->Npair));
     fprintf(fp,"#select_dis  %f\n",mn->select_dis);
     fprintf(fp,"#Ncomponent  %d\n",mn->Ncomponent);
     if ((molA->topodismap.malloced==1)&&(molB->topodismap.malloced==1)){
       fprintf(fp,"#Maxdiff_topodis  %d\n",MaxDiff_Topological_Distance(mn,molA,molB));
     } 

/* 
           #9146  9003  1     1      PG   PG  P    20 20  0  0
           #9147  9004  2     2      O1G  O1G O1    5  5  0  0
           #9148  9005  3     3      O2G  O2G O1    5  5  0  0
           #9149  9006  4     4      O3G  O3G O1    5  5  0  0
           #9150  9010  5     8      O3B  O3B O    11 11  0  0
           #9151  9007  6     5      PB   PB  P    24 24  0  0
*/

    fprintf(fp,"#[num_in_fileA:1] [num_in_fileB:2]  [numA:3]  [numB:4]\n");
    fprintf(fp,"#[atomnameA:5] [atomnameB:6] [atomtype:7] [EC_A:8] [EC_B:9] [ECdiff:10] [Nnei_diff:11]\n");   

    for (i=0;i<mn->Npair;++i){
       a = mn->anumA[i];
       b = mn->anumB[i];
       numA = molA->atoms[a].num_in_file;
       numB = molB->atoms[b].num_in_file;
  
/* 
       fprintf(fp,"%-5d %-5d %4s --- %-5d %-5d %4s %2d %-4s",
         a+1,numA,molA->atoms[a].atomname,b+1,numB,molB->atoms[b].atomname,i+1,molA->atoms[a].atomtype);
*/
       fprintf(fp,"%-5d %-5d %-5d %-5d %4s %4s %-4s",
         numA, numB, a+1, b+1, molA->atoms[a].atomname, molB->atoms[b].atomname,
         molA->atoms[a].atomtype);

 
        if (PAR.levelEC_Dextcon=='0')  fprintf(fp," %2d %2d %2d",molA->atoms[a].EC0,molB->atoms[b].EC0,abs(molA->atoms[a].EC0-molB->atoms[b].EC0));
        else if (PAR.levelEC_Dextcon=='1')  fprintf(fp," %2d %2d %2d",molA->atoms[a].EC1,molB->atoms[b].EC1,abs(molA->atoms[a].EC1-molB->atoms[b].EC1));
        else if (PAR.levelEC_Dextcon=='2')  fprintf(fp," %2d %2d %2d",molA->atoms[a].EC2,molB->atoms[b].EC2,abs(molA->atoms[a].EC2-molB->atoms[b].EC2));
        else if (PAR.levelEC_Dextcon=='3')  fprintf(fp," %2d %2d %2d",molA->atoms[a].EC2,molB->atoms[b].EC2,abs(molA->atoms[a].EC2-molB->atoms[b].EC2));
    
      if (PAR.Wneiatm>0.0){
        Nnei_diff = 0; 
        for (j=0;j<PAR.max_atomtype;++j)
          Nnei_diff += abs(molA->atoms[a].Nnei_atomtype[j] - molB->atoms[b].Nnei_atomtype[j]);
        fprintf(fp," %2d",Nnei_diff);
      }
 
       fprintf(fp,"\n");
     }

     fprintf(fp,"//\n");
   }
 }
 fclose(fp);
} /* end of Write_MATCHlist() */







void Read_MATCH(ifname,M,molA,molB)
  char *ifname;
  struct MATCH *M;
  struct MOLECULE *molA, *molB;
{
 FILE *fp;
 char line[MAX_LENGTH_LINE],buffA[MAX_LENGTH_LINE],buffB[MAX_LENGTH_LINE],end,atomtype[4];
 int  Npair_malloc,numA,numB,num_in_fileA,num_in_fileB; 
 int  Nword, Wsta[20],Wend[20];
 printf("#Read_MATCH('%s' M->Npair %d M->Npair_malloc %d )\n",ifname,M->Npair,M->Npair_malloc);

/*
>> FILE FORMAT (NEW) <<
>1
#[num_in_file A/B][numA/B]   [atomnameA/B][atomtype] [EC_A/B] [ECdiff] [Nnei_diff]
2345  2495  5     24     PB   PN  P    24 22  2  0
2346  2496  6     25     O1B  O1N O1    6  6  0  0
2347  2497  7     26     O2B  O2N O1    6  6  0  0
2348  2498  8     27     O3B  O5D O    11 11  0  2
2349  2472  9     1      PA   PA  P    22 22  0  0
2350  2473  10    2      O1A  O1A O1    6  6  0  0
2351  2474  11    3      O2A  O2A O1    6  6  0  0
2352  2494  12    23     O3A  O3  O    12 12  0  0
2353  2475  13    4      O5'  O5B O    11 11  0  0
2354  2476  14    5      C5'  C5B C    13 13  0  0
2355  2477  15    6      C4'  C4B C@   18 18  0  0
2356  2478  16    7      O4'  O4B O@   15 15  0  0

:
37    3984  37    22     O10  O22 O1    6  6  0  0
38    3985  38    23     C28  C23 C    13 13  0  0
39    3987  39    25     O11  O25 O1    5  5  0  0
40    3986  40    24     O12  O24 O1    5  5  0  0
50    3992  50    30     H10  H30 H     0  0  0  0
51    3991  51    29     H11  H29 H     0  0  0  0
52    3990  52    28     H12  H28 H     0  0  0  0
53    3989  53    27     H13  H27 H     0  0  0  0
54    3993  54    31     H14  H31 H     0  0  0  0


This function "Read_MATCH()" reads 
only first two fields ([num_in_file A/B]).  Other fields are optional.

*/

 fp = fopen(ifname,"r");
 if (fp==NULL){
   printf("#ERROR:Can't open atom match file '%s'\n",ifname);
   exit(1);
 }

 /** [1] Just count Npair **/
 Npair_malloc = 0;
 end = 0;
 while(feof(fp)==0){ 
   line[0] = '\0';
   fgets(line,MAX_LENGTH_LINE-1,fp);
   if (line[strlen(line)-1]=='\n') line[strlen(line)-1] = '\0';
   Split_to_Words(line,' ',&Nword,Wsta,Wend,20);
   atomtype[0] = '\0';
   if (Nword>6){
     Get_Part_Of_Line(atomtype,line,Wsta[6],Wend[6]); 
   }

   if ((Nword>=2)&&(line[0]!='#')&&(line[0]!='/')
      &&  ((PAR.includeHinMCS!='T')||(atomtype[0]!='H'))){ 
     Npair_malloc += 1;
   }
   if (line[0]=='/'){ end = 1; }
 }
 fclose(fp);

 printf("#M->Npair_malloc %d Npair_malloc %d\n",M->Npair_malloc,Npair_malloc);

 if (M->Npair_malloc==0){
   printf("#Malloc_MATCH(Npair_malloc %d)\n",Npair_malloc);
   Malloc_MATCH(M,Npair_malloc);
 } 
 printf("#M->Npair_malloc %d Npair_malloc %d\n",M->Npair_malloc,Npair_malloc);
 /** [2] Read the file **/
 fp = fopen(ifname,"r");
 M->Npair = 0;
 end = 0;
 while(feof(fp)==0){ 
   line[0] = '\0';
   fgets(line, MAX_LENGTH_LINE-1, fp);

   if (line[strlen(line)-1]=='\n') line[strlen(line)-1] = '\0';
   Split_to_Words(line,' ',&Nword,Wsta,Wend,20);
   if ((Nword>=2)&&(line[0]!='#')&&(line[0]!='/')){
     Get_Part_Of_Line(buffA,line,Wsta[0],Wend[0]); 
     num_in_fileA = atoi(buffA);
     numA = atom_num_from_num_in_file(molA,num_in_fileA);
     Get_Part_Of_Line(buffB,line,Wsta[1],Wend[1]); 
     num_in_fileB = atoi(buffB); 
     numB = atom_num_from_num_in_file(molB,num_in_fileB);
     if ((numA>=0) && (numA < molA->Natom)&&(numB>=0) && (numB < molB->Natom) && (M->Npair < M->Npair_malloc) 
          &&  ((PAR.includeHinMCS!='T')||(atomtype[0]!='H'))){ 
       printf("#MATCH[%d]  A %d %d %s %s %s", M->Npair,num_in_fileA, numA, molA->atoms[numA].resi, molA->atoms[numA].rnum, molA->atoms[numA].atomname);
       printf(" ----  B %d %d %s %s %s\n",num_in_fileB, numB, molB->atoms[numB].resi, molB->atoms[numB].rnum, molB->atoms[numB].atomname);
       M->anumA[M->Npair] = numA;
       M->anumB[M->Npair] = numB;
       printf("#numB %d ",numB);
       numA = M->anumA[M->Npair];
       numB = M->anumB[M->Npair];
       printf("--> numB %d\n",numB);
       printf("#MATCH[%d]  A %d %d %s %s %s", M->Npair,num_in_fileA, numA, molA->atoms[numA].resi, molA->atoms[numA].rnum, molA->atoms[numA].atomname);
       printf(" ----  B %d %d %s %s %s\n",num_in_fileB, numB, molB->atoms[numB].resi, molB->atoms[numB].rnum, molB->atoms[numB].atomname);



       M->Npair += 1;

     }
     else if ((numA<0) || (numA > molA->Natom)) {printf("#IMPROPER numA '%s' in ATOM_MATCH_FILE line:'%s'\n",buffA,line); exit(1);}
     else if ((numB<0) || (numB > molB->Natom)) {printf("#IMPROPER numB '%s' in ATOM_MATCH_FILE line:'%s'\n",buffB,line); exit(1);}

   }
   if (line[0]=='/'){ end = 1; }
 }

  for (Nword=0;Nword<M->Npair;++Nword){
     numA = M->anumA[Nword];
     numB = M->anumB[Nword];
     printf("#MATCH A %d %d %s %s %s",  molA->atoms[numA].num_in_file, numA, molA->atoms[numA].resi, molA->atoms[numA].rnum, molA->atoms[numA].atomname);
     printf(" ----  B %d %d %s %s %s\n",molB->atoms[numB].num_in_file, numB, molB->atoms[numB].resi, molB->atoms[numB].rnum, molB->atoms[numB].atomname);
  } 


 fclose(fp);


} /* end of Read_MATCH() */










void Paste_Coordinates_into_matched_molA_atoms(m,molA,molB,PasteDimension)
  struct MATCH *m;
  struct MOLECULE *molA, *molB;
  char   PasteDimension;   /* '2'D , '3'D */
{
 int i,j,k,a,b,anumA1,anumA2,anumB1,anumB2,buff;
 int NneighborA,NneighborB,neighborsA[10],neighborsB[10],Npermu;
 
 printf("#Paste_Coordiantes_into_matched_molA_atoms(PasteDimension '%c')\n",PasteDimension);
 for (i=0;i<molA->Natom;++i){
   molA->atoms[i].Pos[0] = molA->atoms[i].Pos[1] = molA->atoms[i].Pos[2] = 0.0; 
 } 
 for (j=0;j<molB->Natom;++j) molB->atoms[j].mark = -1;
 for (i=0;i<m->Npair;++i){
   molB->atoms[m->anumB[i]].mark = m->anumA[i]; /* molB->atoms[].mark : atom number of corresponding molA */
 }

 /** [1] paste ATOM coordinates (XY or XYZ) **/
 for (i=0;i<m->Npair;++i){
   a = m->anumA[i];
   b = m->anumB[i];
   molA->atoms[a].Pos[0]  = molB->atoms[b].Pos[0]; 
   molA->atoms[a].Pos[1]  = molB->atoms[b].Pos[1];
   if (PasteDimension=='2'){
     molA->atoms[a].Pos[2]  = 0.0;
   }
   else{
     molA->atoms[a].Pos[2]  = molB->atoms[b].Pos[2];
   } 
 }

 /** [2] paste ATOM stereo_parity **/
 for (i=0;i<m->Npair;++i){
   a = m->anumA[i];
   b = m->anumB[i];
   molA->atoms[a].stereo_parity = '0';

   if ((molB->atoms[b].stereo_parity=='1') || (molB->atoms[b].stereo_parity=='2')){
     NneighborA = NneighborB = 0;
     for (j=0;j<molA->atoms[a].Nneighbor;++j){
       k = molA->atoms[a].neighbors[j];
       if ((PAR.includeHinMCS=='T') || (molA->atoms[k].one_char_ele != 'H')){
          neighborsA[NneighborA] = k;
          neighborsB[NneighborB] = molA->atoms[a].mark;
          if (neighborsB[NneighborB]>=0){ NneighborB += 1;}
          NneighborA +=1;
         }
     }
     if (NneighborA==NneighborB){
       Npermu = neighbors_permutation_number(NneighborA,neighborsA,neighborsB);
       if ((Npermu%2)==0)  molA->atoms[a].stereo_parity = molB->atoms[b].stereo_parity;
       if ((Npermu%2)==1){ 
         if (molB->atoms[b].stereo_parity == '1') molA->atoms[a].stereo_parity = '2';
         if (molB->atoms[b].stereo_parity == '2') molA->atoms[a].stereo_parity = '1';
       }
     }
   }
 }

 /** [3] paste BOND information (bondtype, bondstereo) **/
 for (j=0;j<molB->Nbond;++j){
   anumB1 = molB->bonds[j].anum1; 
   anumB2 = molB->bonds[j].anum2; 
   if ((molB->atoms[anumB1].mark>=0)&& (molB->atoms[anumB2].mark>=0)){
     anumA1 = molB->atoms[anumB1].mark; 
     anumA2 = molB->atoms[anumB2].mark; 
     /* printf("#BOND_B bondstereo %c B(%d %d) A (%d %d)\n",molB->bonds[j].bondstereo,anumB1,anumB2,anumA1,anumA2); */
     i = 0;
     while (i<molA->Nbond){

       if ((molA->bonds[i].anum1 == anumA1) && (molA->bonds[i].anum2 == anumA2)){
         molA->bonds[i].bondstereo = molB->bonds[j].bondstereo;
         molA->bonds[i].bondtype   = molB->bonds[j].bondtype;
         i = molA->Nbond + 1;
       }

       if ((molA->bonds[i].anum1 == anumA2) && (molA->bonds[i].anum2 == anumA1)){
         molA->bonds[i].bondstereo = molB->bonds[j].bondstereo;
         molA->bonds[i].bondtype   = molB->bonds[j].bondtype;
         buff = molA->bonds[i].anum1;
         molA->bonds[i].anum1 = molA->bonds[i].anum2;
         molA->bonds[i].anum2 = buff;
         i = molA->Nbond + 1;
       }
       i += 1;
     } /* while (i) */
   }
 } /* j */

 /** [4] If '2'D, generate XY for hydrogens **/
 if (PasteDimension=='2'){
   set_hydrogen_XY(molA);
 }

} /* end of Paste_Coordinates_into_matched_molA_atoms() */
 






void Write_MATCHlist_in_Horizontal_NumStr(ofname,M,molA,molB)
  char *ofname;
  struct MATCH *M;
  struct MOLECULE *molA, *molB;
{
 FILE *fp;
 int i,a,b,N;
 struct MATCH *mn;
 int *align_onA,*align_onB;
 char abcstr[100];
 
 mark1_for_comparison_target_atoms(molA);
 mark1_for_comparison_target_atoms(molB);

 sprintf(abcstr,"1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"); 
 if (ofname[0]=='-') fp = stdout;
 else { fp = fopen(ofname,"w");
   printf("#Write_MATCHlist_in_Horizontal()-->'%s'\n",ofname); 
   if (fp==NULL){
     printf("#ERROR:Can't write to '%s'\n",ofname);
     exit(1);
    }
 }
 
 align_onA = (int *)malloc(sizeof(int)*(molA->Natom+1)); 
 align_onB = (int *)malloc(sizeof(int)*(molB->Natom+1)); 

 write_molA_molB_number_head_lines(fp,molA,molB,M,'T');

 mn = M;
 N = 0;
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){ 
   mn = mn->next;
   N += 1;
   for (a=0;a<molA->Natom;++a) align_onA[a] = -1;
   for (b=0;b<molB->Natom;++b) align_onB[b] = -1;

    for (i=0;i<mn->Npair;++i){ 
      align_onA[mn->anumA[i]] = molB->atoms[mn->anumB[i]].num_heavy;
      align_onB[mn->anumB[i]] = molA->atoms[mn->anumA[i]].num_heavy;
    }

   fprintf(fp,"#[%3d]%3d %5.3f %5.1f:",
     N,mn->Npair,Tanimoto_Coefficient(mn->Npair,molA,molB),mn->select_dis);
  
    for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
       if (align_onA[a]>=0) fprintf(fp,"%c",molB->atoms[align_onA[a]].one_char_ele);
                       else  fprintf(fp,"-");
    }
   }
   fprintf(fp,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
       if (align_onB[b]>=0) fprintf(fp,"%c",molA->atoms[align_onB[b]].one_char_ele);
                       else fprintf(fp,"-");
    }
   }
   fprintf(fp,"\n");
   fprintf(fp,"                     "); 
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
       if (align_onA[a]>=0) fprintf(fp,"%c",abcstr[align_onA[a]%62]);
                       else  fprintf(fp,"-");
    }
   }
   fprintf(fp,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
       if (align_onB[b]>=0) fprintf(fp,"%c",abcstr[align_onB[b]%62]);
                       else fprintf(fp,"-");
    }
   }
  fprintf(fp,"\n");


 }

 if (fp!=stdout) fclose(fp);

 free(align_onB);
 free(align_onA);

} /* end of Write_MATCHlist_in_Horizontal_NumStr() */




void Write_MATCHlist_in_Horizontal_Element(ofname,M,molA,molB)
  char *ofname;
  struct MATCH *M;
  struct MOLECULE *molA, *molB;
{
 FILE *fp;
 int i,a,b,N;
 struct MATCH *mn;
 char *align_onA,*align_onB;
 
 if (ofname[0]=='-') fp = stdout;
 else { fp = fopen(ofname,"w");
   printf("#Write_MATCHlist_in_Horizontal()-->'%s'\n",ofname); 
   if (fp==NULL){
     printf("#ERROR:Can't write to '%s'\n",ofname);
     exit(1);
   }
 }
 
 mark1_for_comparison_target_atoms(molA);
 mark1_for_comparison_target_atoms(molB);
 align_onA = (char *)malloc(sizeof(char)*(molA->Natom+1)); 
 align_onB = (char *)malloc(sizeof(char)*(molB->Natom+1)); 

 write_molA_molB_number_head_lines(fp,molA,molB,M,'F');

 mn = M;
 N = 0;
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){ 
   mn = mn->next;
   N += 1;
   for (a=0;a<molA->Natom;++a){ align_onA[a] = -1;}
   for (b=0;b<molB->Natom;++b){ align_onB[b] = -1;}

    for (i=0;i<mn->Npair;++i){ 
      align_onA[mn->anumA[i]] = mn->anumB[i];
      align_onB[mn->anumB[i]] = mn->anumA[i];
    }
   fprintf(fp,"#[%3d]%3d %5.3f %5.1f:",
     N,mn->Npair,Tanimoto_Coefficient(mn->Npair,molA,molB),mn->select_dis);
/*
   fprintf(fp,"[%3d]%3d %5.3f %5.1f %3.0f %3.0f %3.0f:",
     N,mn->Npair,Tanimoto_Coefficient(mn->Npair,molA->Nheavyatom,molB->Nheavyatom),mn->select_dis,mn->Dneiatm,mn->Dextcon,mn->Dtopodis);

*/
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
       if (align_onA[a]>=0){ 
            fprintf(fp,"%c",molB->atoms[(int)align_onA[a]].one_char_ele);
       }
        else  fprintf(fp,"-");
    }
   }
   fprintf(fp,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
       if (align_onB[b]>=0){ 
            fprintf(fp,"%c",molA->atoms[(int)align_onB[b]].one_char_ele);
       }
        else  fprintf(fp,"-");
    }
   }


  fprintf(fp,"\n");
 }


 fprintf(fp,"\n");
 /** Write AtomPair in Number **/ 
 fprintf(fp,"#   Nmcs|tani|seldis:\n"); 
 mn = M;
 N = 0;
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){ 
   mn = mn->next;
   N += 1;
   fprintf(fp,"#[%3d]%3d %5.3f %5.1f:",
     N,mn->Npair,Tanimoto_Coefficient(mn->Npair,molA,molB),mn->select_dis);
  /* for (i=0;i<PAR.max_atomtype;++i){ fprintf(fp," %d",mn->Natomtype[i]);} */
   for (i=0;i<mn->Npair;++i){
     fprintf(fp,"|%d %d",molA->atoms[mn->anumA[i]].num_in_file,molB->atoms[mn->anumB[i]].num_in_file);
   }
   fprintf(fp,"\n");
 }
 
 if (fp!=stdout) fclose(fp);

 free(align_onB);
 free(align_onA);

} /* end of Write_MATCHlist_in_Horizontal_Element() */







void Write_MATCHlist_in_Horizontal_Number(ofname,M,molA,molB)
  char *ofname;
  struct MATCH *M;
  struct MOLECULE *molA, *molB;
{
 FILE *fp;
 int i,a,b,N;
 struct MATCH *mn;
 int *align_onA,*align_onB;
 
 if (ofname[0]=='-') fp = stdout;
 else { fp = fopen(ofname,"w");
 printf("#Write_MATCHlist_in_Horizontal()-->'%s'\n",ofname);
   if (fp==NULL){
     printf("#ERROR:Can't write to '%s'\n",ofname);
     exit(1);
   }
 }
 mark1_for_comparison_target_atoms(molA);
 mark1_for_comparison_target_atoms(molB);
 
 align_onA = (int *)malloc(sizeof(int)*(molA->Natom+1)); 
 align_onB = (int *)malloc(sizeof(int)*(molB->Natom+1)); 

 write_molA_molB_number_head_lines(fp,molA,molB,M,'F');

 mn = M;
 N = 0;
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){ 
   mn = mn->next;
   N += 1;
   for (a=0;a<molA->Natom;++a) align_onA[a] = -1;
   for (b=0;b<molB->Natom;++b) align_onB[b] = -1;

    for (i=0;i<mn->Npair;++i){ 
      align_onA[mn->anumA[i]] = mn->anumB[i];
      align_onB[mn->anumB[i]] = mn->anumA[i];
    }

   fprintf(fp,"#[%3d]%3d %5.3f %5.2f:",
     N,mn->Npair,Tanimoto_Coefficient(mn->Npair,molA,molB),mn->select_dis);
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
       if (align_onA[a]>=0) fprintf(fp,"%c",molB->atoms[align_onA[a]].one_char_ele);
                      else  fprintf(fp,"-");
    }
   }
   fprintf(fp,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
       if (align_onB[b]>=0) fprintf(fp,"%c",molA->atoms[align_onB[b]].one_char_ele);
                       else fprintf(fp,"-");
    }
   }

  fprintf(fp,"\n");
  fprintf(fp,"     :");
  for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
       if (align_onA[a]>=0){
           if (molB->atoms[align_onA[a]].num_in_file>=10){ 
            fprintf(fp,"%d",(int)floor(molB->atoms[align_onA[a]].num_in_file/10)); } 
            else fprintf(fp," ");
       }  else fprintf(fp,"-");
    }
   }
  fprintf(fp,":");
  for (b=0;b<molA->Natom;++b){
     if (molB->atoms[b].mark==1){
       if (align_onB[b]>=0){
        if (molA->atoms[align_onB[b]].num_in_file>=10){ 
            fprintf(fp,"%d",(int)floor(molA->atoms[align_onB[b]].num_in_file/10)); } 
         else fprintf(fp," ");
       } else fprintf(fp,"-");
     }
   }
  fprintf(fp,"\n");

  fprintf(fp,"     :");
  for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
       if (align_onA[a]>=0){ 
            fprintf(fp,"%d",(int)floor(molB->atoms[align_onA[a]].num_in_file%10)); } else fprintf(fp,"-");
    }
   }
  fprintf(fp,":");
  for (b=0;b<molA->Natom;++b){
     if (molB->atoms[b].mark==1){
       if (align_onB[b]>=0){ 
            fprintf(fp,"%d",(int)floor(molA->atoms[align_onB[b]].num_in_file%10)); } else fprintf(fp,"-");
    }
   }
  fprintf(fp,"\n");
 }




 if (fp!=stdout) fclose(fp);

 free(align_onB);
 free(align_onA);

} /* end of Write_MATCHlist_in_Horizontal_Number() */



void write_molA_molB_number_head_lines(fpo,molA,molB,M,AbcNumOut)
  FILE *fpo;  /* file pointer for output */
  struct MOLECULE *molA, *molB; 
  struct MATCH *M;
  char AbcNumOut; /* T or F */
{
 int a,b,L,max_num_in_fileA, max_num_in_fileB,Npair_mcs;
 char abcstr[100],numstr[8];
 float tanimoto;

 mark1_for_comparison_target_atoms(molA);
 mark1_for_comparison_target_atoms(molB);

 
 sprintf(abcstr,"1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"); 
 fprintf(fpo,"#>> Atom_Number MATCHing file <<\n");
 fprintf(fpo,"#COMMAND '%s'\n",PAR.COMMAND);
 Set_END_DATE();
 fprintf(fpo,"#DATE       %s -- %s\n",PAR.START_DATE,PAR.END_DATE);
 fprintf(fpo,"#COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);
 fprintf(fpo,"#AlgoType %c\n",PAR.AlgoType);
 fprintf(fpo,"#ConnectGraphType %c\n",PAR.ConnectGraphType);
 /* 
 fprintf(fpo,"#Wneiatm %.2f Wextcon %.2f Wtopodis %.2f\n",PAR.Wneiatm, PAR.Wextcon, PAR.Wtopodis);
 fprintf(fpo,"#CalcFinished %c\n",PAR.calc_finish);
 */ 
 fprintf(fpo,"#MoleculeA %s\n",molA->filename);
 fprintf(fpo,"#MoleculeB %s\n",molB->filename);
 fprintf(fpo,"#NatomA %d\n",molA->Natom);
 fprintf(fpo,"#NatomB %d\n",molB->Natom);
 fprintf(fpo,"#NheavyatomA %d\n",molA->Nheavyatom);
 fprintf(fpo,"#NheavyatomB %d\n",molB->Nheavyatom);
 Npair_mcs = M->Npair;

 tanimoto = 0.0;
 if ((M->next != NULL) && (M->nodetype != 'E')){
    Npair_mcs = M->next->Npair;
    tanimoto = Tanimoto_Coefficient(Npair_mcs,molA,molB);
 }
 /*
 fprintf(fpo,"#SUMMARY molA %s Natm %d molB %s NatomB %d Nmcs %d tanimoto %f\n",molA->filename,molA->Nheavyatom,molB->filename,molB->Nheavyatom,Npair_mcs,tanimoto);
 */ 
 fprintf(fpo,"#Nmatched_atom %d\n",Npair_mcs);
 fprintf(fpo,"#tanimoto %f\n",tanimoto);
 fprintf(fpo,"#                    "); 
 fprintf(fpo,"#[Molecule A]");
 for (a=0;a<molA->Nheavyatom-12;++a) fprintf(fpo," ");
 fprintf(fpo,":");
 fprintf(fpo,"#[Molecule B]");
 for (b=0;b<molB->Nheavyatom-12;++b) fprintf(fpo," ");
 fprintf(fpo,"\n");

 /** ring **/ 
 fprintf(fpo,"#[block_string]      :"); 
 for (a=0;a<molA->Natom;++a){
   if (molA->atoms[a].mark==1) fprintf(fpo,"%c",molA->atoms[a].ringblock);
 }
 fprintf(fpo,":");
 for (b=0;b<molB->Natom;++b){
   if (molB->atoms[b].mark==1) fprintf(fpo,"%c",molB->atoms[b].ringblock);
 }
 fprintf(fpo,"\n");

 max_num_in_fileA = molA->atoms[molA->Natom-1].num_in_file;
 max_num_in_fileB = molB->atoms[molB->Natom-1].num_in_file;
 /** number(1000) **/ 
 if ((max_num_in_fileA>=1000)||(max_num_in_fileB>=1000)){
   fprintf(fpo,"#                     "); 
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
      sprintf(numstr,"%d",molA->atoms[a].num_in_file); L = strlen(numstr);
      if (molA->atoms[a].num_in_file>=1000) fprintf(fpo,"%c",numstr[L-1-3]); else fprintf(fpo," ");
     }
   }
   fprintf(fpo,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
       sprintf(numstr,"%d",molB->atoms[b].num_in_file); L = strlen(numstr);
       if (molB->atoms[b].num_in_file>=1000) fprintf(fpo,"%c",numstr[L-1-3]); else fprintf(fpo," ");
     }
   }
  fprintf(fpo,"\n");
 }
 
 /** number(100) **/ 
 if ((max_num_in_fileA>=100)||(max_num_in_fileB>=100)){
   fprintf(fpo,"#                     "); 
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
      sprintf(numstr,"%d",molA->atoms[a].num_in_file); L = strlen(numstr);
      if (molA->atoms[a].num_in_file>=100) fprintf(fpo,"%c",numstr[L-1-2]); else fprintf(fpo," ");
     }
   }
   fprintf(fpo,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
       sprintf(numstr,"%d",molB->atoms[b].num_in_file); L = strlen(numstr);
       if (molB->atoms[b].num_in_file>=100) fprintf(fpo,"%c",numstr[L-1-2]); else fprintf(fpo," ");
     }
   }
  fprintf(fpo,"\n");
 }
 
 /** number(10) **/ 
 if ((max_num_in_fileA>=10)||(max_num_in_fileB>=10)){
   fprintf(fpo,"#                     "); 
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1){
      sprintf(numstr,"%d",molA->atoms[a].num_in_file); L = strlen(numstr);
      if (molA->atoms[a].num_in_file>=10) fprintf(fpo,"%c",numstr[L-1-1]); else fprintf(fpo," ");
     }
   }
   fprintf(fpo,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1){
      sprintf(numstr,"%d",molB->atoms[b].num_in_file); L = strlen(numstr);
      if (molB->atoms[b].num_in_file>=10) fprintf(fpo,"%c",numstr[L-1-1]); else fprintf(fpo," ");
     }
   }
   fprintf(fpo,"\n");
 }
 
 /** number(1) **/ 
 fprintf(fpo,"#[atom_number]       :"); 
 for (a=0;a<molA->Natom;++a){
   if (molA->atoms[a].mark==1){ fprintf(fpo,"%d",molA->atoms[a].num_in_file%10); }
 }
 fprintf(fpo,":");
 for (b=0;b<molB->Natom;++b){
   if (molB->atoms[b].mark==1){ fprintf(fpo,"%d",molB->atoms[b].num_in_file%10); }
 }
 fprintf(fpo,"\n");
 
 /** one_char_element **/ 
 /*
 fprintf(fpo,"      ");
 for (a=0;a<molA->Natom;++a){
   if (molA->atoms[a].one_char_ele!='H') fprintf(fpo,"%c",molA->atoms[a].one_char_ele);
 }
 fprintf(fpo,":");
 for (b=0;b<molB->Natom;++b){
   if (molB->atoms[b].one_char_ele!='H') fprintf(fpo,"%c",molB->atoms[b].one_char_ele);
 }
 fprintf(fpo,"\n");
 */
 /** one_char_element (ring:cno) **/ 
 fprintf(fpo,"#[element]           :"); 
 for (a=0;a<molA->Natom;++a){
   if (molA->atoms[a].mark==1){ fprintf(fpo,"%c",molA->atoms[a].one_char_ele); }
 }
 fprintf(fpo,":");
 for (b=0;b<molB->Natom;++b){
   if (molB->atoms[b].mark==1){ fprintf(fpo,"%c",molB->atoms[b].one_char_ele); }
 }
 fprintf(fpo,"\n");


 /** abcstr_heavy_char **/ 
 if (AbcNumOut=='T'){
   fprintf(fpo,"                     "); 
   for (a=0;a<molA->Natom;++a){
     if (molA->atoms[a].mark==1) fprintf(fpo,"%c",abcstr[molA->atoms[a].num_heavy%62]);
   }
   fprintf(fpo,":");
   for (b=0;b<molB->Natom;++b){
     if (molB->atoms[b].mark==1) fprintf(fpo,"%c",abcstr[molB->atoms[b].num_heavy%62]);
  }
  fprintf(fpo,"\n");
 }
 /** ---line-- **/ 
 fprintf(fpo,"#    Nmcs|tani|seldis:"); 
 for (a=0;a<molA->Natom;++a){
   if (molA->atoms[a].mark==1) fprintf(fpo,"-");
 }
 fprintf(fpo,":");
 for (b=0;b<molB->Natom;++b){
   if (molB->atoms[b].mark==1) fprintf(fpo,"-");
 } 
 fprintf(fpo,"\n");

} /* end of write_molA_molB_number_head_lines() */












void Assign_Aligned_AtomNumber_to_tFactor(m,mol,moltypeAB)
  struct MATCH *m;
  struct MOLECULE *mol;
  char moltypeAB; /* 'A' or 'B' */
{
 int a,i,hit;

 for (a=0;a<mol->Natom;++a){

   hit = -1;
   i = 0;
   while ((i<m->Npair) && (hit==-1)){
    if (moltypeAB=='A') {if (a==m->anumA[i]) hit = i;}
    if (moltypeAB=='B') {if (a==m->anumB[i]) hit = i;}
    i += 1;
   } 

   if (hit>=0) mol->atoms[a].tFactor = (float)(hit+1);
          else mol->atoms[a].tFactor = 0.0;
 } 


} /* end of Assign_Aligned_AtomNumber_to_tFactor() */


void mark1_for_comparison_target_atoms(mol)
  struct MOLECULE *mol;
{
 int a;

 for (a=0;a<mol->Natom;++a){
   mol->atoms[a].mark = 0;
   if ((mol->atoms[a].one_char_ele!='H')||(PAR.includeHinMCS=='T')) mol->atoms[a].mark = 1;
 }
}



int atom_num_from_num_in_file(mol,num_in_file)
  struct MOLECULE *mol;
  int num_in_file;
{
  int a;

  for (a=0;a<mol->Natom;++a){
    if (mol->atoms[a].num_in_file == num_in_file){ return(a);}
  }
  return(-1);

} /* end of atom_num_from_num_in_file() */


