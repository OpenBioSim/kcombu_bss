/*
 
 <ioKCF.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
  input/output functions for KCF files


<example of KCF file : alanine>
---------------------------------------------------------------------
ENTRY       alanin                      Compound
ATOM        6
            1   O6a O     5.1350   -0.2500
            2   O6a O     4.2690    1.2500
            3   N1a N     2.5369    0.2500
            4   C1c C     3.4030   -0.2500
            5   C1a C     3.4030   -1.2500
            6   C6a C     4.2690    0.2500
BOND        5
            1     1   6 1
            2     2   6 2
            3     4   3 1 #Down
            4     4   5 1
            5     4   6 1
///
---------------------------------------------------------------------
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
int Translate_Lines_into_KCF_Molecule();
void Read_and_Append_KCF_Annotations();
void Write_KCF_Molecule();

/** FUNCTIONS (LOCAL) **/




int Translate_Lines_into_KCF_Molecule(HeadLineNode,mol)
 struct LINENODE *HeadLineNode; 
 struct MOLECULE *mol;
{
 struct LINENODE *xn; 
 char line[MAX_LENGTH_LINE],end,buff[128],status; 
 int Nline,Lline,natom,nbond,anum1,anum2;
 int Wsta[50],Wend[50],Nword;
 
 /* printf("#Translate_Lines_into_KCF_Molecule(HeadLineNode,mol)\n"); */
 Nline = natom = nbond = 0; 
 mol->name[0] = mol->comment[0] = '\0';
 mol->Natom = mol->Nheavyatom = 0;
 mol->Nbond = 0; 
 mol->conmap.N = mol->dismap.N = 0; 
 mol->conmap.malloced = mol->dismap.malloced = 0; 
 mol->Nmatch = 0;
 mol->chiral_flag = '0';
 mol->Npermu = 0; 
 mol->permuhead.next = NULL;
 status = ' ';
 end = 0;
 xn = HeadLineNode;
 while ((xn->next!=NULL)&&(end==0)){
   xn = xn->next;

   sprintf(line,"%s",xn->line);  
   /* printf("%s status %c Natom %d Nbond %d\n",line,status,mol->Natom, mol->Nbond); */
   Lline = strlen(line);
 
         if (strncmp("///",line,3)==0){ end = 1;}
   /** (1) ENTRY line ***/
   else  if (strncmp("ENTRY",line,5)==0){ 
/*01234567890123456789012345678901234567890123456789*/
/*ENTRY       alanin                      Compound*/
     Split_to_Words(line,' ',&Nword,Wsta,Wend,10);
     Get_Part_Of_Line(mol->name,line,Wsta[1],Wend[1]);
   }
   /** (2) ATOM line ***/
   else  if (strncmp("ATOM",line,4)==0){ 
/*01234567890123456789012345678901234567890123456789*/
/* ATOM        6 */
     Split_to_Words(line,' ',&Nword,Wsta,Wend,10);
     Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
     mol->Natom = atoi(buff); 
     Malloc_MOLECULE(mol, mol->Natom,0);
     status = 'A';
  }

   /** (3) BOND line ***/
  else  if (strncmp("BOND",line,4)==0){ 
/*01234567890123456789012345678901234567890123456789*/
/*BOND        5*/
     Split_to_Words(line,' ',&Nword,Wsta,Wend,10);
     Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
     mol->Nbond = atoi(buff); 
     /* Malloc_MOLECULE(mol, 0,  mol->Nbond);  */
     mol->bonds = (struct BOND*)malloc(sizeof(struct BOND)*mol->Nbond);
     status = 'B';
  }
  else  if (strncmp("    ",line,4)!=0){ status = ' ';} 

   /** (4) atom or bond information ***/
  else  if ((strncmp("    ",line,4)==0) && (status=='A')){
/*
            1   X   Cl   12.3739   -3.5272
            2   X   F    14.1174   -2.5473
            3   O2x O     2.0000   -0.4728
            4   O2a O     7.1962    0.5272
            5   O2a O     7.1962    2.5272
            6   N1y N     3.7320    0.5272
            7   N1b N    10.6766   -0.5074
            8   N5x N    10.6882    2.5619
            9   N5x N    11.5942    1.0064
            10  C1b C     4.5981    1.0272

*/
      Split_to_Words(line,' ',&Nword,Wsta,Wend,10);
      Get_Part_Of_Line(buff,line,Wsta[0],Wend[0]);
      natom = atoi(buff)-1; 
      if (natom < mol->Natom){
        mol->atoms[natom].num = natom; 
        mol->atoms[natom].num_in_file = natom + 1; 
        Get_Part_Of_Line(mol->atoms[natom].atomname,line,Wsta[1],Wend[1]);
        Get_Part_Of_Line(mol->atoms[natom].element,line,Wsta[2],Wend[2]);
        Get_Part_Of_Line(buff,line,Wsta[3],Wend[3]);
        mol->atoms[natom].Pos[0] = atof(buff); 
        Get_Part_Of_Line(buff,line,Wsta[4],Wend[4]);
        mol->atoms[natom].Pos[1] = atof(buff); 
        mol->atoms[natom].Pos[2] = 0.0;
        sprintf(mol->atoms[natom].resi,"MOL");
        sprintf(mol->atoms[natom].rnum,"   1 ");
        sprintf(mol->atoms[natom].mol2atomtype,"%s",mol->atoms[natom].element);
        mol->atoms[natom].charge = 0.0;
        mol->atoms[natom].chain = 'A';
      } 
     }
     else  if ((strncmp("    ",line,4)==0) && (status=='B')){
/*
            1     1   6 1
            2     2   6 2
            3     4   3 1 #Down
 */
      Split_to_Words(line,' ',&Nword,Wsta,Wend,10);
      Get_Part_Of_Line(buff,line,Wsta[0],Wend[0]);
      nbond = atoi(buff)-1; 
      if (nbond < mol->Nbond){ 
        Get_Part_Of_Line(buff,line,Wsta[1],Wend[1]);
        anum1 = atoi(buff)-1;
        mol->bonds[nbond].anum1 = anum1;
        Get_Part_Of_Line(buff,line,Wsta[2],Wend[2]);
        anum2 = atoi(buff)-1;
        mol->bonds[nbond].anum2 = anum2;
        Get_Part_Of_Line(buff,line,Wsta[3],Wend[3]);
        mol->bonds[nbond].bondtype = buff[0];
        mol->bonds[nbond].bondstereo = '0';
        if ((anum1<mol->conmap.N) && (anum2<mol->conmap.N)){
          mol->conmap.map[anum1][anum2] = mol->bonds[nbond].bondtype;
          mol->conmap.map[anum2][anum1] = mol->bonds[nbond].bondtype;
        }
      }
   }

   Nline += 1; 

  } /* while */

 /* printf("#name '%s' Nline %d\n",mol->name,Nline); */


 if (PAR.Wdistan>0.0){
   Malloc_FLOAT2DMAP(&(mol->dismap),mol->Natom);
   Cal_Distance_Map_for_Heavy_atoms(mol);
 }
 
 /** ending procedure **/
 if (PAR.OutCalProcess=='T')
    printf("#Read_KCF_Moleclue():Natom %3d Nbond %3d '%s'\n",mol->Natom,mol->Nbond,mol->molform);

 return(mol->Natom);

} /* end of Translate_Lines_into_KCF_Molecule() */







void Write_KCF_Molecule(ofname,mol,comment)
 char *ofname;
 struct MOLECULE *mol;
 char   *comment;
{
 FILE *fpo;
 int i,j;

 printf("#Write_KCF_Molecule()-->'%s'\n",ofname);
 fpo = fopen(ofname,"w");
 if (fpo==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }
 fprintf(fpo,"%s\n",mol->name);
 fprintf(fpo,"  %s\n",mol->comment);
 fprintf(fpo,"\n");
 /*  13 12  0     1  0  0  0  0  0999 V2000 */
 /*  39 41  0  0  1  0  0  0  0  0999 V2000 */
 /*  10  9  0     0  0  0  0  0  0999 V2000 */
 /* 107110  0     0  0  0  0  0  0999 V2000 */
 /*  30 31  0  0  1  0  0  0  0  0999 V2000 */
 fprintf(fpo,"%3d%3d  0  0  %c  0  0  0  0  0999 V2000\n",mol->Natom,mol->Nbond,mol->chiral_flag);
 /*    5.6720    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0 */
 
 for (i=0;i<mol->Natom;++i){
   fprintf(fpo,"%10.4f%10.4f%10.4f %-2s %2d  %d  %d  0  0  0  0  0  0  0  0  0\n",
     mol->atoms[i].Pos[0],mol->atoms[i].Pos[1],mol->atoms[i].Pos[2],mol->atoms[i].element,
     mol->atoms[i].mass_diff,mol->atoms[i].char_charge,mol->atoms[i].stereo_parity); 
 }
/*  1  6  1  0  0  0  0 */
 for (j=0;j<mol->Nbond;++j){
  fprintf(fpo,"%3d%3d  %c  %d  0  0  0\n",
        mol->bonds[j].anum1+1,mol->bonds[j].anum2+1,mol->bonds[j].bondtype,mol->bonds[j].bondstereo);
 }

 fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
 if (strlen(comment)>0){ 
   fprintf(fpo,"#%s\n",comment);
 }
 fprintf(fpo,"M  END\n");
 fclose(fpo);

} /* end of Write_KCF_Molecule() */





void Read_and_Append_KCF_Annotations(isdffile,osdffile)
  char *isdffile,*osdffile;
{
 FILE *fpi,*fpo;
 char line[MAX_LENGTH_LINE],status;

 /* printf("#Read_and_Append_KCF_Annotations('%s')-->'%s'\n",isdffile,osdffile); */
 fpi = fopen(isdffile,"r");
 if (fpi==NULL){printf("#ERROR:Can't open sdffile '%s'\n",isdffile); exit(1);}
 fpo = fopen(osdffile,"a");

 status = 'M';

 while (feof(fpi)==0){
   line[0] = '\0';
   fgets(line,MAX_LENGTH_LINE,fpi);
   if (status == 'A') fprintf(fpo,"%s",line); 
   if (strncmp(line,"M  END",6)==0) status = 'A';
 }

 fclose(fpi);
 fclose(fpo);
} /* end of Read_and_Append_KCF_Annotations() */
