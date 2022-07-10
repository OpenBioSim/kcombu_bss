/*
 
 <ioSDF.c> 


==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
  input/output functions for SDF files


<example of SDF file from PubChem : glycine >
---------------------------------------------------------------------
750
  -OEChem-09230802072D

 10  9  0     0  0  0  0  0  0999 V2000
    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1350    0.2500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.6720    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1350   -0.3700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1 10  1  0  0  0  0
  2  5  2  0  0  0  0
  3  4  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
  4  5  1  0  0  0  0
  4  6  1  0  0  0  0
  4  7  1  0  0  0  0
M  END
---------------------------------------------------------------------

<example of SDF from PubChem  '2-chlorophenol'>
7245
  -OEChem-01141020002D

 13 13  0     0  0  0  0  0  0999 V2000
    2.0000    0.3450    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    3.7321    1.3450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321    0.3450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -0.1550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -0.1550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.1550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -1.1550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -1.6550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1350    0.1550    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3291   -1.4650    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1350   -1.4650    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -2.2750    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690    1.6550    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0  0  0  0
  2  3  1  0  0  0  0
  2 13  1  0  0  0  0
  3  4  1  0  0  0  0
  3  5  2  0  0  0  0
  4  6  2  0  0  0  0
  5  7  1  0  0  0  0
  5  9  1  0  0  0  0
  6  8  1  0  0  0  0
  6 10  1  0  0  0  0
  7  8  2  0  0  0  0
  7 11  1  0  0  0  0
  8 12  1  0  0  0  0
M  END

<example of SDF file from PDB : 2CH, '2-chlorophenol'>
2CH
  -ISIS-            3D

 13 13  0  0  0  0  0  0  0  0  0
    0.9260    0.0170    0.2410 C   0  0  0  0  0
    1.0460    0.0180    1.6230 C   0  0  0  0  0
   -0.0860    0.0020    2.4140 C   0  0  0  0  0
   -1.3390   -0.0130    1.8300 C   0  0  0  0  0
   -1.4620   -0.0150    0.4530 C   0  0  0  0  0
   -0.3320    0.0060   -0.3430 C   0  0  0  0  0
    2.0390    0.0320   -0.5380 O   0  0  0  0  0
   -0.4870    0.0050   -2.0720 CL  0  0  0  0  0
    2.0240    0.0300    2.0800 H   0  0  0  0  0
    0.0070    0.0020    3.4900 H   0  0  0  0  0
   -2.2230   -0.0260    2.4500 H   0  0  0  0  0
   -2.4420   -0.0280   -0.0000 H   0  0  0  0  0
    2.2800   -0.8890   -0.6990 H   0  0  0  0  0
  1  2  2  0  0  0
  1  6  1  0  0  0
  1  7  1  0  0  0
  2  3  1  0  0  0
  2  9  1  0  0  0
  3  4  2  0  0  0
  3 10  1  0  0  0
  4  5  1  0  0  0
  4 11  1  0  0  0
  5  6  2  0  0  0
  5 12  1  0  0  0
  6  8  1  0  0  0
  7 13  1  0  0  0
M  END



<example of SDF file from PDB : GLY >
ALA
  -ISIS-            3D

 13 12  0  0  0  0  0  0  0  0  0
   -0.9660    0.4930    1.5000 N   0  0  0  0  0
    0.2570    0.4180    0.6920 C   0  0  0  0  0
   -0.0940    0.0170   -0.7160 C   0  0  0  0  0
   -1.0560   -0.6820   -0.9230 O   0  0  0  0  0
    1.2040   -0.6200    1.2960 C   0  0  0  0  0
:
:
    2.1130   -0.6760    0.6970 H   0  0  0  0  0
    0.4350    0.1820   -2.6470 H   0  0  0  0  0
01234567890123456789012345678901234567890123456789012345678901234567890123456789
          1         2         3         4
  1  2  1  0  0  0
  1  7  1  0  0  0
:
  6 13  1  0  0  0
M  END

--------------------------------------------------------------
<example of SDF file from PDB : FS2 >
FS2
  -ISIS-            3D

  9 14  0  0  0  0  0  0  0  0  0
   20.5910   19.5580   23.9170 FE  0  0  0  0  0
   20.7450   22.2140   24.4580 FE  0  0  0  0  0
   17.2590   18.6350   24.9940 FE  0  0  0  0  0
   19.5360   20.5940   26.8540 FE  0  0  0  0  0
   20.0350   21.3270   22.6020 S   0  0  0  0  0
   21.6780   20.5750   25.6970 S   0  0  0  0  0
   18.9650   22.3420   25.5060 O   0  0  0  0  0
   17.7170   19.9900   26.5730 O   0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0
  1  2  1  0  0  0
  1  5  1  0  0  0
  1  6  1  0  0  0
  1  9  1  0  0  0
  2  4  1  0  0  0
  2  5  1  0  0  0
  2  6  1  0  0  0
  2  7  1  0  0  0
  3  8  1  0  0  0
  3  9  1  0  0  0
  4  6  1  0  0  0
  4  7  1  0  0  0
  4  8  1  0  0  0
  4  9  1  0  0  0
M  END
*/ 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "globalvar.h"
#include "2DMAP.h" 
#include "molecule.h" 
#include "molprop.h" 
#include "ioLINE.h" 
#include "vector3.h" 

/** FUNCTIONS (GLOBAL) **/
int Translate_Lines_into_SDF_Molecule();
void Read_and_Append_SDF_Annotations();
void Write_SDF_Molecule();
int  Set_Atomic_Stereo_Parity_from_2D_BondStereo();

/** FUNCTIONS (LOCAL) **/
static int check_all_Zs_are_zero();
static void set_neighbors_by_order_of_atom_number();




int Translate_Lines_into_SDF_Molecule(HeadLineNode,mol)
 struct LINENODE *HeadLineNode; 
 struct MOLECULE *mol;
{
 struct LINENODE *xn; 
 char line[MAX_LENGTH_LINE],status,end,buff[MAX_LENGTH_LINE],a[MAX_LENGTH_LINE],b[MAX_LENGTH_LINE],natom_ok,nbond_ok; 
 int Nline,Lline,natom,nbond,anum1,anum2;
 int Nword,Nsta[10],Nend[10];
 char key[MAX_LENGTH_LINE],value[MAX_LENGTH_LINE];
 
/*
  status 'I': initial header line
  status 'A': atoms
  status 'B': bonds
  status 'C': comments
  status 'E': END
 
*/
/*
 printf("#Translate_Lines_into_SDF_Molecule(HeadLineNode,mol)\n");  fflush(stdout);
 */ 
 Nline = natom = nbond = 0; 
 mol->name[0] = mol->comment[0] = '\0';
 mol->resi[0] = '\0';
 mol->chain   = ' ';
 mol->Natom = mol->Nheavyatom = 0;
 mol->Nbond = 0; 
 mol->conmap.malloced = mol->dismap.malloced = 0; 
 mol->Nmatch = 0;
 mol->chiral_flag = '0';
 mol->Npermu = 0; 
 mol->permuhead.next = NULL;
 key[0] =  value[0] = '\0';
 status = 'I'; 
 end = 0;
 xn = HeadLineNode;
 while ((xn->next!=NULL)&&(end==0)){
   xn = xn->next;
   sprintf(line,"%s",xn->line);  
   Lline = strlen(line);
/*
   printf("'%s' Natom %d status '%c' filename '%s'\n",line,mol->Natom,status,mol->filename); fflush(stdout);
   printf("%d:'%s' name '%s' Natom %d Nbond %d\n",Nline,line,mol->name,mol->Natom,mol->Natom); 
*/ 
   if (strncmp("M  ",line,3)==0) {status = 'C';} 

   else if (strncmp("$$$",line,3)==0){ status = 'E'; end = 1;} 
   /** (1) molecular name line ***/
   else if ((Nline<3)&&(strlen(line)>0)&&(mol->name[0]=='\0')){
     Remove_Symbol_from_String(buff,line,' ');
     if ((PAR.molnamekey[0] == '\0')&&(strlen(buff)<30)){
       sprintf(mol->name,"%s",buff);
     }
     printf("#line '%s' mol->name '%s'\n",line,mol->name); fflush(stdout); 
   }
/*
 
<PubChem>
6022
 -OEChem-03071001442D

42 44  0     1  0  0  0  0  0999 V2000

<gvk_kinase>

  MOEE2007           2D

 39 41  0  0  1  0  0  0  0  0999 V2000

<LIGAND-EXPO>
ADP
  -ISIS-            3D

 42 44  0  0  0  0  0  0  0  0  0
*/
   /** (2) Natom,Nbond line ***/

   /**012345678901234567890123456789012345678901234567890 **/
   /**           1         2         3         4           **/
   /** 13 12  0  0  0  0  0  0  0  0  0 **/
   /** 10  9  0     0  0  0  0  0  0999 V2000 **/
   /**107110  0     0  0  0  0  0  0999 V2000 **/
   /*  30 31  0  0  1  0  0  0  0  0999 V2000 */
   else if ((Nline<5) && (strlen(line)>30) && (mol->Natom==0) && (mol->Nbond==0)) { 
      natom_ok = nbond_ok = 0;
      Get_Part_Of_Line(a,line,0,2);
      if (((isdigit(a[0])!=0)||(a[0]==' '))&&((isdigit(a[1])!=0)||(a[1]==' '))&&((isdigit(a[2])!=0)||(a[2]==' '))) natom_ok = 1;
      Get_Part_Of_Line(b,line,3,5);
      if (((isdigit(b[0])!=0)||(b[0]==' '))&&((isdigit(b[1])!=0)||(b[1]==' '))&&((isdigit(b[2])!=0)||(b[2]==' '))) nbond_ok = 1;
      if ((natom_ok==1)&&(nbond_ok==1)){ 
        mol->Natom = atoi(a); 
        mol->Nbond = atoi(b); 
        Malloc_MOLECULE(mol,mol->Natom,mol->Nbond);
        mol->chiral_flag = line[14];
        /* printf("#'%s' Natom %d Nbond %d\n",mol->filename,mol->Natom,mol->Nbond);  fflush(stdout); */
      } 
      /* printf("#line '%s' Natom %d Nbond %d natom_ok '%s' %d nbond_ok '%s' %d\n",line,mol->Natom,mol->Nbond,a,natom_ok,b,nbond_ok);  */
    }
   /** (3) ATOM lines ***/
   else if ((Nline>3) && (Lline>30) && (line[0]!='#') && (mol->Natom>0) && (status!='C')){ 
/*
    0.4350    0.1820   -2.6470 H   0  0  0  0  0
   20.5910   19.5580   23.9170 FE  0  0  0  0  0
    2.0000    0.3450    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0  <-- PubChem
   -0.4870    0.0050   -2.0720 CL  0  0  0  0  0                       <-- PDB LIGAND_EXPO
    2.8660   -3.5551    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.4026    1.0102    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    4.9917    1.8182    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    8.4752    2.2720    0.0000 P   0  0  3  0  0  0  0  0  0  0  0  0
   10.1980    2.4500    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
012345678901234567890123456789012345678901234567890123
          1         2         3         4         5 
*/
    mol->atoms[natom].num         = natom;
    mol->atoms[natom].num_in_file = natom + 1;
    Get_Part_Of_Line(buff,line, 0, 9); 
    mol->atoms[natom].Pos[0] = atof(buff);
    Get_Part_Of_Line(buff,line,10,19); 
    mol->atoms[natom].Pos[1] = atof(buff);
    Get_Part_Of_Line(buff,line,20,29); 
    mol->atoms[natom].Pos[2] = atof(buff);
    Get_Part_Of_Line(mol->atoms[natom].element,line,31,32); 
    if (line[32]==' ') mol->atoms[natom].element[1] = '\0';
    if (islower(line[32])!=0) mol->atoms[natom].element[1] = toupper(line[32]);

    mol->atoms[natom].mass_diff = 0;
    if (Lline>=35){
      Get_Part_Of_Line(buff,line,33,35); 
      mol->atoms[natom].mass_diff = atoi(buff);
    }

    mol->atoms[natom].char_charge   = '0';
    if (Lline>=38){
      Get_Part_Of_Line(buff,line,36,38); 
      mol->atoms[natom].char_charge   = line[38];
    }

    mol->atoms[natom].stereo_parity  = '0';
    if (Lline>=41){
      mol->atoms[natom].stereo_parity = line[41]; 
    } 
    /*   
    sprintf(mol->atoms[natom].atomname,"%s%d",mol->atoms[natom].element,mol->atoms[natom].num);
    */ 
    sprintf(mol->atoms[natom].atomname,"%s",mol->atoms[natom].element);
    sprintf(mol->atoms[natom].resi,"MOL");
    sprintf(mol->atoms[natom].rnum,"   1 ");
    mol->atoms[natom].chain = 'A';
 
    sprintf(mol->atoms[natom].mol2atomtype,"%s",mol->atoms[natom].element);
    mol->atoms[natom].Nmatch = 0;
    mol->atoms[natom].charge = 0.0;
   

    ++natom;
    if (natom>mol->Natom){ 
      printf("#ERROR(Translate_Lines_into_SDF_Molecule '%s':'%s'):natom %d is over the header-defined Natom %d.\n",
              mol->name,mol->filename,natom,mol->Natom); 
      exit(1); }
   }
   /** (4) BOND lines ***/
   else if ((Nline>3)&&(Lline<30)&&(line[0]!='#')&&(mol->Nbond>0)&&(status !='C')&&((isdigit(line[0]))||(line[0]==' '))){ 
/*
01234567890123456789012345678901234567890123456789012345678901234567890123456789
          1         2         3         4
  1  2  1  0  0  0
  1  7  1  0  0  0
  3 21  1  0  0  0  0
 18  4  1  6  0  0  0
  4 34  1  0  0  0  0
 19  5  1  6  0  0  0

*/
     mol->bonds[nbond].num = nbond;
     Get_Part_Of_Line(buff,line, 0, 2); 
     anum1 = atoi(buff)-1;
     mol->bonds[nbond].anum1 = anum1;
     Get_Part_Of_Line(buff,line, 3, 5); 
     anum2 = atoi(buff)-1;
     mol->bonds[nbond].anum2 = anum2;
     mol->bonds[nbond].bondtype   = line[8];
     mol->bonds[nbond].bondstereo = line[11];
 
    /* 
    printf("[%d] %d %d '%c'\n", 
      mol->bonds[nbond].num, mol->bonds[nbond].anum1, mol->bonds[nbond].anum2,
      mol->bonds[nbond].bondtype);
     printf("nbond %d Nbond %d\n",nbond,mol->Nbond); fflush(stdout); 
    */
     if ((anum1>=0)&&(anum1<mol->conmap.N)&&(anum2>=0)&&(anum2<mol->conmap.N)&&(nbond<mol->Nbond)){ 
 /*
       printf("line '%s'\n",line); fflush(stdout);
       printf("anum1 %d anum2 %d conmap.N %d nbond %d Nbond %d\n",
       anum1,anum2,mol->conmap.N,nbond,mol->Nbond); fflush(stdout);
 */
       mol->conmap.map[anum1][anum2] = mol->bonds[nbond].bondtype;
       mol->conmap.map[anum2][anum1] = mol->bonds[nbond].bondtype;
     }
     ++nbond; 
     if (nbond > mol->Nbond){ 
      printf("#ERROR(Translate_Lines_into_SDF_Molecule '%s':'%s'):nbond %d is over the header-defined Nbond %d.return(0);\n",
        mol->name,mol->filename,nbond,mol->Nbond); 
      return(0); 
     }
   }
   /** (5) Comment/Annotation lines ***/
/*
> <NAMIKI_ID> (NS-01548881)
NS-01548881

> <SUPPLIERNAME_1> (Vitas-M)
Vitas-M

> <SUPPLIERID_1> (STK838097)
STK838097
*/ 
  else if ( (Nline>3) && (Lline>0) && (line[0]!='#') && (status =='C')){ 
     if (line[0]=='>'){
       Split_to_Words(line,'>',&Nword,Nsta,Nend,10);
       Get_Part_Of_Line(buff,line,Nsta[0]+1,Nend[0]);
       Remove_Symbol_from_String(key,buff,'<');
     }
     else if (key[0]!='\0'){
       sprintf(value,"%s",line);
       if ((PAR.molnamekey[0] != '\0') && (strcmp(key,PAR.molnamekey)==0)){ 
         sprintf(mol->name,"%s",value);
       }
      if (PAR.READ_MOLECULAR_PROPERTY=='T') Add_key_value_string_to_KEYVALUEs(key,value,&(mol->HeadProperty));
     }
   }
   Nline += 1; 

  } /* while */

 /* printf("#name '%s' Nline %d\n",mol->name,Nline); */

  Set_one_char_ele(mol);

  

 if (PAR.Wdistan>0.0){
   Malloc_FLOAT2DMAP(&(mol->dismap),mol->Natom);
   Cal_Distance_Map_for_Heavy_atoms(mol);
 }

 if (mol->Natom==0){ 
   /*
   printf("#WARNING(Translate_Lines_into_SDF_Molecule):the file '%s' contains no atom.\n",mol->filename); 
   */ 
   Set_Molform(mol);
   sprintf(mol->molform,"-");
   return(0);
 }


 /** ending procedure **/
  if (PAR.OutCalProcess=='T')
    printf("#Translate_Lines_into_SDF_Moleclue():Natom %3d Nbond %3d '%s'\n",mol->Natom,mol->Nbond,mol->molform);

  /*
  Show_KEYVALUEs(&(mol->HeadProperty));
  */

  return(mol->Natom);

} /* end of Translate_Lines_into_SDF_Molecule() */







void Write_SDF_Molecule(ofname,mol,comment)
 char *ofname;
 struct MOLECULE *mol;
 char   *comment;
{
 FILE *fpo;
 int i,j;
 struct KEYVALUE *kn;

 printf("#Write_SDF_Molecule()-->'%s'\n",ofname);
 if (ofname[0]=='-') fpo = stdout;
 else{ 
   fpo = fopen(ofname,"w");
   if (fpo==NULL){
     printf("#ERROR:Can't write to '%s'\n",ofname);
     exit(1);
   }
 }
 printf("#mol->name len %d '%s'\n",strlen(mol->name), mol->name);
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
   fprintf(fpo,"%10.4f%10.4f%10.4f %-2s %2d  %c  %c  0  0  0  0  0  0  0  0  0\n",
     mol->atoms[i].Pos[0],mol->atoms[i].Pos[1],mol->atoms[i].Pos[2],mol->atoms[i].element,
     mol->atoms[i].mass_diff,mol->atoms[i].char_charge,mol->atoms[i].stereo_parity); 
 }
/*  1  6  1  0  0  0  0 */
 for (j=0;j<mol->Nbond;++j){
  fprintf(fpo,"%3d%3d  %c  %c  0  0  0\n",
        mol->bonds[j].anum1+1,mol->bonds[j].anum2+1,mol->bonds[j].bondtype,mol->bonds[j].bondstereo);
 }

 /*
 fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
 if (strlen(comment)>0){ 
   fprintf(fpo,"#%s\n",comment);
 }
 */
 fprintf(fpo,"M  END\n");
 kn = &(mol->HeadProperty);
 while (kn->next != NULL){ 
   kn = kn->next;
   fprintf(fpo,"> <%s>\n",kn->key);
   fprintf(fpo,"%s\n",kn->value);
   fprintf(fpo,"\n");
 } 

 fprintf(fpo,"$$$$\n");
 if (ofname[0]!='-') {fclose(fpo);}

} /* end of Write_SDF_Molecule() */




void Read_and_Append_SDF_Annotations(isdffile,osdffile)
  char *isdffile,*osdffile;
{
 FILE *fpi,*fpo;
 char line[MAX_LENGTH_LINE],status;

 /* printf("#Read_and_Append_SDF_Annotations('%s')-->'%s'\n",isdffile,osdffile); */
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
} /* end of Read_and_Append_SDF_Annotations() */








int check_all_Zs_are_zero(mol)
  struct MOLECULE *mol;
{
  int i;
  float eps;

  eps = 0.0001; 
  for (i=0;i<mol->Natom;++i){
    if (fabs(mol->atoms[i].Pos[2])>eps){
       return(0);
    }
  }
  return(1);
} /* end of check_all_Zs_are_zero() */


void set_neighbors_by_order_of_atom_number(Nneighbor,neighbors)
  int Nneighbor;
  int *neighbors;
{
  int Nchange,i,buff;
  do{
    Nchange = 0;
    for (i=1;i<Nneighbor;++i){
      if (neighbors[i-1]>neighbors[i]){
        buff = neighbors[i-1];
        neighbors[i-1] = neighbors[i];
        neighbors[i] = buff;
        Nchange += 1;
      }
    }
  } while (Nchange>0);   

} /* end of set_neighbors_by_order_of_atom_number() */



int Set_Atomic_Stereo_Parity_from_2D_BondStereo(mol,SetUpStrategy)
  struct MOLECULE *mol;
  char   SetUpStrategy; 
  /*
  'C':hecking. not reassign stereo parity.
  'P':rimary use for stereo_parity_byBond, even if original stereo_parity is different.
  'S':econdary use for stereo_parity_byBond. Replace withstereo_parity_byBond 
     only if original stereo_parity = '0' or '3'.
  'B':ond only. Ignore atomic stereo parity.
  */ 
{
 int b,s,e,i,k,ni,Nchiral,Nagree, Ndisagree, Nonlybond, Nonlyatom,  Nbondchiral,Nup,Ndown;
 float neipos[4][3];   /* 3D position of neighboring atoms neipos[0,1,2,3][x,y,z] */
 float blen,nvec[3],v01[3],v12[3],vs3[3];
 char  stereo_parity_byBond;
 /*
 
  char bondstereo;  '0' not stereo (default), '1':Up, '4':Either, '6'Down 

    s =====> e  : Up  ('1'), 'e' has positive Z coordinate
    s - - -> e  : Down('6'), 'e' has negative Z coordinate

     0
     |
  3--s--1 
     |
     2

  (01 x 12) dprod s3 >=0, stereo_parity_of_s = 1
        otherwise       , stereo_parity_of_s = 2

  if (mol->atoms[s].Nneighbor==3):
    make dummy 3 = 4s - (0+1+2) 

 */

/*
  printf("#Set_Atomic_Stereo_Parity__from_2D_BondStereo(filename '%s' Natom %d SetUpStrategy:'%c')\n",mol->filename,mol->Natom,SetUpStrategy);
*/
 if  (check_all_Zs_are_zero(mol)==0){
  /*
   printf("#WARNING: Set_Stereo_Parity_of_Molecule_from_2D_BondStereo() cannot be performed, because some atoms has non-zero Z.\n");
  */
   return(0); 
 }

  Nchiral = Nagree = Ndisagree = Nonlybond = Nonlyatom = 0;
  if (SetUpStrategy=='C'){
    printf(">%s\n",mol->filename);
  }

 /** [1] mark=='0' for BondStereo annotated atoms **/
  for (i=0;i<mol->Natom;++i){ mol->atoms[i].mark = ' ';}

  for (b=0;b<mol->Nbond;++b){
    if (   ((mol->bonds[b].bondstereo=='1') || (mol->bonds[b].bondstereo=='6'))
         &&((mol->atoms[mol->bonds[b].anum1].Nneighbor==3)||(mol->atoms[mol->bonds[b].anum1].Nneighbor==4)) ){
      s = mol->bonds[b].anum1;  /* start */
      mol->atoms[s].mark = '0';}
  }

 /** [2] scan marked atoms as s **/ 
  for (s=0;s<mol->Natom;++s){
    if (mol->atoms[s].mark == '0'){
      Nbondchiral = Nup = Ndown = 0;
      /* (2-1) setup neipos[4][3] */
      for (i=0;i<mol->atoms[s].Nneighbor;++i){
        ni = mol->atoms[s].neighbors[i];
        for (k=0;k<3;++k){ neipos[i][k] = mol->atoms[ni].Pos[k];}
      }

      /* (2-2) setup Z coordinates of e in neipos[4][3] */
      for (b=0;b<mol->Nbond;++b){
        if ( (mol->bonds[b].anum1==s) &&  ((mol->bonds[b].bondstereo=='1') || (mol->bonds[b].bondstereo=='6'))){
          Nbondchiral += 1;
          e = mol->bonds[b].anum2; 
          blen = distance_vec3(mol->atoms[s].Pos, mol->atoms[e].Pos); 
          for (i=0;i<mol->atoms[s].Nneighbor;++i){
            ni = mol->atoms[s].neighbors[i];
            if (e==ni){
                if (mol->bonds[b].bondstereo=='1'){ neipos[i][2] =  blen;}
                if (mol->bonds[b].bondstereo=='6'){ neipos[i][2] = -blen;}

                if ((mol->atoms[s].Nneighbor ==3)&&((Nup>0)||(Ndown>0))){
                   neipos[i][2] = 0.0;
                   /*
                    >> EXCEPTIONAL CASE !! << 
                    If the atom has both Nup and Ndown atoms with Nneighbor=3,
                    then Z-coordinate of dummy is always zero; then, the Stereo parity cannot be decided.
                    See the KEGG_DRUG D00926, the atom C59.
                    */ 
                }
                if (mol->bonds[b].bondstereo=='1'){ Nup   += 1;}
                if (mol->bonds[b].bondstereo=='6'){ Ndown += 1;}
            }
          }
        }
      }
 
       /** (2-3) make dummy vec[3],if Nneighbor == 3 **/      
      if (mol->atoms[s].Nneighbor==3){
        for (k=0;k<3;++k) {
          neipos[3][k] = 4*mol->atoms[s].Pos[k] - (neipos[0][k] + neipos[1][k] + neipos[2][k]);
        }
      }


/*
      if (Nbondchiral > 1){
       printf("#NBONDCHIRAL %s %d Nup %d Ndown %d %s%-2d stereo by_atom %c  Nneighbor %d\n",
         mol->filename,Nbondchiral, Nup, Ndown,mol->atoms[s].element, mol->atoms[s].num_in_file,
         mol->atoms[s].stereo_parity, mol->atoms[s].Nneighbor);
      }
*/
      /** (2-4) Set up stereo_parity_byBond **/      
      stereo_parity_byBond = '0';
      sub_vec3(v01,neipos[1],neipos[0]);  normalize_vec3(v01);
      sub_vec3(v12,neipos[2],neipos[1]);  normalize_vec3(v12);
      cross_vec3(nvec,v01,v12);    normalize_vec3(nvec);
      sub_vec3(vs3,neipos[3],mol->atoms[s].Pos);  normalize_vec3(vs3);
      if (dot_prod3(nvec,vs3)>=0.0){ stereo_parity_byBond = '1';}
                             else  { stereo_parity_byBond = '2';}
  
      mol->atoms[s].mark = stereo_parity_byBond;

 /*
      printf("#%s %d Nneighbor %d dprod %lf stereo_parity_byBond %c\n",
        mol->atoms[s].atomname, mol->atoms[s].num_in_file,mol->atoms[s].Nneighbor,
        dot_prod3(nvec,vs3),stereo_parity_byBond); 
      for (i=0;i<mol->atoms[s].Nneighbor;++i){
            ni = mol->atoms[s].neighbors[i];
        printf("#neighbor[%d] %s %d %f %f %f \n",i,mol->atoms[ni].atomname, mol->atoms[ni].num_in_file
            ,neipos[i][0], neipos[i][1], neipos[i][2]);
      }
*/


/*
       printf("#STEREO_PARITY %s%-2d stereo by_atom %c  by_bond %c Nneighbor %d %f %f %f\n",
         mol->atoms[s].element, mol->atoms[s].num_in_file,
         mol->atoms[s].stereo_parity, stereo_parity_byBond,mol->atoms[s].Nneighbor, 
         mol->atoms[s].Pos[0], mol->atoms[s].Pos[1],mol->atoms[s].Pos[2]);
       printf("HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f\n",mol->atoms[s].num_in_file,mol->atoms[s].element,"CEN",' '," 1",
          mol->atoms[s].Pos[0], mol->atoms[s].Pos[1],mol->atoms[s].Pos[2]);
 
      for (i=0;i<mol->atoms[s].Nneighbor;++i){
        ni = mol->atoms[s].neighbors[i];
        printf("HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f\n",mol->atoms[ni].num_in_file,mol->atoms[ni].element,"NEI",' '," 1",
                neipos[i][0], neipos[i][1], neipos[i][2]);
      }
      if (mol->atoms[s].Nneighbor==3){
        printf("HETATM%5d %4s %3s %c%5s   %8.3f%8.3f%8.3f\n",4,"H","DUM",' '," 1",neipos[3][0], neipos[3][1], neipos[3][2]);
      } 
*/ 
      /** (2-5) Set up stereo_prity by stereo_parity_byBond  **/      
      if (SetUpStrategy=='P'){
         if (stereo_parity_byBond != '0') {mol->atoms[s].stereo_parity = stereo_parity_byBond;}
         /*
         printf("#%s %d stereo_byBond %c\n",mol->atoms[s].atomname, mol->atoms[s].num_in_file, stereo_parity_byBond); 
         */
      }
  
      else if (SetUpStrategy=='S'){
         if ( ((mol->atoms[s].stereo_parity=='0')|| (mol->atoms[s].stereo_parity=='3'))&&
                 (stereo_parity_byBond != '0') ){  mol->atoms[s].stereo_parity = stereo_parity_byBond;}
      }
  
      else if (SetUpStrategy=='B'){
         mol->atoms[s].stereo_parity = stereo_parity_byBond;
      }

    } 
  }  

   /* [3] reset bond stereo */
   /*
    if ((SetUpStrategy=='P')||  (SetUpStrategy=='S')|| (SetUpStrategy=='B')){
      for (b=0;b<mol->Nbond;++b){
         mol->bonds[b].bondstereo = '0';
      } 
   }
   */

   /* [4] Show agree and disagree bwn atom stereo and bond sterero  */
  
  if (SetUpStrategy=='C'){
    for (i=0;i<mol->Natom;++i){
      if( (mol->atoms[i].stereo_parity=='1')||(mol->atoms[i].stereo_parity=='2')||(mol->atoms[i].mark!=' ')){
        printf("#STEREO_PARITY %s%-2d stereo by_atom %c by_bond %c Nneighbor %d ",
          mol->atoms[i].element, mol->atoms[i].num_in_file,
          mol->atoms[i].stereo_parity, mol->atoms[i].mark,mol->atoms[i].Nneighbor);
        Nchiral += 1; 

       if (mol->atoms[i].stereo_parity != mol->atoms[i].mark){ 

          if (((mol->atoms[i].stereo_parity=='1')||(mol->atoms[i].stereo_parity=='2')) &&
               (mol->atoms[i].mark != ' ')){ 
                printf(" DISAGREE"); 
               Ndisagree += 1; 
             }

          if ( ((mol->atoms[i].stereo_parity=='0')|| (mol->atoms[i].stereo_parity=='3'))  &&
               (mol->atoms[i].mark != ' ')){ 
               printf(" ONLY_BOND");
               Nonlybond += 1; 
             }
          if ( ((mol->atoms[i].stereo_parity=='1')|| (mol->atoms[i].stereo_parity=='2'))  &&
               (mol->atoms[i].mark == ' ')){ 
               printf(" ONLY_ATOM");
               Nonlyatom += 1; 
            }
        }
  
       if ((mol->atoms[i].stereo_parity == mol->atoms[i].mark)&&(mol->atoms[i].mark !=' ')){ 
         printf(" AGREE");
         Nagree += 1; 
       }
       printf("\n");
     }
   }
   
   printf("#NCHIRAL_ATOM %d\n",Nchiral); 
   if (Nonlyatom>0){ printf("#NCHIRAL_ATOM_ONLY_ATOM %d\n",Nonlyatom); }
   if (Nonlybond>0){ printf("#NCHIRAL_ATOM_ONLY_BOND %d\n",Nonlybond); }
   if (Nagree>0){ printf("#NCHIRAL_ATOM_AGREE_ATOM_BOND %d\n",Nagree); }
   if (Ndisagree>0){ printf("#NCHIRAL_ATOM_DISAGREE_ATOM_BOND %d\n",Ndisagree); }
   
 }

  return(1);

} /* end of Set_Atomic_Stereo_Parity_from_2D_BondStereo() */


