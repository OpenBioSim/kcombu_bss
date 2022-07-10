/*
 
 <ioLib.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <dirent.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "ioLINE.h"
#include "ioLib.h"
#include "ioSDF.h"
#include "ioPDB.h"
#include "ioMOL2.h"
#include "ioKCF.h"
#include "match.h"
#include "molprop.h"
#include "options.h"
#include "molring.h"
#include "moltopodis.h"
#include "ringblock.h"
#include "ioSMILES.h"
#include "stereo_check.h"


/*** FUNCTIONS (GLOBAL) ***/
int Make_LibFileList_from_Colon_Separated_String();
int Get_List_of_LibFiles_Under_Dir();
int Get_List_of_LIBMOLs_from_one_library_file();
int Read_One_Molecule_from_one_Library_file();
int Read_One_Entry_as_LINENODEs();
int Read_List_of_LIBMOLs();
int Get_List_of_LIBMOLs_Under_Dir();
void Copy_MATCH_Info_to_LIBMOL();
void Write_Search_Result();
int  Write_Similar_Molecules_To_One_File();
int  Write_Similar_Molecules_Under_Dir();
void Write_Clustering_Res_List_Of_NODEs();
void Cal_Nsta_Nend_from_Bunshi_Bunbo_for_all_vs_all();

void Write_Similarities_in_List_format_Header();
void Write_Two_Groups_Similarities_in_List_format_Header();
void Write_Similarities_in_List_format_Tail();
int Number_of_LIBMOL();
struct LIBMOL *Get_LIBMOL_by_number();
struct LIBMOL* Add_LIBMOL_to_LIBMOLlist();
void Increment_Nmatch_for_Matched_Atoms();
void Write_Nmatch_rasmol_Script();
void insert_LIBMOL_into_score_decreasing_proper_position();
void delete_the_first_LIBMOL_node();
void reverse_LIBMOL_list();
void Show_LIST_of_LIBMOL();



/*** FUNCTIONS (LOCAL) ***/




int Make_LibFileList_from_Colon_Separated_String(HeadLibFile,string)
  struct LINENODE *HeadLibFile; 
  char *string; /* colon-separated string, such as '0000.sdf:0001:sdf:0002.sdf' */
{
  int Nsta[100],Nend[100],Nword;
  int i,N;
  char space2[100];

  HeadLibFile->next = NULL;
  N = 0;
  Split_to_Words(string,':',&Nword,Nsta,Nend,100);
  for (i=0;i<Nword;++i){ 
    Get_Part_Of_Line(space2,string,Nsta[i],Nend[i]); 
    if ((Nend[i]-Nsta[i])>1){
      Add_string_to_LINENODEs(space2,HeadLibFile);
      N += 1;
    }
  } 

 return(N);

} /* end of Make_LibFileList_from_Colon_Separated_String() */




int Get_List_of_LibFiles_Under_Dir(idirmlib,HeadLibFile,ftailstr)
 char   *idirmlib;    /* input directory for library molecules */
 struct  LINENODE *HeadLibFile;
 char   *ftailstr;   /* file ftail string. If ftailstr[0]=='\0', accept everything. */
{
  DIR *dp0, *dp1;
  struct dirent *dir0;
  char file[1024];
  int Nlist,len0,Lftailstr;

  Lftailstr = strlen(ftailstr);  
  Nlist = 0;
  HeadLibFile->next = NULL; 
  HeadLibFile->num  = -1; 
  printf("#Get_List_of_LibFiles_Under_Dir('%s', ftailstr '%s')\n",idirmlib,ftailstr);

 
  dp0 = opendir(idirmlib);
  if (dp0==NULL){printf("#ERROR:directory '%s' is not found.\n",idirmlib); exit(1);} 

  /* check under the given directory (idirlib) */
  while((dir0=readdir(dp0)) != NULL){
    len0 = strlen(dir0->d_name);
    if ((len0>0)&&(strcmp(dir0->d_name,".")!=0) && (strcmp(dir0->d_name,"..")!=0)){
        sprintf(file,"%s/%s",idirmlib,dir0->d_name); 
        dp1 = opendir(file);
        if (dp1==NULL){
          Add_string_to_LINENODEs(file,HeadLibFile);
          Nlist += 1;
        }
    }
  } /* dir0 */

  closedir(dp0); 
  printf("#Nlist %d\n",Nlist); 
  return(Nlist);
} /* end of Get_List_of_LibFiles_Under_Dir() */









int Get_List_of_LIBMOLs_from_one_library_file(ifname,HeadLib,filetype,num_file)
 char   *ifname;   /* input file name */ 
 struct LIBMOL *HeadLib;
 char   *filetype; /* 'S':DF, 'P':DB  (to be assigned) */
 int    num_file;  /* number of file for the filename ifname[]. The default should be '0'. */
{
 FILE   *fp;
 char   s[256],name[128],buff[128];
 int    len,Nlibmol,Nline, nline_entry,end,minNline; 
 struct LIBMOL *nn;
 long offset0;
 DIR  *dir;

 minNline = 5;

 if (*filetype=='-'){
    *filetype =  molecular_file_type_from_tail_of_file(ifname);
  }

 if ((*filetype!='S')&&(*filetype != 'P')&&(*filetype != '2')){
   printf("#ERROR:filetype '%c' is not supported.\n",*filetype);
   exit(1);
 }

 printf("#Get_List_of_LIBMOLs_from_one_library_file('%s' filetype '%c')\n",ifname,*filetype);

 fp = fopen(ifname,"r"); 
 if (fp==NULL){printf("#ERROR:Can't open library file '%s'\n",ifname); exit(1);}
 dir = opendir(ifname);
 if (dir!=NULL){printf("#ERROR:the library file '%s' is a directory.\n",ifname); exit(1);}
 closedir(dir); 

 nn = HeadLib;
 while (nn->next != NULL) nn = nn->next;

 Nlibmol =  Nline = nline_entry = 0;
 end = 0;
 offset0 = ftell(fp);
 name[0] = '\0';

 while (feof(fp)==0){
   s[0] = '\0';
   fgets(s,256,fp);
   len = strlen(s);
   if (feof(fp)!=0) end = 1;
/*
   printf("before %d '%s'\n",len,s); 
*/
   if (s[len-1]=='\n') {s[len-1] = '\0'; len -= 1;}
   if (s[len-1]=='\r') {s[len-1] = '\0'; len -= 1;}
   if (s[len-1]=='\n') {s[len-1] = '\0'; len -= 1;}
   if (s[len-1]=='\r') {s[len-1] = '\0'; len -= 1;}
/*
   printf("after  %d '%s'\n",len,s); 
   printf("'%s' Nlibmol %d\n",s,Nlibmol); 
*/ 
   ++Nline;
   ++nline_entry;

   /*** for SDF file ***/
   if (*filetype=='S'){
     if ((nline_entry<3)&&(strlen(s)>0)&&(name[0]=='\0')){
       Remove_Symbol_from_String(buff,s,' ');
       if (strlen(buff)<30) sprintf(name,"%s",buff); 
     }

     if ((s[0]=='$')&&(s[1]=='$')&&(s[2]=='$')&&(s[3]=='$')){ end = 1; }
   }
   /*** for PDB file ***/
   if (*filetype=='P'){
       /* if ((s[0]=='T')&&(s[1]=='E')&&(s[2]=='R')) {end = 1;} */
       if ((s[0]=='E')&&(s[1]=='N')&&(s[2]=='D')) {end = 1;}
   }

   if ((*filetype=='2') && (strncmp(s,"@<TRIPOS>MOLECULE",17)==0)) { end = 1;}

   /* printf("%s:%d:%d:end %d\n",s,Nline,nline_entry,end); */
  
   if (end==1){
     if (nline_entry>=minNline){
       nn->next = (struct LIBMOL*)malloc(sizeof(struct LIBMOL)); 
       nn->next->next = NULL;
       nn->next->prev = nn;
       nn = nn->next;
       nn->num = nn->prev->num + 1;
       nn->name = (char *)malloc(sizeof(char)*(128));
       sprintf(nn->name,"%s",name);
       nn->num_file = num_file;  
       nn->class = (char *)malloc(sizeof(char)*(128));
       nn->class[0] = '\0';
       nn->molform = NULL;
       nn->offset_libfile = offset0;
       /* printf("#offset0 %d\n",(int)offset0); */
       Nlibmol += 1; 
       name[0] = '\0';
    }
    nline_entry = 0; 
    end = 0;
    offset0 = ftell(fp);
   }
 
 }
 fclose(fp);

 printf(" --> %d molecules\n",Nlibmol);
 return(Nlibmol);

} /* end of Get_List_of_LIBMOLs_from_one_library_file() */







int Read_One_Molecule_from_one_Library_file(fp,mol,FileType,AtmHetType,BondType)
  FILE    *fp; 
  struct MOLECULE *mol;
  char   FileType;   /* '-' not assigned. 'S':sdf, 'P':pdb  */
  char   AtmHetType; /* 'A'tom, 'H'etatm, 'B'oth (only for PDB) */
  char   BondType;    /* 'B': guess bonds by distances and cal topo-distances */

{
  int ok,Nline;
  struct LINENODE HeadLineNode;

/* 
  printf("#Read_One_Molecule_from_one_Library_file(fp,mol,FileType:%c AtmHetType:%c,BondType:%c)\n",FileType,AtmHetType,BondType);
 */

  ok = 0;
  Nline = Read_One_Entry_as_LINENODEs(fp,&HeadLineNode,FileType);

  /* printf("#Nline %d\n",Nline); */ 
  if (Nline>0){

  /*
   printf("--begin--\n");
    ln = &HeadLineNode;
    while (ln->next != NULL){
     ln = ln->next;
     printf("%d:'%s'\n",ln->num,ln->line);
    }
    printf("--end--\n");
  */

    mol->filetype = FileType;
    mol->filename[0] = '\0';
         if (mol->filetype=='P'){ok = Translate_Lines_into_PDB_Molecule(&HeadLineNode,mol,AtmHetType,BondType);}
    else if (mol->filetype=='S'){ok = Translate_Lines_into_SDF_Molecule(&HeadLineNode,mol); }
    else if (mol->filetype=='2'){ok = Translate_Lines_into_MOL2_Molecule(&HeadLineNode,mol); }
    else if (mol->filetype=='K'){ok = Translate_Lines_into_KCF_Molecule(&HeadLineNode,mol); }
  
    if (ok==0){
      Free_LINENODEs(&HeadLineNode);
      return(0);
     }
 

   Set_Center_of_Gravity(mol);
 
   if (mol->filetype != 'P'){Set_Unique_atomname_to_Molecule(mol);}

   Set_radius(mol);
   Set_Nneighbor_CHNOSP(mol);
   Set_Extended_Connectivity(mol);
   Set_Molform(mol);

   if (BondType == 'B'){
     Set_Topological_Distance(mol,'-');
     Set_Ring_Block_and_SSSR_Ring(mol);
   }

   Set_atomtype_of_Atoms(mol,PAR.atomtype_class);
   Count_atomtype_in_Natomtype(mol);
   Set_Nneighbor_atomtype(mol);

        if (PAR.bondtype_class == 'A'){ Set_Aromatic_Bond_Atoms(mol);}
   else if (PAR.bondtype_class == 'R'){ Set_Rotatable_Bond(mol);}

   if ((PAR.AtomStereoFromBondStereo!='F') && (mol->filetype=='S')){
     Set_Atomic_Stereo_Parity_from_2D_BondStereo(mol,PAR.AtomStereoFromBondStereo);
   }

   if (mol->SetStereo3D == 'T'){ Set_stereo_parity_of_3D_molecule(mol); }

   if (PAR.Wrank > 0.0){
     Set_Aromatic_Bond_Atoms(mol);
     Set_Connected_Structure(mol);
     Set_Unique_Extended_Connectivity(mol,-1,0);
   }


 
   sprintf(mol->core_molname,"%s",mol->name);

  }

  Free_LINENODEs(&HeadLineNode);
  /* printf("#Natom %d Nbond %d\n",mol->Natom,mol->Nbond); */
  return(ok);

} /* end of Read_One_Molecule_from_one_Library_file() */





int Read_One_Entry_as_LINENODEs(fp,HeadLineNode,filetype)
 FILE *fp; 
 struct LINENODE *HeadLineNode;
 char filetype; /* 'S':DF, 'P':DB '2':mol2 */
{
 char s[150];
 struct LINENODE *ln;
 int L,nline,end_entry,minNline,Nmol2_comment_obs,Nmol2_molecule_obs;
 char include_this_line;
 long file_offset;
/*
<Split string for SDF file>
M  END
$$$$

<Split string for PDB file>
TER
END
<Split string for MOL2 file>
@<TRIPOS>MOLECULE
*/
   if ((filetype!='S')&&(filetype!='P')&&(filetype!='2')){ 
     printf("#ERROR:improper molecular filetype '%c'. Assign 'S'df or 'P'db' or '2':mol2.\n",filetype); 
     exit(1);
   }

   minNline = 5;
   /* printf("#Read_One_Entry_as_LINENODEs(filetype '%c')\n",filetype);  */
   HeadLineNode->next = NULL;
   HeadLineNode->prev = NULL;
   ln = HeadLineNode;
   Nmol2_comment_obs = Nmol2_molecule_obs = 0;

   end_entry = 0;
   nline = 0;

   while ((feof(fp)==0)){
     s[0] = '\0';
     file_offset = ftell(fp);
     fgets(s,149,fp);
     
     L = strlen(s);
     if (s[L-1]=='\n') {s[L-1] = '\0'; L = strlen(s);}
     if (s[L-1]=='\r') {s[L-1] = '\0'; L = strlen(s);}
     if (s[L-1]=='\n') {s[L-1] = '\0'; L = strlen(s);}
     if (s[L-1]=='\r') {s[L-1] = '\0'; L = strlen(s);}

/*
     printf("##file_offset %d %d:%s:end_entry %d\n",file_offset,nline,s,end_entry);  
*/
     include_this_line = 1;

     if (filetype=='2'){
       if (strncmp(s,"@<TRIPOS>COMMENT",16)==0){ Nmol2_comment_obs += 1;}
       if (Nmol2_comment_obs > 1) { 
         end_entry = 1; 
         fseek(fp,file_offset,SEEK_SET); 
         include_this_line = 0;
       }

       if (strncmp(s,"@<TRIPOS>MOLECULE",17)==0){ Nmol2_molecule_obs += 1;}
       if (Nmol2_molecule_obs > 1) { 
         end_entry = 1; 
         fseek(fp,file_offset,SEEK_SET); 
         include_this_line = 0;
       }
     }

     if (include_this_line==1){
       ln->next = (struct LINENODE *)malloc(sizeof(struct LINENODE));
       ln->next->prev = ln;
       ln->next->next = NULL; 
       ln = ln->next;
       /* printf("s '%s'\n",s); fflush(stdout); */
       ln->line = (char *)malloc(sizeof(char)*(strlen(s)+1));
       sprintf(ln->line,"%s",s);
       ln->num = nline;
       ++nline;
     }
/*
     if ((filetype=='S') &&
       (((s[0]=='M')&&(s[1]==' ')&&(s[2]==' ')&&(s[3]=='E')&&(s[4]=='N')&&(s[5]=='D')) ||
        ((s[0]=='$')&&(s[1]=='$')&&(s[2]=='$')&&(s[3]=='$')) ) ){ end_entry = 1;}
*/
     if ((filetype=='S') && (s[0]=='$')&&(s[1]=='$')&&(s[2]=='$')&&(s[3]=='$')){ end_entry = 1;}

     if ((filetype=='P') &&
       (((s[0]=='T')&&(s[1]=='E')&&(s[2]=='R')) ||
        ((s[0]=='E')&&(s[1]=='N')&&(s[2]=='D')))   ){ end_entry = 1;}
 
     if (feof(fp)==1) end_entry = 1;
 
     if (end_entry==1){
       if (nline>=minNline){ return(nline);} 
       else {
         Free_LINENODEs(HeadLineNode);
         HeadLineNode->next = NULL;
         HeadLineNode->prev = NULL;
         ln = HeadLineNode;
         end_entry = 0;
       }
      nline = 0;  
     }
 
   } /* while */
 
 return(0);

} /* end of Read_One_Entry_as_LINENODEs() */







int Read_List_of_LIBMOLs(fname,HeadList)
 char   *fname;
 struct LIBMOL *HeadList;
{
 FILE *fp;
 char line[1024];
 int  Nword,Wsta[100],Wend[100];
 struct LIBMOL *nn;
 int Nlist;

 /**
 >>> FILE FORMAT <<
 [name0] [class0] 
 [name1] [class1] 
 [name2] [class2] 
 [name3] [class3] 
 :       : 
 **/
 HeadList->next = NULL; 
 nn = HeadList;
 
 fp = fopen(fname,"r"); 
 if (fp==NULL){printf("#ERROR:Can't open listfile '%s'\n",fname); exit(1);}
 Nlist = 0;
 while (feof(fp)==0){
   line[0] = '\0';
   fgets(line,1024,fp);
   if (line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';
   if ((line[0]!='#') && (strlen(line)>0)){

     Split_to_Words(line,' ',&Nword,Wsta,Wend,100);
     /* printf("'%s' -> Nword %d\n",line,Nword);  */
     
     nn->next = (struct LIBMOL*)malloc(sizeof(struct LIBMOL)); 
     nn->next->next = NULL;
     nn->next->prev = nn;
     nn = nn->next;
     nn->num = Nlist;
     nn->offset_libfile = 0; 
     nn->name = (char *)malloc(sizeof(char)*(Wend[0]-Wsta[0]+2));
     Get_Part_Of_Line(nn->name,line,Wsta[0],Wend[0]);
     if (Nword>1){
       nn->class = (char *)malloc(sizeof(char)*(strlen(line)-Wsta[1]+2));
       Get_Part_Of_Line(nn->class,line,Wsta[1],strlen(line)-1);
     }
     else {
       nn->class = (char *)malloc(sizeof(char));
       nn->class[0] = '\0';
     } 
     Nlist += 1;
     /* if ((Nlist%100)==0) printf("#%d '%s' '%s'\n",Nlist,nn->name,nn->class);   */
   } 
 }
 fclose(fp);
 printf("#Read_List_of_NODEs('%s') : Nlist %d\n",fname,Nlist);
 return(Nlist);
} /* end of Read_List_of_LIBMOLs() */







int Get_List_of_LIBMOLs_Under_Dir(idirlib,HeadList,fheadstr)
 char   *idirlib;    /* input directory for library molecules */
 struct LIBMOL *HeadList;
 char   *fheadstr;   /* file head string. If fheadstr[0]=='\0', accept everything. */
{
  DIR *dp0, *dp1, *dp2;
  struct dirent *dir0,*dir1;
  char file[1024];
  struct LIBMOL *nn;
  int Nlist,len0,len1,Lfheadstr;
  int Nfile_under_dir;

  Lfheadstr = strlen(fheadstr);  
  Nlist = 0;
  HeadList->next = NULL; 
  nn = HeadList;
  printf("#Get_List_of_LIBMOLs_Under_Dir('%s', fheadstr '%s')\n",idirlib,fheadstr);

 
  dp0 = opendir(idirlib);
  if (dp0==NULL){printf("#ERROR:directory '%s' is not found.\n",idirlib); exit(1);} 

  /* check under the given directory (idirlib) */
  while((dir0=readdir(dp0)) != NULL){
    len0 = strlen(dir0->d_name);
    if ((len0>0)&&(strcmp(dir0->d_name,".")!=0) && (strcmp(dir0->d_name,"..")!=0)){
      if ((fheadstr[0]=='\0') || (strncmp(dir0->d_name,fheadstr,Lfheadstr)==0)){
        sprintf(file,"%s/%s",idirlib,dir0->d_name); 
        dp1 = opendir(file);
        if (dp1==NULL){
          /* (1) make nodes for the file under the given directory (idirlib) */
          nn->next = (struct LIBMOL*)malloc(sizeof(struct LIBMOL)); 
          nn->next->next = NULL;
          nn->next->prev = nn;
          nn = nn->next;
          nn->num = Nlist;
          nn->name = (char *)malloc(sizeof(char)*(len0+1));
          sprintf(nn->name,"%s",dir0->d_name);
          /* printf("nn->name '%s' len0 %d\n",nn->name,len0);  */
          nn->class = (char *)malloc(sizeof(char)*1);
          nn->class = NULL;
          nn->tanimoto_mcs = nn->tanimoto_oneatom =  nn->tanimoto_atompair =  0.0; 
          nn->num_file = 0;
          nn->offset_libfile =  nn->offset_descfile = 0;
          /* printf("#1st_layer:%s\n",nn->name);  */
          Nlist += 1;
        }
        else{
         /* check under the directory under the given directory  (idirlib) */
          printf("#go into the directory '%s' ",file);
          Nfile_under_dir = 0;
          while((dir1=readdir(dp1)) != NULL){
            len1 = strlen(dir1->d_name);
            if ((len1>0)&&(strcmp(dir1->d_name,".")!=0) && (strcmp(dir1->d_name,"..")!=0)){
              sprintf(file,"%s/%s/%s",idirlib,dir0->d_name,dir1->d_name); 
              dp2 = opendir(file);
              if (dp2==NULL){
                /* (2) make nodes for the file under the directory under the given directory (idirlib) */
                nn->next = (struct LIBMOL*)malloc(sizeof(struct LIBMOL)); 
                nn->next->next = NULL;
                nn->next->prev = nn;
                nn = nn->next;
                nn->num = Nlist;
                nn->name = (char *)malloc(sizeof(char)*(len0+len1+2));
                sprintf(nn->name,"%s/%s",dir0->d_name,dir1->d_name);
                nn->class = (char *)malloc(sizeof(char)*1);
                nn->class = NULL;
                nn->tanimoto_mcs = nn->tanimoto_oneatom =  nn->tanimoto_atompair =  0.0; 
                nn->num_file = 0;
                nn->offset_libfile =  nn->offset_descfile = 0;
                /* printf("2nd_layer:%s\n",nn->name); */
                Nfile_under_dir += 1;
                Nlist += 1;
              }
              else { closedir(dp2); }
            } 
          } /* dir1 */
          printf(" Nfile_under_dir %d\n",Nfile_under_dir);
          closedir(dp1); 
        }
      }
    } 
  } /* dir0 */

  closedir(dp0); 
  printf("#Nlist %d\n",Nlist); 
  return(Nlist);
} /* end of Get_List_of_LIBMOLs_Under_Dir() */






void Copy_MATCH_Info_to_LIBMOL(ln,mn,molQ,molL)
  struct LIBMOL *ln;
  struct MATCH   *mn;
  struct MOLECULE *molQ, *molL;
{
  double log128;
  int i;

  ln->tanimoto_mcs   = Tanimoto_Coefficient(mn->Npair,molL,molQ);
  ln->select_dis = mn->select_dis;
  ln->score_for_sort = ln->tanimoto_mcs;
 /*
  ln->Nmcs     =  mn->Npair;
  Malloc_MATCH(&(ln->match), mn->Npair); 
  Copy_MATCH(&(ln->match), mn); 
 */
  ln->Npair = mn->Npair;
  ln->anumA = (unsigned char *)malloc(sizeof(unsigned char)*mn->Npair);
  ln->anumB = (unsigned char *)malloc(sizeof(unsigned char)*mn->Npair);
  ln->num_in_fileA = (int *)malloc(sizeof(int)*mn->Npair);
  ln->num_in_fileB = (int *)malloc(sizeof(int)*mn->Npair);
  for (i=0;i<mn->Npair;++i){
    ln->anumA[i] = mn->anumA[i];
    ln->anumB[i] = mn->anumB[i];
    ln->num_in_fileA[i] = molQ->atoms[mn->anumA[i]].num_in_file;
    ln->num_in_fileB[i] = molL->atoms[mn->anumB[i]].num_in_file;
  }
 /* Set up alph_score : Alphabet score converted from name[] (0>1>2..>A>B>C..>a>b>c) 
     The alph_score "ABC" is larger than "DEF". 
 */
  ln->alph_score = 0.0;
  log128 = log(128.0);
  for (i=0;i<strlen(ln->name);++i){
    ln->alph_score += (float)(128-ln->name[i])/128.0 * exp(-(float)i*log128);
  }
  /* printf("'%s' %.10f\n",ln->name,ln->alph_score); */

} /* end of Copy_MATCH_Info_to_LIBMOL() */





void Write_Search_Result(ofname,molQ,idirlib,ideslib,HeadList,Nmollib,Nmol_sele_desc,thre_tanimoto_mcs,thre_tanimoto_oneatom,thre_tanimoto_atompair,HeadLibFile)
 char   *ofname;
 struct MOLECULE *molQ;
 char   *idirlib;              /* input directory for library */ 
 char   *ideslib;              /* input descriptor file for library */ 
 struct LIBMOL *HeadList;
 int    Nmollib;               /* Number of compound in library. Nmollib may be larger than length of HeadList */
 int    Nmol_sele_desc;        /* Number of compound in selected by descriptor (one-atom, atom-pair) */
 float  thre_tanimoto_mcs;
 float  thre_tanimoto_oneatom;
 float  thre_tanimoto_atompair;
 struct LINENODE *HeadLibFile;
{
 FILE *fpo;
 struct LIBMOL *xn;
 int Nrank,i,q,*align_onQ,c,column, Nmol_over_thre;
 char str[100],spacestr[100],molform[256],class[128];
 char buff1[100],buff2[100],buff3[100],buff4[100],buff5[100],buff6[100];
 int k,LenMolName,LenNumMol,LenNumLib,LenOffset,Nlibfile,max_num,max_num_file,NdoneMCS;
 long max_offset_libfile;
 struct LINENODE *ln; 
 align_onQ = (int*)malloc(sizeof(int)*(molQ->Natom+1)); 
 printf("#Write_Search_Result(thre_tanimoto mcs %f atompair %f oneatom %f)-->'%s'\n",thre_tanimoto_mcs, thre_tanimoto_atompair, thre_tanimoto_oneatom,ofname);
 fpo = fopen(ofname,"w");
 if (fpo==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }

 /*** [0] Check Output Molecules **/ 
 LenMolName = 9; 
 Nmol_over_thre = 0;
 max_offset_libfile = 0; max_num = max_num_file = 0;
 NdoneMCS = 0;
 xn = HeadList;
 while (xn->next != NULL){
   xn = xn->next;
   xn->subnum = 0;
   ++Nrank;
   if (xn->doneMCS == 1){
     NdoneMCS += 1;
   }

   /* check output conditions */
   if (((thre_tanimoto_mcs<0.0)     ||(xn->tanimoto_mcs      >= thre_tanimoto_mcs)) && 
       ((thre_tanimoto_atompair<0.0)||(xn->tanimoto_atompair >= thre_tanimoto_atompair)) && 
       ((thre_tanimoto_oneatom<0.0) ||(xn->tanimoto_oneatom  >= thre_tanimoto_oneatom)) ) {
          xn->subnum = 1;
          Nmol_over_thre += 1; 
          if ((strlen(xn->name)+1)>LenMolName) LenMolName = strlen(xn->name)+1;
          if (xn->offset_libfile>max_offset_libfile) max_offset_libfile = xn->offset_libfile;
          if (xn->num      > max_num) max_num = xn->num;
          if (xn->num_file > max_num_file) max_num_file = xn->num_file;
   }
 }

 Nlibfile = 0; 
 ln = HeadLibFile;
 while (ln->next != NULL){
   ln = ln->next;
   Nlibfile += 1;
 }

 if (LenMolName>=100) LenMolName = 99;
 if (LenMolName<=9)   LenMolName = 10;

 sprintf(buff1,"%d",max_num);
 LenNumMol = strlen(buff1);
 sprintf(buff1,"%d",max_num_file);
 LenNumLib = strlen(buff1);
 sprintf(buff1,"%ld",max_offset_libfile);
 LenOffset = strlen(buff1);

 for (k=0;k<(LenMolName+LenNumMol+LenNumLib+LenOffset+3);++k) spacestr[k] = ' ';
 spacestr[LenMolName+LenNumMol+LenNumLib+LenOffset+3] = '\0';
 

 fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
 Set_END_DATE();
 fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
 fprintf(fpo,"#DATE_END   %s\n",PAR.END_DATE);
 fprintf(fpo,"#COMP_TIME_PROGRESS  %lf seconds\n",PAR.COMP_TIME_SEC_PROGRESS);
 fprintf(fpo,"#COMP_TIME           %lf seconds\n",PAR.COMP_TIME_SEC);
 fprintf(fpo,"#QUERY_FILENAME    %s\n",molQ->filename);
 fprintf(fpo,"#QUERY_NAME        %s\n",molQ->name);
 fprintf(fpo,"#QUERY_FILETYPE    %c\n",molQ->filetype);
 fprintf(fpo,"#QUERY_MOLFORMULA  %s\n",molQ->molform);
 fprintf(fpo,"#QUERY_NATOM       %d\n",molQ->Natom);
 fprintf(fpo,"#QUERY_NHEAVYATOM  %d\n",molQ->Nheavyatom);
 fprintf(fpo,"#ALGORITHM            %c\n",PAR.AlgoType);
 fprintf(fpo,"#CONNECT_GRAPH        %c\n",PAR.ConnectGraphType);
 fprintf(fpo,"#MAX_DIF_TOPODIS      %d\n",PAR.maxDIFtopodis);
 fprintf(fpo,"#THRE_TANIMOTO_MCS      %f\n",thre_tanimoto_mcs);
 fprintf(fpo,"#THRE_TANIMOTO_ATOMPAIR %f\n",thre_tanimoto_atompair);
 fprintf(fpo,"#THRE_TANIMOTO_ONEATOM  %f\n",thre_tanimoto_oneatom);
 fprintf(fpo,"#LIBRARY_DIRECTORY  %s\n",idirlib);
 fprintf(fpo,"#LIBRARY_DESCRIPTOR %s\n",ideslib);
 fprintf(fpo,"#NMOL_IN_LIBRARY               %d\n",Nmollib);
 fprintf(fpo,"#NMOL_SELECTED_BY_DESCRIPTORS  %d\n",Nmol_sele_desc);
 fprintf(fpo,"#NMOL_CAL_MCS                  %d\n",NdoneMCS);
 fprintf(fpo,"#NMOL_OVER_THRE                %d\n",Nmol_over_thre);
 fprintf(fpo,"#N_LIBRARY_FILE %d\n",Nlibfile);
 ln = HeadLibFile;
 while (ln->next != NULL){
   ln = ln->next;
   fprintf(fpo,"#LIBRARY_FILE_%d %s\n",ln->num,ln->line); 
 }
 fprintf(fpo,"[SIMILAR_MOLECULE_LIST]\n");
 fprintf(fpo,"#[rank(1)] [molecular_num(2)] [file_num(3)] [file_offset(4)] [molecular_name(5)]\n");
 fprintf(fpo,"#[similarity MCS(6)] [similarity ONEATOM(7)] [similarity ATOMPAIR(8)]\n");
 fprintf(fpo,"#[Natompair(9)] [Nheavyatom(10)] [select_dis(11)] [molecular formula(12)]\n");
 
 /** output atomname and atomnumbe for query molecul e**/ 
 fprintf(fpo,"#%s                                            [element_for_query]:",spacestr);
 for (q=0;q<molQ->Natom;++q){ 
   if (molQ->atoms[q].one_char_ele!='H') fprintf(fpo,"%c",molQ->atoms[q].one_char_ele);
 }
 fprintf(fpo,"\n#%s                                               [ring_for_query]:",spacestr);
 for (q=0;q<molQ->Natom;++q){ 
   if (molQ->atoms[q].one_char_ele!='H') fprintf(fpo,"%c",molQ->atoms[q].ringblock);
 }
 fprintf(fpo,"\n");
 fprintf(fpo,"#%s                                        [atom_number_for_query]:",spacestr);

 column = 1;
 if (molQ->Natom >0){
        if (molQ->atoms[molQ->Natom-1].num_in_file>=10000) column = 5;
   else if (molQ->atoms[molQ->Natom-1].num_in_file>=1000)  column = 4;
   else if (molQ->atoms[molQ->Natom-1].num_in_file>=100)   column = 3;
   else if (molQ->atoms[molQ->Natom-1].num_in_file>=10)    column = 2;
 }

 /* printf("# %d column %d\n",molQ->atoms[molQ->Natom-1].num_in_file,column); */
 for (c=0;c<column;++c){
   for (q=0;q<molQ->Natom;++q){ 
     if (molQ->atoms[q].one_char_ele!='H'){
        if (column==5) sprintf(str,"%5d",molQ->atoms[q].num_in_file);
        if (column==4) sprintf(str,"%4d",molQ->atoms[q].num_in_file);
        if (column==3) sprintf(str,"%3d",molQ->atoms[q].num_in_file);
        if (column==2) sprintf(str,"%2d",molQ->atoms[q].num_in_file);
        if (column==1) sprintf(str,"%1d",molQ->atoms[q].num_in_file);
         fprintf(fpo,"%c",str[c]); 
     }
   }
   if (c<(column-1)) fprintf(fpo,"\n#%s                                                               :",spacestr);
 }
 fprintf(fpo,"\n");

 spacestr[strlen(spacestr)-10] = '\0';
 fprintf(fpo,"#[rk][num][name]%s[sMCS][sONE][sPAIR][Npair][Nhvyatm][seldis][molform]      :",spacestr);

 for (q=0;q<molQ->Natom;++q){ 
   if (molQ->atoms[q].one_char_ele!='H') fprintf(fpo,"-");
 }
 fprintf(fpo,"\n"); 
 /*** [1] output searched library molecules in one line **/
 Nrank = 0;
 xn = HeadList;
 while (xn->next != NULL){
   xn = xn->next;
   if (xn->subnum==1){
     ++Nrank;
     for (k=0;k<LenMolName;++k){
       if (k<strlen(xn->name)) str[k] = xn->name[k];
                          else str[k] = ' ';
     } 
     str[LenMolName] = '\0';
     
     if (xn->molform==NULL){ molform[0] = '\0';} else { sprintf(molform,"%s",xn->molform); }
     if (xn->class==NULL)  { class[0] = '\0';}   else {sprintf(class,"%s",xn->class);}
     sprintf(buff1,"%d",xn->num);
     /* printf("LenNumMol %d buff1 '%s'\n",LenNumMol,buff1); fflush(stdout); */
     for (k=0;k<(LenNumMol-strlen(buff1));++k) buff2[k] = ' ';  buff2[LenNumMol-strlen(buff1)] = '\0';
     
     sprintf(buff3,"%d",xn->num_file);
     for (k=0;k<(LenNumLib-strlen(buff3));++k) buff4[k] = ' ';  
     buff4[LenNumLib-strlen(buff3)] = '\0';
    
     sprintf(buff5,"%ld",xn->offset_libfile);
     for (k=0;k<(LenOffset-strlen(buff5));++k) buff6[k] = ' ';  buff6[LenOffset-strlen(buff5)] = '\0';
  
     fprintf(fpo,"%-4d %s%s %s%s %s%s %s %5.3f %5.3f %5.3f %3d %3d %8.3f %-23s ",
      Nrank,buff1,buff2,buff3,buff4,buff5,buff6,str,xn->tanimoto_mcs,xn->tanimoto_oneatom,xn->tanimoto_atompair,xn->Npair,xn->Nheavyatom,xn->select_dis,molform);

     if (xn->doneMCS == 1){ 
       for (q=0;q<molQ->Natom;++q) align_onQ[q] = 0;
       for (i=0;i<xn->Npair;++i){
           if (xn->anumA[i] >= molQ->Natom)
            {printf("#ERROR: i-th anumA %d is over molQ->Natom %d.\n",xn->anumA[i],molQ->Natom);  exit(1); }
          align_onQ[xn->anumA[i]] = 1;
       } 
       for (q=0;q<molQ->Natom;++q){
         if (molQ->atoms[q].one_char_ele!='H'){ 
           if (align_onQ[q]==1) {fprintf(fpo,"%c",molQ->atoms[q].one_char_ele);}
                          else  {fprintf(fpo,"-");}
         }
       }
     /*
      fprintf(fpo,"%8.3f %8.3f %s",xn->rmsd,xn->rmsd_frm_init,xn->prop);
     */
     } 
      fprintf(fpo," %s",class);
      fprintf(fpo,"\n");
   } 
 }
 /** [2] output atom matching in one line **/ 
 fprintf(fpo,"[ATOM_MATCHING]\n");
 Nrank = 0;
 xn = HeadList;
 while (xn->next != NULL){
   xn = xn->next;
   if (xn->subnum==1){
     ++Nrank;
     if (xn->doneMCS == 1){ 
       fprintf(fpo,"%-4d ",Nrank);
       for (i=0;i<xn->Npair;++i){fprintf(fpo,"|%d %d", xn->num_in_fileA[i],xn->num_in_fileB[i]);} 
       fprintf(fpo,"\n"); 
    }
  }
 }
 /** [3] output Raxisang and Tvec to superpose on the query molecule**/ 
 fprintf(fpo,"[SUPERPOSITION]\n");
 fprintf(fpo,"#[rank(1)] [Natompair(2)] [RMSD(3)]\n");
 fprintf(fpo,"#[rot_axisX(4)] [rot_axisY(5)] [rot_axisZ(6)] [rot_ang(7)]\n");
 fprintf(fpo,"#[trans_vecX(8)] [trans_vecY(9)] [trans_vecZ(10)]\n");
 Nrank = 0;
 xn = HeadList;
 while (xn->next != NULL){
   xn = xn->next;
   if (xn->subnum==1){
     ++Nrank;
     if (xn->doneMCS == 1){ 
       fprintf(fpo,"%-4d %3d",Nrank,xn->Npair);
       fprintf(fpo," %+9.5f %+9.6lf %+9.6lf %+9.6lf %+9.6lf %+9.5lf %+9.5lf %+9.5lf\n",xn->RMSD,
          xn->Raxisang[0], xn->Raxisang[1], xn->Raxisang[2], xn->Raxisang[3],
          xn->Tvec[0], xn->Tvec[1], xn->Tvec[2]);
    }
  }
 }

 Nrank = 0;




 free(align_onQ); 
 fclose(fpo);
} /* end of Write_Search_Resultk() */




int Write_Similar_Molecules_To_One_File(ofname,ilibfile,idirlib,HeadList,thre_tanimoto_mcs,thre_tanimoto_oneatom,thre_tanimoto_atompair,filetype,Nlibfile,HeadLibFile,Nmax_outmol)
 char   *ofname;
 char   *ilibfile;             /* input library file */
 char   *idirlib;              /* input directory for library */
 struct LIBMOL *HeadList;
 float  thre_tanimoto_mcs;
 float  thre_tanimoto_oneatom;
 float  thre_tanimoto_atompair;
 char   filetype;  /* 'S'df, '2':mol2, 'P'db */ 
 int    Nlibfile;
 struct LINENODE *HeadLibFile;
 int    Nmax_outmol; /* maximum number of output molecules */
{
  FILE *fpi,*fpo;
  struct LIBMOL *xn;
  char line[500],ifname[500],ifname0[500];
  struct LINENODE HeadLineNode,*ln;
  int Noutmol;

  printf("#Write_Similar_Molecules_To_One_File(filetype %c) --> '%s'\n",filetype,ofname);
  fpo = fopen(ofname,"w");
  if (fpo==NULL){
    printf("#WARNING:Can't write to '%s'\n",ofname);
    return(0);
  }     
  ifname[0] = ifname0[0] = '\0';
  Noutmol = 0;
  fpi = NULL;

  xn = HeadList;
  while (xn->next != NULL){
    xn = xn->next;
    /* check output conditions */
    if (((thre_tanimoto_mcs<0.0)     ||(xn->tanimoto_mcs      >= thre_tanimoto_mcs)) && 
        ((thre_tanimoto_atompair<0.0)||(xn->tanimoto_atompair >= thre_tanimoto_atompair)) && 
        ((thre_tanimoto_oneatom<0.0) ||(xn->tanimoto_oneatom  >= thre_tanimoto_oneatom)) && 
        ((Nmax_outmol<0) || (Noutmol < Nmax_outmol)) 
       ) {

      Noutmol += 1;      

      if (Nlibfile==0){
        sprintf(ifname,"%s/%s",idirlib,xn->name);
        fpi = fopen(ifname,"r");
        while (feof(fpi)==0){
          line[0] = '\0';
          fgets(line,500,fpi);
          fprintf(fpo,"%s",line); 
        } 
        fclose(fpi);
      }
      else {
       String_from_LINENODE_num(ifname,HeadLibFile,xn->num_file);
       if (strcmp(ifname,ifname0)!=0){
         if (fpi != NULL) fclose(fpi);
         fpi = fopen(ifname,"r");
       }
       fseek(fpi,xn->offset_libfile,SEEK_SET);
       Read_One_Entry_as_LINENODEs(fpi,&HeadLineNode,filetype);
       ln = &HeadLineNode;
       while (ln->next != NULL){
         ln = ln->next;
         fprintf(fpo,"%s\n",ln->line);
       }
       Free_LINENODEs(&HeadLineNode);
       sprintf(ifname0,"%s",ifname);
     }
    }

  } /* xn */

  fclose(fpo);

  printf("  write %d molecules.\n",Noutmol); 
  return(Noutmol);

} /* end of Write_Similar_Molecules_To_One_File() */




int Write_Similar_Molecules_Under_Dir(outdir,ilibfile,idirlib,HeadList,thre_tanimoto_mcs,thre_tanimoto_oneatom,thre_tanimoto_atompair,filetype,Nlibfile,HeadLibFile,Nmax_outmol)
 char   *outdir;
 char   *ilibfile;         /* input library file */
 char   *idirlib;          /* input directory for library */
 struct LIBMOL *HeadList;
 float  thre_tanimoto_mcs;
 float  thre_tanimoto_oneatom;
 float  thre_tanimoto_atompair;
 char   filetype;  /* 'S'df, '2':mol2, 'P'db */ 
 int    Nlibfile;
 struct LINENODE *HeadLibFile;
 int    Nmax_outmol; /* maximum number of output molecules */
{
  FILE *fpi,*fpo;
  struct LIBMOL *xn;
  char line[500],ifname[500],ifname0[500],ofname[500],buff[128];
  struct LINENODE HeadLineNode,*ln;
  int Noutmol;
  int Nword, Wsta[100], Wend[100];

  printf("#Write_Similar_Molecules_Under_Dir(filetype %c) --> '%s/'\n",filetype,outdir);
  ifname[0] = ifname0[0] = '\0';
  Noutmol = 0;
  fpi = NULL;

  xn = HeadList;
  while (xn->next != NULL){
    xn = xn->next;
    /* check output conditions */
    if (((thre_tanimoto_mcs<0.0)     ||(xn->tanimoto_mcs      >= thre_tanimoto_mcs)) && 
        ((thre_tanimoto_atompair<0.0)||(xn->tanimoto_atompair >= thre_tanimoto_atompair)) && 
        ((thre_tanimoto_oneatom<0.0) ||(xn->tanimoto_oneatom  >= thre_tanimoto_oneatom))  &&
        ((Nmax_outmol<0) || (Noutmol < Nmax_outmol)) 
       ) {

      Noutmol += 1;
     
      sprintf(ofname,"%s/%s",outdir,xn->name); 
      /* printf("-->'%s'\n",ofname); */
      fpo = fopen(ofname,"w");
      if (fpo==NULL){
        printf("#WARNING:Can't write to '%s'\n",ofname);
        if (First_Appear_Pos_in_String(xn->name,'/')>=0){
          Split_to_Words(xn->name,'/',&Nword,Wsta,Wend,100);
          Get_Part_Of_Line(buff,xn->name,Wsta[Nword-1],Wend[Nword-1]);  
          sprintf(ofname,"%s/%s",outdir,buff); 
          fpo = fopen(ofname,"w");
          if (fpo==NULL){ 
            printf("#WARNING:Can't write to '%s'\n",ofname);
            sprintf(ofname,"%s/%d",outdir,xn->num); 
            fpo = fopen(ofname,"w");
            if (fpo==NULL){ 
              printf("#WARNING:Can't write to '%s'\n",ofname);
              return(0);
            }
          }
        }
        else{
          return(0);
        }
      }     
      printf("#write to -> '%s'\n",ofname);
 
      if (Nlibfile==0){
        sprintf(ifname,"%s/%s",idirlib,xn->name);
        fpi = fopen(ifname,"r");
        while (feof(fpi)==0){
          line[0] = '\0';
          fgets(line,500,fpi);
          fprintf(fpo,"%s",line); 
        } 
        fclose(fpi);
      }
      else {
       String_from_LINENODE_num(ifname,HeadLibFile,xn->num_file);
       if (strcmp(ifname,ifname0)!=0){
         if (fpi != NULL) fclose(fpi);
         fpi = fopen(ifname,"r");
       }
       fseek(fpi,xn->offset_libfile,SEEK_SET);
       Read_One_Entry_as_LINENODEs(fpi,&HeadLineNode,filetype);
       ln = &HeadLineNode;
       while (ln->next != NULL){
         ln = ln->next;
         fprintf(fpo,"%s\n",ln->line);
       }
       Free_LINENODEs(&HeadLineNode);
       sprintf(ifname0,"%s",ifname);
     }

     fclose(fpo);

    }

  } /* xn */

  printf("  write %d molecules.\n",Noutmol); 
  return(Noutmol);

} /* end of Write_Similar_Molecules_Under_Dir() */







void Write_Clustering_Res_List_Of_NODEs(ofname,HeadList,tanimoto_thre)
 char *ofname;
 struct LIBMOL *HeadList;
 float  tanimoto_thre;
{
 FILE *fpo;
 struct LIBMOL *xn;
 int Nrank;

 printf("#Write_Clustering_Res_List_Of_NODEs()-->'%s'\n",ofname);
 fpo = fopen(ofname,"w");
 if (fpo==NULL){
   printf("#ERROR:Can't write to '%s'\n",ofname);
   exit(1);
 }
 fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
 fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
 fprintf(fpo,"#DATE_END   %s\n",Get_Date_String());
 fprintf(fpo,"#TANIMOTO_THRE    %f\n",tanimoto_thre);
 fprintf(fpo,"#NMOL_IN_LIBRARY  %d\n",Number_of_LIBMOL(HeadList));

 xn = HeadList;
 Nrank = 0;
 while (xn->next != NULL){
   xn = xn->next;
   if (xn->prev->score_for_sort != xn->score_for_sort)
     fprintf(fpo,"%s cluster_num %.0f Nmember %d Nheavyatom %d %s\n",xn->name,xn->score_for_sort,xn->subnum,xn->Nheavyatom,xn->molform);
    else
     fprintf(fpo,"#%s Nheavyatom %d %s\n",xn->name,xn->Nheavyatom,xn->molform);
 }
 fclose(fpo);
} /* end of Write_Clustering_Res_List_Of_NODEs() */




void Cal_Nsta_Nend_from_Bunshi_Bunbo_for_all_vs_all(Nlist,bunshi,bunbo,Nsta,Nend)
  int Nlist, bunshi, bunbo; /* (input) */
  int *Nsta,*Nend;          /* (to be calculated) */
{
 /*
   Assumed all-vs-all comparison, such as:
     for (i=Nsta;i<Nend;++i)
       for (j=i+1;j<Nlist;++j) Compare(i,j)
 */
 float f0,f1,r0,r1;

 r0 = (float)bunshi/(float)bunbo;
 r1 = (float)(bunshi+1)/(float)bunbo;
 f0 = 1.0 - sqrt(1-r0);
 f1 = 1.0 - sqrt(1-r1);
 /* printf("r0 %f r1 %f f0 %f f1 %f\n",r0,r1,f0,f1); */
 *Nsta = (int)floor(Nlist*f0);
 *Nend = (int)floor(Nlist*f1);
 if (bunshi==0)         *Nsta = 0;
 if (bunshi==(bunbo-1)) *Nend = Nlist;

} /* end of Cal_Nsta_Nend_from_Bunshi_Bunbo_for_all_vs_all() */



int Number_of_LIBMOL(head)
 struct LIBMOL *head;
{
 struct LIBMOL *xn;
 int N;
 N = 0; 
 xn = head;
 while (xn->next !=NULL){
  xn = xn->next;
  N += 1;
 }
 return(N);
} /* end of Number_of_LIBMOL() */



struct LIBMOL *Get_LIBMOL_by_number(head,num)
 struct LIBMOL *head;
 int num;
{
  struct LIBMOL *xn;
  xn = head;
  while (xn->next !=NULL){
    xn = xn->next;
    if (xn->num == num){
      return(xn);
    } 
  }
  return(NULL);
} /* end of Get_LIBMOL_by_number() */



struct LIBMOL* Add_LIBMOL_to_LIBMOLlist(add,HeadLibMol)
 struct LIBMOL *add;
 struct LIBMOL *HeadLibMol;
{
  struct LIBMOL *new,*last;

  last = HeadLibMol;
  while (last->next != NULL) { last = last->next; }

  new = (struct LIBMOL *)malloc(sizeof(struct LIBMOL));
  last->next = new;
  new->prev = last;
  new->next = NULL;
  new->num = last->num + 1;
  
  new->Nheavyatom = add->Nheavyatom;
  new->select_dis = add->select_dis;
  new->score_for_sort = add->score_for_sort;
  new->doneMCS = add->doneMCS;
  new->Npair = add->Npair;
  new->num_file       = add->num_file;
  new->offset_libfile = add->offset_libfile;
  new->offset_descfile = add->offset_descfile;
  new->tanimoto_mcs  = add->tanimoto_mcs;
  new->tanimoto_oneatom  = add->tanimoto_oneatom;
  new->tanimoto_atompair  = add->tanimoto_atompair;
  new->subnum  = add->subnum;

  new->mol = add->mol;

  if ((add->name != NULL)&&(strlen(add->name)>0)){
    new->name = (char *)malloc(sizeof(char)*(strlen(add->name)+1));
    sprintf(new->name,"%s",add->name);
  }  else new->name = NULL;

  if ((add->class != NULL)&&(strlen(add->class)>0)){
    new->class = (char *)malloc(sizeof(char)*(strlen(add->class)+1));
    sprintf(new->class,"%s",add->class);
  }  else new->class = NULL;

  if ((add->molform != NULL)&&(strlen(add->molform)>0)){
    new->molform = (char *)malloc(sizeof(char)*(strlen(add->molform)+1));
    sprintf(new->molform,"%s",add->molform);
  }  else new->molform = NULL;

  new->anumA = NULL;
  new->anumB = NULL;
  new->num_in_fileA = NULL;
  new->num_in_fileB = NULL;
  new->mol = NULL;
  new->oneatom_descriptor  = NULL;
  new->atompair_descriptor = NULL;
  new->ring_descriptor     = NULL;
  
  return(new);

} /* end of Add_LIBMOL_to_LIBMOLlist() */






void Write_Similarities_in_List_format_Header(fpo,HeadOfList,idirlib)
   FILE *fpo; 
   struct LIBMOL *HeadOfList;  
   char *idirlib;
{
  struct LIBMOL *nn;
  int Nlist,i;

  Nlist = Number_of_LIBMOL(HeadOfList);

  fprintf(fpo,"#>>Similarities in List format<<\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
  if (idirlib[0] !='\0') fprintf(fpo,"#LIBRARY_DIRECTORY %s\n",idirlib);
  fprintf(fpo,"#NMOL_IN_LIBRARY %d\n",Nlist);
  nn = HeadOfList;
  for (i=0;i<Nlist;++i){
    nn = nn->next;
    fprintf(fpo,"#DATA %d %s %s",i,nn->mol->filename,nn->name);
    if (nn->class!=NULL) fprintf(fpo," %s",nn->class);
    fprintf(fpo,"\n");
  }
} /* end of Write_Similarities_in_List_format_Header() */




void Write_Two_Groups_Similarities_in_List_format_Header(fpo,HeadOfListQ,HeadOfListL)
   FILE *fpo; 
   struct LIBMOL *HeadOfListQ; 
   struct LIBMOL *HeadOfListL; 
{
  struct LIBMOL *nn;
  int NlistL,NlistQ;

  NlistQ = Number_of_LIBMOL(HeadOfListQ);
  NlistL = Number_of_LIBMOL(HeadOfListL);

  fprintf(fpo,"#>>Similarities in List format<<\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
  fprintf(fpo,"#NMOL_IN_QUERY   %d\n",NlistQ);
  fprintf(fpo,"#NMOL_IN_LIBRARY %d\n",NlistL);
  nn = HeadOfListQ;
  while (nn->next !=NULL){
    nn = nn->next;
    fprintf(fpo,"#DATAQ %d %s\n",nn->num,nn->name);
  }
  nn = HeadOfListL;
  while (nn->next !=NULL){
    nn = nn->next;
    fprintf(fpo,"#DATAL %d %s\n",nn->num,nn->name);
  }

} /* end of Write_Two_Groups_Similarities_in_List_format_Header() */










void Write_Similarities_in_List_format_Tail(fpo)
  FILE *fpo;
{
  Set_END_DATE();
  fprintf(fpo,"#COMMAND    %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
  fprintf(fpo,"#DATE_END   %s\n",PAR.END_DATE);
  fprintf(fpo,"#COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);
}





void Increment_Nmatch_for_Matched_Atoms(M,molA,molB)
 struct MATCH *M;
 struct MOLECULE *molA, *molB;
{
 int m,a,b;

  for (m=0;m<M->Npair;++m){
    a = M->anumA[m];
    b = M->anumB[m];
    if ((a>=0) && (b>=0)){
      molA->atoms[a].Nmatch += 1;
      molA->Nmatch += 1;
      molB->atoms[b].Nmatch += 1;
      molB->Nmatch += 1;
    }
  }
} /* end of Increment_Nmatch_for_Matched_Atoms() */






void Write_Nmatch_rasmol_Script(ofname,mol,Nmollib,comment)
 char *ofname;
 struct MOLECULE *mol;
 int Nmollib;
 char *comment;
{
 FILE *fp;
 char color[10][16];
 int Ncolor,i,cnum,maxNmatch,size;
 float percent;

 printf("#Write_Nmatch_rasmol_Script(Nmollib:%d)-->'%s'\n",Nmollib,ofname);

 Ncolor = 9;
 sprintf(color[8],"red");    sprintf(color[7],"redorange");
 sprintf(color[6],"orange"); sprintf(color[5],"pink");
 sprintf(color[4],"yellow"); sprintf(color[3],"green");
 sprintf(color[2],"cyan");   sprintf(color[1],"purple");
 sprintf(color[0],"blue");
 maxNmatch = 0;
 for (i=0;i<mol->Natom;++i){
   if (mol->atoms[i].Nmatch > maxNmatch) maxNmatch = mol->atoms[i].Nmatch;
   /* printf("%d %d %d\n",i,mol->atoms[i].Nmatch,maxNmatch);  */
 }
 fp = fopen(ofname,"w");
 fprintf(fp,"# Write_Nmatch_rasmol_Script(Nmollib%d)-->'%s' #\n",Nmollib, ofname);
 fprintf(fp,"#>>Procedures to show Nmatch-colored and -sized molecular view by RasMol\n");  
      if (mol->filetype=='S'){ fprintf(fp,"#rasmol -mdl "); } 
 else if (mol->filetype=='2'){ fprintf(fp,"#rasmol -mol2 "); } 
 else { fprintf(fp,"#rasmol "); } 
 fprintf(fp,"%s\n",mol->filename);
 fprintf(fp,"#RasMol> script \"%s\"\n",ofname);
 if (comment[0] != '\0'){ fprintf(fp,"#%s\n",comment); }
 fprintf(fp,"echo \">> Write_Nmatch_rasmol_Script(Nmollib%d)-->'%s' <<\"\n",Nmollib, ofname);
 fprintf(fp,"echo \"COMMAND %s\"\n",PAR.COMMAND);
 fprintf(fp,"echo \"%s %-20s Nmatch %4d maxNmatch %4d\"\n",mol->filename,mol->name,mol->Nmatch,maxNmatch);
 if (comment[0] != '\0'){ fprintf(fp,"echo \"%s\"\n",comment); }
 

 fprintf(fp,"select all\n");
 fprintf(fp,"wireframe\n");
 fprintf(fp,"spacefill false\n");
 fprintf(fp,"color gray\n");
 for (i=0;i<mol->Natom;++i){
   if (mol->atoms[i].Nmatch>0){
     cnum = (int)((float)Ncolor*mol->atoms[i].Nmatch/(float)(maxNmatch+1));
     percent = 100.0*mol->atoms[i].Nmatch/Nmollib;
     size = (int)(200.0*percent/100.0);

     fprintf(fp,"echo \"atom %3d %3s Nmatch %3d/%3d %6.1f %% %-9s size %d\"\n"
     ,mol->atoms[i].num_in_file,mol->atoms[i].element,mol->atoms[i].Nmatch,Nmollib, percent,color[cnum],size);
     fprintf(fp,"select atomno=%d\n",mol->atoms[i].num_in_file);
     fprintf(fp,"color  %s\n",color[cnum]);
     fprintf(fp,"spacefill %d\n",size);
   }
 }
 fprintf(fp,"select all\n");
 fclose(fp);
} /* end of Write_Nmatch_rasmol_Script() */






void insert_LIBMOL_into_score_decreasing_proper_position(n,head)
  struct LIBMOL *n;     /* newly adding node */
  struct LIBMOL *head;
/*
 
 Sort LIBMOL node in the Increasing order.
 
 >> example for n = [3] << 

 
 [head]->[1] ==> [head]->[1]->[3]
 [head]->[5] ==> [head]->[3]->[5]
 [head]->[1]->[5] ==> [head]->[1]->[3]->[5]
 [head]->[1]->[2]->[6] ==> [head]->[1]->[2]->[3]->[6]

     [head]->[ ]->[o]->[p]->[q]
 ==> [head]->[ ]->[o]->[p]->[n]->[q]

*/
{
  struct LIBMOL *p; /* insert position */
  struct LIBMOL *x,*o,*q;
  char n_is_min;
/*
  printf("#insert_LIBMOL_into_score_decreasing_proper_position(n:%f head)\n",n->score_for_sort);
*/
  x = head;
  p = NULL;
  n_is_min = 1;
  while ((x->next != NULL)&&(p==NULL)){
    x = x->next;
    if (    (x->score_for_sort <= n->score_for_sort ) 
         && ((x->next ==NULL) || (n->score_for_sort <=x->next->score_for_sort)) ){
       p = x;
    } 
    if (x->score_for_sort <  n->score_for_sort ) { n_is_min = 0;}
  }

  if (p!=NULL){
    q = p->next; 
    o = p->prev; 
    p->next = n;
    n->next = q;
    n->prev = p;
    if (q != NULL){ q->prev = n;}
  }

  /* put in just after the head */

  /* printf("#new_value %f n_is_min %d\n",n->score_for_sort,n_is_min); */

  if ((p==NULL) && (n_is_min == 1)) {
    x = head->next; 
    head->next = n;
    n->prev = head;
    n->next = x;
    if (x !=NULL){ x->prev = n;}
  }

} /* end of insert_LIBMOL_into_score_decreasing_proper_position() */


void delete_the_first_LIBMOL_node(head)
  struct LIBMOL *head;
{
  struct LIBMOL *x;

  /* printf("#delete_the_first_LIBMOL_node(head)\n"); */
  x = head->next;
  if (x != NULL){
    head->next = x->next;
    if (x->next != NULL){
      x->next->prev = head;
    }
    free(x->name);
    free(x->molform);
    free(x);
  }

} /* end of delete_the_first_LIBMOL_node() */



void reverse_LIBMOL_list(head)
  struct LIBMOL *head;
{
/*

 [head]->[1]->[2]->[3]->[4]->NULL
 ==>
 [head]->[4]->[3]->[2]->[1]


 [head]->[1]->[2]->[3]->[4]->[5]->NULL
 ==>
 [head]->[5]->[4]->[3]->[2]->[1]


 [head]->[]->[sta]->[]->[]->[end]->[]-> NULL
 [head]->[]->[sta]->[]->[end]->[]-> NULL
 [head]->[]->[sta]->[end]->[]-> NULL

*/
  struct LIBMOL *end,*sta;
  struct LIBMOL *enext,*eprev,*snext,*sprev;
  int n,finish,neighbor;
  printf("#reverse_LIBMOL_list(head)\n");

  end = head;
  while (end->next != NULL){
    end = end->next;
  }
  printf("#end %f\n",end->score_for_sort);
  sta = head->next;
 
  n = 0;
  finish = 0;
  while ((finish==0)&&(end != NULL) && (sta !=NULL)){
     /* printf("sta %f end %f\n",sta->score_for_sort,end->score_for_sort); */
     enext = end->next;
     eprev = end->prev;
     snext = sta->next;
     sprev = sta->prev;

     neighbor = 0;
 /*    [head]->[]->[sta]->[end]->[]-> NULL */
     if ((snext == end)&&(eprev==sta)){ finish = 1; neighbor = 1;}
 /*    [head]->[]->[sta]->[]->[end]->[]-> NULL */
     else if (snext == eprev){ finish = 1;}

     if (neighbor ==0){
       sta->next = enext;
       sta->prev = eprev;
       end->next = snext;
       end->prev = sprev;
       if ((eprev!=sta)&&(eprev != NULL)){ eprev->next = sta;}
       if (enext != NULL){ enext->prev = sta;}
       if (sprev != NULL){ sprev->next = end;}
       if ((snext!=end)&&(snext != NULL)){ snext->prev = end;}
     }
     else if (neighbor ==1){
 /*    [head]->[]->[sta]->[end]->[]-> NULL */
 /* => [head]->[]->[end]->[sta]->[]-> NULL */
       sta->next = enext;
       sta->prev = end;
       end->next = sta;
       end->prev = sprev;
       if (enext != NULL){ enext->prev = sta;}
       if (sprev != NULL){ sprev->next = end;}
     }
  
     /*
     printf(">n %d\n",n);
     Show_LIST_of_LIBMOL(head);
     */

     if (finish != 1){
       end = eprev; 
       sta = snext; 
     }
  }

} /* end of reverse_LIBMOL_list() */



void Show_LIST_of_LIBMOL(head)
  struct LIBMOL *head;
{
  struct LIBMOL *x;
  int n;
  x = head;
  n = 0;
  while (x->next != NULL){
    x = x->next;
    printf("[%d] %f\n",n,x->score_for_sort);
    n += 1;
  } 
 
} /* end of Show_LIST_of_LIBMOL() */
