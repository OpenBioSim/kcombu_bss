/*

 <ioLINE.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 functions for "struct LINENODE" and "struct KEYVALUE" 

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <string.h>
#include <sys/stat.h>
#include "globalvar.h"
#include "molecule.h"

#if defined(_WIN32) || defined(WIN32)
#define popen _popen
#endif

/** FUNCTIONS (GLOBAL) **/
int Read_LINENODEs();
void Add_string_to_LINENODEs();
void Free_LINENODEs();
void Show_LINENODEs();
int  Number_of_LINENODEs();
void String_from_LINENODE_num();

void Add_key_value_string_to_KEYVALUEs();
void Show_KEYVALUEs();
void Free_KEYVALUEs();
int  Get_Value_from_KEYVALUE();
void Write_Property_for_One_Molecule();
void Write_Headers_of_Property_File();
int wildcard_match();
int substring_match();
int match_properties_and_patterns();


int Read_LINENODEs(fname,HeadLineNode)
 char *fname;
 struct LINENODE *HeadLineNode;
{
 FILE *fp;
 struct LINENODE *ln;
 char str[MAX_LENGTH_LINE], command[MAX_LENGTH_LINE];
 int L,nline,max_length_line;
 struct stat stat_Buf;
/*
 printf("#int Read_LINENODEs('%s')\n",fname); fflush(stdout); 
 */
 

 L = strlen(fname);
 fp = NULL; 
 if ((L>2)&&(fname[L-2]=='.')&&(fname[L-1]=='Z')){
   sprintf(command,"uncompress -c %s",fname);
   fp = popen(command,"r");
 }
 else if ((L>2)&&(fname[L-3]=='.')&&(fname[L-2]=='g')&&(fname[L-1]=='z')){
   sprintf(command,"zcat %s",fname);
   fp = popen(command,"r");
 }
 else{
   fp = fopen(fname,"r");
 }
 
 if (fp==NULL) {fprintf(stderr,"#ERROR(Read_LINENODEs):Can't open \"%s\"\n",fname); return(0);}
 if (stat(fname,&stat_Buf)!=0){
   fprintf(stderr,"#ERROR(Read_LINENODEs):Can't open  \"%s\"\n",fname); 
   return(0);
 }

 /** If fname is directory, return(0) : added by T.Kawabata (2019/10/24) **/ 
 if ((stat_Buf.st_mode &  S_IFMT)==S_IFDIR){
   fprintf(stderr,"#ERROR(Read_LINENODEs): \"%s\" is directory.\n",fname); 
   return(0);
 }


 HeadLineNode->next = NULL;
 HeadLineNode->prev = NULL;
 ln = HeadLineNode;
 max_length_line = MAX_LENGTH_LINE - 1;
 nline = 0;
 while (feof(fp)==0){
   str[0] = '\0';
   fgets(str,max_length_line,fp);
   /* printf("str '%s'\n",str);  */
   L = strlen(str);
   if (L>0){
     if (str[L-1]=='\n') {str[L-1] = '\0'; L -= 1;}
     if (str[L-1]=='\r') {str[L-1] = '\0'; L -= 1;}
     if (str[L-1]=='\n') {str[L-1] = '\0'; L -= 1;}
     if (str[L-1]=='\r') {str[L-1] = '\0'; L -= 1;}

     if (ln->next==NULL){
       ln->next = (struct LINENODE *)malloc(sizeof(struct LINENODE));
       ln->next->prev = ln;
       ln->next->next = NULL; }
     ln = ln->next;
     ln->line = (char *)malloc(sizeof(char)*(L+1));
     sprintf(ln->line,"%s",str);
     ln->num = nline;
     nline += 1;
   }
 }

 fclose(fp);
 return(nline);
} /* end of Read_LINENODEs() */



void Add_string_to_LINENODEs(newstr,HeadLineNode)
 char *newstr;
 struct LINENODE *HeadLineNode;
{
  struct LINENODE *newN,*lastN;

  lastN = HeadLineNode;
  HeadLineNode->num = -1;

  while (lastN->next != NULL){
    lastN = lastN->next;
  }

  lastN->next = (struct LINENODE *)malloc(sizeof(struct LINENODE));
  newN = lastN->next;
  newN->prev = lastN;
  newN->next = NULL;
  newN->line = (char *)malloc(sizeof(char)*(strlen(newstr)+1)); 
  sprintf(newN->line,"%s",newstr);
  newN->num = lastN->num + 1;

} /* end of Add_string_to_LINENODEs() */




void Free_LINENODEs(HeadLineNode)
 struct LINENODE *HeadLineNode;
{
 struct LINENODE *ln;

 ln = HeadLineNode;
 while (ln->next != NULL) ln = ln->next;

 while (ln->prev!=NULL){
   ln = ln->prev;
   free(ln->next->line); 
   free(ln->next);
 }

 HeadLineNode->next = NULL;

} /* end of Free_LINENODEs() */


void Show_LINENODEs(HeadLineNode)
 struct LINENODE *HeadLineNode;
{
  struct LINENODE *bn;
  bn = HeadLineNode;
  while (bn->next != NULL){
    bn = bn->next;
    printf("%d '%s'\n",bn->num,bn->line);
  }
} /* end of Show_LINENODEs() */


int  Number_of_LINENODEs(HeadLineNode)
 struct LINENODE *HeadLineNode;
{
  struct LINENODE *bn;
  int n;
  n = 0;
  bn = HeadLineNode;
  while (bn->next != NULL){
    bn = bn->next;
    n += 1;
  }
  return(n);
} /* end of Show_LINENODEs() */


void String_from_LINENODE_num(string,HeadLineNode,num)
 char   *string;
 struct LINENODE *HeadLineNode;
 int    num;
{
  struct LINENODE *bn;
/*
 *  printf("#String_from_LINENODE_num(string,HeadLineNode,num:%d)\n",num);
 *  */
  string[0] = '\0';
  bn = HeadLineNode;
  while (bn->next != NULL){
    bn = bn->next;
    if (bn->num==num)
    sprintf(string,"%s",bn->line);
  }

} /* end of Number_of_LINENODEs() */




void Add_key_value_string_to_KEYVALUEs(newkey,newvalue,HeadKeyValue)
 char *newkey,*newvalue;
 struct KEYVALUE *HeadKeyValue;
{
  struct KEYVALUE *newN,*lastN;
  char keyhit; 
  /*
  printf("#Add_key_value_string_to_KEYVALUEs('%s':'%s',HeadKeyValue)\n",newkey,newvalue);
 */
  lastN = HeadKeyValue;
  HeadKeyValue->num = -1;
  keyhit = 0;
  while ((lastN->next != NULL)&&(keyhit==0)){
     lastN = lastN->next;
     if (strcmp(lastN->key,newkey)==0) {keyhit = 1; }
  }

  if (keyhit==0){
    lastN->next = (struct KEYVALUE *)malloc(sizeof(struct KEYVALUE));
    newN = lastN->next;
    newN->prev = lastN;
    newN->next = NULL;
    newN->key   = (char *)malloc(sizeof(char)*(strlen(newkey)+1));
    newN->value = (char *)malloc(sizeof(char)*(strlen(newvalue)+1));
    sprintf(newN->key,  "%s",newkey);
    sprintf(newN->value,"%s",newvalue);
    newN->num = lastN->num + 1;
  }
  else if (keyhit==1){
    free(lastN->value);
    lastN->value = (char *)malloc(sizeof(char)*(strlen(newvalue)+1));
    sprintf(lastN->value,"%s",newvalue);
  }
} /* end of Add_key_value_string_to_KEYVALUEs() */





void Show_KEYVALUEs(HeadKeyValue)
 struct KEYVALUE *HeadKeyValue;
{
  struct KEYVALUE *bn;
/*
  printf("#void Show_KEYVALUEs(HeadKeyValue)\n");
*/
  bn = HeadKeyValue;
  while (bn->next != NULL){
    bn = bn->next;
    printf("%d '%s':'%s'\n",bn->num,bn->key,bn->value);
  }
} /* end of Show_KEYVALUEs() */



void Free_KEYVALUEs(HeadKeyValue)
 struct KEYVALUE *HeadKeyValue;
{
 struct KEYVALUE *ln;

 ln = HeadKeyValue;
 while (ln->next != NULL) ln = ln->next;

 while (ln->prev!=NULL){
   ln = ln->prev;
   free(ln->next->key);
   free(ln->next->value);
   free(ln->next);
 }

 HeadKeyValue->next = NULL;
} /* end of Free_KEYVALUEs() */


int Get_Value_from_KEYVALUE(valstr,HeadKeyValue,key)
  char  *valstr;
  struct KEYVALUE *HeadKeyValue;
  char *key;
{
 struct KEYVALUE *ln;
 int i;
 char wildcard;

 wildcard = 0;
 for (i=0;i<strlen(key);++i){
   if ((key[i]=='*')||(key[i]=='?')) wildcard = 1;
 }
 valstr[0] = '\0';

 if (wildcard==0){
   ln = HeadKeyValue;
   while (ln->next != NULL){
     ln = ln->next;
     if (strcmp(key,ln->key)==0){
       sprintf(valstr,"%s",ln->value); return(1);
     }
   }
   sprintf(valstr,"---"); return(1);
 }
 else if (wildcard==1){
   ln = HeadKeyValue;
   while (ln->next != NULL){
     ln = ln->next;
     if (wildcard_match(key,ln->key)==1){
       if (valstr[0]!='\0') strcat(valstr,":");
       strcat(valstr,ln->value);
     }
   }
   if (valstr[0]=='\0') sprintf(valstr,"---"); 
   return(1); 
 }
  
 return(0);

} /* end of Get_Value_from_KEYVALUE() */



int wildcard_match(pattern,string)
  char *pattern;
  char *string;
{
  int i,Lpat,Lstr;

  Lpat = strlen(pattern);
  Lstr = strlen(string);

  if (Lstr > Lpat) return(0);

  for (i=0;i<Lstr;++i){
    if ((pattern[i]!='?') && (pattern[i]!=string[i])) return(0);
  }

  if (Lpat>Lstr){
    for (i=Lstr;i<Lpat;++i){ if (pattern[i]!='?') return(0);}
  }

  return(1);
} /* end of wildcard_match() */



int substring_match(pattern,string)
  char *pattern;
  char *string;
{
  int i,o,Lpat,Lstr,match;

  /* printf("#pat '%s' str '%s' \n",pattern,string); */
  Lpat = strlen(pattern);
  Lstr = strlen(string);

  if (Lpat > Lstr) return(0);

  for (o=0;o<=(Lstr-Lpat);++o){
    match = 1;
    for (i=0;i<Lpat;++i){
      if ((pattern[i]!='?') && (pattern[i]!=string[i+o])) match = 0;
    }
    if (match==1) { /* printf("-->1 \n"); */ return(1);}
  }

  /* printf("-->0\n"); */
  return(0);
} /* end of substring_match() */








void Write_Property_for_One_Molecule(fpo,mol,num_mol,molname,num_libfile,offset_libfile)
  FILE *fpo;
  struct MOLECULE *mol;
  int  num_mol;    /* number of molecule */
  char *molname;
  int  num_libfile;    /* number of libfile */
  long offset_libfile; /* file offset for library file */
{
  int i;
  char value[MAX_LENGTH_LINE],spl;

  spl = '\t';
  
  /*
  printf("#Write_Property_for_One_Molecule(num_mol %d molname '%s')\n",num_mol,molname); fflush(stdout);
  */


  fprintf(fpo,"%d%c%s%c%d%c%ld%c%d%c%s"
    ,num_mol,spl,molname,spl,num_libfile,spl,offset_libfile,spl,mol->Nheavyatom,spl,mol->molform);
  for (i=0;i<PAR.Nproperty;++i){
    Get_Value_from_KEYVALUE(value,&(mol->HeadProperty),PAR.PROP_NAME[i]);
    fprintf(fpo,"%c%s",spl,value);
  }
  fprintf(fpo,"\n");

} /* end of Write_Property_for_One_Molecule() */


void Write_Headers_of_Property_File(fpo, Nlist,idirlib,HeadLibFile,lib_filetype)
  FILE *fpo; 
  int  Nlist;
  char *idirlib;
  struct LINENODE *HeadLibFile;
  char  lib_filetype; /* 'S'df,'2':mol2, 'P'db */
{
  struct LINENODE *ln;
  int n;
  fprintf(fpo,"#>>File_for_Molecular_Properties<<\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE    %s\n",PAR.START_DATE);
  if (idirlib[0] !='\0') fprintf(fpo,"#LIBRARY_DIRECTORY %s\n",idirlib);
  fprintf(fpo,"#LIBRARY_FILETYPE %c\n",lib_filetype);
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
/*
  fprintf(fpo,"#NPROPERTY %d\n",PAR.Nproperty);
  for (n=0;n<PAR.Nproperty;++n){
    fprintf(fpo,"#PROPERTY_%d %s\n",n,PAR.PROP_NAME[n]);
  }
 */
  fprintf(fpo,"#NCOLUMN %d\n",PAR.Nproperty+6);
  fprintf(fpo,"#COLUMN_1 num_mol\n");
  fprintf(fpo,"#COLUMN_2 compound_filename\n");
  fprintf(fpo,"#COLUMN_3 num_file\n");
  fprintf(fpo,"#COLUMN_4 file_offset\n");
  fprintf(fpo,"#COLUMN_5 Nheavyatom\n");
  fprintf(fpo,"#COLUMN_6 molecular_formula\n");
  for (n=0;n<PAR.Nproperty;++n){
    fprintf(fpo,"#COLUMN_%d %s\n",n+7,PAR.PROP_NAME[n]);
  }
/*
  fprintf(fpo,"#COLUMN [1:num_mol] [2:compound_filename] [3:num_file] [4:file_offset] [5:Nheavyatom] [6:molecular_formula]\n");
  fprintf(fpo,"#COLUMN ");
  for (n=0;n<PAR.Nproperty;++n){
    fprintf(fpo," [%d:%s]",7+n,PAR.PROP_NAME[n]);
  }
  fprintf(fpo,"\n");
 */

} /* end of Write_Headers_of_Property_File() */



int match_properties_and_patterns(HeadProperty,Nproperty_focus,prop_name,prop_min,prop_max,prop_str)
  struct KEYVALUE *HeadProperty;
  int    Nproperty_focus;
  char   prop_name[MAX_PROPERTY][64];  /* name of i-th property */
  float  prop_min[MAX_PROPERTY];       /* min value for i-th property */
  float  prop_max[MAX_PROPERTY];       /* max value for i-th property */
  char   prop_str[MAX_PROPERTY][64];   /* substring for i-th property */
{
  int p;
  char valstr[MAX_LENGTH_LINE],min_ok,max_ok,str_ok;
  float value;

  value = 0.0;
  for (p=0;p<Nproperty_focus;++p){
    Get_Value_from_KEYVALUE(valstr,HeadProperty,prop_name[p]);
    min_ok = max_ok = str_ok = 1;
    if ((prop_min[p]!=12345678.0)||(prop_max[p]!=12345678.0)) value = atof(valstr);
    if ((prop_min[p]!=12345678.0)&& (value < prop_min[p])) min_ok = 0;
    if ((prop_max[p]!=12345678.0)&& (value > prop_max[p])) max_ok = 0;
    if ((prop_str[p]!='\0')&& (substring_match(prop_str[p],valstr)==0)) str_ok = 0;
    if ((min_ok==0) || (max_ok==0) || (str_ok==0)) return(0);
  }

  return(1);

} /* end of match_properties_and_patterns() */

