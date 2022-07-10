/*

 <ReadOpt.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 for reading options from arguments or the option file.

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "ReadOpt.h"

#define MAX_LEN 1024

/*** FUNCTIONS (GLOBAL) ***/
void Show_Options();
void Read_Options_From_Arguments();
void Read_Options_From_File();
void Make_COMMAND_String_From_OptVal_List();

/*** FUNCTIONS (LOCAL) ***/
static struct OPTION_VALUE *Add_OPTION_VALUE_ListEnd();


void Show_Options(OptValHead)
 struct OPTION_VALUE *OptValHead;
{
 struct OPTION_VALUE *ov;
 printf("#Show_Options()\n");
 ov = OptValHead;
 while (ov->next != NULL){
   ov = ov->next;
   printf("#OPTION '%s' VALUE '%s'\n",ov->opt,ov->val);
 }
 printf("#End_Show_Options()\n");
} /* end of Show_Options() */


void Read_Options_From_Arguments(argc,argv,command,OptValHead)
 int argc;   /* (input) */
 char **argv; /* (input) */
 char   *command;   /* (to be assigned) */
 struct OPTION_VALUE *OptValHead; /* (to be assigned) */
{
 int k,L,i;
 char opt_str[MAX_LEN],val_str[MAX_LEN];

 sprintf(command,"%s",argv[0]); 
 for (k=1;k<argc;++k){
   strcat(command," ");
   strcat(command,argv[k]);
 }

 k = 1;
 while (k<argc){

    if ((argv[k][0]=='-')&&(k<(argc-1))){
      L = strlen(argv[k]);
      for (i=1;i<L;++i){ opt_str[i-1] = argv[k][i];}
      opt_str[L-1] = '\0'; 
      sprintf(val_str,"%s",argv[k+1]);
      /* printf("Option '%s' Value '%s'\n",opt_str,val_str);  */
      Add_OPTION_VALUE_ListEnd(OptValHead,opt_str,val_str);
      k += 2; 
   }
   else{
     k += 1;
   }

 } 

  /* printf("#COMMAND '%s'\n",command); */

} /* end of Read_Options_From_Arguments() */







void Read_Options_From_File(ifname,OptValHead)
 char *ifname;
 struct OPTION_VALUE *OptValHead;
{
 int L,i;
 FILE *fp;
 char line[MAX_LEN],buff[MAX_LEN],opt_str[MAX_LEN],val_str[MAX_LEN];

/*
>> FILE EXAMPLE <<
#option file for "gmfit"
-igi 1finAB_2.gau
-ig0 1finA_4.gau
-ig1 1finB_4.gau
-nm 3
-MI T
-Nmi 100
-hom F
*/

 printf("#Read_Options_From_File('%s')\n",ifname);
 fp = fopen(ifname,"r");
 if (fp==NULL){printf("#ERROR:Can't open option_script file '%s'\n",ifname); exit(1);}

/* Add_OPTION_VALUE_ListEnd(OptValHead,"script",ifname); */
 
 while (feof(fp)==0){
   line[0] = '\0'; 
   fgets(line,511,fp);
   /* printf("%s\n",line);  */

   if ((line[0]!='#')&&(strlen(line)>3)){
      sscanf(line,"%s %s",buff,val_str); 
      L = strlen(buff);
      for (i=1;i<L;++i) opt_str[i-1] = buff[i];
      opt_str[L-1] = '\0'; 
      printf("#Option '%s' Value '%s'\n",opt_str,val_str); 
      Add_OPTION_VALUE_ListEnd(OptValHead,opt_str,val_str);
   }
 } 

 fclose(fp);

} /* end of Read_Options_From_File() */





void Make_COMMAND_String_From_OptVal_List(command,program_name,OptValHead)
 char   *command;   /* to be assigned */
 char   *program_name;  /* (input) */ 
 struct OPTION_VALUE *OptValHead; /* (input) */
{
 struct OPTION_VALUE *ov;
 int L;
 /* 
 printf("#Make_COMMAND_String_From_OptVal_List(OptValHead)\n");
 */

 /* (1) Check the length of COMMAND string */
 L = 20;
 ov = OptValHead;
 while (ov->next != NULL){
  ov = ov->next;
  L += (strlen(ov->opt) + strlen(ov->val) + 2); 
 } 
  printf("#Length of COMMAND string : %d\n",L); 
 /* 
 PAR.COMMAND = (char *)malloc(sizeof(char)*L);
 */ 
 /* (2) Make COMMAND string */
 sprintf(command,"%s", program_name);
 ov = OptValHead;
 while (ov->next != NULL){
   ov = ov->next;
   printf("#option '%s' value '%s\n",ov->opt,ov->val); 
   strcat(command," -");
   strcat(command,ov->opt);
   strcat(command," ");
   strcat(command,ov->val);
 }

/*
 while (ov->prev != OptValHead){
  ov = ov->prev;
  printf("#option '%s' value '%s\n",ov->opt,ov->val); 
  strcat(PAR.COMMAND," -");
  strcat(PAR.COMMAND,ov->opt);
  strcat(PAR.COMMAND," ");
  strcat(PAR.COMMAND,ov->val);
 }
 */
 
 printf("#COMMAND %s\n",command);
} /* end of Make_COMMAND_String_From_OptVal_List() */


struct OPTION_VALUE *Add_OPTION_VALUE_ListEnd(OptValHead,opt_str,val_str)
 struct OPTION_VALUE *OptValHead;
 char *opt_str,*val_str;
{
 struct OPTION_VALUE *new,*last;

 last = OptValHead;
 while (last->next !=NULL){ last = last->next; }

 last->next = (struct OPTION_VALUE*)malloc(sizeof(struct OPTION_VALUE));
 new = last->next;
 new->next = NULL;
 new->prev = last;
 new->opt  = (char *)malloc(sizeof(char)*(strlen(opt_str)+1));
 new->val  = (char *)malloc(sizeof(char)*(strlen(val_str)+1));
 sprintf(new->opt,"%s",opt_str);
 sprintf(new->val,"%s",val_str);
 new->Lopt = strlen(opt_str);
 new->accept = 0;
 return(new);
} /* end of Add_OPTION_VALUE_ListEnd() */







