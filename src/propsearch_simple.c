/*
  <prop_search_simple.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================



*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <sys/time.h>
#include <dirent.h>

#define MAX_LINE     2048 
#define MAX_FILESIZE 512
#define MAX_WORD     100 
#define MAX_QUERY 10 

static char LastModDate[] = "2017/08/24";

char *Get_Date_String()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec"};
 static char string[64];
 now_t = time(NULL);
 loc_t = localtime(&now_t);

 sprintf(string,"%s %d,%d %d:%d:%d",
  Mon[loc_t->tm_mon],loc_t->tm_mday,loc_t->tm_year+1900,
  loc_t->tm_hour,loc_t->tm_min,loc_t->tm_sec);
 return(string);

} /* end of Get_Date_String() */




void Split_to_Words(str,splsym,Nword,Wsta,Wend,Nwordmax,MultiSplitSym2Single)
 char *str;          /* Input String */
 char splsym;        /* Symbol for split */
 int *Nword;         /* Number of words  */
 int Wsta[];         /* Start point of str for a word */
 int Wend[];         /* End point of str for a word */
 int Nwordmax;       /* Maxium number of Nword  */
 char MultiSplitSym2Single;  /* 'T':multiple split sym --> regarded as single. */
{
 /* [Example]
                      11 
            012345678901 
     str = "//abc/d//ef/"
     splsym = '/'

   (1) MultiSplitSym2Single == 'S':
      Nword 3 
     (Wsta,Wend) = {(2,4),(6,6),(9,10)}

   (2) MultiSplitSym2Single != 'S':
      Nword 7 
     (Wsta,Wend) = {(0,-1),(1,0),(2,4),(6,6) (8,7),(9,10),(12,11)}

    **NOTE** if (Wsta>Wend) word[] = ''.

 */
 int i,L,sta;
 L = strlen(str);

 if (MultiSplitSym2Single=='S'){
   *Nword = 0; i = 0;
   while ((i<L)&&(*Nword < Nwordmax)){
    if (str[i]!=splsym){
       Wsta[*Nword] = i;
       while ((str[i]!=splsym)&&(i<=(L-1))) { ++i; }
       Wend[*Nword] = i-1;
       *Nword += 1;
     }
    ++i;
   }
 }
 else if (MultiSplitSym2Single!='S'){
   *Nword = 0; 
   sta = 0;
   for (i=0;i<=L;++i){
     if ((str[i]==splsym)||(str[i]=='\0')){
       if (*Nword<Nwordmax){
         Wsta[*Nword] = sta;
         Wend[*Nword] = i-1;
 /*
         if ((i-1)>sta){ Wend[*Nword] = i-1;}
                   else{ Wend[*Nword] = sta;}
  */
         sta = i + 1;
         *Nword += 1;
       }
     }
   }
 }
 
} /* end of Split_to_Words() */


void Get_Part_Of_Line(part,line,s,e)
  char *part;
  char *line;
  int  s,e;
{
 int i,E,L;
 L = strlen(line)-1;
 if (s>e) part[0] = '\0';
 if (line[L] == '\n') L -= 1;
 if (s<0) s = 0;
 if (e>L) E = L; else E = e;
 for (i=s;i<=E;++i) part[i-s] = line[i];
 part[E-s+1] = '\0';
} /* end of Get_Part_of_Line() */



int findsym_in_str(sym,str)
  char sym;
  char *str;
{
  int L,i;
  L = strlen(str);
  for (i=0;i<L;++i){
    if (str[i]==sym) return(1);
  }
  return(0);
} /* end of findsym_in_str() */

int substring_match(pat,str)
  char *pat;
  char *str; 
{
  int Lpat,Lstr,i,j;
  char hit;
  /* printf("#pattern '%s' str '%s' ",pat,str); */
  Lpat = strlen(pat);
  Lstr = strlen(str);
  if (Lpat>Lstr) return(0);
  for (i=0;i<=(Lstr-Lpat);++i){
    hit = 1;
    j = 0;
    while ((j<Lpat)&&(hit==1)){
      if (pat[j]!=str[i+j]) hit = 0;
      j += 1;
    }
    if (hit==1) {  return(1);}
  }
  return(0);
} /* end of substring_match() */



int main(argc,argv)
  int argc;
  char **argv;
{
  char line[MAX_LINE],valstr[MAX_LINE];
  char COMMAND[MAX_FILESIZE],buff[MAX_FILESIZE];
  char ipropfile[MAX_FILESIZE],ofname[MAX_FILESIZE];
  int    NQUERY;                    /* NUMBER OF QUERY */
  int    QUERY_NUM[MAX_QUERY];      /* property number of i-th query */
  char   QUERY_STR[MAX_QUERY][64];  /* string pattern for i-th query */
  float  QUERY_VAL[MAX_QUERY];      /* value  for i-th query */
  char   QUERY_OPR[MAX_QUERY];      /* operator for i-th query. 'U'pper, 'L'ower, 'S'ubstring */
  char   SplitSymType,MultiSplitSym2Single,hit; 
  int Nword, Wsta[MAX_WORD],Wend[MAX_WORD];
  int i,k,n,L;
  FILE *fpi,*fpo;
  int Ndata,Ndatahit,Ndatashow, Ndatashow_max;
  float valfloat;

  SplitSymType = 'T';
  MultiSplitSym2Single = 'F';
  sprintf(ofname,"-");
  NQUERY = 0;
  for (i=0;i<MAX_QUERY;++i){  
    QUERY_NUM[i] = -1;
    QUERY_VAL[i] = 0.0;
    QUERY_OPR[i] = 'S';
    QUERY_STR[i][0] = '\0';
  }

  Ndatashow_max = 1000;

  if (argc<2){

   /* 
    sprintf(line,"//abc/d//ef/");
    Split_to_Words(line,'/',&Nword,Wsta,Wend,MAX_WORD,'S');
    printf("line '%s' Nword %d\n",line,Nword);
    for (i=0;i<Nword;++i){ Get_Part_Of_Line(buff,line,Wsta[i],Wend[i]); printf("%d %d '%s'\n",Wsta[i],Wend[i],buff);} 
    
    sprintf(line,"//abc/d//ef/");
    Split_to_Words(line,'/',&Nword,Wsta,Wend,MAX_WORD,'M');
    printf("line '%s' Nword %d\n",line,Nword);
    for (i=0;i<Nword;++i){ Get_Part_Of_Line(buff,line,Wsta[i],Wend[i]); printf("%d %d '%s'\n",Wsta[i],Wend[i],buff);} 
    exit(1); 
   */

    printf("propsearch [input property_file] <options>\n");
    printf(" for searching tab-splited file.\n");
    printf(" coded by T.Kawabata. LastModified:%s\n",LastModDate);
    printf("<options>\n");
    printf(" -qn0  : property number of i-th query.[%d]\n",QUERY_NUM[0]);
    printf(" -qv0  : value/string for i-th query.[%s]\n",QUERY_STR[0]);
    printf(" -qo0  : operator for i-th query. 'U'pper, 'L'ower, 'S'ubstring[%c]\n",QUERY_OPR[0]);
    printf(" -of   : outputfile name [%s]\n",ofname);
    printf(" -mxsh : maximum number of showing hits(<0:not assign the limit) [%d]\n",Ndatashow_max); 
    printf(" -spl  : spliting symbol. 'T'ab,'S'pace [%c]\n",SplitSymType); 
    printf(" -mspl : Multi split symbols are considered as a single symbol (T or F)[%c]\n",MultiSplitSym2Single); 
    exit(1);
  }
  sprintf(ipropfile,"%s",argv[1]);

  /*** [1] Read options ***/

  for (k=0;k<argc;++k){ 
    if (k>0) strcat(COMMAND," ");
    strcat(COMMAND,argv[k]);
   }

  k = 0;
  while (k<argc){
    if ((argv[k][0]=='-') && ((k+1)<argc)){
       L = strlen(argv[k]);
        if ((L==4)&&(argv[k][1]=='q')&&(argv[k][2]=='n')) {
       sprintf(buff,"%c",argv[k][3]); n = atoi(buff); 
      if (n>=NQUERY) NQUERY = n + 1; 
       ++k; QUERY_NUM[n] = atoi(argv[k]); }
   else if ((L==4)&&(argv[k][1]=='q')&&(argv[k][2]=='o')) {
       sprintf(buff,"%c",argv[k][3]); n = atoi(buff); ++k; QUERY_OPR[n] = argv[k][0]; }
   else if ((L==4)&&(argv[k][1]=='q')&&(argv[k][2]=='v')) {
       sprintf(buff,"%c",argv[k][3]); n = atoi(buff); ++k; sprintf(QUERY_STR[n],"%s",argv[k]); }
   else if ((L==3)&&(argv[k][1]=='o')&&(argv[k][2]=='f')) {
       ++k; sprintf(ofname,"%s",argv[k]); }
   else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='x')&&(argv[k][3]=='s')&&(argv[k][4]=='h')) {
       ++k; Ndatashow_max = atoi(argv[k]); }
   else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='p')&&(argv[k][3]=='l')) {
       ++k; SplitSymType = argv[k][0]; }
   else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='s')&&(argv[k][3]=='p')&&(argv[k][4]=='l')) {
       ++k; MultiSplitSym2Single = argv[k][0]; }
   else { printf("#ERROR:Can't understand option %s\n",argv[k]); exit(1);}
   }
   ++k;
  } /* while k */

  if (ofname[0]=='-') fpo = stdout;
  else{
    fpo = fopen(ofname,"w");
    printf("#write_search_result() -->'%s'\n",ofname);
  } 
  fprintf(fpo,"#COMMAND %s\n",COMMAND);
  fprintf(fpo,"#DATE    %s\n",Get_Date_String());

  fprintf(fpo,"#NQUERY %d\n",NQUERY);
  for (n=0;n<NQUERY;++n){
    fprintf(fpo,"#QUERY  %2d COLUMN %2d OPR %c VALUE '%s'\n",n,QUERY_NUM[n],QUERY_OPR[n],QUERY_STR[n]);
    if ((QUERY_OPR[n]=='U')||(QUERY_OPR[n]=='L')) QUERY_VAL[n] = atof(QUERY_STR[n]);
    QUERY_NUM[n] = QUERY_NUM[n] -1 ;
  }


  /*** [2] Searching ***/
  fpi = fopen(ipropfile,"r");
  if (fpi == NULL){
    printf("#ERROR:Can't open ipropfile '%s'\n",ipropfile);
    exit(1);
  }
  printf("#INPUT_FILENAME %s\n",ipropfile);

  Ndata = Ndatahit = Ndatashow = 0;
  while (feof(fpi)==0){ 
    fgets(line,MAX_LINE-1,fpi);
    L = strlen(line);
    if ((L>0)&&(line[L-1] == '\n')) {line[L-1] = '\0'; L = L -1;}
    if (strncmp(line,"#COLUMN",7)==0) fprintf(fpo,"%s\n",line);
    if ((line[0]!='#') && (strlen(line)>10)){
       Ndata += 1;
       /* printf("[%d]%s\n",Ndata,line); fflush(stdout); */
            if (SplitSymType=='T') Split_to_Words(line,'\t',&Nword,Wsta,Wend,MAX_WORD,MultiSplitSym2Single);
       else if (SplitSymType=='S') Split_to_Words(line,' ', &Nword,Wsta,Wend,MAX_WORD,MultiSplitSym2Single);
       hit = 1;
       for (n=0;n<NQUERY;++n){
         if (QUERY_NUM[n]>=Nword){ hit = 0;}
         else{
           Get_Part_Of_Line(valstr,line,Wsta[QUERY_NUM[n]],Wend[QUERY_NUM[n]]);
           if (QUERY_OPR[n] != 'S') { valfloat = atof(valstr); }
           if ((QUERY_OPR[n]=='L') && (valfloat<QUERY_VAL[n])) { hit = 0;}
           if ((QUERY_OPR[n]=='U') && (valfloat>QUERY_VAL[n])) { hit = 0;}
           if ((QUERY_OPR[n]=='S') && (substring_match(QUERY_STR[n],valstr)==0)) {hit = 0;}
         }
       }
       if (hit==1){
         Ndatahit += 1; 
         if ((Ndatashow_max<0)||(Ndatashow<Ndatashow_max)){
           fprintf(fpo,"%s\n",line);
           Ndatashow += 1; 
         }
       }
    }
  }
  fclose(fpi);
  fprintf(fpo,"#Ndata     %d\n",Ndata);
  fprintf(fpo,"#Ndatahit  %d\n",Ndatahit);
  fprintf(fpo,"#Ndatashow %d\n",Ndatashow);

  if ((Ndatashow_max>0)&&(Ndatahit>Ndatashow_max)){
    fprintf(fpo,"#Ndatahit %d is over Ndatashow_max %d\n",Ndatahit,Ndatashow_max);
  }
  if (ofname[0]!='-') fclose(fpo);

 return(1);

} /* end of main() */

