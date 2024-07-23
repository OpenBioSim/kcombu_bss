/*
 
 <propsearch.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#if defined(_WIN32) || defined(WIN32)
  #include "windows.h"
  #include "dirent_windows.h"
  #include "time_windows.h"
#else
  #include <sys/time.h>
  #include <dirent.h>
#endif

#define MAX_LINE     2500 
#define MAX_FILESIZE 512
#define MAX_WORD     100 
#define MAX_QUERY 20 
#define MAX_OR    10 

static char LastModDate[] = "2017/08/24";

struct OFFSET_NODE{
  long  offset;
  char  mark;         /* if mark==1, output this line */
  int   score;        /* such as number of Nhit_OR */
  struct OFFSET_NODE *next,*prev;
};



void add_offset_to_OFFSET_NODEs(offset,score,HeadOffset)
  long   offset; 
  int    score;
  struct OFFSET_NODE *HeadOffset; 
{
  struct OFFSET_NODE *nextN;
  
  nextN = HeadOffset->next;
  HeadOffset->next = (struct OFFSET_NODE *)malloc(sizeof(struct OFFSET_NODE));
  HeadOffset->next->next = nextN;
  HeadOffset->next->prev = HeadOffset;
  HeadOffset->next->offset = offset; 
  HeadOffset->next->mark   = 1;
  HeadOffset->next->score  = score;
  if (nextN!=NULL){
     /* HeadOffset->next->num = nextN->num + 1; */
     nextN->prev = HeadOffset->next;
  }
} /* end of Add_string_to_LINENODEs() */


void bubble_sort_OFFSET_NODEs_by_score(offset_head)
  struct OFFSET_NODE *offset_head;
{
  struct OFFSET_NODE *an,*bn;
  int Nchange;
  long  offset;
  char  mark; 
  int   score;

  /* printf("#bubble_sort_OFFSET_NODEs_by_score()\n"); */

  do{
    Nchange = 0;
    an = offset_head;
    while (an->next != NULL){
      an = an->next;
      if (an->next != NULL){
        bn = an->next;
        if (an->score < bn->score){
          Nchange += 1;
          offset = an->offset;
          mark   = an->mark;
          score  = an->score;
          an->offset = bn->offset;
          an->mark   = bn->mark;
          an->score  = bn->score;
          bn->offset = offset;
          bn->mark   = mark;
          bn->score  = score;
        }
      }
    }
  }while (Nchange>0);
} /* end of bubble_sort_OFFSET_NODEs_by_score(offset_head) */


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





int Get_Time_Second()
{
 time_t      now_t;
 struct tm  *loc_t;
 now_t = time(NULL);
 loc_t = localtime(&now_t);
 return(loc_t->tm_sec);
} /* end of Get_Time_Second() */


double Get_Time_in_Second_by_Double()
{
 struct timeval tv;
 gettimeofday(&tv,NULL);
 return(tv.tv_sec + (double)tv.tv_usec*1e-6);
} /* end of Get_Time_in_Second_by_Double() */



double Set_Comp_Time(double start_time_sec){
 double end_sec;
 end_sec = Get_Time_in_Second_by_Double();
 return(end_sec - start_time_sec);
}





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



int substring_match(pat,str,LowUppIgnore)
  char *pat;   /* pattern */
  char *str;   /* target string */
  char LowUppIgnore; /* If 'T', ignore lower and upper of alphabet.  */
      /* If (LowUppIgnore=='T'), pat[] should be  lower. For example, pat[] should be 'ABC', not 'aBc'. */
{
  int Lpat,Lstr,i,j,ij;
  char hit,s;
  /* printf("#pattern '%s' str '%s' ",pat,str); */
  Lpat = strlen(pat);
  Lstr = strlen(str);
  if (Lpat>Lstr) return(0);

  if (LowUppIgnore != 'T'){
    for (i=0;i<=(Lstr-Lpat);++i){
      hit = 1;
      j = 0;
      while ((j<Lpat)&&(hit==1)){
        if (pat[j]!=str[i+j]) hit = 0;
        j += 1;
      }
      if (hit==1) {  return(1);}
    }
  }
  else if (LowUppIgnore == 'T'){
    for (i=0;i<=(Lstr-Lpat);++i){
      hit = 1;
      j = 0;
      while ((j<Lpat)&&(hit==1)){
        ij = i + j;
        if ((str[ij] >='A') && (str[ij] <= 'Z')){ s = str[ij] + 'a' - 'A'; } else { s = str[ij]; }
        if (pat[j] != s) hit = 0;
        j += 1;
      }
      if (hit==1) {  return(1);}
    }
  }

  return(0);
} /* end of substring_match() */




int isfloat(str)
  char *str;
{
  int i,L;
  L = strlen(str);
  for (i=0;i<L;++i){
    if (isdigit(str[i])!=0) return(1); 
  }
  return(0);
} /* end of isfloat() */



int main(argc,argv)
  int argc;
  char **argv;
{
  char line[MAX_LINE],valstr[MAX_LINE];
  char COMMAND[MAX_FILESIZE],buff[MAX_FILESIZE];
  char ipropfile[MAX_FILESIZE],ofname[MAX_FILESIZE];
  int    NQUERY;                    /* NUMBER OF QUERY */
  int    QUERY_NUM[MAX_QUERY][MAX_OR];      /* property number of i-th query */
  char   QUERY_STR[MAX_QUERY][MAX_OR][64];  /* string pattern for i-th query */
  float  QUERY_VAL[MAX_QUERY][MAX_OR];      /* value  for i-th query */
  char   QUERY_OPR[MAX_QUERY][MAX_OR];      /* operator for i-th query. 'U'pper, 'L'ower, 'S'ubstring,'N'on-substring */
  char   QUERY_NOR[MAX_QUERY];              /* Number of OR queries (0,1,2,..) */
  int    max_QUERY_NUM;                    /* maximum property number  */
  
  char   SplitSymType,MultiSplitSym2Single,LowerUpperIgnore, SortByNhitOR,hit; 
  int  NDATASHOW_MAX,NRANDSELECT, RANDSEED;
  int Nword, Wsta[MAX_WORD],Wend[MAX_WORD];
  int i,k,n,L,o,Nhitquery,nhitor,NhitOR;
  FILE *fpi,*fpo;
  int Ndata,Ndatahit,Ndatashow, Ndata_select;
  float valfloat, Paccept, Prand;
  struct OFFSET_NODE HeadOffset,*node;
  long offset;
  double start_time_sec;

  SplitSymType = 'T';
  MultiSplitSym2Single = 'F';
  LowerUpperIgnore = 'F';
  SortByNhitOR = 'T';

  sprintf(ofname,"-");
  NQUERY = 0;
  max_QUERY_NUM = 0;
  for (n=0;n<MAX_QUERY;++n){  
    QUERY_NOR[n] = -1; 
    for (o=0;o<MAX_OR;++o){  
      QUERY_NUM[n][o] = -1;
      QUERY_VAL[n][o] = 0.0;
      QUERY_OPR[n][o] = '-';
      QUERY_STR[n][o][0] = '\0';
    }
  }

  NDATASHOW_MAX = 1000;
  HeadOffset.next = NULL;
  HeadOffset.prev = NULL;
  NRANDSELECT = -1;
  RANDSEED    = 0;


  if (argc<2){
    printf("propsearch [input property_file] <options>\n");
    printf(" for searching tab-splited file, and random selection\n");
    printf(" coded by T.Kawabata. LastModified:%s\n",LastModDate);
    printf("<options for query>\n");
    printf("*LOGIC: (q0||r0||s0||..||x0) && (q1||r1||s1||..||x1) && .. && (q%d||r%d||s%d||..||x%d).\n",MAX_QUERY-1,MAX_QUERY-1,MAX_QUERY-1,MAX_QUERY-1); 
    printf(" -qn0,-rn0,-sn0,..,-xn0: property number of 0-th query (1,2,...).[%d]\n",QUERY_NUM[0][0]);
    printf(" -qv0,-rv0,-sv0,..,-xv0: value/string for 0-th query.[%s]\n",QUERY_STR[0][0]);
    printf(" -qo0,-ro0,-so0,..,-xo0: operator for 0-th query. 'U'pper, 'L'ower, 'S'ubstring, 'N'on-substring[%c]\n",QUERY_OPR[0][0]);
    printf("<options for output>\n");
    printf(" -of   : output file name [%s]\n",ofname);
    printf(" -mxsh : maximum number of showing hits(<0:not assign the limit) [%d]\n",NDATASHOW_MAX); 
    printf(" -sort:  sort by NhitOR (T or F) [%c]\n",SortByNhitOR); 
    printf("<options for random selection>\n");
    printf(" -nrnd : number of random selections. <0:no random selection. [%d]\n",NRANDSELECT); 
    printf(" -rseed: seed number for random selection [%d]\n",RANDSEED); 
    printf("<options for input file parsing>\n");
    printf(" -spl  : spliting symbol. 'T'ab,'S'pace [%c]\n",SplitSymType); 
    printf(" -mspl : Multi split symbols are considered as a single symbol ([spl][spl]-->[spl]) (T or F)[%c]\n",MultiSplitSym2Single); 
    printf("<options for pattern matching>\n");
    printf(" -luig : Lower and Upper types of alphabet are ignored. (T or F) [%c]\n",LowerUpperIgnore); 
    exit(1);
  }
  sprintf(ipropfile,"%s",argv[1]);

/************************/
/*** [1] Read options ***/
/************************/

  for (k=0;k<argc;++k){ 
    if (k>0) strcat(COMMAND," ");
    strcat(COMMAND,argv[k]);
   }

  k = 0;
  while (k<argc){
    if ((argv[k][0]=='-') && ((k+1)<argc)){
       L = strlen(argv[k]);
   if (((L==4)||(L==5))&&((argv[k][2]=='n')||(argv[k][2]=='v')||(argv[k][2]=='o'))&&(isdigit(argv[k][3]))){
          if (argv[k][1]=='q') o = 0;
     else if (argv[k][1]=='r') o = 1;
     else if (argv[k][1]=='s') o = 2;
     else if (argv[k][1]=='t') o = 3;
     else if (argv[k][1]=='u') o = 4;
     else if (argv[k][1]=='v') o = 5;
     else if (argv[k][1]=='w') o = 6;
     else if (argv[k][1]=='x') o = 7;
     else if (argv[k][1]=='y') o = 8;
     else if (argv[k][1]=='z') o = 9;
     else {printf("#ERROR:Can't understand the option '%s'.\n",argv[k]); exit(1); }
    
     if (o>=MAX_OR){printf("#ERROR:%d is over MAX_OR %d.\n",o,MAX_OR); exit(1); }
 
     if (L==4) {sprintf(buff,"%c"  ,argv[k][3]); n = atoi(buff);}
     if (L==5) {sprintf(buff,"%c%c",argv[k][3],argv[k][4]); n = atoi(buff);}

     if (n>=MAX_QUERY){
       printf("#ERROR:NumOfQuery %d is over MAX_QUERY %d.\n",n,MAX_QUERY);
       exit(1);
     }

     if ((n+1)>NQUERY) NQUERY = n+1; 
 
     if (QUERY_NOR[n]<(o+1)) QUERY_NOR[n] = o+1;

     /* printf("#argv '%s' value '%s'\n",argv[k],argv[k+1]); */
          if (argv[k][2]=='n'){ QUERY_NUM[n][o] = atoi(argv[k+1]); } 
     else if (argv[k][2]=='v'){ sprintf(QUERY_STR[n][o],"%s",argv[k+1]); } 
     else if (argv[k][2]=='o'){ QUERY_OPR[n][o] = argv[k+1][0]; } 
     else printf("woops!!\n");
     ++k;
     /* printf("#n %d o %d NUM '%d' STR '%s' OPR '%c'\n",n,o,QUERY_NUM[n][o],QUERY_STR[n][o],QUERY_OPR[n][o]); */
   }
   else if ((L==3)&&(argv[k][1]=='o')&&(argv[k][2]=='f')) {
       ++k; sprintf(ofname,"%s",argv[k]); }
   else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='x')&&(argv[k][3]=='s')&&(argv[k][4]=='h')) {
       ++k; NDATASHOW_MAX = atoi(argv[k]); }
   else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='p')&&(argv[k][3]=='l')) {
       ++k; SplitSymType = argv[k][0]; }
   else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='s')&&(argv[k][3]=='p')&&(argv[k][4]=='l')) {
       ++k; MultiSplitSym2Single = argv[k][0]; }
   else if ((L==5)&&(argv[k][1]=='n')&&(argv[k][2]=='r')&&(argv[k][3]=='n')&&(argv[k][4]=='d')) {
       ++k; NRANDSELECT = atoi(argv[k]); }
   else if ((L==6)&&(argv[k][1]=='r')&&(argv[k][2]=='s')&&(argv[k][3]=='e')&&(argv[k][4]=='e')&&(argv[k][5]=='d')) {
       ++k; RANDSEED = atoi(argv[k]); }
   else if ((L==5)&&(argv[k][1]=='l')&&(argv[k][2]=='u')&&(argv[k][3]=='i')&&(argv[k][4]=='g')) {
       ++k; LowerUpperIgnore = argv[k][0]; }
   else if ((L==5)&&(argv[k][1]=='s')&&(argv[k][2]=='o')&&(argv[k][3]=='r')&&(argv[k][4]=='t')) {
       ++k; SortByNhitOR = argv[k][0]; }
   else { printf("#ERROR:Can't understand option %s\n",argv[k]); exit(1);}
   }
   ++k;
  } /* while k */

  if (RANDSEED<0) srand(Get_Time_Second());
            else  srand(RANDSEED);

  if (ofname[0]=='-') fpo = stdout;
  else{
    fpo = fopen(ofname,"w");
    printf("#write_search_result() -->'%s'\n",ofname);
  } 
  fprintf(fpo,"#COMMAND %s\n",COMMAND);
  start_time_sec = Get_Time_in_Second_by_Double();
  fprintf(fpo,"#START_DATE %s\n",Get_Date_String());

  fprintf(fpo,"#NQUERY %d\n",NQUERY);
  for (n=0;n<NQUERY;++n){
    for (o=0;o<QUERY_NOR[n];++o){
      if ((QUERY_NUM[n][o]<=0)||(QUERY_OPR[n][o]=='-')){
         printf("#ERROR:QUERY[%d][%d] :QUERY_NUM:%d QUERY_OPR:'%c'\n",n,o,QUERY_NUM[n][o],QUERY_OPR[n][o]);
         exit(1);
      }
      if ((QUERY_OPR[n][o]=='U')||(QUERY_OPR[n][o]=='L')) QUERY_VAL[n][o] = atof(QUERY_STR[n][o]);

      if ((LowerUpperIgnore=='T') && ((QUERY_OPR[n][o]=='S')||(QUERY_OPR[n][o]=='N'))){
        for (i=0;i<strlen(QUERY_STR[n][o]);++i){
          if ((QUERY_STR[n][o][i] >='A') && (QUERY_STR[n][o][i] <= 'Z')){ 
              QUERY_STR[n][o][i] = QUERY_STR[n][o][i] + 'a' - 'A'; 
          }
        }
      }

     fprintf(fpo,"#QUERY %2d OR %2d COLUMN %2d OPR %c VALUE '%s'\n",n,o,QUERY_NUM[n][o],QUERY_OPR[n][o],QUERY_STR[n][o]);
     QUERY_NUM[n][o] = QUERY_NUM[n][o] -1 ;
     if (QUERY_NUM[n][o]>max_QUERY_NUM) max_QUERY_NUM = QUERY_NUM[n][o];
    }
  }
  fprintf(fpo,"#NDATASHOW_MAX %d\n",NDATASHOW_MAX);
  fprintf(fpo,"#NRANDSELECT   %d\n",NRANDSELECT);
  fprintf(fpo,"#RANDSEED      %d\n",RANDSEED);

/*********************/
/*** [2] Searching ***/
/*********************/
  fpi = fopen(ipropfile,"r");
  if (fpi == NULL){
    printf("#ERROR:Can't open ipropfile '%s'\n",ipropfile);
    exit(1);
  }
  printf("#INPUT_FILENAME %s\n",ipropfile);

  Ndata = Ndatahit = Ndatashow = 0;
  while (feof(fpi)==0){ 
    offset = ftell(fpi);
    line[0] = '\0';
    fgets(line,MAX_LINE-1,fpi);
    L = strlen(line);
    if ((L>0)&&(line[L-1] == '\n')) {line[L-1] = '\0'; L = L -1;}
    if (strncmp(line,"#COLUMN",7)==0) fprintf(fpo,"%s\n",line);
    if ((line[0]!='#') && (strlen(line)>10)){
       /* printf("[%d]%s\n",Ndata,line); fflush(stdout); */
            if (SplitSymType=='T') Split_to_Words(line,'\t',&Nword,Wsta,Wend,MAX_WORD,MultiSplitSym2Single);
       else if (SplitSymType=='S') Split_to_Words(line,' ', &Nword,Wsta,Wend,MAX_WORD,MultiSplitSym2Single);
       if (Nword>=max_QUERY_NUM){
         Ndata += 1;
         Nhitquery = 0;
         NhitOR = 0;
         for (n=0;n<NQUERY;++n){
           hit = 0;
           nhitor = 0;
           for (o=0;o<QUERY_NOR[n];++o){
             if (QUERY_NUM[n][o]<Nword){
               Get_Part_Of_Line(valstr,line,Wsta[QUERY_NUM[n][o]],Wend[QUERY_NUM[n][o]]);
               if (isfloat(valstr)==1){
                  valfloat = atof(valstr); 
                 if ((QUERY_OPR[n][o]=='L') && (valfloat >= QUERY_VAL[n][o])) { nhitor += 1;}
                 if ((QUERY_OPR[n][o]=='U') && (valfloat <= QUERY_VAL[n][o])) { nhitor += 1;}
               }
               if ((QUERY_OPR[n][o]=='S') && (substring_match(QUERY_STR[n][o],valstr,LowerUpperIgnore)==1)) {nhitor += 1;}
               if ((QUERY_OPR[n][o]=='N') && (substring_match(QUERY_STR[n][o],valstr,LowerUpperIgnore)==0)) {nhitor += 1;}
             }
           } /* o */
           NhitOR += nhitor;
           if (nhitor>0) {Nhitquery += 1;}
         } /* n */

         if (Nhitquery==NQUERY){
           Ndatahit += 1; 
           add_offset_to_OFFSET_NODEs(offset,NhitOR,&HeadOffset);
         }
      }

    } /* line != '^#' */

  } /* while */

/*****************************************/
/** [3] Mark lines randomly (optional) ***/
/*****************************************/
  if (NRANDSELECT>0){
    Ndata_select = 0;
    node = &HeadOffset;
    k = 0;
    while (node->next != NULL){ 
      node = node->next;
      node->mark = 0;
      if (Ndata_select < NRANDSELECT){ 
        Paccept = (float)(NRANDSELECT-Ndata_select)/(float)(Ndatahit-k);
        Prand   = (float)rand()/(float)RAND_MAX;
        if (Prand<=Paccept){ 
           node->mark = 1;
           Ndata_select += 1;
        }
        k += 1;
      }
    }
  }

/**************************/
/** [4] Output Hit Line ***/
/**************************/
  if (SortByNhitOR=='T') bubble_sort_OFFSET_NODEs_by_score(&HeadOffset); 

  node = &HeadOffset;
  while (node->next != NULL){
    node = node->next;
    if ((node->mark==1)&&((NDATASHOW_MAX<0)||(Ndatashow<NDATASHOW_MAX))){
      fseek(fpi,node->offset,SEEK_SET);
      line[0] = '\0';
      fgets(line,MAX_LINE-1,fpi);
      if (line[strlen(line)-1]=='\n') line[strlen(line)-1] = '\0';
      /* fprintf(fpo,"%s",line); */
      if (SplitSymType=='T') fprintf(fpo,"%s\t%d\n",line,node->score); 
      if (SplitSymType=='S') fprintf(fpo,"%s %d\n",line,node->score); 
      Ndatashow += 1; 
    }
  }

  fprintf(fpo,"#Ndata       %d\n",Ndata);
  fprintf(fpo,"#Ndatahit    %d\n",Ndatahit);
  fprintf(fpo,"#Ndatashow   %d\n",Ndatashow);
  if ((NRANDSELECT<0)&&(NDATASHOW_MAX>0)&&(Ndatahit>NDATASHOW_MAX)){
    fprintf(fpo,"#Ndatahit %d is over NDATASHOW_MAX %d\n",Ndatahit,NDATASHOW_MAX);
  }
  
  fprintf(fpo,"#END_DATE %s\n",Get_Date_String());
  fprintf(fpo,"#COMP_TIME_SEC %lf\n",Get_Time_in_Second_by_Double() - start_time_sec);


  fclose(fpi);
  if (ofname[0]!='-') fclose(fpo);
  return(1);

} /* end of main() */

