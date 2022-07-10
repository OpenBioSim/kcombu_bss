/*
 
 <SmiComp.c>

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

#define HASH_SIZE_MALLOC  10000
#define MAX_LINE          2500
#define MAX_FILE_LEN       512 

static char LastModDate[] = "2016/11/17";
static char COMMAND[MAX_FILE_LEN];
static char START_DATE[128];

struct SIM_MOLECULE{
  struct MOLECULE     *mol;  /* pointer of MOLECULE for the similar molecule */
  struct SIM_MOLECULE *next; /* pointer of next SIM_MOLECULE */ 
};


struct MOLECULE{
  char *id;           /* ideintifier string of MOLECULE (malloc later)*/ 
  char *smiles;       /* SMILES string (malloc later)*/ 
  int  Lsmiles;       /* length of smiles */ 
  char mark;          /* mark for various purposes */
  struct MOLECULE *next; /* pointer of next MOLECULE */
  struct MOLECULE *prev; /* pointer of previous MOLECULE */
  struct SIM_MOLECULE simhead; /* head of SIM_MOLECULE */
  int    Nsimmol;     /* number of similar molecule */
};


int    HASH_SIZE;
char   HASH_TYPE;
struct MOLECULE  HEAD_MOL_HASH[HASH_SIZE_MALLOC];
int               NUM_MOL_HASH[HASH_SIZE_MALLOC];

/*** FUNCTIONS ***/
int  hash_func();
void add_new_MOLECULE_after_the_head();
void add_new_SIM_MOLECULE_after_the_head();
void Free_MOLECULE();
int  read_SMILES_file();
void split_string_by_first_appearing_symbol();
char *Get_Date_String();
int element_number();



int hash_func(key)
  char *key;
{
  int i,s;
  s = 0;
  /** 'L':Length of key **/
  if (HASH_TYPE=='L'){
    s = strlen(key); 
  }
  /** 'S': Sum_of_each_char */ 
  else if (HASH_TYPE=='S'){
    for (i=0;i<strlen(key);++i){
      s += key[i];  
    }
  }
  else if (HASH_TYPE=='E'){
    for (i=0;i<strlen(key);++i){
      s += element_number(key[i]);
    }
  }
  /** 'C':carbon-binary **/
  else if (HASH_TYPE=='C'){
    for (i=0;i<strlen(key);++i){
      if ((key[i]=='C')||(key[i]=='c')) {s += 2*s + 1;}
                  else {s += 2*s;}
      if (s>HASH_SIZE){
        s = s%HASH_SIZE; 
      }
    }
  }
  else{
    s = strlen(key); 
  }
 
  /* if ((s%HASH_SIZE)==42){printf("'%s' %d %d\n",key,s,s%HASH_SIZE);} */

  return(s % HASH_SIZE);
}

int element_number(x)
  char x;
{
       if (x=='C') return(1);
  else if (x=='c') return(2);
  else if (x=='N') return(3);
  else if (x=='n') return(4);
  else if (x=='O') return(5);
  else if (x=='o') return(6);
  else if (x=='P') return(7);
  else if (x=='S') return(8);
  else if (x=='(') return(9);
  else if (x==')') return(10);
  else if (x=='[') return(11);
  else if (x==']') return(12);
  else if (x=='1') return(13);
  else if (x=='2') return(14);
  else if (x=='=') return(15);
  else if (x=='#') return(16);
  else return(17);
}




void add_new_MOLECULE_after_the_head(head,new_smiles,new_id)
  struct MOLECULE *head;
  char *new_smiles;
  char *new_id;
{
  struct MOLECULE *next,*new;

  next = head->next;
  head->next = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));
  new = head->next;
  head->next->prev = head;
  if (next != NULL) next->prev = head->next;
  new->next = next;
  new->smiles = (char *)malloc(sizeof(char)*(strlen(new_smiles)+1));
  new->id     = (char *)malloc(sizeof(char)*(strlen(new_id)+1));
  sprintf(new->smiles,"%s",new_smiles); 
  sprintf(new->id,"%s",new_id); 
  new->simhead.mol  = NULL;
  new->simhead.next = NULL;
  new->Lsmiles = strlen(new_smiles);
  new->mark = 0;
  new->Nsimmol = 0;
} /* end of add_new_MOLECULE_after_the_head() */


void add_new_SIM_MOLECULE_after_the_head(head,new_mol)
  struct SIM_MOLECULE *head;
  struct MOLECULE *new_mol;
{
  struct SIM_MOLECULE *next,*new;

  next = head->next;
  head->next = (struct SIM_MOLECULE*)malloc(sizeof(struct SIM_MOLECULE));
  new = head->next;
  new->next = next;
  new->mol  = new_mol;

} /* end of add_new_SIM_MOLECULE_after_the_head() */


void Free_MOLECULE(head)
  struct MOLECULE *head;
{
  struct MOLECULE *m;
  m = head; 
  while (m->next != NULL) m = m->next;
  while (m->prev != NULL){
    m = m->prev;
    if (m->next != NULL){
      free(m->next->smiles); 
      free(m->next->id); 
      free(m->next); 
    }
 }
 head->next = NULL;

} /* end of Free_MOLECULE() */


int read_SMILES_file(ifname,HeadMolSmilesLen,DIVbunshi,DIVbunbo)
  char *ifname;
  struct MOLECULE HeadMolSmilesLen[HASH_SIZE];
  int  DIVbunshi,DIVbunbo;
  /*
    read only if (Lsmiles%DIVbunbo)==DIVbunshi.
  */
{
  FILE *fp;
  char line[MAX_LINE],smiles[MAX_LINE],tail[MAX_LINE],tail2[MAX_LINE],id[100];
  int i,Lline,hash_value,Nmolecule,Ncum;

  printf("#read_SMILES_file('%s')\n",ifname);
  fp = fopen(ifname,"r");
  if (fp==NULL){ printf("#ERROR:Can't open smiles file '%s'\n",ifname); exit(1);}
  for (i=0;i<HASH_SIZE;++i){
    NUM_MOL_HASH[i] = 0;
    HEAD_MOL_HASH[i].next    = NULL;
    HEAD_MOL_HASH[i].Lsmiles = 0;
  }

  Nmolecule = 0; 
  while (feof(fp)==0){
    line[0] = '\0';
    fgets(line,MAX_LINE-1,fp);
    Lline = strlen(line);
    if (line[Lline-1]== '\n') {line[Lline-1] = '\0'; Lline -= 1;} 
    if ((line[0]!='#') && (Lline>5)){
      split_string_by_first_appearing_symbol(smiles,tail,line,' ');
      split_string_by_first_appearing_symbol(id,tail2,tail,' ');
      /* printf("line '%s' smiles '%s' tail '%s' id '%s'\n",line,smiles,tail,id); */
      /* hash_value = strlen(smiles); */
      hash_value = hash_func(smiles); 
      if (hash_value>=HASH_SIZE){ hash_value = HASH_SIZE-1;} 

      if ((hash_value%DIVbunbo)==DIVbunshi){
        add_new_MOLECULE_after_the_head(&HeadMolSmilesLen[hash_value],smiles,id);
        NUM_MOL_HASH[hash_value] += 1;
        Nmolecule += 1; 
      }
    }
  }

  fclose(fp);

  Ncum = 0;
  /*
  for (i=0;i<HASH_SIZE;++i){
    printf("%d %d %d %f\n",i,NUM_MOL_HASH[i],Ncum,(float)Ncum/(float)Nmolecule);
    Ncum += NUM_MOL_HASH[i];
  } 
  */

  return(Nmolecule);
} /* end of read_SMILES_file() */




void split_string_by_first_appearing_symbol(head,tail,line,split_symbol)
  char *head,*tail,*line;
  char split_symbol;
/*
 >> example <<
  line = "NC(=O)c1ncc(n(c1)) D00144" -> head = "NC(=O)c1ncc(n(c1))", tail = "D00144"
  line = 'AB CD E'-> head = 'AB', tail = 'CD E'
  line = 'ABCDE'  -> head = 'ABCDE', tail = ''
*/
{
  int i,isp,Lline;

  Lline = strlen(line);
  head[0] = '\0';
  tail[0] = '\0';
 
  isp = -1;
  i = 0;
  while ((i<Lline)&&(isp==-1)){
    if (line[i]==split_symbol) isp = i;
    i += 1;
  } 
  if (isp>=0){
    for (i=0;i<isp;++i) head[i] = line[i];
    head[i] = '\0';
    for (i=isp+1;i<Lline;++i) tail[i-isp-1] = line[i];
    tail[Lline-isp-1] = '\0';
  }
  else{
    sprintf(head,"%s",line);
  }
} /* end of split_string_by_first_appearing_symbol() */



char *Get_Date_String()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec"};
 static char string[64];
 now_t = time(NULL);
 loc_t = localtime(&now_t);

 sprintf(string,"%s %d,%d %d:%d:%d", Mon[loc_t->tm_mon],loc_t->tm_mday,loc_t->tm_year+1900,
  loc_t->tm_hour,loc_t->tm_min,loc_t->tm_sec);
 return(string);

} /* end of Get_Date_String() */



void mk_string_check_hash_distribution(string)
  char *string;
{
 int k,minNhash,maxNhash,val_minNhash,val_maxNhash;
 int min_hash_appear, max_hash_appear;
 double S,SS,M,SD;


 minNhash = maxNhash = -1;
 min_hash_appear = max_hash_appear = -1;
 S = SS = 0.0;

 for (k=0;k<HASH_SIZE;++k){
   S  += NUM_MOL_HASH[k];
   SS += NUM_MOL_HASH[k] * NUM_MOL_HASH[k];
   if (NUM_MOL_HASH[k]>0){
     if (min_hash_appear<0){ min_hash_appear = k;}
     max_hash_appear = k;
   }
   if ((minNhash <0) || (NUM_MOL_HASH[k]<minNhash)){
     minNhash = NUM_MOL_HASH[k];
     val_minNhash = k;
   }
   if ((maxNhash <0) || (NUM_MOL_HASH[k]>maxNhash)){
     maxNhash = NUM_MOL_HASH[k];
     val_maxNhash = k;
   }
 }

  M = S/HASH_SIZE;
  SD = sqrt(SS/HASH_SIZE - M*M);
  sprintf(string,"Ndata %.0f HASH_RANGE %d...%d NUM_IN_EACH_HASH_BIN M %lf SD %lf min %d [%d] max %d [%d]",S,min_hash_appear, max_hash_appear,M,SD,minNhash,val_minNhash,maxNhash,val_maxNhash);

} /* end of mk_string_check_hash_distribution() */





int main(argc,argv)
  int argc;
  char **argv;
{
  int k,L,Nmolecule,NmoleculeA,NmoleculeB,Nidpair,Lline,hash_val;
  int DIVbunshi,DIVbunbo; 
  char MODE, ismifile[MAX_FILE_LEN], ismifileA[MAX_FILE_LEN], ismifileB[MAX_FILE_LEN],init,hit;
  char oidfile[MAX_FILE_LEN],oidpairfile[MAX_FILE_LEN],oclusfile[MAX_FILE_LEN],ouniqfile[MAX_FILE_LEN];
  char line[MAX_LINE],smiles[MAX_LINE],tail[MAX_LINE],tail2[MAX_LINE],id[MAX_FILE_LEN];
  char str_hash_dist[MAX_LINE];
  struct MOLECULE *mA,*mB;
  struct SIM_MOLECULE *sm; 
  FILE *fpi,*fpo,*fpop,*fpou;

  MODE = 'A';
  ismifile[0] = ismifileA[0] = ismifileB[0]  = '\0';
  oidfile[0] = oidpairfile[0] = oclusfile[0] = ouniqfile[0] = '\0';
  sprintf(oidfile,"same_smiles.out");
  DIVbunshi = 0; DIVbunbo = 1;
  HASH_TYPE = 'C';
  HASH_SIZE = HASH_SIZE_MALLOC;
 
  sprintf(START_DATE,"%s",Get_Date_String());
  if (argc<2){ 
    printf("SmiComp <options>\n");
    printf(" to compare SMILES strings by the exact match, using hash table.\n");
    printf(" coded by T.Kawabata. LastModified: %s\n",LastModDate);
    printf("<simple usage>\n");
    printf("  SmiComp -M A -if  [SMILES_file]\n");
    printf("  SmiComp -M 2 -ifA [SMILES_fileA(smaller)] -ifB [SMILES_fileB(larger)]\n");
    printf("<options>\n");
    printf(" -M   : Mode. 'A'll-vs-all, '2':libA vs libB. [%c]\n",MODE);
    printf("<for both MODE 'A' and MODE '2'>\n");
    printf(" -oid : out file for same-smiles IDs for each ID[%s]\n",oidfile);
    printf(" -opr : out file for same-smiles ID pairs.[%s]\n",oidpairfile);
    printf(" -hashtype: hash_function type. 'L'ength, 'S'um_of_each_char, 'C':arbon-binary [%c]\n",HASH_TYPE);
    printf("          : (Nsmi=1 million and hashsize=10000, 'C' shows the best performance.)\n");
    printf(" -hashsize: size of hash_table (<=HASH_SIZE_MALLOC %d) [%d]\n",HASH_SIZE_MALLOC,HASH_SIZE);
    printf(" -div : Job division. {0..Ndivision-1}/(Ndivision). [%d/%d]\n",DIVbunshi,DIVbunbo);
    printf("<for MODE 'A'>\n");
    printf(" -if  : input SMILES file for library (stored in Memory)[%s]\n",ismifile);
    printf(" -ocl : out file for clusters for same SMILES [%s]\n",oclusfile);
    printf(" -ouq : out file for unique SMILES [%s]\n",ouniqfile);
    printf("<for MODE '2'>\n");
    printf(" -ifA : input SMILES file for library A[%s]\n",ismifileA);
    printf(" -ifB : input SMILES file for library B (stored in Memory)[%s]\n",ismifileB);
    exit(1);
  }

 COMMAND[0] = '\0';
 for (k=0;k<argc;++k){ 
   if (k>0) strcat(COMMAND," "); strcat(COMMAND,argv[k]);
 }
 k = 1;
 while (k<argc){
   if (argv[k][0]=='-'){
     L = strlen(argv[k]);
          if ((L==3)&&(argv[k][1]=='i')&&(argv[k][2]=='f')){ ++k; sprintf(ismifile,"%s",argv[k]);}
     else if ((L==2)&&(argv[k][1]=='M')){ ++k; MODE = argv[k][0];}
     else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='f')&&(argv[k][3]=='A')){ ++k; sprintf(ismifileA,"%s",argv[k]);}
     else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='f')&&(argv[k][3]=='B')){ ++k; sprintf(ismifileB,"%s",argv[k]);}
     else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='i')&&(argv[k][3]=='d')){ ++k; sprintf(oidfile,"%s",argv[k]);}
     else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='p')&&(argv[k][3]=='r')){ ++k; sprintf(oidpairfile,"%s",argv[k]);}
     else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='c')&&(argv[k][3]=='l')){ ++k; sprintf(oclusfile,"%s",argv[k]);}
     else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='u')&&(argv[k][3]=='q')){ ++k; sprintf(ouniqfile,"%s",argv[k]);}
     else if ((L==9)&&(argv[k][1]=='h')&&(argv[k][2]=='a')&&(argv[k][3]=='s')&&(argv[k][4]=='h')&&
                      (argv[k][5]=='t')&&(argv[k][6]=='y')&&(argv[k][7]=='p')&&(argv[k][8]=='e')) { ++k; HASH_TYPE = argv[k][0];}
     else if ((L==9)&&(argv[k][1]=='h')&&(argv[k][2]=='a')&&(argv[k][3]=='s')&&(argv[k][4]=='h')&&
                      (argv[k][5]=='s')&&(argv[k][6]=='i')&&(argv[k][7]=='z')&&(argv[k][8]=='e')) { ++k; HASH_SIZE = atoi(argv[k]);}
     else if ((L==4)&&(argv[k][1]=='d')&&(argv[k][2]=='i')&&(argv[k][3]=='v')){
       ++k;
       split_string_by_first_appearing_symbol(line,smiles,argv[k],'/');
       DIVbunshi = atoi(line);
       DIVbunbo  = atoi(smiles);
       printf("#DIVbunshi %d DIVbunbo %d\n",DIVbunshi,DIVbunbo);
       }
     else { printf("#Can't understand the option '%s'\n",argv[k]); exit(1);}
   }
   k += 1;
 }

 

 /***************************************/
 /*** MODE 'A': all-vs-all comparison ***/
 /***************************************/
 if (MODE == 'A'){

   Nmolecule = read_SMILES_file(ismifile,HEAD_MOL_HASH,DIVbunshi,DIVbunbo);

   printf("#Nmolecule %d HASH_SIZE %d DIV %d/%d Nmol_cal %d\n",Nmolecule,HASH_SIZE,DIVbunshi,DIVbunbo,L); fflush(stdout);
   mk_string_check_hash_distribution(str_hash_dist);
   printf("#HASH_DISTRIBUTION %s\n",str_hash_dist);
   /*** ALL-VS-ALL COMPARISON ***/
   fpop = NULL; fpo = NULL;
   if (oidpairfile[0] != '\0'){
     if (oidpairfile[0]=='-') fpop = stdout;
     else { fpop = fopen(oidpairfile,"w"); 
     }
   }

 
   for (k=0;k<HASH_SIZE;++k){
     if ((k%DIVbunbo)==DIVbunshi){
     mA = &(HEAD_MOL_HASH[k]);
       while (mA->next != NULL){
         mA = mA->next;
         if (mA->next != NULL){
           mB = mA;
           while (mB->next != NULL){
             mB = mB->next;
             if (strcmp(mA->smiles,mB->smiles)==0){
               if (fpop!=NULL){
                 fprintf(fpop,"%s %s %s\n",mA->id,mB->id,mA->smiles);
               }
               add_new_SIM_MOLECULE_after_the_head(&(mA->simhead),mB);
               add_new_SIM_MOLECULE_after_the_head(&(mB->simhead),mA);
               mA->Nsimmol = mA->Nsimmol + 1;
               mB->Nsimmol = mB->Nsimmol + 1;
               Nidpair += 1;
               if ((Nidpair%1000)==0){ printf("%s %s %s\n",mA->id,mB->id,mA->smiles); fflush(stdout);}
             }
           } /* mB */
         }
       } /* mA */
     }
   }
   if (fpop !=NULL){
     printf("#identical pairs -->'%s'\n",oidpairfile);
     fclose(fpop);
   }
 


 
   /*** output same-smiles IDs list for each ID ***/
    if (oidfile[0] != '\0'){
      fpo = fopen(oidfile,"w");
      printf("#output_same_ids_file()-->'%s'\n",oidfile);
      fprintf(fpo,"#>>> SAME_SMILES_IDS_LIST -->'%s' <<<\n",oidfile);
      fprintf(fpo,"#COMMAND '%s'\n",COMMAND);
      fprintf(fpo,"#HASH_TYPE %c\n",HASH_TYPE);
      fprintf(fpo,"#HASH_SIZE %d\n",HASH_SIZE);
      fprintf(fpo,"#HASH_DISTRIBUTION %s\n",str_hash_dist);
      fprintf(fpo,"#START_DATE '%s'\n",START_DATE);
      fprintf(fpo,"#END_DATE   '%s'\n",Get_Date_String());
      for (k=0;k<HASH_SIZE;++k){
         mA = &(HEAD_MOL_HASH[k]);
         while (mA->next != NULL){
           mA = mA->next;
           if (mA->simhead.next != NULL){
             fprintf(fpo,">%s\n",mA->id);
             sm = &(mA->simhead);
             while (sm->next != NULL){
               sm = sm->next;
               fprintf(fpo,"%s ",sm->mol->id);
             }         
             fprintf(fpo,"\n");
           }
        } /* mA */
     } /* k */
    fclose(fpo); 
   } 


   /*** output cluster of same-SMILES  ***/
    if (oclusfile[0] != '\0'){
      fpo = fopen(oclusfile,"w");
      printf("#output_cluster_with_same_SMILES()-->'%s'\n",oclusfile);
      fprintf(fpo,"#COMMAND '%s'\n",COMMAND);
      fprintf(fpo,"#START_DATE '%s'\n",START_DATE);
      fprintf(fpo,"#END_DATE   '%s'\n",Get_Date_String());
      for (k=0;k<HASH_SIZE;++k){
         mA = &(HEAD_MOL_HASH[k]);
         while (mA->next != NULL){
           mA = mA->next;
           if (mA->mark == 0){
             mA->mark = 1;
             fprintf(fpo,">%s %d\n",mA->id,mA->Nsimmol+1);
             fprintf(fpo,"%s",mA->id);
             if (mA->simhead.next != NULL){
               sm = &(mA->simhead);
               while (sm->next != NULL){
                 sm = sm->next;
                 sm->mol->mark = 1;
                 fprintf(fpo," %s",sm->mol->id);
               }       
             }
             fprintf(fpo,"\n");
          }
        } /* mA */
     } /* k */
    fclose(fpo); 
   } 

   /*** output unique SMILES  ***/
    if (ouniqfile[0] != '\0'){
      fpo = fopen(ouniqfile,"w");
      printf("#output_unique_SMILES()-->'%s'\n",ouniqfile);
      fprintf(fpo,"#COMMAND '%s'\n",COMMAND);
      fprintf(fpo,"#START_DATE '%s'\n",START_DATE);
      fprintf(fpo,"#END_DATE   '%s'\n",Get_Date_String());

      for (k=0;k<HASH_SIZE;++k){
         mA = &(HEAD_MOL_HASH[k]);
         while (mA->next != NULL){
           mA = mA->next;
           mA->mark = 0;
         }
       }
      for (k=0;k<HASH_SIZE;++k){
         mA = &(HEAD_MOL_HASH[k]);
         while (mA->next != NULL){
           mA = mA->next;
           if (mA->mark == 0){
             fprintf(fpo,"%s %s",mA->smiles,mA->id);
             mA->mark = 1;
             if (mA->simhead.next != NULL){
               sm = &(mA->simhead);
               while (sm->next != NULL){
                 sm = sm->next;
                 sm->mol->mark = 1;
                 fprintf(fpo," %s",sm->mol->id);
               }       
             }
             fprintf(fpo,"\n"); 
           }
         }
      }
      fclose(fpo); 
    }
 }

 /***************************************/
 /*  MODE '2':  libA-vs-libB comparison */
 /***************************************/
 if (MODE == '2'){
   NmoleculeB = read_SMILES_file(ismifileB,HEAD_MOL_HASH,DIVbunshi,DIVbunbo);
   printf("#NmoleculeB %d\n",NmoleculeB);

  /*
   for (k=0;k<HASH_SIZE;++k){
     if (HEAD_MOL_HASH[k].next != NULL){
       printf(">Lsmiles %d\n",k);
       mB = &(HEAD_MOL_HASH[k]);
       while (mB->next != NULL){
         mB = mB->next;
         printf("%s    [%s]\n",mB->smiles,mB->id);
       }
     } 
   }
    */
  
   fpi = fopen(ismifileA,"r");
   if (fpi==NULL){
     printf("#ERROR: Can't open '%s'.\n",ismifileA);
     exit(1); }

   fpo = fopen(oidfile,"w");
   if (fpo==NULL){
     printf("#ERROR: Can't write to '%s'.\n",oidfile);
     exit(1); }
   fprintf(fpo,"#COMMAND    '%s'\n",COMMAND);
   fprintf(fpo,"#START_DATE '%s'\n",START_DATE);

   fpop = NULL; 
   if (oidpairfile[0] != '\0'){
     fpop = fopen(oidpairfile,"w"); 
     fprintf(fpop,"#COMMAND    '%s'\n",COMMAND);
     fprintf(fpop,"#START_DATE '%s'\n",START_DATE);
   }

   fpou = NULL;
   if (ouniqfile[0] != '\0'){
      fpou = fopen(ouniqfile,"w");
      fprintf(fpou,"#COMMAND    '%s'\n",COMMAND);
      fprintf(fpou,"#START_DATE '%s'\n",START_DATE);
   }

   NmoleculeA = 0;
   while (feof(fpi)==0){
     line[0] = '\0';
     fgets(line,MAX_LINE-1,fpi);
     Lline = strlen(line);
     if (line[Lline-1]== '\n') {line[Lline-1] = '\0'; Lline -= 1;} 
     if ((line[0]!='#') && (Lline>5)){
       split_string_by_first_appearing_symbol(smiles,tail,line,' ');
       split_string_by_first_appearing_symbol(id,tail2,tail,' ');
       /* Lsmiles = strlen(smiles); */
       hash_val  = hash_func(smiles);

       if ((hash_val % DIVbunbo)==DIVbunshi){
         init = 1;
         hit = 0;
         mB = &(HEAD_MOL_HASH[hash_val]);
         while (mB->next != NULL){
            mB = mB->next;
           if (strcmp(smiles,mB->smiles)==0){
             if (init==1){ fprintf(fpo,">%s\n%s",id,mB->id); init = 0;}
             else fprintf(fpo," %s",mB->id); 
             if (fpop!=NULL){fprintf(fpop,"%s %s\n",id,mB->id);}
             hit = 1;
           }
         }
         if (init==0) fprintf(fpo,"\n");
         if ((hit==0) && (fpou!=NULL)){
           fprintf(fpou,"%s %s\n",smiles,id); 
         }
       }
       NmoleculeA += 1;
     } 
   } 
   
   fclose(fpi);

   if (fpo!=NULL){
     fprintf(fpo,"#END_DATE '%s'\n",Get_Date_String());
     printf("#write_identical_smiles_id() -->'%s'\n",oidfile);
     fclose(fpo);
   }

   if (fpop!=NULL){
     fprintf(fpop,"#END_DATE '%s'\n",Get_Date_String());
     printf("#write_identical_smiles_pairs() -->'%s'\n",oidpairfile);
     fclose(fpop);
   }

   if (fpou!=NULL){
     fprintf(fpou,"#END_DATE '%s'\n",Get_Date_String());
     printf("#write_Unique_smiles_id() -->'%s'\n",ouniqfile);
     fclose(fpou);
   }


 }
 return(1);
} /* end of main() */

