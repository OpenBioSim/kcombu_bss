/**

<TreeUN.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 Programs for distance-based phylogenetic tree 
 ( for UPGMA and Neighbor-Joining )
 
 coded by Takeshi Kawabata.

**/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*************************
** STRUCTURE DEFINITION ** 
**************************/

struct NODE { char *name;                 /* Only for Leaf (malloc later) */
              int  num;                   /* Node Number   */
              double  height;             /* Height (for UPGMA). heights for leaves are 0.*/
              struct NODE *ch1,*ch2,*ch3; /* Child Pointer  (for TREE graph )*/
              struct NODE *par;           /* Parent Pointer (for TREE graph ) */
              double len1,len2,len3;      /* Length from parent to child */
              struct NODE *back,*next;    /* Back and next Pointer (for LINEAR graph) */
              char   type;                /* 'L'eaf, 'A'ncestor, '3'-tribranched node  */
              int    nmem;                /* Number of leaves under thie node (for UPGMA) */
              double  *Lleaf;             /* [NLEAF] (malloc later) Lengths for leaves  */
              double alph;                /* Alphabetical Value of "name" */    
             };


/**********************
** GLOBAL VARIABLES  ** 
***********************/

/** Tree-style data  **/
struct NODE   LStaNode;     /* Start Node for LINEAR graph */
struct NODE   TRootNode;    /* Root Node for TREE graph    */
int NLEAF;    /* Number of LEAF */
int NNODE;    /* Number of NODE (in Llist or Tree ) */

/** Array-style data **/
double **DMAT; /* [Nnode][Nnode] : (malloc later) */
struct NODE    **NodeArray;  /* [Nnode] : Array of leaf pointers (malloc later) */

char AlphType;  /* [T or F] : Alphabetical Order Output for names of leaves*/

static char LastModDate[] = "Apr 24, 2013";

/***************
*** FUNCTIONS ** 
****************/

void Get_First_Word_from_File();
void Read_Distance_Matrix();
void Read_Pair_Distance_List();
void Neighbor_Joining();
void Hierarchical_Clustering();
void Set_RSUM();
int Number_Node_Llist();
void Output_NewHamp();
void Set_Lleaf_BottomUp();
void Set_Lleaf_TopDown();
int  Find_MidPoint_Branch();
void Test_Scan();
void Make_New_Root();
void Reverse_Relation();
void Malloc_Lleaf();
void Copy_Node();
void Show_Nodes();
void Add_Node_Llist();
void Delete_Node_Llist();
double Alphabet_Value();
void Output_Tree_Distance();
void Output_Tree_Distance_Vertical();
void Output_Leaves();
double Calculate_Stress();
void Output_Cluster();



int main(argc,argv)
 int argc;
 char **argv;
{
 FILE *fp;
 char dmatfile[128],ipairfile[128],outfile[128],rootfile[128],outdmatfile[128],oleaffile[128],oclusfile[128];
 char odmatvfile[128];
 char Gtype;
 int L,i,j,k,max_i,max_j,Ncluster;
 double len,maxL,mlen,HthreCluster;
 int stanode,endnode;

 /** Set Default **/
 Gtype = 'N'; sprintf(outfile,"-"); rootfile[0] = outdmatfile[0] = ipairfile[0] = '\0';
 odmatvfile[0] = oleaffile[0] = oclusfile[0] = '\0';


 AlphType = 'T'; 
 HthreCluster = -1.0;
 
 if (argc<2){
  printf("TreeUN [dmatfile] (-option)\n");
  printf(" coded by Takeshi Kawabata. LastModified:%s\n",LastModDate);
  printf("  -G :'N'eighbor-Joining 'U'PGMA, 'G'roup-average, 'W'ard method[%c]\n",Gtype);
  printf("     :'S'ingle linkage(nearest neighbor) 'C'omplete linkage(furthest neighbor)\n");
  printf("  -ip: input  file for pair distance vertical file[%s]\n",ipairfile);
  printf("  -of: output file for tree in new-hampshire format[%s]\n",outfile);
  printf("  -df: output file for tree-distance matrix[%s]\n",outdmatfile);
  printf("  -vf: output file for tree-distance matrix vertical[%s]\n",odmatvfile);
  printf("  -rf: output file of rooted-tree[%s] (only for NJ)\n",rootfile);
  printf("  -lf: output file of tree-ordered leaves[%s]\n",oleaffile);
  printf("  -A : Alphabetical Order Ourput (T or F)[%c]\n",AlphType);
  printf("  -ct: threshold for clustering [%lf]\n",HthreCluster); 
  printf("  -cf: output file threshold for clustering [%s]\n",oclusfile); 
  exit(1);
 }

 
 /** Read Option **/
 if (argc>1) sprintf(dmatfile,"%s",argv[1]);
 k = 1;
 while (k<argc){
   if (argv[k][0]=='-'){
     L = strlen(argv[k]);
          if ((L==2)&&(argv[k][1]=='G')) { ++k; Gtype = argv[k][0]; printf("#Gtype %c\n",Gtype);}
     else if ((L==3)&&(argv[k][1]=='i')&&(argv[k][2]=='p')) { ++k; sprintf(ipairfile,"%s",argv[k]); }
     else if ((L==3)&&(argv[k][1]=='o')&&(argv[k][2]=='f')) { ++k; sprintf(outfile,"%s",argv[k]); }
     else if ((L==3)&&(argv[k][1]=='r')&&(argv[k][2]=='f')) { ++k; sprintf(rootfile,"%s",argv[k]); }
     else if ((L==3)&&(argv[k][1]=='d')&&(argv[k][2]=='f')) { ++k; sprintf(outdmatfile,"%s",argv[k]); }
     else if ((L==3)&&(argv[k][1]=='v')&&(argv[k][2]=='f')) { ++k; sprintf(odmatvfile,"%s",argv[k]); }
     else if ((L==3)&&(argv[k][1]=='l')&&(argv[k][2]=='f')) { ++k; sprintf(oleaffile,"%s",argv[k]); }
     else if ((L==2)&&(argv[k][1]=='A')) { ++k; AlphType = argv[k][0]; }
     else if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='f')) { ++k; sprintf(oclusfile,"%s",argv[k]); }
     else if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='t')) { ++k; HthreCluster = atof(argv[k]); }
     else {printf("#ERROR:Can't understand option \"%s\"\n",argv[k]); exit(1);}
    }
   ++k;
  } /* while */

 /*** [1] READ DISTANCE MATRIX   ***/
 Read_Distance_Matrix(dmatfile,&LStaNode);

 /* Show_Nodes(&LStaNode); */

 /*** [2] CALCULATE PHYLOGENETIC TREE ***/
 
      if (Gtype == 'N') Neighbor_Joining(&LStaNode,&TRootNode);
 else if ((Gtype == 'U')||(Gtype=='G')||(Gtype=='S')||(Gtype=='C')||(Gtype=='W'))
   Hierarchical_Clustering(&LStaNode,&TRootNode,Gtype);
 else {printf("#ERROR: Gtype=='%c' is not accepted.",Gtype); exit(1);}

 /*** [3] OUTPUT TREES   ***/
 if (outfile[0] == '-') fp = stdout;
 else
 { fp = fopen(outfile,"w"); printf("#Write to '%s'\n",outfile);
   if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",outfile); exit(1); } }
 Output_NewHamp(fp,&TRootNode);
 fprintf(fp,";\n");
 fclose(fp);


 /*** [4] OUTPUT VARIOUS KINDS OF TREE-RELATED VARIABLES ***/

 /***  (4-1) Output tree-ordered leaf name file ***/
 if (oleaffile[0]!='\0'){
   fp = fopen(oleaffile,"w"); 
   if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",oleaffile); exit(1); } 
   printf("#Output_Leaves --> '%s'\n",oleaffile);
   Output_Leaves(fp,&TRootNode);
   fclose(fp);
 }

 /***  (4-2) Output Cluster file (only for UPGMA) ***/
 if ((oclusfile[0]!='\0') && (HthreCluster>0.0)){
  fp = fopen(oclusfile,"w"); 
  Ncluster = 0;
  if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",oclusfile); exit(1); } 
  printf("#Output_Cluster --> '%s'\n",oclusfile);
  fprintf(fp,"#Clusters by TreeUN.  HthreHeight:%lf\n",HthreCluster);
  Output_Cluster(fp,&TRootNode,HthreCluster,0,&Ncluster);
  fclose(fp);
 }

 /*** (5) TREE-DISTANCE TREATMENT  ***/
  
 if ((rootfile[0] != '\0')||(outdmatfile[0] != '\0')||(odmatvfile[0] != '\0')){
  /*** CALCULATE BRANCH-SUM LENGTH BETWEEN LEAVES ***/ 
  Malloc_Lleaf(&TRootNode); 
  while (TRootNode.nmem == 0) Set_Lleaf_BottomUp(&TRootNode); 
  Set_Lleaf_TopDown(&TRootNode);
  printf("#STRESS %lf\n",Calculate_Stress());
  }

 /*** OUTPUT DISTANCE MATRIX   ***/
 if (outdmatfile[0] != '\0') Output_Tree_Distance(outdmatfile);
 if (odmatvfile[0] != '\0')  Output_Tree_Distance_Vertical(odmatvfile);


 /*** (7) OUTPUT ROOTED NJ  ***/
 if ((rootfile[0] != '\0')&&(Gtype == 'N')){

 maxL = 0; max_i = max_j = 0; 
 for (i=0;i<NLEAF;++i)
  for (j=i+1;j<NLEAF;++j)
  { len = NodeArray[i]->Lleaf[j];
   if (len>maxL) {maxL = len; max_i = i; max_j = j;} }

 /** set midpoint root **/ 
 i = Find_MidPoint_Branch(&TRootNode,max_i,max_j,maxL,&stanode,&endnode,&mlen);

 printf("maxL %lf max_i %d max_j %d\n",maxL,max_i,max_j);
 printf("stanode %d endnode %d mlen %lf case %d\n",stanode,endnode,mlen,i);
 
 /** change root and output **/ 
 Make_New_Root(NodeArray[stanode],NodeArray[endnode],mlen,&TRootNode);

 printf("FINAL ROOT %d ch1 %d  ch2 %d\n",TRootNode.num,TRootNode.ch1->num,TRootNode.ch2->num);

 if (rootfile[0] == '-') fp = stdout;
 else
 { fp = fopen(rootfile,"w"); printf("Write to %s\n",rootfile);
   if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",rootfile); exit(1); } }
 Output_NewHamp(fp,&TRootNode);
 fprintf(fp,";\n");
 fclose(fp);
 }
 return(1);
} /* end of main() */








void Test_Scan(node)
 struct NODE *node;
{
 printf("[%2d] type %c height %8.3f nmem %d",node->num,node->type,node->height,node->nmem);
 if (node->type == 'A'){ 
  printf(" ch1 [%2d] ch2 [%2d] len %8.3lf %8.3lf\n",node->ch1->num,node->ch2->num,node->len1,node->len2);
  Test_Scan(node->ch1);
  Test_Scan(node->ch2);
 } 
 else if (node->type == '3'){ 
  printf(" ch1 [%2d] ch2 [%2d] ch3 [%2d] len %8.3lf %8.3lf %8.3lf\n",
     node->ch1->num,node->ch2->num,node->ch3->num,node->len1,node->len2,node->len3);
  Test_Scan(node->ch1);
  Test_Scan(node->ch2);
  Test_Scan(node->ch3);
  }

  else if (node->type == 'L'){
  printf(" %s\n",node->name);
  }

} /* end of Test_Scan() */





void Make_New_Root(snode,enode,mlen,rootnode)
 struct NODE *snode,*enode,*rootnode;
 double mlen;
{
 struct NODE *newN;

 /* REPLACE OLD NODE TO newN */
 newN = (struct NODE *) malloc(sizeof(struct NODE));
 Copy_Node(newN,rootnode);

 /* ASSIGN NEW NODE */
 rootnode->num = NNODE;
 NodeArray[rootnode->num] = rootnode;
 ++NNODE;
 rootnode->type = 'A';
 rootnode->ch1 = snode;
 rootnode->ch2 = enode;

 rootnode->len1 = mlen; 
 if (snode->ch1->num == enode->num) 
  { snode->ch1 = snode->par;
    rootnode->len2 = snode->len1 - mlen;
    Reverse_Relation(snode->par,snode,&(snode->len1)); } 
 
 if (snode->ch2->num == enode->num) 
  { snode->ch2 = snode->par;
    rootnode->len2 = snode->len2 - mlen;  
    Reverse_Relation(snode->par,snode,&(snode->len2)); } 
 

} /* end of Make_New_Root() */





void Copy_Node(A,B)
 struct NODE *A,*B;
{
 A->type = B->type;
 A->num  = B->num;
 NodeArray[A->num] = A;
 A->par  = B->par;
 A->height = B->height;
 if ((A->type == 'A')||(A->type == '3'))
  { A->ch1 = B->ch1;  
    A->ch2 = B->ch2;  
    B->ch1->par = A;  
    B->ch2->par = A;  
    A->len1 = B->len1;
    A->len2 = B->len2; }
 if (A->type == '3')
 { A->ch3  = B->ch3;  
   B->ch3->par = A;  
   A->len3 = B->len3; } 

} /* end of Copy_Node() */ 




void Reverse_Relation(pnode,cnode,clen)
 struct NODE *pnode,*cnode;
 double *clen;
{
 struct NODE *par,*ch1,*ch2,*ch3;

 par = pnode->par;

 if (pnode->type == 'A'){
  pnode->par = cnode;
   if (pnode->ch1->num == cnode->num) 
  { pnode->ch1 = par; *clen = pnode->len1;
    Reverse_Relation(par,pnode,&(pnode->len1)); }
   if (pnode->ch2->num == cnode->num) 
  { pnode->ch2 = par; *clen = pnode->len2;
    Reverse_Relation(par,pnode,&(pnode->len2)); }
 }

 if (pnode->type == '3'){
   ch1 = pnode->ch1; ch2 = pnode->ch2; ch3 = pnode->ch3;
 
   pnode->type = 'A';
   pnode->par = cnode;
  if (pnode->ch1->num == cnode->num) 
   { pnode->ch1 = ch2; pnode->ch2 = ch3;
     *clen = pnode->len1;
     pnode->len1 = pnode->len2;
     pnode->len2 = pnode->len3; }
  else if (pnode->ch2->num == cnode->num) 
    { pnode->ch2 = ch3; 
      *clen = pnode->len2;
      pnode->len2 = pnode->len3; }
  else if (pnode->ch3->num == cnode->num) 
    { *clen = pnode->len3; 
      printf("ch1 %d %lf\n",pnode->ch1->num,pnode->len1);
      printf("ch2 %d %lf\n",pnode->ch2->num,pnode->len2);
     }
 }

} /* end of Reverse_Relation() */





int Find_MidPoint_Branch(node,i,j,Lmax,sta,end,mlen)
 struct NODE *node;
 int i,j;    /* Leaf number */ 
 double Lmax; /* Half Length between i and j */
 int *sta,*end;
 double *mlen;
{
 double lenP,len1,len2,len3,Lhalf;
 double Li,Li1,Li2,Li3, Lj,Lj1,Lj2,Lj3;
 double eps;

 if (node->type == 'L') return(0);
 
 eps = 0.0000001; Lhalf = Lmax * 0.5;

 /*** Minimum Distance Path Condition  ***/ 
 lenP = node->Lleaf[i] + node->Lleaf[j]; 
 Li = node->Lleaf[i]; Lj = node->Lleaf[j];

 if (fabs(lenP-Lmax) < eps){

  len1 = node->ch1->Lleaf[i] + node->ch1->Lleaf[j]; 
  len2 = node->ch2->Lleaf[i] + node->ch2->Lleaf[j]; 
  Li1  = node->ch1->Lleaf[i]; Lj1  = node->ch1->Lleaf[j]; 
  Li2  = node->ch2->Lleaf[i]; Lj2  = node->ch2->Lleaf[j];

  if (fabs(len1-Lmax)<eps){
        if((Li >= Lhalf)&&(Li1 < Lhalf)) 
     {  *mlen = 0.5 * (Li1 - Lj + node->len1);    
        *sta = node->num; *end = node->ch1->num; return(1);} 
   else if  ((Li <= Lhalf)&&(Li1 > Lhalf))  
     {  *mlen = 0.5 * (Lj1 - Li + node->len1);    
        *sta = node->num; *end = node->ch1->num; return(2);}  
   }
  
  if (fabs(len2-Lmax)<eps){     
        if((Li >= Lhalf)&&(Li2 < Lhalf)) 
     {  *mlen = 0.5 * (Li2 - Lj + node->len2);    
        *sta = node->num; *end = node->ch2->num; return(3);} 
 
   else if  ((Li <= Lhalf)&&(Li2 > Lhalf))  
     {  *mlen = 0.5 * (Lj2 - Li + node->len2);    
        *sta = node->num; *end = node->ch2->num; return(4);}  
   }

  if (node->type == '3'){ 

   len3 = node->ch3->Lleaf[i] + node->ch3->Lleaf[j]; 
   Li3  = node->ch3->Lleaf[i]; Lj3  = node->ch3->Lleaf[j];

  if (fabs(len3-Lmax)<eps){
          if((Li >= Lhalf)&&(Li3 < Lhalf)) 
     {  *mlen = 0.5 * (Li3 - Lj + node->len3);    
        *sta = node->num; *end = node->ch2->num; return(5);} 
   else if  ((Li <= Lhalf)&&(Li3 > Lhalf))  
     {  *mlen = 0.5 * (Lj3 - Li + node->len3);    
        *sta = node->num; *end = node->ch3->num; return(6);}  
   }
 
  } /* if '3' */ 
  
 } /* if lenP = Lmax */

 
  Find_MidPoint_Branch(node->ch1,i,j,Lmax,sta,end,mlen); 
  Find_MidPoint_Branch(node->ch2,i,j,Lmax,sta,end,mlen); 
  if (node->type == '3') Find_MidPoint_Branch(node->ch3,i,j,Lmax,sta,end,mlen); 

 return(1);
} /* end of Find_MidPoint_Branch() */





void Malloc_Lleaf(node)
 struct NODE *node;
{ 
 int i;

 /* Malloc double Lleaf[NLEAF].      */
 /* Initial value of Lleaf[] is -1.0 */
 node->nmem = 0; 
 node->Lleaf = (double *)malloc(sizeof(double)*NLEAF);
 for (i=0;i<NLEAF;++i) node->Lleaf[i] = -1.0; 
 
 if ((node->type=='A')||(node->type=='3'))
 { Malloc_Lleaf(node->ch1); Malloc_Lleaf(node->ch2);}
 
 if (node->type=='3') Malloc_Lleaf(node->ch3);

} /* end of Malloc_Lleaf() */



void Set_Lleaf_BottomUp(node)
 struct NODE *node;
{
 int i;

 if (node->nmem == 0){
 
 if (node->type == 'L')
  { node->nmem = 1;
    node->Lleaf[node->num] = 0.0; 
  } /* if 'L' */ 

  if (node->type == 'A'){
   if ((node->ch1->nmem > 0)&&(node->ch2->nmem >0)){ 
     node->nmem = node->ch1->nmem + node->ch2->nmem; 

     for (i=0;i<NLEAF;++i){ 
      if (node->ch1->Lleaf[i] >= 0.0)
        node->Lleaf[i] = node->ch1->Lleaf[i] + node->len1; 
 
      if (node->ch2->Lleaf[i] >= 0.0)
        node->Lleaf[i] = node->ch2->Lleaf[i] + node->len2;
      } 
    }
   
    else 
    { Set_Lleaf_BottomUp(node->ch1); Set_Lleaf_BottomUp(node->ch2); }

   } /* if 'A' */

  if (node->type == '3'){
   if ( (node->ch1->nmem >0 )&&(node->ch2->nmem>0)&&(node->ch3->nmem >0)){ 
     node->nmem = node->ch1->nmem + node->ch2->nmem + node->ch3->nmem; 
     for (i=0;i<NLEAF;++i){ 
       if (node->ch1->Lleaf[i] >= 0.0)
         node->Lleaf[i] = node->ch1->Lleaf[i] + node->len1; 

       if (node->ch2->Lleaf[i] >= 0.0)
         node->Lleaf[i] = node->ch2->Lleaf[i] + node->len2; 

       if (node->ch3->Lleaf[i] >= 0.0)
         node->Lleaf[i] = node->ch3->Lleaf[i] + node->len3; 
       } 
    
    }
    else 
    { Set_Lleaf_BottomUp(node->ch1); 
      Set_Lleaf_BottomUp(node->ch2); 
      Set_Lleaf_BottomUp(node->ch3); }

   } /* if '3' */

 } /* if nmem == 0 */ 

} /* end of Set_Lleaf_BottomUp() */





void Set_Lleaf_TopDown(node)
 struct NODE *node;
{
 int i;

 if ((node->type == 'A')||(node->type == '3')){ 
  for (i=0;i<NLEAF;++i)
   if ((node->ch1->Lleaf[i]==-1.0)&&(node->Lleaf[i]>=0.0))
    node->ch1->Lleaf[i] = node->Lleaf[i] + node->len1;
  
  for (i=0;i<NLEAF;++i)
   if ((node->ch2->Lleaf[i]==-1.0)&&(node->Lleaf[i]>=0.0))
    node->ch2->Lleaf[i] = node->Lleaf[i] + node->len2;

  if (node->type == '3')
  { for (i=0;i<NLEAF;++i)
     if ((node->ch3->Lleaf[i]==-1.0)&&(node->Lleaf[i]>=0.0))
      node->ch3->Lleaf[i] = node->Lleaf[i] + node->len3; }

   Set_Lleaf_TopDown(node->ch1);
   Set_Lleaf_TopDown(node->ch2);
  
  if (node->type == '3') Set_Lleaf_TopDown(node->ch3);

 }

} /* end of Set_Lleaf_TopDown() */




void Hierarchical_Clustering(lstanode,trootnode,Gtype)
 struct NODE *lstanode,*trootnode;
 char  Gtype; /* type of Clustering 'U'PGMA,'S'ingle-linkage, 'C'omplete-linkage, 'W'ard*/
{
 int end,NLnode;
 double minDij;
 struct NODE *I,*J,*minI,*minJ,*newN,*M;

 end = 0; 
 while (end == 0){

  /** (1) Finding the closest node pairs (minI, minJ) **/ 
   minDij = 1000000.0;
   minI = minJ = NULL;
   /** Find Mimimum Dij **/
   I = lstanode;
   while (I->next !=NULL){
     I = I->next;
     J = I;
     while (J->next != NULL) { 
       J = J->next;
       if (DMAT[I->num][J->num] < minDij) 
        {  minDij = DMAT[I->num][J->num]; minI = I; minJ = J; } 
      } /* while (J) */ 
   
    } /* while (I) */ 
 
   /* 
   printf("min I %d J %d minDij %lf\n",minI->num,minJ->num,minDij);
   */


   /** (2) Making New Node **/
   newN = (struct NODE *) malloc(sizeof(struct NODE));
   newN->num = NNODE; 
   NodeArray[newN->num] = newN;
   newN->type = 'A';
   ++NNODE;

   /** (3) Set new DMAT for Linear List **/
   M = lstanode;
   while (M->next !=NULL){
     M = M->next;
     if ((Gtype=='U')||(Gtype=='G')){
       DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] 
         = (minI->nmem * DMAT[minI->num][M->num] + minJ->nmem*DMAT[minJ->num][M->num])/(minI->nmem + minJ->nmem);
      }

     if (Gtype=='S'){
       if (DMAT[minI->num][M->num]<DMAT[minJ->num][M->num]) 
         DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] = DMAT[minI->num][M->num];
       else
         DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] = DMAT[minJ->num][M->num];
      }

     if (Gtype=='C'){
       if (DMAT[minI->num][M->num]>DMAT[minJ->num][M->num]) 
         DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] = DMAT[minI->num][M->num];
       else
         DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] = DMAT[minJ->num][M->num];
      }

     if (Gtype=='W'){
       DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] 
         = ((minI->nmem + M->nmem)*DMAT[minI->num][M->num] + (minJ->nmem+M->nmem)*DMAT[minJ->num][M->num] 
            - (M->nmem)*DMAT[minI->num][minJ->num])/(minI->nmem + minJ->nmem + M->nmem);
     }

   } /* while (M) */
 

   /** (4) Considering AlphabetOrder **/
   if ((AlphType == 'T')&&(minI->alph > minJ->alph))
     { I = minI; J = minJ; 
       minI = J; minJ = I; }

   /** (5) Include newN to TREE **/
   newN->ch1 = minI;
   newN->ch2 = minJ;
   minI->par = newN; 
   minJ->par = newN; 
   if (Gtype=='U')
    newN->height = DMAT[minI->num][minJ->num]/2.0;
   else
    newN->height = DMAT[minI->num][minJ->num];

   /* if you consider sum of branches as "distance", 
      you must define newN->height = DMAT[minI->num][minJ->num]/2.0; */

   newN->len1 = newN->height - minI->height;
   newN->len2 = newN->height - minJ->height;
   newN->nmem = minI->nmem + minJ->nmem;
   newN->alph = minI->alph;

   /* printf("len1 %lf len2 %lf\n",newN->len1,newN->len2); */
 

   /** (6) Add newN and delete minI and minJ **/
   Add_Node_Llist(lstanode,newN);
   Delete_Node_Llist(minI);
   Delete_Node_Llist(minJ);
 
   NLnode  = Number_Node_Llist(lstanode);

   /** Terminal Operation **/
   if (NLnode <=1 ) 
    { end = 1;
      I = lstanode->next; 
      Copy_Node(trootnode,I); 
      trootnode->type = 'A';
     } 
 
  } /* while (end==0) */


} /* end of Hierarchical_Clustering() */






void Neighbor_Joining(lstanode,trootnode)
 struct NODE *lstanode,*trootnode;
{
 int end,NLnode,nmax_node;
 double Dij,minDij;
 struct NODE *I,*J,*minI,*minJ,*newN,*M;
 double *RSUM;  /* [Nnode]        : (malloc later, for NJ)  */

 nmax_node = NLEAF * 2 + 1;
 RSUM = (double *)malloc(sizeof(double)*nmax_node);
 minI = minJ = NULL; 
 end = 0; 
 while (end == 0)
 {
   /** Set RSUM **/
   Set_RSUM(RSUM,lstanode);
   
  /** Find Mimimum Dij **/
  I = lstanode;
  minDij = 1000000.0;
 
   while (I->next !=NULL){
    I = I->next;
    J = I;
    while (J->next != NULL) {
      J = J->next;
      /* printf("%d %d\n",I->num,J->num);  */
      Dij = DMAT[I->num][J->num] - RSUM[I->num] - RSUM[J->num];
      if (Dij < minDij) {  minDij = Dij; minI = I; minJ = J; } 
     } /* while (J) */ 
  
   } /* while (I) */ 
/*
    printf("min I %d J %d minDij %lf\n",minI->num,minJ->num,minDij); 
*/

  /** Making New Node **/
  newN = (struct NODE *) malloc(sizeof(struct NODE));
  newN->num = NNODE;
  NodeArray[newN->num] = newN;
  newN->type = 'A';
  ++NNODE;

  /** Set DMAT for Linear List**/
  M = lstanode;
  while (M->next !=NULL){
   M = M->next;
   DMAT[newN->num][M->num]  =  DMAT[M->num][newN->num] 
     = (DMAT[minI->num][M->num] + DMAT[minJ->num][M->num] - DMAT[minI->num][minJ->num])/2.0;  
  } 

  /** Considering AlphabetOrder **/
  if ((AlphType == 'T')&&(minI->alph > minJ->alph))
    { I = minI; J = minJ; 
      minI = J; minJ = I; }
  
 /** Include newN to TREE **/
  newN->ch1 = minI;
  newN->ch2 = minJ;
  minI->par = newN; 
  minJ->par = newN; 
  newN->len1 = (DMAT[minI->num][minJ->num] + RSUM[minI->num] - RSUM[minJ->num])/2.0;
  newN->len2 = (DMAT[minI->num][minJ->num] + RSUM[minJ->num] - RSUM[minI->num])/2.0;
  newN->alph = minI->alph;
  /* printf("len1 %lf len2 %lf\n",newN->len1,newN->len2); */
 

  /** Add newN and delete minI and minJ **/
  Add_Node_Llist(lstanode,newN);
  Delete_Node_Llist(minI);
  Delete_Node_Llist(minJ);
 
  NLnode  = Number_Node_Llist(lstanode);

  /** Terminal Operation **/
  if (NLnode <=2 ) 
   { 
     end = 1;
     I = lstanode->next;  
     J = I->next;
 
     if (I->type == 'A')
     { Copy_Node(trootnode,I);
       trootnode->ch3 = J;
       J->par = trootnode;  trootnode->len3 = DMAT[I->num][J->num]; }
     else  
     { Copy_Node(trootnode,J);
       trootnode->ch3 = I;
       I->par = trootnode;  trootnode->len3 = DMAT[I->num][J->num]; }
     
      trootnode->type = '3';
    } 
 
 } /* while (end==0) */


 free(RSUM);

} /* end of Neighbor_Joininig() */






void Set_RSUM(RSUM,lstanode)
 double *RSUM;
 struct NODE *lstanode;
{
 struct NODE *I,*J;
 int nnode;

 I = lstanode;
 while (I->next != NULL) 
 { I = I->next;
  RSUM[I->num] = 0.0;
  J = lstanode; nnode = 0; 

  while (J->next != NULL) 
  { J = J->next;
    RSUM[I->num] += DMAT[I->num][J->num];
    ++nnode;
   } /* while (J) */

  RSUM[I->num] /= (nnode-2);
  /* printf("RSUM %d %lf\n",I->num,RSUM[I->num]); */

 } /* while (I) */


} /* end of Set_RSUM() */




void Read_Distance_Matrix(fname,lstanode)
 char   *fname;
 struct NODE *lstanode;
{
 /* 
 <FORMAT EXAMPLE>
 4
 1mbd-        0.0  34.9  29.3  28.0
 1ecd-       34.9   0.0  43.0  40.9
 4hhbA       29.3  43.0   0.0  23.1
 4hhbB       28.0  40.9  23.1   0.0
 */
 FILE *fp;
 int i,j,nmax_node;
 char word[128],line[512];
 struct NODE *bnode;

 fp = fopen(fname,"r");
 if (fp==NULL) { printf("#ERROR:Can't open dismatfile \"%s\"\n",fname); exit(1);}

 /** Read NLEAF (skipping "#'-headed line) **/ 
 NLEAF = 0;
 while ((NLEAF==0)&&(feof(fp)==0)){ 
   line[0] = '\0'; 
   fgets(line,511,fp);
   if ((line[0]!='#')&&(strlen(line)>0)){
     NLEAF = atoi(line); 
     if (NLEAF<3){
       printf("#ERROR(TreeUN:NLEAF=%d):at least three leaves are necessary.'%s'\n",NLEAF,line); exit(1);
     }
   } 
 }

 /*  printf("NLEAF %d\n",NLEAF);  */

 /** malloc DMAT,RSUM,NodeArray **/ 
 nmax_node = NLEAF * 2 + 1;
 DMAT = (double **)malloc(sizeof(double *)*nmax_node);
 for (i=0;i<nmax_node;++i)
  DMAT[i] = (double *)malloc(sizeof(double )*nmax_node);

 NodeArray = (struct NODE **)malloc(sizeof(struct NODE*)*nmax_node);

 /** Reading file using line and buffline **/ 
 bnode = lstanode;

 for (i=0;i<NLEAF;++i){
   bnode->next = (struct NODE *)malloc(sizeof(struct NODE));
   bnode->next->back = bnode;
   bnode = bnode->next; 
   bnode->type = 'L'; bnode->nmem = 1;
   bnode->height = 0.0; 
   bnode->num = i; 
   NodeArray[bnode->num] = bnode;
   bnode->next = NULL;
   
   Get_First_Word_from_File(fp,word,127);
  
   bnode->name = (char *)malloc(sizeof(char)*(strlen(word)+1));
   sprintf(bnode->name,"%s",word);
   bnode->alph = Alphabet_Value(bnode->name);

   /* printf("i %d %s %lf\n",i,bnode->name,bnode->alph); */

   for (j=0;j<NLEAF;++j){ 
     Get_First_Word_from_File(fp,word,127);
     DMAT[i][j] = atof(word);
     if (DMAT[i][j]<0.0)
     { printf("#ERROR:Distance %d and %d is negative (%lf)\n",i,j,DMAT[i][j]); 
       exit(1); }

     /* printf("[%d %d] %lf\n",i,j,DMAT[i][j]);  */
     
    } /* j */ 

  } /* i loop */ 

 fclose(fp);

 NNODE = NLEAF;

} /* end of Read_Distance_Matrix() */









double Alphabet_Value(name)
 char *name;
{
 int i,L;
 double base,val; 

 L = strlen(name);
 val = 0.0; base = 1.0;

 for (i=0;i<L;++i){
   val += (double)name[i] * base;
   base /= 128.0;
  }

 return(val);
}  /* end of Alphabet_Value() */






void Get_First_Word_from_File(fp,word,LwordMax)
 FILE *fp;
 char *word;
 int LwordMax;
{
 /*
  Return a [!isgraph]-splitted word, which appears first.
  ex)
   if fp contains   = " ABCDE  FGH ";
   then return word = "ABCDE".

  isgraph(sym)==0 means sym = ' ', '\t' or '\n'.
*/

 char sym;
 int  end,wstart,Lword;

 word[0] = '\0';
 end = Lword =  wstart = 0;

 while (end==0){
  sym = fgetc(fp);

  if ((wstart==0)&&(isgraph(sym)!=0)) {wstart = 1;}

  if ((wstart==1)
      &&((isgraph(sym)==0)||(feof(fp)==1)||(Lword>=LwordMax)))
    { wstart = 0; word[Lword] = '\0'; end = 1;}

  if (wstart==1) {word[Lword] = sym; ++Lword;}

  if (feof(fp)==1) end = 1;
 }

} /* end of Get_First_Word_from_File() */



void Output_NewHamp(fp,node)
 FILE *fp;
 struct NODE *node;
{
 if ((node->type == 'A')||(node->type == 'R')){
   fprintf(fp,"(\n");
   if (node->ch1->type == 'L')
     fprintf(fp,"%s", node->ch1->name); else Output_NewHamp(fp,node->ch1); 

  fprintf(fp,":%lf,\n",node->len1);

  if (node->ch2->type == 'L')
   fprintf(fp,"%s", node->ch2->name); else Output_NewHamp(fp,node->ch2);

  fprintf(fp,":%lf)\n",node->len2);
  }

 if (node->type == '3'){
  fprintf(fp,"(\n");

  if (node->ch1->type == 'L')
   fprintf(fp,"%s", node->ch1->name); else Output_NewHamp(fp,node->ch1);

  fprintf(fp,":%lf,\n",node->len1);

  if (node->ch2->type == 'L')
  fprintf(fp,"%s", node->ch2->name); else Output_NewHamp(fp,node->ch2);

  fprintf(fp,":%lf,\n",node->len2);

  if (node->ch3->type == 'L')
  fprintf(fp,"%s", node->ch3->name); else Output_NewHamp(fp,node->ch3);

  fprintf(fp,":%lf)\n",node->len3);
 }

} /* end of Output_NewHamp() */






void Show_Nodes(head)
 struct NODE *head;
{
 struct NODE *bn;
 char pt;

 bn = head;
 while (bn->next != NULL){
  bn = bn->next;
  if (bn->next == NULL) pt = 'N'; else pt = '*';
  printf("%s %c\n",bn->name,pt);
  }

} /* end of Show_Nodes() */



void Add_Node_Llist(stanode,newnode)
 struct NODE *stanode,*newnode;
{
 struct NODE *bn;

 bn = stanode;
 while (bn->next != NULL) { bn = bn->next; }

 bn->next = newnode;
 newnode->back = bn;
 newnode->next = NULL;

} /* end of Add_Node_Llist() */


int Number_Node_Llist(stanode)
 struct NODE *stanode;
{
 struct NODE *bn;
 int n;
 n = 0;
 bn = stanode;
 while (bn->next != NULL) { bn = bn->next; ++n;}
 return(n);

} /* end of Number_Node_Llist() */



void Delete_Node_Llist(delnode)
 struct NODE *delnode;
{
 delnode->back->next = delnode->next;
 delnode->next->back = delnode->back; 

} /* end of Delete_Node_Llist() */


void Output_Tree_Distance(fname)
 char *fname;
{
 FILE *fp;
 int i,j;

 if (fname[0] == '-') fp = stdout;
 else
  { fp = fopen(fname,"w"); printf("#Output_Tree_Distance() -->\"%s\"\n",fname);
    if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1); } }

  fprintf(fp,"%d\n",NLEAF);
  for (i=0;i<NLEAF;++i)
  { fprintf(fp,"%-10s ",NodeArray[i]->name);
    for (j=0;j<NLEAF;++j) fprintf(fp,"%5.1lf ",NodeArray[i]->Lleaf[j]);
    fprintf(fp,"\n"); }

 fprintf(fp,"#STRESS %lf\n",Calculate_Stress());
 if (fp!=stdout) fclose(fp);

} /* end of Output_Tree_Distance() */




void Output_Tree_Distance_Vertical(fname)
 char *fname;
{
 FILE *fp;
 int i,j;

 if (fname[0] == '-') fp = stdout;
 else
  { fp = fopen(fname,"w"); printf("#Output_Tree_Distance_Vertical()-->\"%s\"\n",fname);
    if (fp==NULL) {printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1); } }

 for (i=0;i<NLEAF;++i)
  for (j=i;j<NLEAF;++j) 
  fprintf(fp,"%s %s %lf\n",NodeArray[i]->name,NodeArray[j]->name,NodeArray[i]->Lleaf[j]);
 if (fp!=stdout) fclose(fp);

} /* end of Output_Tree_Distance_Vertical() */



void Output_Leaves(fp,node)
 FILE *fp;
 struct NODE *node;
{
 if (node->type == 'L') fprintf(fp,"%s\n", node->name);
 
 else if ((node->type == 'A')||(node->type == 'R')){
  Output_Leaves(fp,node->ch1); 
  Output_Leaves(fp,node->ch2); 
  }
 else if (node->type == '3'){
  Output_Leaves(fp,node->ch1); 
  Output_Leaves(fp,node->ch2); 
  Output_Leaves(fp,node->ch3); 
  }

} /* end of Output_Leaves() */








double Calculate_Stress()
{
 int i,j;
 double SSdiff, SSreal,stress;

 SSdiff = SSreal = 0.0; 
  for (i=0;i<NLEAF;++i){ 
    for (j=i+1;j<NLEAF;++j){
     SSdiff += (DMAT[i][j] - NodeArray[i]->Lleaf[j])* (DMAT[i][j]-NodeArray[i]->Lleaf[j]);
     SSreal +=  DMAT[i][j]*DMAT[i][j];
     printf("%s %s %lf -> %lf\n",NodeArray[i]->name,NodeArray[j]->name,DMAT[i][j],NodeArray[i]->Lleaf[j]);
    }
 }

 if ((SSdiff>0.0)&&(SSreal>0.0)) stress = sqrt(SSdiff/SSreal); else stress = 0.0;

 return(stress);

} /* end of Calculate_Stress() */





void Output_Cluster(fp,node,height_thre,out_cluster,Ncluster)
 FILE   *fp;
 struct NODE *node;
 double height_thre;
 char   out_cluster;  /* 0 or 1 */
 int    *Ncluster;    /* number of cluster */
{
  char new_out_cluster;

  new_out_cluster = out_cluster;
  if ((out_cluster==0) && (node->height <= height_thre)){
   (*Ncluster) += 1;
   fprintf(fp,">CLUSTER %d Nmember %d\n",*Ncluster,node->nmem);
   new_out_cluster = 1;
  }

   if (node->type == 'A'){ 
    Output_Cluster(fp,node->ch1,height_thre,new_out_cluster,Ncluster);
    Output_Cluster(fp,node->ch2,height_thre,new_out_cluster,Ncluster);
   }

   if (node->type == '3'){ 
    Output_Cluster(fp,node->ch1,height_thre,new_out_cluster,Ncluster);
    Output_Cluster(fp,node->ch2,height_thre,new_out_cluster,Ncluster);
    Output_Cluster(fp,node->ch3,height_thre,new_out_cluster,Ncluster);
   }

   if (node->type == 'L'){ 
    fprintf(fp,"%s\n",node->name); 
   }

} /* end of Output_Cluster() */
