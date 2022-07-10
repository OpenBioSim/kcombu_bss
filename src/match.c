/*

 <match.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 functions for dealing with "structure MATCH"

 
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
#include "clique.h"
#include "c_clique.h"
#include "buildup.h"
#include "MrgSrtMATCH.h"
#include "PCAfit.h"
#include "molprop.h"
#include "stereo_check.h"
#include "io_match.h"



/** FUNCTIONS (GLOBAL) **/
int  Match_Two_Molecules();
int  Make_All_Atom_Pair_MATCH();
int  Make_All_Atom_Pair_MATCH_for_Substructure_or_Isomorphic();
float cal_select_dis_for_MATCH();
void  cal_select_dis_for_MATCHlist();
void Malloc_MATCH();
void Free_MATCH();
void Malloc_MATCHlist();
void Add_MATCH_to_MATCHlist();
int  Add_MATCH_to_MATCHlist_Keep_Only_Best();
void Free_MATCHlist();
void Copy_MATCH();
void Copy_MATCHlist();
int Length_of_MATCHlist();
int Nmark_in_MATCHlist();
float Tanimoto_Coefficient();
int Smaller_Int();
int Larger_Int();
void Keep_Only_Nonredundant_MATCH();
void Keep_Only_Substructure_or_Isomorphic_MATCH();
void add_hydrogen_pairs_to_MATCH();
float Compare_Electric_Charges();
void Sort_MATCHlist_by_Npair_and_select_dis();

/** FUNCTIONS (LOCAL) **/
static int Common_Pair_Bwn_Two_MATCHes();
static void Bubble_Sort_anumAB_in_MATCHlist_by_anumA();
static int atomtype_match();



int Match_Two_Molecules(molA, molB,optMlist,AlgoType,ConGraphType,maxDIFtopodis,Nkeep,iatom_matchfile)
 /*
   Match Molecule using various type of Algorithm 
 */
  struct MOLECULE *molA, *molB;
  struct MATCH *optMlist;
  char   AlgoType;         /* Algorithm type. 'B'uild-up, 'X':exact */
  char   ConGraphType;     /* Connection of MCS. 'C'onnected, 'D'isconnected, 'T'D-MCS,'t'C-MCS */
  int    maxDIFtopodis;    /* accepted maximum difference of topological distance (number of bonds in the shortest path) */
  int    Nkeep; 
  char   *iatom_matchfile;
{
 struct MATCH pairMlist,CS,*m; 
 int  Npair_max,Npair,ret;
 char  *CA_NOTstring;

/*
 printf("#Match_Two_Molecules(AlgoType:%c ConGraphType:%c maxDIFtopodis:%d Nkeep:%d)\n",AlgoType,ConGraphType,maxDIFtopodis,Nkeep);
*/

 Npair_max =  Smaller_Int(molA->Nheavyatom,molB->Nheavyatom);
 PAR.calc_finish = 'N';
 ret = 0;

 /* If iatom_matchfile != NULL, read MATCH from the file **/ 
 if ((iatom_matchfile[0]!='\0')&&(AlgoType!='B')){
   optMlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
   optMlist->next->Npair_malloc = 0;
   optMlist->next->Npair        = 0;
   optMlist->next->next = NULL;
   Read_MATCH(iatom_matchfile,optMlist->next,molA,molB);
 }
 /** ['B'] : Build-up method **/
  else if (AlgoType=='B'){
    /** (1) Preparing "optMlist" and "pairMlist" **/ 
    Malloc_MATCHlist(optMlist, Npair_max,Nkeep);


    if ((ConGraphType=='S') || (ConGraphType=='I')){
      Npair = Make_All_Atom_Pair_MATCH_for_Substructure_or_Isomorphic(&pairMlist,molA,molB,ConGraphType);
    }
    else{ 
      Npair = Make_All_Atom_Pair_MATCH(&pairMlist,molA,molB);
    }
 
    PAR.TotalNatompair = Npair;
    Merge_Sort_Double_Linked_List_MATCH(&pairMlist,'I'); 
  
   if (iatom_matchfile[0]!='\0'){
     Read_MATCH(iatom_matchfile,optMlist->next,molA,molB);
     printf("#input Npair %d\n",optMlist->next->Npair);
     optMlist->next->select_dis =  cal_select_dis_for_MATCH(optMlist->next,molA,molB);
     optMlist->nodetype = 'I';
     optMlist->next->nodetype = 'N';
     optMlist->next->next->nodetype = 'E';
   }
   else if (PAR.Wdistan<=0.0) {
     Copy_MATCHlist(optMlist,&pairMlist);
   }



    /** (2) Do Build-up **/
       if (((ConGraphType=='C')||(ConGraphType=='D')||(ConGraphType=='T')
           ||(ConGraphType=='t')||(ConGraphType=='S')||(ConGraphType=='I')))
      BuildUp(optMlist,&pairMlist,molA,molB,Nkeep,ConGraphType,maxDIFtopodis,PAR.maxNcomponent);
    else if (ConGraphType=='n')
      Incremental_Disconnected_BuildUp(optMlist,&pairMlist,molA,molB,Nkeep,maxDIFtopodis,PAR.maxNcomponent);
    else if (ConGraphType=='i'){
      printf("#BuildUp(optMlist,&pairMlist,molA,molB,Nkeep,'C',PAR.maxNcomponent);\n");
      BuildUp(optMlist,&pairMlist,molA,molB,Nkeep,'C',maxDIFtopodis,PAR.maxNcomponent);
      printf("#BuildUp(optMlist,&pairMlist,molA,molB,Nkeep,'D',PAR.maxNcomponent);\n");
      BuildUp(optMlist,&pairMlist,molA,molB,Nkeep,'D',maxDIFtopodis,PAR.maxNcomponent);
    }
    else {printf("#ERROR:Can't understand the option -con '%c'.\n",ConGraphType); exit(1);}
    if ((ConGraphType=='S')||(ConGraphType=='I'))  
      Keep_Only_Substructure_or_Isomorphic_MATCH(optMlist,molA, molB,ConGraphType);
    cal_select_dis_for_MATCHlist(optMlist,molA,molB);
    Set_Number_of_Connected_Components_for_MATCHlist(optMlist,molA);
    Free_MATCHlist(&pairMlist);
    PAR.calc_finish = 'F';
  }




 
 /* ['X'] or ['x'] && ConGraphType=='D': Bron-Kerbosch method (disconnected MCS) */
 else if (((AlgoType=='X') || (AlgoType=='x')) && ((ConGraphType=='D')||(ConGraphType=='T'))){
    PAR.M = 0;
    PAR.N = 0;
    PAR.O = 0;
    PAR.P = 0;
    /** Preparing "pairMlist" **/ 
    Npair = Make_All_Atom_Pair_MATCH(&pairMlist,molA,molB);
    PAR.TotalNatompair = Npair;
    Evaluate_Degree_of_Atom_Pair(&pairMlist,molA,molB,ConGraphType,maxDIFtopodis); 
    Npair = Delete_Small_Degree_of_Atom_Pair(&pairMlist,1);
    Merge_Sort_Double_Linked_List_MATCH(&pairMlist,'D'); 

    if ((PAR.BKmaxNpair>0) && (Npair>PAR.BKmaxNpair)){
     printf("#WARNING:too much Npair(%d). Reduce into Npair %d.\n",Npair,PAR.BKmaxNpair); 
     Npair = Keep_Only_First_K_MATCH(&pairMlist,PAR.BKmaxNpair);
    }

    Renumbering_MATCHlist(&pairMlist);
    /** Preparing "CS" **/ 
    Malloc_MATCH(&CS,1);
    CS.Npair = 0;
    optMlist->next = NULL; optMlist->nodetype = 'E';
  
    CA_NOTstring  = (char *)malloc(sizeof(char)*Npair); 
    Init_Bit_String(CA_NOTstring,'c',Npair);
    
    /* printf("#Npair %d\n",Npair); */

         if (AlgoType=='X') 
      ret = BronKerbosch_Pivot_BitString(optMlist,&CS,&pairMlist,Npair,CA_NOTstring,molA,molB,ConGraphType,maxDIFtopodis);
    else if (AlgoType=='x') 
      ret = BronKerbosch_Version1_BitString(optMlist,&CS,&pairMlist,Npair,CA_NOTstring,molA,molB,ConGraphType,maxDIFtopodis);
    printf("#FINAL(Ntry %d Nfree %d Nkept %d Nstop %d Nout %d)\n",PAR.N,PAR.M,PAR.N-PAR.M,PAR.O,PAR.P); fflush(stdout);
 
    cal_select_dis_for_MATCHlist(optMlist,molA,molB);
    Merge_Sort_Double_Linked_List_MATCH(optMlist,'I');
    Bubble_Sort_anumAB_in_MATCHlist_by_anumA(optMlist);
    Set_Number_of_Connected_Components_for_MATCHlist(optMlist,molA);
    Free_MATCHlist(&pairMlist);
    if (ret==-1) PAR.calc_finish = 'B'; else PAR.calc_finish = 'F';
  }




 /* ['X'] or ['x'] && ConGraphType=='C': Koch-Cazals-Karande method (connected MCS) */
 else if (((AlgoType=='X') || (AlgoType=='x')) && ((ConGraphType=='C')||(ConGraphType=='t'))){

    /** Preparing "pairMlist" **/ 
    Npair = Make_All_Atom_Pair_MATCH(&pairMlist,molA,molB);
    PAR.TotalNatompair = Npair;
    Evaluate_Connected_Degree_of_Atom_Pair(&pairMlist,molA,molB);
    Npair = Delete_Small_Degree_of_Atom_Pair(&pairMlist,1); 
    Merge_Sort_Double_Linked_List_MATCH(&pairMlist,'D'); 
    /* Print_MATCHlist(&pairMlist,"pairMlist"); */
    if ((PAR.BKmaxNpair>0) && (Npair>PAR.BKmaxNpair)){
     printf("#WARNING:too much Npair(%d). Reduce into Npair %d.\n",Npair,PAR.BKmaxNpair); 
     Npair = Keep_Only_First_K_MATCH(&pairMlist,PAR.BKmaxNpair);
    }

    Renumbering_MATCHlist(&pairMlist);
    
    /** Preparing "CS" **/ 
    Malloc_MATCH(&CS,1);
    CS.Npair = 0;
    optMlist->next = NULL; optMlist->nodetype = 'E';
    /* printf("#Npair %d\n",Npair); */

    ret = C_clique_Cazals_Karande(optMlist,&pairMlist,Npair,molA,molB,ConGraphType,maxDIFtopodis); 
   
    cal_select_dis_for_MATCHlist(optMlist,molA,molB);
    Merge_Sort_Double_Linked_List_MATCH(optMlist,'I');
    Bubble_Sort_anumAB_in_MATCHlist_by_anumA(optMlist); 
    Set_Number_of_Connected_Components_for_MATCHlist(optMlist,molA);
    Free_MATCHlist(&pairMlist);
    if (ret==-1) PAR.calc_finish = 'B'; else PAR.calc_finish = 'F';
  }





 
 /** ['P']: PCA-based 3D matching **/
 else if (AlgoType=='P'){
    Cal_Molecule_EigVec_Covariance_Matrix(molA);
    Cal_Molecule_EigVec_Covariance_Matrix(molB);
    optMlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
    optMlist->next->next = NULL;
    optMlist->next->ConnectGraphType  = PAR.ConnectGraphType;
    optMlist->next->maxDIFtopodis     = PAR.maxDIFtopodis;
    optMlist->next->atomtype_class    = PAR.atomtype_class;
    Cal_MATCH_for_PCAfit_Superimposed_3D_Molecule_Pair(optMlist->next,molA,molB,ConGraphType);
    Bubble_Sort_anumAB_in_MATCHlist_by_anumA(optMlist);
    optMlist->next->Ncomponent = Number_of_Connected_Components_for_MATCH(optMlist->next->anumA,optMlist->next->Npair,molA);
    PAR.calc_finish = 'F';
  }

 /** ['3']: Raw-3D-matching **/
  else if (AlgoType=='3'){
    optMlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
    optMlist->next->next = NULL;
    optMlist->next->ConnectGraphType  = PAR.ConnectGraphType;
    optMlist->next->maxDIFtopodis     = PAR.maxDIFtopodis;
    optMlist->next->atomtype_class    = PAR.atomtype_class;
    Cal_MATCH_for_Superimposed_3D_Molecule_Pair(optMlist->next,molA,molB,ConGraphType);
    optMlist->next->Ncomponent = Number_of_Connected_Components_for_MATCH(optMlist->next->anumA,optMlist->next->Npair,molA);
    Bubble_Sort_anumAB_in_MATCHlist_by_anumA(optMlist);
    PAR.calc_finish = 'F';
  }
  else {
    printf("#ERROR: -alg '%c', -con '%c' is not type valid.\n",AlgoType,ConGraphType);
    exit(1);
  }

  /** Check redundancy against best one */
 /*
  m = optMlist;
  while ((m->next != NULL)&&(m->next->nodetype!='E')){
    m = m->next; 
    i = Common_Pair_Bwn_Two_MATCHes(optMlist->next,m);
  }
 */

 if (PAR.includeHinMCS=='T'){
   add_hydrogen_pairs_to_MATCH(optMlist,molA,molB);
 }



 m = optMlist;
 while ((m->next != NULL)&&(m->next->nodetype!='E')){
    m = m->next; 
    m->ConnectGraphType = ConGraphType;   
    m->maxDIFtopodis    = maxDIFtopodis;
 }
 return(ret);

} /* end of Match_Two_Molecules() */






int Make_All_Atom_Pair_MATCH(Mlist,molA,molB)
 struct MATCH    *Mlist; 
 struct MOLECULE *molA,*molB;
{
 int a,b,Nmatch;
 struct MATCH *mn;
 int memory;

 /*
  making all the heavy atom pairs having the same atomtype. 
 */ 
 memory = 0;
 
 Mlist->num  = -1;
 Mlist->next = NULL;
 Nmatch = 0; 
 for (a=0;a<molA->Natom;++a){
  if ((molA->atoms[a].atomtype[0]!='\0')&&(molA->atoms[a].one_char_ele!='H')){
   for (b=0;b<molB->Natom;++b){
    if ((molB->atoms[b].atomtype[0]!='\0')&&(molB->atoms[b].one_char_ele!='H')){
      if (atomtype_match(molA->atoms[a].atomtype,molB->atoms[b].atomtype)==1){  
     /* if (strcmp(molA->atoms[a].atomtype,molB->atoms[b].atomtype)==0){  */
        /* printf("'%s' vs '%s'\n",molA->atoms[a].atomtype,molB->atoms[b].atomtype);  */
        mn = Mlist->next; 
        Mlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
        memory += sizeof(struct MATCH);
        Mlist->next->next = mn;
        if (mn!=NULL) mn->prev = Mlist->next; 
        Mlist->next->prev = Mlist;
        mn = Mlist->next;
        mn->num = Nmatch;
        Malloc_MATCH(mn,1);
        mn->Npair = 1;
        mn->anumA[0] = a;
        mn->anumB[0] = b;
        mn->select_dis =  cal_select_dis_for_MATCH(mn,molA,molB);
        /* printf(" (%d %d) select_dis %f\n",a+1,b+1,mn->select_dis);  */
        mn->nodetype = 'N';
        mn->Ncomponent = 1; 
        mn->ConnectGraphType  = PAR.ConnectGraphType;
        mn->maxDIFtopodis     = PAR.maxDIFtopodis;
        mn->atomtype_class    = PAR.atomtype_class;

        ++Nmatch;
      }
     }
   }
  } 
 }
 /*
 printf("#Make_All_Atom_Pair_MATCH(Mlist,molA,molB) Nmatch %d memory %d byte (MATCH %d byte)\n",Nmatch,memory,sizeof(struct MATCH));
 */
 return(Nmatch);

} /* end of Make_All_Atom_Pair_MATCH() */





int Make_All_Atom_Pair_MATCH_for_Substructure_or_Isomorphic(Mlist,molA,molB,SItype)
  struct MATCH    *Mlist; 
  struct MOLECULE *molA,*molB;
  char SItype;  /* 'S'ubstructure, 'I'somorphic */
{
 int a,b,j,Nmatch;
 struct MATCH *mn;
 int memory;
 char neighbor_ok,match_a,ec_ok;
 /*
  making all the heavy atom pairs having the same atomtype. 
 */ 
 memory = 0;
 
 Mlist->num  = -1;
 Mlist->next = NULL;
 Nmatch = 0; 
 for (a=0;a<molA->Natom;++a){
  if ((molA->atoms[a].atomtype[0]!='\0')&&(molA->atoms[a].one_char_ele!='H')){
   match_a = 0;
   for (b=0;b<molB->Natom;++b){
    if ((molB->atoms[b].atomtype[0]!='\0')&&(molB->atoms[b].one_char_ele!='H')){
      if (atomtype_match(molA->atoms[a].atomtype,molB->atoms[b].atomtype)==1){  
      /* if (strcmp(molA->atoms[a].atomtype,molB->atoms[b].atomtype)==0){ */
        neighbor_ok = 1;
        ec_ok = 1;
        j = 0;
        while ((neighbor_ok==1)&&(j<PAR.max_atomtype)){
          /*
          printf("a %d b %d atomtype %d A %d B %d\n",a,b,j,molA->atoms[a].Nnei_atomtype[j],molB->atoms[b].Nnei_atomtype[j]);
          */ 
          if ((SItype=='S')&&(molA->atoms[a].Nnei_atomtype[j] >  molB->atoms[b].Nnei_atomtype[j])) neighbor_ok = 0;
          if ((SItype=='I')&&(molA->atoms[a].Nnei_atomtype[j] != molB->atoms[b].Nnei_atomtype[j])) neighbor_ok = 0;
          j += 1;
        }
 
        if (SItype=='I'){
           if (PAR.levelEC_Dextcon=='0'){
            if ((molA->atoms[a].EC0 != molB->atoms[b].EC0) || (molA->atoms[a].EC1 != molB->atoms[b].EC1)) ec_ok = 0;
           }

           else if (PAR.levelEC_Dextcon=='1'){
            if ((molA->atoms[a].EC0 != molB->atoms[b].EC0) || (molA->atoms[a].EC1 != molB->atoms[b].EC1)) ec_ok = 0;
           }
           else if ((PAR.levelEC_Dextcon=='2')||(PAR.levelEC_Dextcon=='3')){
            if ((molA->atoms[a].EC0 != molB->atoms[b].EC0) || (molA->atoms[a].EC1 != molB->atoms[b].EC1) ||
                (molA->atoms[a].EC2 != molB->atoms[b].EC2)) ec_ok = 0;
           }
       }

        if ((neighbor_ok==1) && (ec_ok==1)){
          /* printf("'%s' vs '%s'\n",molA->atoms[a].atomtype,molB->atoms[b].atomtype);  */
          mn = Mlist->next; 
          Mlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
          memory += sizeof(struct MATCH);
          Mlist->next->next = mn;
          if (mn!=NULL) mn->prev = Mlist->next; 
          Mlist->next->prev = Mlist;
          mn = Mlist->next;
          mn->num = Nmatch;
          Malloc_MATCH(mn,1);
          mn->Npair = 1;
          mn->anumA[0] = a;
          mn->anumB[0] = b;
          mn->select_dis =  cal_select_dis_for_MATCH(mn,molA,molB);
          /* printf(" (%d %d) select_dis %f\n",a+1,b+1,mn->select_dis);  */
          mn->nodetype = 'N';
          mn->Ncomponent = 1;  
          mn->ConnectGraphType  = PAR.ConnectGraphType;
          mn->maxDIFtopodis     = PAR.maxDIFtopodis;
          mn->atomtype_class    = PAR.atomtype_class;
          ++Nmatch;
          match_a += 1;
        }
      }
     }
    }
    if (match_a==0){
      Free_MATCHlist(Mlist);
      return(0);
    }
  } 
 }
/*
 printf("#Make_All_Atom_Pair_MATCH_for_Substructure_or_Isomorphic(Mlist,molA,molB) Nmatch %d memory %d byte (MATCH %d byte)\n",Nmatch,memory,sizeof(struct MATCH));
 */ 
 return(Nmatch);

} /* end of Make_All_Atom_Pair_MATCH_for_Substructure_or_Isomorphic() */









float cal_select_dis_for_MATCH(M,molA,molB)
  struct MATCH    *M;
  struct MOLECULE *molA,*molB;
{
 int p,q,a,b,x,y,j,Nnei_diff,EC_diff,Dtopo_diff,rank_diff;
 int min;

 Nnei_diff = EC_diff = Dtopo_diff = rank_diff = 0;
 for (a=0;a<PAR.max_atomtype;++a) M->Natomtype[a] = 0;


 for (p=0;p<M->Npair;++p){
  a = M->anumA[p];
  b = M->anumB[p];
  M->Natomtype[Number_of_atomtype(molA->atoms[a].atomtype)] += 1;

  if (PAR.Wneiatm>0.0){
    for (j=0;j<PAR.max_atomtype;++j){
      Nnei_diff += abs(molA->atoms[a].Nnei_atomtype[j] - molB->atoms[b].Nnei_atomtype[j]);
    }
    if (PAR.includeHinMCS=='T'){
      Nnei_diff += abs ( (molA->atoms[a].Nneighbor - molA->atoms[a].Nnei_heavy)
                        -(molB->atoms[b].Nneighbor - molB->atoms[b].Nnei_heavy) );
    }
  }

  if (PAR.Wextcon>0.0){
         if (PAR.levelEC_Dextcon=='0')
      EC_diff   += abs(molA->atoms[a].EC0 - molB->atoms[b].EC0);
    else if (PAR.levelEC_Dextcon=='1')
      EC_diff   += abs(molA->atoms[a].EC1 - molB->atoms[b].EC1);
    else if (PAR.levelEC_Dextcon=='2')
      EC_diff   += abs(molA->atoms[a].EC2 - molB->atoms[b].EC2);
    else if (PAR.levelEC_Dextcon=='3')
      EC_diff   += abs(molA->atoms[a].EC3 - molB->atoms[b].EC3);
  }

  if (PAR.Wrank>0.0){
      rank_diff   += abs(molA->atoms[a].rank - molB->atoms[b].rank);
   }
 }

/*
 for (a=0;a<PAR.max_atomtype;++a) printf("%d ",M->Natomtype[a]);
 printf("\n");
 */

 if (PAR.Wtopodis>0.0){
   for (p=0;p<M->Npair;++p){
     a = M->anumA[p];
     b = M->anumB[p];
     for (q=p+1;q<M->Npair;++q){
       x = M->anumA[q];
       y = M->anumB[q];
       /*
       Dtopo_diff += abs(molA->topodismap.map[a][x] - molB->topodismap.map[b][y]);
       Dtopo_diff += 2.0 * abs(molA->topodismap.map[a][x] - molB->topodismap.map[b][y])/(molA->topodismap.map[a][x] + molB->topodismap.map[b][y]);
       Dtopo_diff += abs(molA->topodismap.map[a][x] - molB->topodismap.map[b][y])/sqrt(molA->topodismap.map[a][x] * molB->topodismap.map[b][y]);
       */
       if (molA->topodismap.map[a][x]<molB->topodismap.map[b][y])
         min = molA->topodismap.map[a][x];
       else
         min = molB->topodismap.map[b][y];
       if (min<=PAR.maxD_Dtopodis) { Dtopo_diff += abs(molA->topodismap.map[a][x] - molB->topodismap.map[b][y]); }
     }
   }
     
 }

 M->Dneiatm  =  (float)Nnei_diff;
 M->Dextcon  =  (float)EC_diff;
 M->Dtopodis =  (float)Dtopo_diff;

 /*
 printf("Nnei_diff %d EC_diff %d Dtopo_diff %d Ncomponent %d\n", Nnei_diff,EC_diff,Dtopo_diff,M->Ncomponen);
 */
 return(PAR.Wneiatm * Nnei_diff + PAR.Wextcon * EC_diff  + PAR.Wtopodis * Dtopo_diff + PAR.Wncompo * M->Ncomponent + PAR.Wrank * rank_diff);
} /* end of cal_select_dis_for_MATCH() */



void cal_select_dis_for_MATCHlist(Mlist,molA,molB)
  struct MATCH    *Mlist;
  struct MOLECULE *molA,*molB;
{
 struct MATCH *mn;
 mn = Mlist;
 /* printf("#void cal_select_dis_for_MATCHlist(Mlist,molA,molB)\n"); */
 while ((mn->next != NULL) && (mn->next->nodetype!='E')){
   mn = mn->next;
   mn->select_dis =  cal_select_dis_for_MATCH(mn,molA,molB); 
 }
} /* end of cal_select_dis_for_MATCHlist() */






void Malloc_MATCHlist(Mlist,Npair_max,lenMATCHlist)
  struct MATCH *Mlist;
  int    Npair_max;  /* length of anumA and anumB */
  int    lenMATCHlist;     /* number of 'struct MATCH' in linked list 'Mlist' */
{
 int m,n;
 struct MATCH *mn;
 /*
 printf("#Malloc_MATCHlist(M,Nmax_pair %d Nmatch %d)\n",Nmax_pair,Nmatch);
 */ 
 Mlist->next = NULL;
 Mlist->nodetype = 'I';
 for (m=0;m<lenMATCHlist;++m){
   mn = Mlist->next;
   Mlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
   Mlist->next->next = mn;
   Mlist->next->prev = Mlist;
   mn = Mlist->next;
   Malloc_MATCH(mn,Npair_max);
   for (n=0;n<Npair_max;++n) mn->anumA[n] = mn->anumB[n] = -1; 
   mn->Npair = 0;
   mn->select_dis = Npair_max * 1000.0;
   mn->nodetype   = 'E';
   mn->ConnectGraphType  = PAR.ConnectGraphType;
   mn->maxDIFtopodis     = PAR.maxDIFtopodis;
   mn->atomtype_class    = PAR.atomtype_class;
 } 
} /* end of Malloc_MATCHlist() */




void Add_MATCH_to_MATCHlist(newM,Mlist)
  struct MATCH *newM;
  struct MATCH *Mlist;
{
  struct MATCH *nnext;  

 /*
  printf("#Add_MATCH_to_MATCHlist(newM,Length of Mlist:%d)\n",Length_of_MATCHlist(Mlist));
 */
  nnext = Mlist->next;
  Mlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
  Malloc_MATCH(Mlist->next,newM->Npair);
  Copy_MATCH(Mlist->next,newM);
  Mlist->next->nodetype = 'N'; 
  Mlist->next->next = nnext; 
  Mlist->next->prev = Mlist; 
  Mlist->next->ConnectGraphType  = PAR.ConnectGraphType;
  Mlist->next->maxDIFtopodis     = PAR.maxDIFtopodis;
  Mlist->next->atomtype_class    = PAR.atomtype_class;
  if (nnext != NULL){
    nnext->prev = Mlist->next;
  }

 /*
  printf("#--> Mlist:%d\n",Length_of_MATCHlist(Mlist));
 */

} /* end of Add_MATCH_to_MATCHlist() */



int Add_MATCH_to_MATCHlist_Keep_Only_Best(newM,Mlist)
  struct MATCH *newM;
  struct MATCH *Mlist;
{
  struct MATCH *nnext;  
  int Nbest_orig;

  if (Mlist->next==NULL) Nbest_orig = 0;
                   else  Nbest_orig = Mlist->next->Npair;

  if ((newM->Npair > Nbest_orig) && (Nbest_orig>0)){
     /* printf("#Free_MATCHlist\n"); fflush(stdout); */
     Free_MATCHlist(Mlist);
  }

  if ((Mlist->next==NULL)||(newM->Npair>=Nbest_orig)){   
    nnext = Mlist->next;
    Mlist->next = (struct MATCH*)malloc(sizeof(struct MATCH));
    Malloc_MATCH(Mlist->next,newM->Npair);
    Copy_MATCH(Mlist->next,newM);
    Mlist->next->nodetype = 'N'; 
    Mlist->next->next = nnext; 
    Mlist->next->prev = Mlist; 
    Mlist->next->ConnectGraphType  = PAR.ConnectGraphType;
    Mlist->next->maxDIFtopodis     = PAR.maxDIFtopodis;
    Mlist->next->atomtype_class    = PAR.atomtype_class;
    if (nnext != NULL){
      nnext->prev = Mlist->next;
    }
    return(1); 
  }
  return(0);

} /* end of Add_MATCH_to_MATCHlist_Keep_Only_Best() */





void Free_MATCHlist(Mlist)
  struct MATCH *Mlist;
{
 struct MATCH *mn,*mnext;
 int Nmatch_free;

 Nmatch_free = 0;
 mn = Mlist->next;
 do {
    if (mn!=NULL){
     mnext = mn->next;
     Free_MATCH(mn);
     free(mn);
     ++Nmatch_free;
     mn = mnext; 
   }
 }while (mn != NULL);
 
 Mlist->next = NULL;
 /* printf("#Free_MATCHlist(Nmatch_free %d)\n",Nmatch_free); */
} /* end of Free_MATCHlist() */





void Copy_MATCH(Mnew,Morig)
 struct MATCH *Mnew,*Morig;
{
  int n;


  Mnew->Npair   = Morig->Npair;
  Mnew->select_dis   = Morig->select_dis;
  Mnew->Ncomponent = Morig->Ncomponent;
  Mnew->nodetype = Morig->nodetype;
  Mnew->select_dis = Morig->select_dis;
  Mnew->Dneiatm  = Morig->Dneiatm;
  Mnew->Dextcon  = Morig->Dextcon;
  Mnew->Dtopodis = Morig->Dtopodis;
  Mnew->mark     = Morig->mark;
  Mnew->subfloat = Morig->subfloat;

  Mnew->ConnectGraphType = Morig->ConnectGraphType;
  Mnew->maxDIFtopodis    = Morig->maxDIFtopodis;
  Mnew->atomtype_class   = Morig->atomtype_class;

  for (n=0;n<Morig->Npair;++n){
    Mnew->anumA[n] = Morig->anumA[n]; 
    Mnew->anumB[n] = Morig->anumB[n]; 
   }

  for (n=0;n<PAR.max_atomtype;++n){
    Mnew->Natomtype[n] = Morig->Natomtype[n]; 
  }

} /* end of Copy_MATCH() */



void Copy_MATCHlist(MlistNew,MlistOrig)
  struct MATCH *MlistNew;
  struct MATCH *MlistOrig;
{
 struct MATCH *m0,*m1;
 int n;

 m0 = MlistOrig;
 m1 = MlistNew;

 while ((m0->next!=NULL) && (m1->next!=NULL)){
   m0 = m0->next;
   m1 = m1->next;
   m1->Npair  = m0->Npair;
   m1->select_dis = m0->select_dis;
   m1->Ncomponent   = m0->Ncomponent;
   m1->nodetype     = m0->nodetype;
   m1->select_dis   = m0->select_dis;
   m1->Dneiatm      = m0->Dneiatm;
   m1->Dextcon      = m0->Dextcon;
   m1->Dtopodis     = m0->Dtopodis;
   for (n=0;n<m0->Npair;++n){
     m1->anumA[n] = m0->anumA[n]; 
     m1->anumB[n] = m0->anumB[n]; 
    }
   for (n=0;n<PAR.max_atomtype;++n){
     m1->Natomtype[n] = m0->Natomtype[n]; 
   }
 }

} /* end of Copy_MATCHlist() */








void Malloc_MATCH(m,Npair)
  struct MATCH *m;
  int Npair;
{ int n;

 /* printf("#void Malloc_MATCH(m,Npair %d)\n",Npair); */
 m->Npair        = Npair;
 m->Npair_malloc = Npair;
 m->anumA = (int *)malloc(sizeof(int)*Npair);
 m->anumB = (int *)malloc(sizeof(int)*Npair);
 /* 
 m->anumA = (unsigned char *)malloc(sizeof(unsigned char)*Npair);
 m->anumB = (unsigned char *)malloc(sizeof(unsigned char)*Npair);
 */ 

 for (n=0;n<m->Npair_malloc;++n){ 
   m->anumA[n] = m->anumB[n] = -1; 
 }
 m->mark = 0;
 m->ConnectGraphType = PAR.ConnectGraphType;
 m->maxDIFtopodis    = PAR.maxDIFtopodis;
 m->atomtype_class   = PAR.atomtype_class;

} /* end of  Malloc_MATCH() */


void Free_MATCH(m)
  struct MATCH *m;
{
 m->Npair_malloc = m->Npair = 0;
 m->next = m->prev = NULL; 
 free(m->anumA);
 free(m->anumB);
} /* end of  Free_MATCH() */



int Length_of_MATCHlist(Mlist)
 struct MATCH *Mlist;
{
 struct MATCH *mn;
 int N;
 N = 0;
 mn = Mlist; 
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){
   mn = mn->next;
   N += 1; 
 }
 return(N);
} /* end of Length_of_MATCHlist() */


int Nmark_in_MATCHlist(Mlist)
 struct MATCH *Mlist;
{
 struct MATCH *mn;
 int N;
 N = 0;
 mn = Mlist; 
 while ((mn->next != NULL)&&(mn->next->nodetype!='E')){
   mn = mn->next;
   if (mn->mark==1) N += 1; 
 }
 return(N);
} /* end of Nmark_in_MATCHlist() */


float Tanimoto_Coefficient(Ncommon,molA,molB)
 int Ncommon;
 struct MOLECULE *molA,*molB;
{
  int NA,NB;

  if (PAR.includeHinMCS=='T'){
    NA = molA->Natom; NB = molB->Natom;
  } 
  else{
    NA = molA->Nheavyatom; NB = molB->Nheavyatom;
  }
  return((float)(Ncommon)/(float)(NA+NB-Ncommon));
}


int Smaller_Int(A,B)
 int A,B;
{ if (A<B) return(A); else return(B); }

int Larger_Int(A,B)
 int A,B;
{ if (A>B) return(A); else return(B); }



int Common_Pair_Bwn_Two_MATCHes(mI,mJ)
  struct MATCH *mI;
  struct MATCH *mJ;
{
  int i,j,NIandJ;
  char hit;

  NIandJ = 0;
  for (i=0;i<mI->Npair;++i){
    j = 0;
    hit = 0;
    while ((j<mJ->Npair)&&(hit==0)){
      if ((mI->anumA[i]==mJ->anumA[j]) && (mI->anumB[i]==mJ->anumB[j])) hit = 1;
      j += 1;
    }
    if (hit==1) NIandJ += 1;
  }
  return(NIandJ);
} /* end of Common_Pair_Bwn_Two_MATCHes() */




void Keep_Only_Nonredundant_MATCH(Mlist)
  struct MATCH *Mlist;
{
  struct MATCH *mI,*mJ;
  int Ncommon_pair;

  /** [1] Initialize mark **/ 
  mI = Mlist;
  while ((mI->next!=NULL)&&(mI->next->nodetype!='E')){
    mI = mI->next;
    mI->mark = 0; 
  }

  /** [2] Cal Ncommon_pair for all-vs-all **/ 
  mI = Mlist;
  while ((mI->next!=NULL)&&(mI->next->nodetype!='E')){
    mI = mI->next;
    mJ = mI;
    while ((mJ->next!=NULL)&&(mJ->next->nodetype!='E')){
      mJ = mJ->next;
      if (mI->Npair == mJ->Npair){
        Ncommon_pair = Common_Pair_Bwn_Two_MATCHes(mI,mJ);
        if (Ncommon_pair==mI->Npair) mJ->mark += 1; 
      }
    }
  }

  /** [3] Remove non-zero marked MATCH **/ 
  mI = Mlist;
  while ((mI->next!=NULL)&&(mI->next->nodetype!='E')){
    mI = mI->next;
    if (mI->mark > 0){
      mJ = mI;
      if (mI->prev != NULL) mI->prev->next = mI->next;
      if (mI->next != NULL) mI->next->prev = mI->prev;
      mI = mI->prev; 
      Free_MATCH(mJ);
    }
  }
} /* end of Keep_Only_Nonredundant_MATCH() */



void Keep_Only_Substructure_or_Isomorphic_MATCH(Mlist,molA, molB,ConnectType)
  struct MATCH *Mlist;
  struct MOLECULE *molA, *molB;
  char ConnectType;   /* 'S'ubstructure, 'I'somorphic */
{
  struct MATCH *mI, *mJ;

  /** [1] Initialize mark **/ 
  mI = Mlist;
  while ((mI->next!=NULL)&&(mI->next->nodetype!='E')){
    mI = mI->next;
    mI->mark = 0;
    if (ConnectType=='S'){ 
     if (mI->Npair<molA->Nheavyatom) mI->mark = 1;
    }
    if (ConnectType=='I'){ 
       if ((mI->Npair<molA->Nheavyatom)||(mI->Npair<molB->Nheavyatom)) mI->mark = 1; 
    } 
  }

  /** [2] Remove non-zero marked MATCH **/ 
  mI = Mlist;
  while ((mI->next!=NULL)&&(mI->next->nodetype!='E')){
    mI = mI->next;
    if (mI->mark == 1){
      mJ = mI;
      if (mI->prev != NULL) mI->prev->next = mI->next;
      if (mI->next != NULL) mI->next->prev = mI->prev;
      mI = mI->prev; 
      Free_MATCH(mJ); 
    }
  }
} /* end of Keep_Only_Substructure_or_Isomorphic_MATCH() */











void Bubble_Sort_anumAB_in_MATCHlist_by_anumA(Mlist)
  struct MATCH    *Mlist;
{
 struct MATCH *mn;
 int i,change,buff;

 mn = Mlist;
 while ((mn->next != NULL) && (mn->next->nodetype!='E')){
   mn = mn->next;
   do{
     change = 0;
     for (i=1;i<mn->Npair;++i){
       if (mn->anumA[i-1]>mn->anumA[i]){
         change += 1;
         buff           = mn->anumA[i-1];
         mn->anumA[i-1] = mn->anumA[i];
         mn->anumA[i]   = buff;
         buff           = mn->anumB[i-1];
         mn->anumB[i-1] = mn->anumB[i];
         mn->anumB[i]   = buff;
       }
     }
   } while (change>0);
 
 }

} /* end of Bubble_Sort_anumAB_in_MATCHlist_by_anumA() */




int atomtype_match(atomtypeA,atomtypeB)
  char *atomtypeA,*atomtypeB;
/*

  match     --> return(1)
  not_match --> return(0)

*/
{
  if (PAR.atomtype_class=='D'){
         if (atomtypeA[0]==atomtypeB[0]) return(1);
    else if ((atomtypeA[0]=='B') && ((atomtypeB[0]=='D')||(atomtypeB[0]=='A'))) return(1);
    else if ((atomtypeB[0]=='B') && ((atomtypeA[0]=='D')||(atomtypeA[0]=='A'))) return(1);
    else return(0);
  }
  else{
    if (strcmp(atomtypeA,atomtypeB)==0) return(1);
    else return(0);
  }
}




void add_hydrogen_pairs_to_MATCH(Mlist,mA,mB)
  struct MATCH    *Mlist; 
  struct MOLECULE *mA,*mB;
{
  struct MATCH *m;
  int i,j,a,b;
  int neihydA[10],neihydB[10],NneihydA,NneihydB; /* hydrogen neighbors */
  int Npair_new;
  unsigned char *anumA_new,*anumB_new;          /* new matchings */
  printf("#add_hydrogen_pairs_to_MATCH(NatomA %d NatomB %d)\n",mA->Natom,mB->Natom);

  anumA_new = (unsigned char *)malloc(sizeof(unsigned char)*mA->Natom);
  anumB_new = (unsigned char *)malloc(sizeof(unsigned char)*mA->Natom);

  m = Mlist; 
  while ((m->next != NULL) && (m->next->nodetype!='E')){
    m = m->next;
    /* printf("#Npair %d\n",m->Npair); */
    Npair_new = m->Npair;
    for (i=0;i<m->Npair;++i){ 
      anumA_new[i] = m->anumA[i]; 
      anumB_new[i] = m->anumB[i];
    }

    for (i=0;i<m->Npair;++i){
      a = m->anumA[i];
      b = m->anumB[i];
      /* printf("a %d b %d\n",m->anumA[i],m->anumB[i]);  */
      /** (1) set up hydrogen neighbors for (a,b)**/ 
      NneihydA = NneihydB = 0;
      for (j=0;j<mA->Natom;++j){
         /* printf("i %d/%d a %d j %d Natom %d\n",i,m->Npair,a,j,mA->Natom); fflush(stdout); */
        if ((j!=a)&&(mA->conmap.map[a][j] != '0') && (mA->atoms[j].one_char_ele == 'H')){
          neihydA[NneihydA] = j;  NneihydA += 1;
        }
      } 

      for (j=0;j<mB->Natom;++j){
        if ((j!=b)&&(mB->conmap.map[b][j] != '0') && (mB->atoms[j].one_char_ele == 'H')){
          neihydB[NneihydB] = j;  NneihydB += 1;
        }
      } 

      /** (2) set up hydrogen pairs around the pair (a,b)**/ 
      if (NneihydA==NneihydB){
         for (j=0;j<NneihydA;++j){
           anumA_new[Npair_new] = neihydA[j];
           anumB_new[Npair_new] = neihydB[j];
           Npair_new += 1;
         }
      }
    } /* i */

    if (Npair_new > m->Npair){ 
      free(m->anumA);
      free(m->anumB);
/*
      m->anumA = (unsigned char *)malloc(sizeof(unsigned char)*Npair_new);
      m->anumB = (unsigned char *)malloc(sizeof(unsigned char)*Npair_new);
 */
      m->anumA = (int *)malloc(sizeof(unsigned char)*Npair_new);
      m->anumB = (int *)malloc(sizeof(unsigned char)*Npair_new);
      for (i=0;i<Npair_new;++i){ 
        m->anumA[i] = anumA_new[i]; 
        m->anumB[i] = anumB_new[i];
      }
      m->Npair = Npair_new;
    }
  } /* while (m) */


  free(anumA_new);
  free(anumB_new);

} /* end of add_hydrogen_pairs_to_MATCH() */


float Compare_Electric_Charges(m,molA,molB)
  struct MATCH *m;
  struct MOLECULE *molA, *molB;
{
  int i,a,b,N;
  float A,B,AB,AA,BB,qA,qB,hodgkin,CC;
  float mA,vA,mB,vB,cAB;
  printf("#Compare_Electric_Charges(m,molA,molB)\n");
  N = 0;
  A = B = AB = AA = BB = 0.0;
  printf("#[atmA] -- [atmB]  [chargeA] [chargeB]\n");
  for (i=0;i<m->Npair;++i){ 
    a = m->anumA[i]; 
    b = m->anumB[i];
    if ((a>=0) && (b>=0)){
      qA = molA->atoms[a].charge;
      qB = molB->atoms[b].charge;
      printf("#%-3d %-2s -- %-3d %-2s %9.5f %9.5f\n",
        molA->atoms[a].num_in_file, molA->atoms[a].element, 
        molB->atoms[b].num_in_file, molB->atoms[b].element, qA,qB);
      A += qA;
      B += qB;
      AB += qA * qB;
      AA += qA * qA; 
      BB += qB * qB;
      N += 1;
    }
  }
  mA = A/N;
  mB = B/N;
  vA = AA/N - mA*mA;
  vB = BB/N - mB*mB;
  cAB = AB/N - mA*mB;
  CC = cAB/sqrt(vA)/sqrt(vB);

  hodgkin = 0.0; 
  if ((AA+BB)>0.0) hodgkin = 2.0 * AB /(AA+BB);

  printf("#sumQA %f sumQB %f CorrCoeff %f Hodgkin %f\n",A,B,CC,hodgkin);

  return(hodgkin);
} /* end of Compare_Electric_Charges() */




void Sort_MATCHlist_by_Npair_and_select_dis(Mlist)
  struct MATCH    *Mlist; 
{
  struct MATCH *m0,*m1,*m00,*m11;
  int Nchange;

/*
  printf("#void Sort_MATCHlist_by_Npair_and_select_dis(Mlist)\n");
 */

  /*
  m = Mlist;
  while ((m->next != NULL) && (m->next->nodetype != 'E')){
    m = m->next;
    printf("#PRE_MATCH Npair %d select_dis %f\n",m->Npair,m->select_dis);
  }
  */


  /*
  ## ORDER for m0 -> m1 ##
   Npair(m0) > Napir(m1)
   if (Npair(m0) == Napir(m1)):
     select_dis(m0) < select_dis(m1)

  #### BUBBLE SORT #### 
 
  (before) 
   -m00-> m0 -> m1 -> m11
 
  (after) 
   -m00-> m1 -> m0 -> m11

    

  */


  do{
    Nchange = 0;
    m0 = Mlist;
    while ((m0->next !=NULL) && (m0->next->nodetype != 'E') && (Nchange==0)){
      m0  = m0->next;
      m1  = m0->next;
      if ((m1 !=NULL) && (m1->nodetype != 'E')&&
          ( (m0->Npair < m1->Npair) ||
           ((m0->Npair == m1->Npair)&&(m0->select_dis > m1->select_dis))) ){
       /*  
        printf("# Npair %d %d select_dis %f %f\n",m0->Npair,m1->Npair,m0->select_dis,m1->select_dis);
        */
        m00 = m0->prev;
        m11 = m1->next;

        m1->next = m0;
        m1->prev = m00;
        m0->prev = m1;
        m0->next = m11;

        if (m00 != NULL){ m00->next = m1;}
        if (m11 != NULL){ m11->prev = m0;}

        Nchange += 1;

      }
    }
  } while (Nchange>0);

  /*
  m0 = Mlist;
  while ((m0->next != NULL) && (m0->next->nodetype != 'E')){
    m0 = m0->next;
    printf("#SORT_MATCH Npair %d select_dis %f\n",m0->Npair,m0->select_dis);
  }
  */

} /* end of Sort_MATCHlist_by_Npair_and_select_dis() */


