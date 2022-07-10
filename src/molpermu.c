/**
 
 <molpermu.c>

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


 Functions for finding symmetry in the molecule.
 by generating all the permutaions of heavy atoms

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "match.h"
#include "buildup.h"
#include "io_match.h"
#include "MrgSrtMATCH.h"
#include "qRMS.h"

/** FUNCTIONS (GLOBAL) **/
int Make_All_Permutation_of_Equivalent_Heavyatoms();
int Make_EquivMap_for_Heavyatoms();
int Number_of_EquivMap_for_Heavyatoms();
void show_atom_permutation();
void show_all_the_atom_permutations();
int Number_of_PERMUTATION();
int Make_Permutated_MATCHlist_from_Top_MATCH();
int Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation();
int Make_NonRedundant_MATCHlist_by_Permutation();
int Check_Equivalence_bwn_Two_MATCHes_by_Permutation();
void Free_PERMUTATION();
void Rerank_MATCHlist_by_transformed_RMSD();

/** FUNCTIONS (LOCAL) **/
static void swap();
static int heavy_atom_number_permutation();
static int check_conmap_equivalence();
static int check_Nnei_atomtype();
static void Add_PERMUTATION();
static int number_of_permutations();
 

/*** SAMPLE OF THE SIMPLEST PERMUTATION FUNCTIONS ***/

/*
static int permutation();

int permutation(ind,depth,Lind)
  int *ind;
  int depth;
  int Lind;
{
  int i;

  if (depth>=Lind){
    ## here, output ind[] ## 
    return(1);
  }

  i = depth; 
  while (i<Lind){
    swap(ind,depth,i);
    permutation(ind,depth+1,Lind);
    swap(ind,depth,i); 
    i += 1;
  } 

}
*/


void swap(ind,i,j)
  int *ind; /* index */
  int i,j;
{
  int tmp;
  tmp = ind[i];
  ind[i] = ind[j]; 
  ind[j] = tmp; 
}


void show_atom_permutation(mol,ind)
  struct MOLECULE *mol;
  int *ind; /* index */
{
  int i,ii,jj;
  printf("#");
  for (i=0;i<mol->Nheavyatom;++i) {
   /* show only permutated sites */
    if (ind[i]!=i){
      ii = mol->num_frm_heavy[i];
      jj = mol->num_frm_heavy[ind[i]];
      printf("(%c %d-",mol->atoms[ii].one_char_ele,mol->atoms[ii].num_in_file); 
      printf("%c %d)",mol->atoms[jj].one_char_ele,mol->atoms[jj].num_in_file); 
    /*
      printf("(%c %d %d-",mol->atoms[ii].one_char_ele,mol->atoms[ii].num_in_file,ii); 
      printf("%c%d %d)",mol->atoms[jj].one_char_ele,mol->atoms[jj].num_in_file,jj); 
     */
    }
 }
 printf("\n");
}


void show_all_the_atom_permutations(mol)
  struct MOLECULE *mol;
{
  int n;
  struct PERMUTATION *pn;
  n = 0;
  pn = &(mol->permuhead);
  while (pn->next != NULL){
    pn = pn->next;
    printf("#[%3d]",n);
    show_atom_permutation(mol,pn->index_heavy);
    n += 1;
   }
}


int number_of_permutations(permuhead)
  struct PERMUTATION *permuhead;
{
  int n;
  struct PERMUTATION *pn;
  n = 0;
  pn = permuhead;
  while (pn->next != NULL){
    pn = pn->next;
    n += 1;
   }
  return(n);
}



int heavy_atom_number_permutation(mol,Eq,ind,depth,Lind,over_maxPermutations)
  struct MOLECULE *mol;
  struct CHAR2DMAP *Eq; /* EquivMap_for_Heavyatoms */
  int *ind;             /* numbers of heavy atom numbers */
  int depth;
  int Lind; /* length of ind[] */
  char *over_maxPermutations; /* if 'T', over_maxPermutations */
{
  int j;

  if (*over_maxPermutations=='T') return(0);
/*
  printf("#int heavy_atom_number_permutation(mol,index,depth:%d Lind:%d maxPermutations %d)\n",depth,Lind,PAR.maxPermutations); fflush(stdout);
 */
   if (depth>=Lind){
     if (check_conmap_equivalence(mol,ind,mol->Nheavyatom)==1){
       if (number_of_permutations(&(mol->permuhead))<PAR.maxPermutations){
         Add_PERMUTATION(&(mol->permuhead),ind,mol->Nheavyatom);
         /* show_atom_permutation(mol,ind);   */
         return(1);
        }
        else { *over_maxPermutations = 'T';}
     } else {
      /*
      printf("FAIL");
      show_atom_permutation(mol,ind);
      */
    }
    return(0);
  }

  
  j = depth; 
  while (j<Lind){
    if ((Eq->map[depth][ind[j]]==1)&&(Eq->map[j][ind[depth]]==1)){
      swap(ind,depth,j);
      if (check_conmap_equivalence(mol,ind,depth+1)==1){
         heavy_atom_number_permutation(mol,Eq,ind,depth+1,Lind,over_maxPermutations);
      }
      swap(ind,depth,j); 
    }
   j += 1;
  }  
  
  return(0);

}  /* end of  heavy_atom_number_permutation() */



int check_conmap_equivalence(mol,index_heavy,Nupper)
  struct MOLECULE *mol;
  int *index_heavy;
  int Nupper;
{
  int jj,kk,j0,j1,k0,k1; 

  if (Nupper>mol->Nheavyatom){
    printf("#ERROR:Nupper %d is over Nheavyatom %d.\n",Nupper,mol->Nheavyatom);
    exit(1);
  } 

  for (jj=0;jj<Nupper;++jj){
    j0  = mol->num_frm_heavy[jj];
    j1  = mol->num_frm_heavy[index_heavy[jj]];

    if ((j0<0) || (j0>mol->Natom)) {printf("#ERROR:atom num error(%d)\n",j0); exit(1);} 
    if ((j1<0) || (j1>mol->Natom)) {printf("#ERROR:atom num error(%d)\n",j1); exit(1);} 
     
    for (kk=jj+1;kk<Nupper;++kk){
      k0 = mol->num_frm_heavy[kk];
      k1 = mol->num_frm_heavy[index_heavy[kk]];
      if ((k0<0) || (k0>mol->Natom))  {printf("#ERROR:atom num error(%d)\n",k0); exit(1);} 
      if ((k1<0) || (k1>mol->Natom)) {printf("#ERROR:atom num error(%d)\n",k1); exit(1);} 
      if ((mol->conmap.map[j0][k0]=='0')&&(mol->conmap.map[j1][k1]!='0')) return(0);
      if ((mol->conmap.map[j0][k0]!='0')&&(mol->conmap.map[j1][k1]=='0')) return(0);
    }
  }
  return(1);
} /* end of check_conmap_equivalence() */





int check_Nnei_atomtype(mol,i,j)
  struct MOLECULE *mol;
  int i,j;
{
  int m;
  for (m=0;m<PAR.max_atomtype;++m){
    if (mol->atoms[i].Nnei_atomtype[m] != mol->atoms[j].Nnei_atomtype[m]) {return(0);}
  }
  return(1);
}




int Make_EquivMap_for_Heavyatoms(mol,Eq)
  struct MOLECULE *mol;
  struct CHAR2DMAP *Eq; 
{
 int i,j;  /* heavy atom number */
 int ii,jj; /* atom number */
 int Nequivpair;
 char ECok; 

 for (i=0;i<mol->Nheavyatom;++i){
   for (j=0;j<mol->Nheavyatom;++j){ Eq->map[i][j] = 0; }
 }

  /* 
 printf("#Natom %d Nheavyatom %d\n",mol->Natom,mol->Nheavyatom);
 for (i=0;i<mol->Nheavyatom;++i){
   ii = mol->num_frm_heavy[i];
   printf("#heavy_num_i %d atom_num_ii %d element '%s'\n",i,ii,mol->atoms[i].element);
 }
 exit(1);
 */
 Nequivpair = 0;
 for (i=0;i<mol->Nheavyatom;++i){
   ii = mol->num_frm_heavy[i];
   Eq->map[i][i] = 1;
   for (j=i+1;j<mol->Nheavyatom;++j){
     jj = mol->num_frm_heavy[j];
/*
     printf("#element i %d ii %d '%s' j %d jj %d '%s'\n",i,ii,mol->atoms[ii].element,j,jj,mol->atoms[jj].element);
*/
     if ((strcmp(mol->atoms[ii].element,mol->atoms[jj].element)==0) &&
        (strncmp(mol->atoms[ii].atomtype,mol->atoms[jj].atomtype,4)==0) &&
        (check_Nnei_atomtype(mol,ii,jj)==1) ){ 

        ECok = 1;
        if ((PAR.levelEC_Dextcon=='1')&&((mol->atoms[ii].EC0!=mol->atoms[jj].EC0)||
            (mol->atoms[ii].EC1!=mol->atoms[jj].EC1))) ECok = 0; 

        if ((PAR.levelEC_Dextcon=='2')&&((mol->atoms[ii].EC0!=mol->atoms[jj].EC0)||
            (mol->atoms[ii].EC1!=mol->atoms[jj].EC1)||(mol->atoms[ii].EC2!=mol->atoms[jj].EC2))) ECok = 0; 

        if ((PAR.levelEC_Dextcon=='3')&&((mol->atoms[ii].EC0!=mol->atoms[jj].EC0)||
            (mol->atoms[ii].EC1 != mol->atoms[jj].EC1)||(mol->atoms[ii].EC2 != mol->atoms[jj].EC2) ||  
            (mol->atoms[ii].EC3 != mol->atoms[jj].EC3))) ECok = 0; 

        if ((PAR.levelEC_Dextcon=='4')&&((mol->atoms[ii].EC0!=mol->atoms[jj].EC0)||
            (mol->atoms[ii].EC1 != mol->atoms[jj].EC1)||(mol->atoms[ii].EC2 != mol->atoms[jj].EC2) ||  
            (mol->atoms[ii].EC3 != mol->atoms[jj].EC3)||(mol->atoms[ii].EC4 != mol->atoms[jj].EC4))) ECok = 0; 

        if (PAR.Wrank>0.0){
          if (mol->atoms[ii].rank != mol->atoms[jj].rank) ECok = 0;
        }

        if (ECok==1){
        /*
         printf("#equiv pair %d %d %s %s %s %s EC %d %d\n",mol->atoms[ii].num_in_file, mol->atoms[jj].num_in_file,
           mol->atoms[ii].atomname, mol->atoms[jj].atomname,
           mol->atoms[ii].atomtype, mol->atoms[jj].atomtype,
           mol->atoms[ii].EC2, mol->atoms[jj].EC2
          );
         */
         Eq->map[i][j] = Eq->map[j][i] = 1;
         Nequivpair += 1;
       }
    }
  }
 }

 return(Nequivpair);
} /* end of Make_EquivMap_for_Heavyatoms() */




int Number_of_EquivMap_for_Heavyatoms(mol)
  struct MOLECULE *mol;
{
 int i,j,ii,jj,Nequivpair;
 
 Nequivpair = 0;
 for (i=0;i<mol->Nheavyatom;++i){
   ii = mol->num_frm_heavy[i];
   for (j=i+1;j<mol->Nheavyatom;++j){
     jj = mol->num_frm_heavy[j];
     if ((strcmp(mol->atoms[ii].element,mol->atoms[jj].element)==0) &&
        (strncmp(mol->atoms[ii].atomtype,mol->atoms[jj].atomtype,4)==0) &&
        (strncmp(mol->atoms[ii].atomtype, mol->atoms[jj].atomtype,4)==0) &&
        (mol->atoms[ii].EC1==mol->atoms[jj].EC1) &&
        (mol->atoms[ii].EC2==mol->atoms[jj].EC2) && (check_Nnei_atomtype(mol,ii,jj)==1) ){ 
       Nequivpair += 1;
    }
  }
 }

 return(Nequivpair);
} /* end of Number_of_EquivMap_for_Heavyatoms() */



void Add_PERMUTATION(Phead,new_index_heavy,newN)
  struct PERMUTATION *Phead;  /* head of list of PERMUTATION */
  int    *new_index_heavy;
  int    newN;
{
  struct PERMUTATION *nnext;
  int i;

  nnext = Phead->next;
  Phead->next = (struct PERMUTATION*)malloc(sizeof(struct PERMUTATION));
  Phead->next->N = newN;
  Phead->next->index_heavy = (int *)malloc(sizeof(int)*newN);
  for (i=0;i<newN;++i){ 
    Phead->next->index_heavy[i] = new_index_heavy[i]; 
  }
  Phead->next->next = nnext;
  Phead->next->prev = Phead;
  if (nnext != NULL){
     nnext->prev = Phead->next;
  }

} /* end of Add_PERMUTATION() */


void Free_PERMUTATION(HeadPermu)
 struct PERMUTATION *HeadPermu;
{
 struct PERMUTATION *pn;

 pn = HeadPermu;
 printf("HeadPermu %d\n",pn->N); fflush(stdout);
 if (pn->next !=NULL){

   while (pn->next != NULL){ 
     pn = pn->next;
   }

   while (pn->prev!=NULL){
     pn = pn->prev;
     free(pn->next->index_heavy);
     free(pn->next);
   }
   HeadPermu->next = NULL;
 } 
} /* end of Free_PERMUTATION() */



int Number_of_PERMUTATION(Phead)
  struct PERMUTATION *Phead;  /* head of list of PERMUTATION */
{
  struct PERMUTATION *pn;
  int n;
  n = 0;
  pn = Phead;
  while (pn->next != NULL){
    pn = pn->next;
    n += 1;
  }
  return(n);
} /* end of Number_of_PERMUTATION() */



int Make_All_Permutation_of_Equivalent_Heavyatoms(mol)
  struct MOLECULE *mol;
{
  int *ind,i,Nequivpair;
  struct CHAR2DMAP Eq;
  char over_maxPermutations;


  printf("#Make_All_Permutation_of_Equivalent_Heavyatoms('%s')\n",mol->filename);
  mol->permuhead.next = mol->permuhead.prev = NULL;
  mol->permuhead.N = 0;
  mol->permuhead.index_heavy = NULL;
  ind = (int *)malloc(sizeof(int)*mol->Nheavyatom);
  for (i=0;i<mol->Nheavyatom;++i){ind[i] = i;}
  Malloc_CHAR2DMAP(&Eq,mol->Nheavyatom);
  Nequivpair = Make_EquivMap_for_Heavyatoms(mol,&Eq);
  printf("#>Nheavyatom %d Nequivpair %d %6.2f%%\n",
    mol->Nheavyatom,Nequivpair,100.0*2.0*Nequivpair/(mol->Nheavyatom)/(mol->Nheavyatom-1));
  over_maxPermutations = 'F';
  heavy_atom_number_permutation(mol,&Eq,ind,0,mol->Nheavyatom,&over_maxPermutations);
  show_all_the_atom_permutations(mol);
  Free_CHAR2DMAP(&Eq);
  free(ind);
  mol->Npermu = Number_of_PERMUTATION(&(mol->permuhead));
  
  return(mol->Npermu);
} /* end of Make_All_Permutation_of_Equivalent_Heavyatoms() */


int Make_Permutated_MATCHlist_from_Top_MATCH(Mlist,molA,molB)
  struct MATCH *Mlist;
  struct MOLECULE *molA, *molB;
{
  struct PERMUTATION *pA,*pB;
  struct MATCH newM,origM,*mn; 
  int i,Npair,Nmatch,nA,nB,equiv;
 
  if (Mlist->next==NULL) return(0); 


  Npair = Mlist->next->Npair;
  Malloc_MATCH(&origM,Npair);
  Copy_MATCH(&origM,Mlist->next);
  Malloc_MATCH(&newM,Npair);
  Copy_MATCH(&newM,Mlist->next);

  Free_MATCHlist(Mlist);
  Mlist->next = NULL;
 
  Nmatch = 0; 
  nA = 0;
  pA = &(molA->permuhead);
  while (pA->next != NULL){
    pA = pA->next;
    nA += 1;
    nB = 0;
    pB = &(molB->permuhead);
    while (pB->next != NULL){
      pB = pB->next;
      nB += 1;
      for (i=0;i<Npair;++i){
        newM.anumA[i] = pA->index_heavy[molA->atoms[origM.anumA[i]].num_heavy];
        newM.anumB[i] = pB->index_heavy[molB->atoms[origM.anumB[i]].num_heavy];
      } 

      Bubble_Sort_anumA_anumB_array_by_anumA(&newM);
      mn = Mlist;
      equiv = 0;
      while ((mn->next != NULL)&&(equiv==0)){
        mn = mn->next;
        if (Check_Equivalence_bwn_Two_MATCHes(&newM,mn)==1) equiv = 1;
      }

      if (equiv==0){
        Add_MATCH_to_MATCHlist(&newM,Mlist);
        Nmatch += 1;   
      }

    }
 }
 Free_MATCH(&newM);
 Free_MATCH(&origM);
 printf("#newNmatch %d\n",Nmatch);  
 /*
 Keep_Only_Nonredundant_MATCH(Mlist); 
 */ 
 cal_select_dis_for_MATCHlist(Mlist,molA,molB);
 return(Nmatch);

} /* end of Make_Permutated_MATCHlist_from_Top_MATCH(Mlist) */



int Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation(Mlist,molA,molB)
  struct MATCH *Mlist;
  struct MOLECULE *molA, *molB;
{
  struct PERMUTATION *pA;
  struct MATCH newM,origM,*mn; 
  int i,Npair,Nmatch,equiv;
 
  /* show_all_the_atom_permutations(molA); */
  
  if (Mlist->next==NULL) return(0); 

  Npair = Mlist->next->Npair;
  Malloc_MATCH(&origM,Npair);
  Copy_MATCH(&origM,Mlist->next);
  Malloc_MATCH(&newM,Npair);
  Copy_MATCH(&newM,Mlist->next);

  Free_MATCHlist(Mlist);
  Mlist->next = NULL;
 
  Nmatch = 0; 
  pA = &(molA->permuhead);
  while (pA->next != NULL){
    pA = pA->next;
    for (i=0;i<Npair;++i){
      newM.anumA[i] = molA->num_frm_heavy[pA->index_heavy[molA->atoms[origM.anumA[i]].num_heavy]];
      newM.anumB[i] = origM.anumB[i];
    } 

    Bubble_Sort_anumA_anumB_array_by_anumA(&newM);
    mn = Mlist;
    equiv = 0;
    while ((mn->next != NULL)&&(equiv==0)){
      mn = mn->next;
      if (Check_Equivalence_bwn_Two_MATCHes(&newM,mn)==1) equiv = 1;
    }

    if (equiv==0){
      Add_MATCH_to_MATCHlist(&newM,Mlist);
      Nmatch += 1;   
    }
 }
 Free_MATCH(&newM);
 Free_MATCH(&origM);
 printf("#newNmatch %d\n",Nmatch);  
 /*
 Keep_Only_Nonredundant_MATCH(Mlist); 
 */ 
 Print_MATCHlist(Mlist,"Mlist(permutated)");
 cal_select_dis_for_MATCHlist(Mlist,molA,molB);
 return(Nmatch);

} /* end of Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation() */





int Make_NonRedundant_MATCHlist_by_Permutation(Mlist,molA,molB)
  struct MATCH *Mlist;
  struct MOLECULE *molA, *molB;
{
  struct MATCH *mn,*mo; 
  int Npair,Nsame,lenMlist,lenMlist_new;
  struct MATCH origMlist; 
  
 if (Mlist->next==NULL) return(0); 
  Npair = Mlist->next->Npair;
  lenMlist = Length_of_MATCHlist(Mlist);
  Malloc_MATCHlist(&origMlist,Npair,lenMlist); 
  Copy_MATCHlist(&origMlist,Mlist); 
  Free_MATCHlist(Mlist);
 
  Mlist->next = NULL;
  lenMlist_new = 0; 
  mn = &origMlist;
  while (mn->next != NULL){
    mn = mn->next;

    mo = Mlist;
    Nsame = 0;
    while ((mo->next != NULL)&&(Nsame==0)){
      mo = mo->next;
      if (Check_Equivalence_bwn_Two_MATCHes_by_Permutation(mn,mo,molA,molB)==1){ Nsame += 1;}
    }
    if (Nsame==0) {Add_MATCH_to_MATCHlist(mn,Mlist); lenMlist_new += 1; }
  }
 
 Free_MATCHlist(&origMlist);
 return(lenMlist_new);

} /* end of Make_NonRedundant_MATCHlist_by_Permutation() */







int Check_Equivalence_bwn_Two_MATCHes_by_Permutation(mP,mQ,molA,molB)
  struct MATCH *mP,*mQ;   /* two atomic correspondence to be compared */
  struct MOLECULE *molA, *molB;
{
  struct PERMUTATION *pA,*pB;
  struct MATCH newmP; 
  int i,Nmatch,nA,nB;
 
  Malloc_MATCH(&newmP,mP->Npair);
 
  Nmatch = 0; 
  nA = 0;

  pA = &(molA->permuhead);

  while (pA->next != NULL){
    pA = pA->next;
    nA += 1;
    nB = 0;
    pB = &(molB->permuhead);
    while (pB->next != NULL){
      pB = pB->next;
      nB += 1;
      for (i=0;i<mP->Npair;++i){
        newmP.anumA[i] = pA->index_heavy[molA->atoms[mP->anumA[i]].num_heavy];
        newmP.anumB[i] = pB->index_heavy[molB->atoms[mP->anumB[i]].num_heavy];
      } 
      Bubble_Sort_anumA_anumB_array_by_anumA(&newmP);
      if (Check_Equivalence_bwn_Two_MATCHes(&newmP,mQ)==1) { Free_MATCH(&newmP); return(1);}
    }
 }

 Free_MATCH(&newmP);
 return(0);
} /* end of Check_Equivalence_bwn_Two_MATCHes_by_Permutation() */ 




void Rerank_MATCHlist_by_transformed_RMSD(Mlist,molA,molB)
  struct MATCH *Mlist;
  struct MOLECULE *molA,*molB;
{
  struct MATCH *m;
  double g1[3],g2[3],Rmat[3][3];
  int Nmatch;

  g1[0] = g1[1] = g1[2] = 0.0;
  g2[0] = g2[1] = g2[2] = 0.0;
  Rmat[0][0] =  Rmat[1][1] =  Rmat[2][2] = 1.0; 

  Nmatch = 0;
  m = Mlist;
  while ((m->next != NULL) && (m->next->nodetype != 'E')){ 
    m = m->next;  Nmatch += 1;
 }
  printf("#Rerank_MATCHlist_by_transformed_RMSD(Nmatch:%d)\n",Nmatch);

  if (Nmatch>1){ 
    m = Mlist;
    while ((m->next != NULL) && (m->next->nodetype != 'E')){ 
       m = m->next;
       m->subfloat = m->select_dis;
       m->select_dis = Calculate_CRMS_MATCH_Quaternion(m,molA,molB,g1,g2,Rmat);
     }
    Merge_Sort_Double_Linked_List_MATCH(Mlist,'-');
  }

} /* end of Rerank_MATCHlist_by_transformed_RMSD() */




/*
main(argc,argv)
  int argc;
  char **argv;
{
  int N,index[100],i; 
  if (argc<2){
    printf("permu [N]\n");
    exit(1);
  } 
  N = atoi(argv[1]);

  for (i=0;i<N;++i){index[i] = i;}

  permutation(index,0,N);
} 
*/
