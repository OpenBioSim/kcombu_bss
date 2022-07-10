/*

<options.h> 

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


/** GLOBAL FUNCTIONS **/
void Set_Default_Global_PAR();
int Set_Global_PAR_from_OPTION_VALUE();
void Show_Option_Help_MCS();
char *Get_Date_String();
int Get_Time_in_Second();
double Get_Time_in_Second_by_Double();
void Set_END_DATE();



void Set_Default_Global_PAR(){
  int i;

  PAR.MAX_COMP_TIME_SEC = -1.0;
  PAR.M = PAR.N = PAR.O = PAR.P = 0;
  PAR.COMMENT[0] = '\0';
  PAR.Wneiatm    = 1.0;
  PAR.Wextcon    = 1.0;
  PAR.Wtopodis   = 0.0;
  PAR.Wrank      = 0.0;
  PAR.Wncompo    = 0.0;
  PAR.Wdistan = 0.0;
  PAR.Tole_distan = 2.0;
  PAR.atomtype_class = 'K';
  PAR.bondtype_class = 'X';
  PAR.ConnectGraphType = 'C';
  PAR.levelEC_Dextcon = '2';
  PAR.OutCalProcess = 'F';
  PAR.maxNcomponent = -1;
  PAR.AlgoType = 'B';
  PAR.BKmaxNpair = -1;
  PAR.maxDIFtopodis = -1;
  PAR.cal_topodis = 'F';
  PAR.GenEquivAtomPermu = 'F';
  PAR.subchar = ' ';
  PAR.max_separation = 10;
  PAR.maxD_Dtopodis = 4;
  PAR.NonRedSelectDis = 'T';
  PAR.Nkeep = 100;
  PAR.WeightScheme  = 'D';
  PAR.maxNatom_in_ring = 20;
  PAR.maxNatom_in_aromatic_ring = 8;
  PAR.includeHinMCS = 'F';


  sprintf(PAR.START_DATE,"%s",Get_Date_String());
  PAR.START_TIME_SEC = Get_Time_in_Second_by_Double();
  PAR.COMP_TIME_SEC          = 0.0;
  PAR.COMP_TIME_SEC_PROGRESS = 0.0;
 
  PAR.molnamekey[0] = '\0';
  PAR.SeedRand = 1;
  PAR.DAMRLassign = 'H';

  PAR.READ_MOLECULAR_PROPERTY = 'F';
  PAR.Nproperty = 0;
  for (i=0;i<MAX_PROPERTY;++i){
    PAR.PROP_NAME[i][0] = '\0';
    PAR.PROP_NUM[i] = 0;
    PAR.PROP_MIN[i] = 12345678.0;
    PAR.PROP_MAX[i] = 12345678.0;
    PAR.PROP_STR[i][0] =  '\0';
  } 

  PAR.LengthOfBond2D = 1.5;
  PAR.CollisionResolve2D = 'T';
  PAR.CHANGE_KEY_SPACE_TO_UNDERBAR = 'T';
  PAR.CHANGE_NOTE_TO_SUPPLIERS    = 'T';
  PAR.OutRankForSMILES = 'F';

  PAR.tole_dis_clash      =  1.0;
  PAR.tole_dis_atrct      =  2.0;

  PAR.alpha_sigmoid_clash = 10.0;
  PAR.Dmax_clash  = 6.0;
  PAR.minAngleFlipOfFragment = 10.0;

  PAR.AtomStereoFromBondStereo = 'P';
  /*
   The reason why PAR.BondStereo2StereoParity is set to 'P'rimary,
   is that BondStereros are often written in articles, and KEGG anotaters 
   may copy them exactly.
   */

 

  PAR.maxPermutations = 10000;

} /* end of Set_Default_Global_PAR() */





int Set_Global_PAR_from_OPTION_VALUE(ov)
  struct OPTION_VALUE *ov;
{

       if (strcmp(ov->opt,"nk")==0) {PAR.Nkeep = atoi(ov->val); return(1);}
  else if (strcmp(ov->opt,"w")==0)  {PAR.WeightScheme = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"H")==0)  {PAR.includeHinMCS = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"at")==0) {PAR.atomtype_class  = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"ec")==0) {PAR.levelEC_Dextcon = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"bo")==0) {PAR.bondtype_class  = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"con")==0){PAR.ConnectGraphType = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"wna")==0) {PAR.Wneiatm   = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"wec")==0) {PAR.Wextcon   = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"wtd")==0) {PAR.Wtopodis = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"wrk")==0) {PAR.Wrank    = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"wco")==0) {PAR.Wncompo = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"wds")==0) {PAR.Wdistan = atof(ov->val); return(1);} 
  else if (strcmp(ov->opt,"tds")==0) {PAR.Tole_distan = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"mtd")==0) {PAR.maxDIFtopodis = atoi(ov->val); return(1);}
  else if (strcmp(ov->opt,"xtd")==0) {PAR.maxD_Dtopodis = atoi(ov->val); return(1);}
  else if (strcmp(ov->opt,"nrs")==0) {PAR.NonRedSelectDis = ov->val[0]; return(1);} 
  else if (strcmp(ov->opt,"ctd")==0) {PAR.cal_topodis   = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"mco")==0) {PAR.maxNcomponent = atoi(ov->val); return(1);}
  else if (strcmp(ov->opt,"alg")==0) {PAR.AlgoType = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"per")==0) {PAR.GenEquivAtomPermu = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"xsec")==0){PAR.MAX_COMP_TIME_SEC = atof(ov->val); return(1);}
  else if (strcmp(ov->opt,"sbc")==0)  {PAR.subchar = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"sbi")==0)  {PAR.subint  = atoi(ov->val); return(1);}
  else if (strcmp(ov->opt,"bkmx")==0)  {PAR.BKmaxNpair = atoi(ov->val); return(1);}
  else if (strcmp(ov->opt,"nml")==0) {sprintf(PAR.molnamekey,"%s",ov->val); return(1);}
  else if (strcmp(ov->opt,"srd")==0) {PAR.SeedRand  = atoi(ov->val); return(1); }
  else if (strcmp(ov->opt,"prp")==0) {PAR.READ_MOLECULAR_PROPERTY  = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"damrl")==0){PAR.DAMRLassign = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"sp2ub")==0){PAR.CHANGE_KEY_SPACE_TO_UNDERBAR = ov->val[0]; return(1);} 
  else if (strcmp(ov->opt,"no2su")==0){PAR.CHANGE_NOTE_TO_SUPPLIERS = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"frmbdste")==0){PAR.AtomStereoFromBondStereo = ov->val[0]; return(1);}
  else if (strcmp(ov->opt,"mxper")==0){PAR.maxPermutations = atoi(ov->val); return(1);}

  return(0); 
} /* end of Set_Global_PAR_from_OPTION_VALUE() */


void Show_Option_Help_MCS()
{
    printf("<common options for finding MCS>\n");
    printf(" -alg : Algorithm. 'B'uild-up, 'X':eXact,'P'ca-3D-fit '3'D(raw_fit) 'F':from file (-iam option).[%c]\n",PAR.AlgoType);
    printf(" -con : Connectivity.'D'isconnected,'C'onnected, 'S'ubstructure(A<=B), 'I'somorphic (A=B),\n");
    printf("      :              'T'opo_constrained_disconnected  't'opo_constrained_connected)[%c]\n",PAR.ConnectGraphType);
    printf(" -mtd : max difference of topological distance(num of bonds in the shortest path).<0:don't care. [%d]\n",PAR.maxDIFtopodis);
    printf(" -at  : atomtype classification. 'X':don't care 'E'lement(C,N,O) ele+'R'ing (C,C@,N,N@)\n");
    printf("      :   ele+'B'ond_num (C1,C2,C3) 'T':ele+bond+ring (C2,C2@,C3,C3@)\n");
    printf("      :   'K':combu-recommend(C,C@,O,O@,O1,N,N@,N1,S,S@,P,X) 'D'amrl (D,A,M,R,L) 'k'cf\n");
    printf("      :   'F':Halogen-unified kcombu scheme. (F,Cl,Br,I)-->'F'.\n");
    printf("      :   'H':alogen=carbon kcombu scheme. (F,Cl,Br,I)-->'C'. [%c]\n",PAR.atomtype_class);
    printf(" -bo  :bondtype class. 'X':don't care,'C'are SDF/MOL2-defined bondtype\n");
    printf("      :                'A'rommatic_single_double.'R'otatable_or_not. [%c]\n",PAR.bondtype_class);
    printf(" -mco :maximum number of connected component. <=0:don't care. [%d]\n",PAR.maxNcomponent);
    printf(" -H   :include hydrogens in MCS calcuclation (T or F)[%c]\n",PAR.includeHinMCS); 
    printf(" -ctd :calculate topological distance (T or F) [%c]\n",PAR.cal_topodis);
    printf(" -srd  :seed for random number (>0) [%d]\n",PAR.SeedRand);
    printf(" -damrl:Assignment of DAMRL type.  '3':by 3D conf, 'H':by hydrogens, 'X':pdb->'3',other_format->'H' [%c]\n",PAR.DAMRLassign);
    printf("<common options for build-up MCS algorithm>\n");
    printf(" -nk  :Nkeep for build up [%d]\n",PAR.Nkeep);
    printf(" -w   :weighting scheme. 'D'efault weighting. 'S'pecified by user  [%c]\n",PAR.WeightScheme);
    printf(" -wna :weight_select_dis for neighoboring atom type     [%f]\n",PAR.Wneiatm);
    printf(" -wec :weight_select_dis for extended connectivity(EC)  [%f]\n",PAR.Wextcon);
    printf(" -wtd :weight_select_dis for topological distance       [%f]\n",PAR.Wtopodis);
    printf(" -wrk :weight_select_dis for rank by unique EC          [%f]\n",PAR.Wrank);
    printf(" -wco :weight_select_dis for num of connected component [%f]\n",PAR.Wncompo);
    printf(" -nrs :Nonredundancy check using selection distance score ('T' or 'F')[%c]\n",PAR.NonRedSelectDis);
    printf(" -ec  :level of EC for scoring ('0','1','2','3','4') [%c]\n",PAR.levelEC_Dextcon);
    printf(" -xtd :maximum topological distance for considering Dtopodis [%d]\n",PAR.maxD_Dtopodis);
    printf("<common options for exact MCS algorithm>\n");
    printf(" -bkmx: maximum number of atom pairs for Bron-Kerbosch [%d]\n",PAR.BKmaxNpair);
    printf(" -xsec: maximum computational time (in second) [%lf]\n",PAR.MAX_COMP_TIME_SEC);
    printf("<common options for input>\n");
    printf(" -nml : key for molecular name, such as 'NAMIKI_ID' or 'LIGANDBOX_ID' or 'chembl_id' [%s]\n",PAR.molnamekey);
    printf(" -prp : read molecular properties in SDF/MOL2 file. ('T' or 'F') [%c]\n",PAR.READ_MOLECULAR_PROPERTY);
    printf(" -sp2ub :change space to underbar for key string (T or F)[%c]\n",PAR.CHANGE_KEY_SPACE_TO_UNDERBAR);
    printf(" -no2su :change 'NOTE' into 'SUPPLIERS' (T or F)[%c]\n",PAR.CHANGE_NOTE_TO_SUPPLIERS);
    printf(" -frmbdste: Assign Atomic_stereo from Bond_stereo.  'F'alse,'C':just checking, 'P'rimary,'S'econdary,'B':only_bond [%c]\n",
         PAR.AtomStereoFromBondStereo);
    printf("<other options>\n"); 
    printf(" -per :Equivalent Atom Permutation. 'F':false,'N'on-redundant 'B'uild-up_filter\n");
    printf("      :  'A':make all matches by (molA X molB) permu. 'a': make all matches by molA permu (for molA=molB) [%c]\n",PAR.GenEquivAtomPermu);
    printf(" -mxper:maximum permutation number for '-per'. [%d]\n",PAR.maxPermutations);
}


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


int Get_Time_in_Second()
{
 time_t      now_t;
 int sec;
  now_t = time(NULL);
  sec = now_t;
  return(sec);
} /* end of Get_Time_in_Second() */


double Get_Time_in_Second_by_Double()
{
 struct timeval tv;
 gettimeofday(&tv,NULL);
 return(tv.tv_sec + (double)tv.tv_usec*1e-6);
} /* end of Get_Time_in_Second_by_Double() */




void Set_END_DATE(){
 double end_sec;
 sprintf(PAR.END_DATE,"%s",Get_Date_String());
 end_sec = Get_Time_in_Second_by_Double();
 PAR.COMP_TIME_SEC = end_sec - PAR.START_TIME_SEC;
}

