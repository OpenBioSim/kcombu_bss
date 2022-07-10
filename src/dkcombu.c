/*
 
 <dkcombu.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

  finger-print-based chemical comparison

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <dirent.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "ioLINE.h"
#include "ioSDF.h"
#include "ioLib.h"
#include "ioSMILES.h"
#include "match.h"
#include "io_match.h"
#include "buildup.h"
#include "clique.h"
#include "c_clique.h"
#include "PCAfit.h"
#include "MrgSrtMATCH.h"
#include "MrgSrtLIBMOL.h"
#include "qRMS.h"
#include "molprop.h"
#include "molring.h"
#include "molpermu.h"
#include "moltopodis.h"
#include "AtmPairDesc.h"
#include "DualAtmPair.h"
#include "options.h"
#include "ringblock.h"
#include "RingDesc.h"

struct PARAMETERS PAR;

/** FUNCTIONS (LOCAL) **/
static void write_search_result_for_multiple_query_search();
static float fusion_scores();
static void write_progress_message();

int main(argc,argv)
 int argc;
 char **argv;
{
  struct MOLECULE molQ, molL;
  int L,a,n,ok,prop_ok;
  char MODE,DetailOptHelp,CalMCS;
  int NmaxMCS;
  struct MATCH optMlist;
  char oschfile[MAX_FILENAME],osimfile[MAX_FILENAME],oprgfile[MAX_FILENAME],opropfile[MAX_FILENAME];
  char filetype_osim;
  char odeslib[MAX_FILENAME],ideslib[MAX_FILENAME];
  char omolsfile[MAX_FILENAME],omolsdir[MAX_FILENAME];  
  int  Nmax_outmol;

  char fname[MAX_FILENAME],fname0[MAX_FILENAME],buffA[MAX_FILENAME],buffB[MAX_FILENAME];
  
  /*** variables for library */ 
  char   ilistlib[MAX_FILENAME],idirlib[MAX_FILENAME],idirmlib[MAX_FILENAME],iheadlib[MAX_FILENAME];
  char   ltailstr[MAX_FILENAME]; 
  struct LINENODE HeadLibFile;
  struct LIBMOL  HeadLibMol;
  int    Nlibfile; 
  int    Nmollib;
  int    Nmol_sele_desc,Nmol_in_list;

 /** variables for query library */
  char ilistlibQ[MAX_FILENAME],idirlibQ[MAX_FILENAME],iheadlibQ[MAX_FILENAME],idirmlibQ[MAX_FILENAME];
  struct LIBMOL  HeadLibMolQ;
  int    NmollibQ;


  /*** variables for multiple queries */ 
  char idirquery[MAX_FILENAME];
  struct LIBMOL HeadQueryMol;
  int Nmolquery;

  char FusionScDescriptor, FusionScMCS;
 
  float  ThreTanimotoOneAtom,ThreTanimotoAtomPair,ThreTanimotoMCS;
 
  struct LIBMOL *an,*bn,*xn,*qn;
  struct FLOAT2DMAP Dmap;
  int    Nmcs,Nmol;
  double Nmolpair;
  char  OMatWithProperty,des_filemode, accept,AnnotationSmiles,*string; 
  int   DIVbunshi, DIVbunbo, Nsta, Nend;
  FILE *fpl,*fpo,*fpop,*fpL;
  DIR  *dir; 
  int  Ncalmatch;
  unsigned char *descriptorL,*descriptorA,*descriptorB;
  float sim_one,sim_ring,sim_pair,sim_mcs,sim_pair_norm,sim_oneringA, sim_oneringI;
  struct LINENODE *ln;
  long offset;
  struct LUNIT HeadLunit;
  int maxNheavyatom, minNheavyatom;
  char KeepLargestConGraph;
  double Rmat[3][3],gA[3],gB[3];

  struct OPTION_VALUE *ov;

  /****** SET INITIAL PARAMETERS ***********/
  Set_Default_Global_PAR();
  PAR.PROGRAM_TYPE = 'd'; 
  MODE = 'S';
  Initialize_MOLECULE(&molQ);
  Initialize_MOLECULE_string(&molQ);
  molQ.filetype = molL.filetype = '-'; 
  molQ.atmhet = molL.atmhet =  'B';
  molQ.chain = molL.chain =  '-';
  molQ.resi[0] = molL.resi[0] = '\0';
  molQ.SetStereo3D = molL.SetStereo3D = 'F'; 
  molQ.BondType = molL.BondType = 'B'; 

  sprintf(oschfile,"sch.out");
 
  ilistlib[0] = iheadlib[0] = idirlib[0] = idirmlib[0] = '\0';
  osimfile[0] =  oprgfile[0]  =  opropfile[0] = '\0';
  filetype_osim = 'N';

  odeslib[0] = ideslib[0] = omolsfile[0] = omolsdir[0] = '\0'; 
  idirquery[0] = '\0'; 
  ltailstr[0] = '\0'; 
  OMatWithProperty = 'F';
  DIVbunshi = 0; DIVbunbo = 1;

  ilistlibQ[0] = idirlibQ[0] = idirmlibQ[0] = iheadlibQ[0] = '\0';
  NmollibQ = 0;


  ThreTanimotoAtomPair =   0.1;
  ThreTanimotoOneAtom  =   0.1;
  ThreTanimotoMCS      =   0.1;

  PAR.simtype_descriptor  = 'T'; 
  CalMCS = 'T';
  NmaxMCS = 100; 
  Nmax_outmol = -1; 
  des_filemode = 'N';

  PAR.max_ring_size  = 30;
  PAR.max_block_size = 40;
  PAR.max_ring_descriptor  = PAR.max_ring_size + PAR.max_block_size  - 3;
  PAR.type_atompair = 'N';

  PAR.max_separation_cyclic  = 5;
  PAR.max_separation_acyclic = 5;
  PAR.cal_topodis = 'T';
  Nlibfile = 0;
  
  HeadLibFile.next = NULL;
  HeadLibFile.num = -1; 

  HeadLibMol.next = NULL;
  HeadLibMol.num = -1; 
  

  HeadLibMolQ.next = NULL;
  HeadLibMolQ.num  = -1;
  
  Nmollib = 0;
 
  FusionScDescriptor = 'X';
  FusionScMCS        = 'X';

  maxNheavyatom = minNheavyatom = -1;

  AnnotationSmiles = 'F';

  KeepLargestConGraph = 'F';

  ln =NULL;
  fpo = fpl = fpop = NULL;
  n = 0;
  sim_one  = 0;
  accept = 0;
  Nmolquery = 0;

  DetailOptHelp = '-';
  if((argc>1) && ((strcmp(argv[1],"-h")==0) || (strcmp(argv[1],"-H")==0)) )        { DetailOptHelp = 'D';}
  if((argc>1) && ((strncmp(argv[1],"-hm",3)==0) || (strncmp(argv[1],"-HM",3)==0)) ){ DetailOptHelp = 'M';}

  if ((argc<2)&&(DetailOptHelp=='-')){
    printf("dkcombu <options>\n");
    /* printf("sizeof(LIBMOL) %d\n",sizeof(struct LIBMOL)); */
    printf("  %s\n",KCOMBU_STANDS_FOR);
    printf("  'dkcombu' is a hybrid search using both descriptor and MCS comparisons.\n");
    printf("  %s LastModified:%s\n",CODED_BY,LAST_MOD_DATE);
    printf(">'D':for making descriptor file from the molecule library\n");  
    printf(" $dkcombu -M D -ifl  [files of lib multi mol] -odes [output_descriptor_file]\n");
    printf(" $dkcombu -M D -idl  [dir. of  lib mol]       -odes [output_descriptor_file]\n");
    printf(" $dkcombu -M D -idml [dir. of  lib multi mol] -odes [output_descriptor_file]\n");
    printf(">'S': Query vs Library chemical structure search\n");
    printf(" $dkcombu -M S -Q [query_mol] -ides [input_lib_descriptor_file] -osc [result]\n");
    printf(">'A': All-vs-all chemical structure comparison\n");
    printf(" $dkcombu -M A -ides [input_lib_descriptor_file] -osl [output_similarity_listfile]\n");
    printf(" $dkcombu -M A -ides [input_lib_descriptor_file] -osm [output_similarity_matrix]\n");
    printf(">'G': Comparison bwn Query-lib vs Library\n");
    printf(" $dkcombu -M G -idlq [dir. of query-lib mol] -ides [input_lib_descriptor_file] -osm [result]\n");
    printf(">for showing options in detail\n");
    printf(" $dkcombu -h\n");
    printf(">for showing options for MCS\n");
    printf(" $dkcombu -hmcs\n");
   }
  if (DetailOptHelp=='D'){
    printf("dkcombu <options>\n");
    printf("  'dkcombu' is a hybrid search using both descriptor and MCS comparisons.\n");
    printf("  %s\n",KCOMBU_STANDS_FOR); 
    printf("   coded by T.Kawabata. LastModified:%s\n",LAST_MOD_DATE);
    printf("<options>\n");
    printf(" -M   : MODE 'D':make descriptor file 'S':library_search 'M'ultiple_query_vs_library search\n");
    printf("      :      'A'll_vs_all 'a'll_vs_all_for_huge_lib,'G':comparison bwn Query-lib vs Library [%c]\n",MODE);
    printf("<options for MODE 'D'>\n"); 
    printf(" -odes : output file for atompair descriptors   [%s]\n",odeslib);
    printf(" -ides : input  file for atompair descriptors   [%s]\n",ideslib);
    printf(" -ap   : atom pair type use 'R'ingblock 'D'ual 'N'ormal [%c]\n",PAR.type_atompair);
    printf(" -sep  : max separation for atom pair descriptor [%d]\n",PAR.max_separation);  
    printf(" -sepa : max separation_acyclic for atom pair descriptor, only for '-ap D'.[%d]\n",PAR.max_separation_cyclic);  
    printf(" -sepc : max separation_cyclic for atom pair descriptor,  only for '-ap D'.[%d]\n",PAR.max_separation_acyclic);  
    printf(" -fdes : file mode for output/input atompair descriptor file  'T'ext 'R'educed_text 'B'inary [%c]\n",des_filemode);
    printf("<options for input library>\n"); 
    printf(" -ifl   : input files for all the library molecules (file1:file2:..) [%s]\n","");
    printf(" -idl   : input directory for library molecular file          (one-mol-one-file) [%s]\n",idirlib);
    printf(" -idml  : input directory for multi-mol-library molecular file(multi-mol-one-file) [%s]\n",idirmlib);
    printf(" -ihl   : required head string of filename under the directory specied by '-idl' option. [%s]\n",iheadlib);
    printf(" -ill   : input list of library molecules[%s]\n",ilistlib);
    printf(" -ides  : input descriptor file for library molecules[%s]\n",ideslib);
    printf(" -ltail : tail string of library name[%s]\n",ltailstr);
    printf("<options for MODE 'S' (for similarity thredhold)>\n");
    printf(" -tap : Tanimoto threshold value for arom pair descriptor [%f]\n",ThreTanimotoAtomPair);
    printf(" -tao : Tanimoto threshold value for one atom  descriptor [%f]\n",ThreTanimotoOneAtom);
    printf(" -tam : Tanimoto threshold value for MCS                  [%f]\n",ThreTanimotoMCS);
    printf(" -sde : similarity type for descriptor search 'T'animoto, 'C'ahart [%c]\n",PAR.simtype_descriptor); 
    printf(" -nmcs: maximum number of compounds for MCS calculation. No limit if -nmcs -1. [%d]\n",NmaxMCS);    
    printf(" -mcs : Do MCS search after atompair descriptor search (T or F)[%c]\n",CalMCS); 
    printf(" -mxnh: maximum value of Nheavyatom (<0:don't care) [%d]\n",maxNheavyatom); 
    printf(" -minh: minimum value of Nheavyatom (<0:don't care) [%d]\n",minNheavyatom); 
    printf("<options for MODE 'S'>\n"); 
    printf(" -Q   : query molecule Q (*.sdf|*.mol2|*.pdb|*.kcf|*.cif)[%s]\n",molQ.filename);
    printf(" -osc : output search result file   [%s]\n",oschfile);
    printf(" -oprg: output search progress message   [%s]\n",oprgfile);
    printf(" -oms : output file for similar molecules     [%s]\n",omolsfile); 
    printf(" -omsd: output dir  for files for similar molecules [%s]\n",omolsdir); 
    printf(" -nms : max number of output similar molecules. -1:don't care. [%d]\n",Nmax_outmol); 
    printf(" -kplcn : keep only the largest connected subgraph for query Q (T or F)[%c]\n",KeepLargestConGraph);
    printf("<options for MODE 'A'>\n"); 
    printf(" -osm : output file for similarities  [%s]\n",osimfile);
    printf(" -fsm : file formats for osm. 'L':simple pairwise list.'N':list with molecular names\n");
    printf("      :                       'D':distance matrix.  'S'imilarity matrix [%c]\n",filetype_osim);
    printf(" -wpr : ourput dis/sim matrix with property names ('T' or 'F')[%c]\n",OMatWithProperty);
    printf("<options for MODE 'a'>\n"); 
    printf(" -div : Job division. (bunshi)/(bunbo)[%d/%d]\n",DIVbunshi,DIVbunbo);
    printf("<options for '-M G'>\n");
    printf(" -iflq  : input files for the query library molecules (file1:file2:..) [%s]\n","");
    printf(" -idlq  : input directory for query library molecules   [%s]\n",idirlibQ);
    printf(" -idmlq : input directory for multi-mol query library molecular file(multi-mol-one-file) [%s]\n",idirmlibQ);
    printf(" -illq  : input list of query library molecules[%s]\n",ilistlibQ);
    printf("<options for MODE 'M'>\n"); 
    printf(" -idq  : input directory for library molecular file    (one-mol-one-file) [%s]\n",idirquery);
    printf(" -fud  : fusion score scheme for descriptor 'A'verage 'X':maximum [%c]\n",FusionScDescriptor);
    printf(" -fum  : fusion score scheme for descriptor 'A'verage 'X':maximum [%c]\n",FusionScMCS);
    printf("<other options>\n"); 
    printf(" -o   : output calculation process 'T' or 'F' [%c]\n",PAR.OutCalProcess);
    printf(" -fQ,-fL : file formats.'P'db,'S'df,'K'cf, '2':MOL2, 'C':cif [%c%c]\n",molQ.filetype,molL.filetype);    
    printf(" -aQ,-aL : AtomHetero type. 'A'tom 'H'etatm 'B'oth[%c%c]\n",molQ.atmhet,molL.atmhet);    
    printf(" -bQ,-bL : Bond Type. 'B':consider bonds, otherwise, not consider bonds[%c%c]\n",molQ.BondType,molL.BondType);
    printf("<options for properties for MODE 'D' and 'S'>\n");
    printf(" -pn0   : name of 0-th property [%s]\n",PAR.PROP_NAME[0]);
    printf(" -pi0   : minimum value of 0-th property [%f]\n",PAR.PROP_MIN[0]);
    printf(" -px0   : maximum value of 0-th property [%f]\n",PAR.PROP_MAX[0]);
    printf(" -oprp  : output filename for properties [%s]\n",opropfile);
    printf(" -anosmi: add SMILES string as property annotation (T or F)[%c]\n",AnnotationSmiles);
    printf(" -ste3Q : Setup Stereo Parity from 3D for mol Q (T or F)[%c]\n",molQ.SetStereo3D);
    printf(" -ste3L : Setup Stereo Parity from 3D for mol L (T or F)[%c]\n",molL.SetStereo3D);
  }

  if (DetailOptHelp=='M'){
    Show_Option_Help_MCS();
  }

  if ((argc<2)||(DetailOptHelp!='-')) exit(1);
 

 /****** READ ARGUMENTS ***********/

 sprintf(PAR.START_DATE,"%s",Get_Date_String());
 PAR.START_TIME_SEC = Get_Time_in_Second_by_Double(); 
 PAR.OptValHead.next = NULL;
 Read_Options_From_Arguments(argc,argv,PAR.COMMAND,&(PAR.OptValHead));
 ov = &(PAR.OptValHead);
 while (ov->next != NULL){
   ov = ov->next;
   if (Set_Global_PAR_from_OPTION_VALUE(ov)==0){
         if (strcmp(ov->opt,"M")==0) { MODE = ov->val[0];}
    else if (strcmp(ov->opt,"Q")==0) { sprintf(molQ.filename,"%s",ov->val);}
    else if (strcmp(ov->opt,"ap")==0)  { PAR.type_atompair = ov->val[0];}
    else if (strcmp(ov->opt,"ifl")==0) { 
      Nlibfile = Make_LibFileList_from_Colon_Separated_String(&HeadLibFile,ov->val);
      Show_LINENODEs(&HeadLibFile);
    }
    else if (strcmp(ov->opt,"idl")==0) { sprintf(idirlib,"%s",ov->val);}
    else if (strcmp(ov->opt,"idq")==0) { sprintf(idirquery,"%s",ov->val);}
    else if (strcmp(ov->opt,"idml")==0) { sprintf(idirmlib,"%s",ov->val);}
    else if (strcmp(ov->opt,"ihl")==0)  { sprintf(iheadlib,"%s",ov->val);}
    else if (strcmp(ov->opt,"ill")==0)  {sprintf(ilistlib,"%s",ov->val);}
    else if (strcmp(ov->opt,"osm")==0)  {sprintf(osimfile,"%s",ov->val);}
    else if (strcmp(ov->opt,"oprg")==0)  { sprintf(oprgfile,"%s",ov->val);}
    else if (strcmp(ov->opt,"wpr")==0)  { OMatWithProperty = ov->val[0];}
    else if (strcmp(ov->opt,"osc")==0)   { sprintf(oschfile,"%s",ov->val);}
    else if (strcmp(ov->opt,"tao")==0)  { ThreTanimotoOneAtom  = atof(ov->val);}
    else if (strcmp(ov->opt,"tap")==0)  { ThreTanimotoAtomPair = atof(ov->val);}
    else if (strcmp(ov->opt,"tam")==0)  { ThreTanimotoMCS      = atof(ov->val);}
    else if (strcmp(ov->opt,"odes")==0) { sprintf(odeslib,"%s",ov->val);}
    else if (strcmp(ov->opt,"ides")==0) { sprintf(ideslib,"%s",ov->val);}
    else if (strcmp(ov->opt,"fdes")==0) { des_filemode = ov->val[0]; }
    else if (strcmp(ov->opt,"mcs")==0) { CalMCS = ov->val[0]; }
    else if (strcmp(ov->opt,"sep")==0)  { PAR.max_separation = atoi(ov->val); }
    else if (strcmp(ov->opt,"sepc")==0) { PAR.max_separation_cyclic = atoi(ov->val); }
    else if (strcmp(ov->opt,"anosmi")==0) { AnnotationSmiles = ov->val[0]; }
    else if (strcmp(ov->opt,"sepa")==0) { PAR.max_separation_acyclic = atoi(ov->val); }
    else if (strcmp(ov->opt,"nmcs")==0) { NmaxMCS = atoi(ov->val); }
    else if (strcmp(ov->opt,"oms")==0) {  sprintf(omolsfile,"%s",ov->val); }
    else if (strcmp(ov->opt,"omsd")==0) {  sprintf(omolsdir,"%s",ov->val); }
    else if (strcmp(ov->opt,"nms")==0) {   Nmax_outmol = atoi(ov->val); }
    else if (strcmp(ov->opt,"ltail")==0) { sprintf(ltailstr,"%s",ov->val); }
    else if (strcmp(ov->opt,"fQ")==0)  { molQ.filetype = MoleculeFileType(ov->val);}
    else if (strcmp(ov->opt,"fL")==0)  { molL.filetype = MoleculeFileType(ov->val);}
    else if (strcmp(ov->opt,"aQ")==0)  { molQ.atmhet = ov->val[0];}
    else if (strcmp(ov->opt,"aL")==0)  { molL.atmhet = ov->val[0];}
    else if (strcmp(ov->opt,"bQ")==0)  { molQ.BondType = ov->val[0];}
    else if (strcmp(ov->opt,"bL")==0)  { molL.BondType = ov->val[0];}
    else if (strcmp(ov->opt,"ste3Q")==0)  { molQ.SetStereo3D = ov->val[0];}
    else if (strcmp(ov->opt,"ste3L")==0)  { molL.SetStereo3D = ov->val[0];}
    else if (strcmp(ov->opt,"kplcn")==0)  { KeepLargestConGraph = ov->val[0];}
    else if (strcmp(ov->opt,"div")==0) { 
       a = First_Appear_Pos_in_String(ov->val,'/');
       L = strlen(ov->val);
       if ((L>1)&&(a>0)){ 
         Get_Part_Of_Line(buffA,ov->val,a+1,L-1); DIVbunbo  = atoi(buffA);
         Get_Part_Of_Line(buffA,ov->val,0,a-1);   DIVbunshi = atoi(buffA);
       }
    }
    else if (strcmp(ov->opt,"fud")==0) {  FusionScDescriptor = ov->val[0]; }
    else if (strcmp(ov->opt,"fum")==0) {  FusionScMCS        = ov->val[0]; }
    else if (strcmp(ov->opt,"idlq")==0) { sprintf(idirlibQ,"%s",ov->val);}
    else if (strcmp(ov->opt,"mxnh")==0) { maxNheavyatom = atoi(ov->val);}
    else if (strcmp(ov->opt,"minh")==0) { minNheavyatom = atoi(ov->val);}
    else if (strcmp(ov->opt,"oprp")==0) { sprintf(opropfile,"%s",ov->val);}
    else if (strcmp(ov->opt,"oprp")==0) { sprintf(opropfile,"%s",ov->val);}
    else if (((strlen(ov->opt)==3)||(strlen(ov->opt)==4))&&(ov->opt[0]=='p')&&(ov->opt[1]=='n')) {
       Get_Part_Of_Line(buffA,ov->opt,3,strlen(ov->opt)-1); 
       n = atoi(buffA);
       if (n<MAX_PROPERTY){
         sprintf(PAR.PROP_NAME[n],"%s",ov->val);
         if ((n+1)>PAR.Nproperty) PAR.Nproperty = n + 1;
         PAR.READ_MOLECULAR_PROPERTY = 'T';
       } 
      else {printf("#WARNING:-pn%d is ignored. %d is over MAX_PROPERTY %d.\n",n,n,MAX_PROPERTY);}
    }
    else if (((strlen(ov->opt)==3)||(strlen(ov->opt)==4))&&(ov->opt[0]=='p')&&(ov->opt[1]=='i')) {
       Get_Part_Of_Line(buffA,ov->opt,3,strlen(ov->opt)-1); 
       n = atoi(buffA);
       if (n<MAX_PROPERTY){
         PAR.PROP_MIN[n] = atof(ov->val);
         if ((n+1)>PAR.Nproperty) PAR.Nproperty = n + 1;
         PAR.READ_MOLECULAR_PROPERTY = 'T';
       }
      else {printf("#WARNING:-pi%d is ignored. %d is over MAX_PROPERTY %d.\n",n,n,MAX_PROPERTY);}
    }
    else if (((strlen(ov->opt)==3)||(strlen(ov->opt)==4))&&(ov->opt[0]=='p')&&(ov->opt[1]=='x')) {
       Get_Part_Of_Line(buffA,ov->opt,3,strlen(ov->opt)-1); 
       n = atoi(buffA);
       if (n<MAX_PROPERTY){
         PAR.PROP_MAX[n] = atof(ov->val);
         if ((n+1)>PAR.Nproperty) PAR.Nproperty = n + 1;
         PAR.READ_MOLECULAR_PROPERTY = 'T';
       }
      else {printf("#WARNING:-px%d is ignored. %d is over MAX_PROPERTY %d.\n",n,n,MAX_PROPERTY);}
    }
    else if (((strlen(ov->opt)==3)||(strlen(ov->opt)==4))&&(ov->opt[0]=='p')&&(ov->opt[1]=='s')) {
       Get_Part_Of_Line(buffA,ov->opt,3,strlen(ov->opt)-1); 
       n = atoi(buffA);
       if (n<MAX_PROPERTY){
         sprintf(PAR.PROP_STR[n],"%s",ov->val);
         if ((n+1)>PAR.Nproperty) PAR.Nproperty = n + 1;
         PAR.READ_MOLECULAR_PROPERTY = 'T';
       }
      else {printf("#WARNING:-ps%d is ignored. %d is over MAX_PROPERTY %d.\n",n,n,MAX_PROPERTY);}
    }
    else { printf("#ERROR:Can't understand option '%s'\n",ov->opt); exit(1);}

   }
 }


 
 if ((PAR.ConnectGraphType=='T') || (PAR.ConnectGraphType=='t')){
   if (PAR.maxDIFtopodis <0) PAR.maxDIFtopodis = 0;
 }


/* Default Weighting */
 if ((PAR.WeightScheme=='D') && ((PAR.AlgoType=='B')||(PAR.AlgoType=='b'))){ 
     PAR.Wneiatm  = 1.0; 
     PAR.Wextcon  = 1.0; 
     PAR.Wdistan  = 0.0; 
     PAR.Wncompo  = 0.0;
     PAR.Wtopodis = 0.0; 
     if ((PAR.ConnectGraphType=='T')||(PAR.ConnectGraphType=='t')||(PAR.ConnectGraphType=='D')){ PAR.Wtopodis = 1.0;} 
 }

  Set_max_atomtype_by_atomtype_class();

  if (PAR.type_atompair == 'R'){ 
    PAR.max_atompair_descriptor = 
       (PAR.max_atomtype  * (PAR.max_atomtype + 1))/2* PAR.max_separation + 
       (PAR.ring_atomtype * (PAR.ring_atomtype + 1))/2* PAR.max_separation;
  }
  else if (PAR.type_atompair=='D'){
    PAR.max_atompair_descriptor   = (PAR.max_atomtype * (PAR.max_atomtype + 1))/2 * 
                                    (PAR.max_separation_cyclic + 1) * (PAR.max_separation_acyclic + 1);
  }
 else if (PAR.type_atompair=='N'){
    PAR.max_atompair_descriptor   = (PAR.max_atomtype * (PAR.max_atomtype + 1))/2* PAR.max_separation;
  }


  if (CalMCS!='T') ThreTanimotoMCS = 0.0;


printf("#COMMAND '%s'\n",PAR.COMMAND);


 /*******************************************************/ 
 /*** MODE 'S': Library Search for the query molecule ***/
 /*******************************************************/ 
 if (MODE=='S'){ 
  write_progress_message(oprgfile,"START_CALCULATION");

  if (PAR.ConnectGraphType=='S'){
   ThreTanimotoMCS = ThreTanimotoOneAtom = ThreTanimotoAtomPair = 0.0;  
  }


  if (ideslib[0]!='\0') {
    HeadLibFile.next = NULL;
    HeadLibFile.num = -1; 
    Nmollib = Read_Headers_of_AtomPair_Descriptor(ideslib,idirlib,&des_filemode,&HeadLibFile,&(molL.filetype));
    /* Show_LINENODEs(&HeadLibFile); */
    Nlibfile = Number_of_LINENODEs(&HeadLibFile);
  }


  descriptorL = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor); 

  /** [1] Read molecule Q **/
  Read_MOLECULE(molQ.filename,&molQ,molQ.filetype,molQ.atmhet,molQ.BondType);
  printf("#molQ.name '%s'\n",molQ.name);
  printf("#mol->topodismap.malloced %d\n",molQ.topodismap.malloced);
  if (KeepLargestConGraph == 'T'){
    /* added by T.Kawabata (2018/10/19) */
    Set_Connected_Structure(&molQ);
    Keep_Only_Largest_Connected_Structure(&molQ);
  }
 
  
  if (PAR.type_atompair=='D'){
    Set_DualAtomPair_Descriptor(&molQ);
  }
  else if (PAR.type_atompair=='N'){
    Set_AtomPair_Descriptor(&molQ);
  }

/*
  for (a=0;a<molQ.Natom;++a){
    printf("%d %s %s\n",molQ.atoms[a].num_in_file,molQ.atoms[a].element,molQ.atoms[a].atomtype);
  }
  exit(1);
 */
  /* Set_Ring_Descriptor(&molQ); */

  printf(">%s Natom %d Nbond %d\n",molQ.name, molQ.Natom, molQ.Nbond); 

  /** [2] Read Library List One by One and compare with molQ by atompair descriptor **/

  printf("## DO ATOMPAIR DESCRIPTOR SEARCH ##\n"); 
  Nmol = 0;
  if (ideslib[0]!='\0'){ fpl = fopen(ideslib,"r");}
  if (fpl==NULL) {printf("#ERROR:Can't open description lib file '%s'\n",ideslib); exit(1);}
  Nmol_sele_desc = 0;
  Nmol_in_list   = 0;

  HeadLibMol.num  = 0;
  HeadLibMol.next = NULL;

  do{
    xn = Read_AtomPair_Descriptor_into_LIBMOL_One_by_One(fpl,descriptorL,des_filemode,'A'); 
    accept = 1; 

    if(xn!=NULL){
    
      /*  
        printf(">%d '%s' %d\n",xn->num,xn->name,xn->offset_libfile);
      */ 
      if ((minNheavyatom>=0) && (xn->Nheavyatom < minNheavyatom)){ accept = 0;}
      if ((maxNheavyatom>=0) && (xn->Nheavyatom > maxNheavyatom)){ accept = 0;}

      if (PAR.ConnectGraphType=='I'){
        if ((Isomorphic_condition_oneatom_descriptors(molQ.oneatom_descriptor,xn->oneatom_descriptor)==0) || 
            (isomorphic_condition_atompair_descriptors(molQ.atompair_descriptor,descriptorL)==0)){accept = 0;}
       }

      if (PAR.ConnectGraphType=='S'){
         if ((Substructure_condition_oneatom_descriptors(molQ.oneatom_descriptor,xn->oneatom_descriptor)==0) || 
            (substructure_condition_atompair_descriptors(molQ.atompair_descriptor,descriptorL)==0)){accept = 0;}
       }

      if (accept==1){
        xn->tanimoto_oneatom = Tanimoto_bwn_oneatom_descriptors(molQ.oneatom_descriptor,xn->oneatom_descriptor);
        xn->tanimoto_atompair = similarity_bwn_atompair_descriptors_oneatom_normalized(molQ.atompair_descriptor, descriptorL);
        xn->score_for_sort = xn->tanimoto_atompair; 

        if (((ThreTanimotoOneAtom >=0.0)&&(xn->tanimoto_oneatom  < ThreTanimotoOneAtom)) || 
            ((ThreTanimotoAtomPair>=0.0)&&(xn->tanimoto_atompair < ThreTanimotoAtomPair)) ){accept = 0;}
        if (accept == 1) {Nmol_sele_desc += 1;}
      }
   
 
      if ((CalMCS=='T')&&(NmaxMCS>=1) && (Nmol_in_list >= NmaxMCS)){
        if ((HeadLibMol.next != NULL) && (xn->score_for_sort < HeadLibMol.next->score_for_sort)){
          accept = 0;
        }
      } 

      if (accept == 1) {
        if ((CalMCS=='T')&& (NmaxMCS>=1) && (Nmol_in_list >= NmaxMCS)){
          delete_the_first_LIBMOL_node(&HeadLibMol);
        } 
        insert_LIBMOL_into_score_decreasing_proper_position(xn,&HeadLibMol);
        Nmol_in_list   += 1;
      }

     /* Show_LIST_of_LIBMOL(&HeadLibMol);  */

      if (accept==0){
        /* printf("#accept=0 %s %s\n",xn->name,xn->molform); */
        free(xn->name);
        free(xn->molform);
        free(xn);
      }

       Nmol += 1; 
       if ((Nmol%5000)==0){ 
          sprintf(buffA,"Step I:descriptor search: %5.1f %% %d/%d",100.0*(float)Nmol/Nmollib,Nmol,Nmollib);
          write_progress_message(oprgfile,buffA);
       }

    }


  } while (xn != NULL); 

  printf("#NMOL_SELECTION_BY_ONEATOM_AND_ATOM_PAIR_DESCRIPTOR  %d\n",Nmol_sele_desc);
 
  reverse_LIBMOL_list(&HeadLibMol);
 /* 
  Merge_Sort_Double_Linked_List_LIBMOL(&HeadLibMol,'D');  
 */

  if (CalMCS != 'T'){
    Set_END_DATE();
    printf("#COMMAND   %s\n",PAR.COMMAND);
    printf("#COMP_TIME  %f seconds\n",PAR.COMP_TIME_SEC);
    printf("## [%d/%d] comounds are over tanimoto_oneatom %f atompair %f##\n",Nmol_sele_desc,Nmollib,ThreTanimotoOneAtom,ThreTanimotoAtomPair); 
    if (oschfile[0] != '\0'){
      Write_Search_Result(oschfile,&molQ,idirlib,ideslib,&HeadLibMol,Nmollib,Nmol_sele_desc,
            ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,&HeadLibFile);  
     }
    write_progress_message(oprgfile,"FINISH_CALCULATION");
    exit(1);
  }

  Set_END_DATE();
  PAR.COMP_TIME_SEC_PROGRESS = PAR.COMP_TIME_SEC;
  printf("#COMP_TIME_FOR_DESCRIPTOR_SEARCH  %f seconds\n",PAR.COMP_TIME_SEC_PROGRESS);

  /** [3] Calculate MCS-based similarity  **/
  sprintf(buffA,"Step I: descriptor search: %5.1f %% %d/%d",100.0*Nmol/Nmollib,Nmol,Nmollib);
  write_progress_message(oprgfile,buffA);

  Nmcs = Nmol_in_list;
  if ((NmaxMCS>0)&&(Nmcs>NmaxMCS)) {Nmcs = NmaxMCS;}

  printf("## DO MCS COMPARISON for [%d/%d] over tanimoto_atompair %f (NmaxMCS %d)##\n",Nmcs,Nmollib,ThreTanimotoAtomPair,NmaxMCS); 


  Nmol = 0; 
  prop_ok = 1;
  xn = &HeadLibMol;
  fpL = NULL;
  fname[0] = fname0[0] = '\0';
  Ncalmatch  =  0;
  while (xn->next !=NULL){
    xn = xn->next;
    xn->score_for_sort = 0.0;
    accept = 0;
    if ((NmaxMCS < 0) || (Ncalmatch < NmaxMCS)){

      sprintf(molL.name,"%s",xn->name);
      sprintf(molL.filename,"%s/%s",idirlib,xn->name);
 
      if (Nlibfile==0){ 
        ok = Read_MOLECULE(molL.filename,&molL,molL.filetype,molL.atmhet,molL.BondType);
        printf("#molL filename '%s' Natom %d\n",molL.filename, molL.Natom);
      }
      else {
       String_from_LINENODE_num(fname,&HeadLibFile,xn->num_file);
       if (strcmp(fname,fname0)!=0){
         if (fpL!= NULL) fclose(fpL);
         fpL = fopen(fname,"r"); 
       }
       fseek(fpL,xn->offset_libfile,SEEK_SET);
       ok = Read_One_Molecule_from_one_Library_file(fpL,&molL,molL.filetype,molL.atmhet,'L');
       if ((xn->name[0]=='\0')&&(molL.name[0]!='\0')) sprintf(xn->name,"%s",molL.name); 
       sprintf(fname0,"%s",fname);
     }


      if (ok>0){ 
        /* printf(">DoMCS '%s' score_for_sort %f\n",xn->name,xn->score_for_sort); */
        xn->Nheavyatom = molL.Nheavyatom; 
        xn->molform = (char *)malloc(sizeof(char)*(strlen(molL.molform)+1));
        sprintf(xn->molform,"%s",molL.molform);
        if (PAR.Nproperty>0)
          prop_ok = match_properties_and_patterns(&(molL.HeadProperty),PAR.Nproperty, PAR.PROP_NAME,PAR.PROP_MIN,PAR.PROP_MAX,PAR.PROP_STR);

        if ((molL.Nheavyatom>0)&&(prop_ok==1)){
           xn->doneMCS = 1;
           Match_Two_Molecules(&molQ,&molL,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
           Ncalmatch += 1;
           if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E') && (optMlist.next->Npair>0)){

              xn->RMSD = Calculate_CRMS_MATCH_Quaternion(optMlist.next,&molQ,&molL,gA,gB,Rmat,"BonA");
              make_rot_axis_angle_from_rot_matrix(xn->Raxisang,Rmat);
              make_tvec_from_rmat_gorig_gtra(xn->Tvec, Rmat, gB, gA);
 /*
              printf("#RMSD %f Raxisang %lf %lf %lf %lf Tvec %lf %lf %lf\n",xn->RMSD,
xn->Raxisang[0], xn->Raxisang[1], xn->Raxisang[2], xn->Raxisang[3], xn->Tvec[0], xn->Tvec[1], xn->Tvec[2]);
  */
              Copy_MATCH_Info_to_LIBMOL(xn,optMlist.next,&molQ,&molL);
              xn->score_for_sort = xn->tanimoto_mcs; 
              accept = 1;
           }
/*
          printf("[%d] %s Natom %d Npair %d tanimoto_atompair %f tanimoto_mcs %f select_dis %f\n",
            xn->num,xn->name,molQ.Nheavyatom,xn->Npair,xn->tanimoto_atompair,xn->tanimoto_mcs,xn->select_dis);
  */
          Free_MOLECULE(&molL);
          Free_MATCHlist(&optMlist);
        }
     }

     Nmol += 1;
     if ((Nmol%50)==0){ 
       sprintf(buffA,"Step II: MCS search: %5.1f %% %d/%d",100.0*(float)Nmol/Nmcs,Nmol,Nmcs);
       write_progress_message(oprgfile,buffA);
     }
   }

     if (accept==0){ 
       /* printf("REMOVE %s\n",xn->name); */
       if (xn->prev != NULL) xn->prev->next = xn->next;
       if (xn->next != NULL) xn->next->prev = xn->prev;
     }

  } /* while (xn) */
 
  Merge_Sort_Double_Linked_List_LIBMOL(&HeadLibMol,'D');
 
  if (oschfile[0] != '\0'){
   Write_Search_Result(oschfile,&molQ,idirlib,ideslib,&HeadLibMol,Nmollib,Nmol_sele_desc,ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,&HeadLibFile);  
  }

  if (omolsfile[0] != '\0'){
   Write_Similar_Molecules_To_One_File(omolsfile,"",idirlib,&HeadLibMol,
         ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,molL.filetype,Nlibfile,&HeadLibFile,Nmax_outmol);  
  }

  if (omolsdir[0] != '\0'){
   Write_Similar_Molecules_Under_Dir(omolsdir,"",idirlib,&HeadLibMol,
     ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,molL.filetype,Nlibfile,&HeadLibFile,Nmax_outmol);  
  }


  write_progress_message(oprgfile,"FINISH_CALCULATION");
}


 /**************************************************************************/
 /*** MODE 'D': make descriptor file by scanning the molecular library   ***/
 /**************************************************************************/
 if (MODE == 'D'){ 
   fpl = NULL;
   /**  Read Library List  **/
   printf("## Read Library List ##\n"); 
   if (Nlibfile==0){
          if (ilistlib[0]!='\0'){ Nmollib  =  Read_List_of_LIBMOLs(ilistlib,&HeadLibMol);}
     else if (idirlib[0]!='\0'){  Nmollib  =  Get_List_of_LIBMOLs_Under_Dir(idirlib,&HeadLibMol,iheadlib);}
     else if (idirmlib[0]!='\0'){ Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlib,&HeadLibFile,"");}
   }
   Show_LINENODEs(&HeadLibFile);

/*
   if (Nlibfile>0 ){  
      ln = &HeadLibFile;
      while (ln->next != NULL){
        ln = ln->next;
        n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMol,&(molL.filetype),ln->num);
        Nmollib += n;
      }
    }
 */

  printf("#Nmollib  %d\n",Nmollib);
  printf("#Nlibfile %d\n",Nlibfile);
  if ((Nmollib==0) && (Nlibfile==0)){
    printf("#ERROR:Nmollib=0 and Nlibfile=0.\n"); exit(1);
  }
   /*
  an = &HeadLibMol;
  while (an->next != NULL){
    an = an->next;
    printf(">%d %s\n",an->num,an->name);
  }
  */ 

  if (odeslib[0]  !='\0'){ 
    fpo  = fopen(odeslib,"w"); 
    if (fpo==NULL){printf("#ERROR:Can't write to odeslib '%s'\n",odeslib); exit(1);}
  }
  if (opropfile[0]!='\0'){ 
    fpop = fopen(opropfile,"w");
    if (fpop==NULL){printf("#ERROR:Can't write to opropfile '%s'\n",opropfile); exit(1);}
  }

  /** Read all molecules **/
   printf("## Read all molecules ##\n"); 

   if (Nlibfile==0){
     if (odeslib[0]!='\0')  { Write_Headers_of_AtomPair_Descriptor(fpo,Nmollib,idirlib,des_filemode,&HeadLibFile,molL.filetype); }
     if (opropfile[0]!='\0'){ Write_Headers_of_Property_File(fpop,Nmollib,idirlib,&HeadLibFile,molL.filetype);}
     printf("#Nmollib %d\n",Nmollib); fflush(stdout); 
     n = 0; 

     an = &HeadLibMol;
     while (an->next !=NULL){
       an = an->next;  
     
        printf(molL.name,"%s",an->name);
        sprintf(molL.filename,"%s/%s",idirlib,an->name);
        ok = Read_MOLECULE(molL.filename,&molL,molL.filetype,molL.atmhet,molL.BondType);
        if (an->name[0]=='\0'){ sprintf(an->name,"%d",an->num); }

       /*
        if (PAR.Nproperty>0){
           printf(">%s\n",molL.name);
          for (i=0;i<PAR.Nproperty;++i){
            printf("[%d] '%s':'%s'\n",i,PAR.PROP_NAME[i],Get_Value_from_KEYVALUE(&(molL.HeadProperty),PAR.PROP_NAME[i]));
          }
        }
        */

        if ((molL.Natom>0) && ((ltailstr[0]=='\0')||(Match_Tail_String(molL.name,ltailstr)==1)) ){ 
         if ((n%100)==0){ printf(">[%d/%d] %s %s Natom %d\n",n,Nmollib,molL.filename,molL.name,molL.Natom); }
       
         /* Set_Connected_Structure(&molL); */
         if (PAR.type_atompair=='D'){Set_DualAtomPair_Descriptor(&molL);}
         else if (PAR.type_atompair == 'N'){Set_AtomPair_Descriptor(&molL);}
         /* Set_Ring_Descriptor(&molL); */

        if ((AnnotationSmiles == 'T') && (molL.Natom>0)){
          make_LUNITs_from_Molecule_by_DFS(&HeadLunit,&molL);
          string = (char *)malloc(sizeof(char)*(length_string_from_LUNITs(&molL,&HeadLunit)));
          make_string_from_LUNITs(&molL,string,&HeadLunit); printf("#string '%s'\n",string);
          Add_key_value_string_to_KEYVALUEs("SMILES",string,&(molL.HeadProperty));
          Free_LUNITs(&HeadLunit);
          free(string);

         }

         if (odeslib[0]!='\0'){
           Write_AtomPair_Descriptor(fpo,&molL,an->num,an->name,an->class,des_filemode, an->num_file,an->offset_libfile);
         }
         if (opropfile[0]!='\0'){ Write_Property_for_One_Molecule(fpop,&molL,an->num,an->name,an->num_file,an->offset_libfile);}
        Free_MOLECULE(&molL);
        n += 1;
       }
     }
  }

   if (Nlibfile>0){
     if ((molL.filetype=='-')&&(ln->next != NULL)){
       printf("#ln->next->line '%s'\n",ln->next->line);
       molL.filetype = molecular_file_type_from_tail_of_file(ln->next->line);
     }
     if (odeslib[0]!='\0'){  Write_Headers_of_AtomPair_Descriptor(fpo,Nmollib,idirlib,des_filemode,&HeadLibFile,molL.filetype); }
     n = 0; 
     ln = &HeadLibFile;
     fname[0] = '\0'; 
     fpl = NULL;
     while (ln->next !=NULL){
       ln = ln->next;  
       printf("'%s'\n",ln->line);  
       sprintf(fname,"%s",ln->line); 
       fpl = fopen(fname,"r"); 
       if (fpl==NULL){
         printf("#ERROR:Can't open a library file '%s'.\n",fname);
         exit(1); }

       dir = opendir(fname);
       if (dir!=NULL){
         printf("#ERROR:a library file '%s' is a directory.\n",fname);
         exit(1); }
       closedir(dir);

       while (feof(fpl)==0){
         offset = ftell(fpl);
         ok = Read_One_Molecule_from_one_Library_file(fpl,&molL,molL.filetype,molL.atmhet,'L');
         if ((molL.Natom>0) && ((ltailstr[0]=='\0')||(Match_Tail_String(molL.name,ltailstr)==1)) ){ 
           if ((n%100)==0) 
             printf(">[%d: %d/%d] %s %s Natom %d\n",n,ln->num,Nlibfile,molL.filename,molL.name,molL.Natom); 
       
           /* Set_Connected_Structure(&molL);  */
           if (PAR.type_atompair=='D')      Set_DualAtomPair_Descriptor(&molL);
           else if (PAR.type_atompair=='N') Set_AtomPair_Descriptor(&molL);
           /* Set_Ring_Descriptor(&molL);  */
           if (odeslib[0]!='\0')  { Write_AtomPair_Descriptor(fpo,&molL,n,molL.name,"",des_filemode, ln->num,offset);}
           if (opropfile[0]!='\0'){ Write_Property_for_One_Molecule(opropfile,&molL,n,molL.name,'a',ln->num,offset);}
          n += 1;
         }
        Free_MOLECULE(&molL);
       }
     fclose(fpl);
     } /* ln */
  }

  Nmollib = n;
  if (odeslib[0]!='\0'){  
    printf("Write_AtomDescriptor(Nmollib %d) -->'%s'\n",Nmollib,odeslib);
    Write_Footers_of_AtomPair_Descriptor(fpo,Nmollib);
    fclose(fpo);
  }

  if (opropfile[0]!='\0'){ 
    printf("#Write_Molecular_Properties() -->'%s'\n",opropfile);
    fclose(fpop);
  }

 }



 /**********************************************/
 /*** MODE 'A': all-vs-all search            ***/
 /** (keep all the molecular data in memory) ***/
 /**********************************************/
 if (MODE=='A'){
   printf("#PAR.max_atompair_descriptor %d\n",PAR.max_atompair_descriptor);
   if (ideslib[0]!='\0') {
     HeadLibFile.next = NULL;
     HeadLibFile.num = -1; 
     Nmollib = Read_Headers_of_AtomPair_Descriptor(ideslib,idirlib,&des_filemode,&HeadLibFile,&(molL.filetype));
     Nlibfile = Number_of_LINENODEs(&HeadLibFile);
   }
   HeadLibMol.num  = 0;
   HeadLibMol.next = NULL;
   sim_one = sim_pair = sim_pair_norm = sim_ring = sim_oneringA = sim_oneringI = 0.0;
   if (ideslib[0]!='\0'){ 
     Nmollib = Read_AtomPair_Descriptor_into_LIBMOL(ideslib,&HeadLibMol,des_filemode,'A');
   }

   Nmolpair = Nmollib * (Nmollib + 1.0)/2.0; 
   printf("#Nmollib %d Nmolpair %.0lf\n",Nmollib,Nmolpair); 
  
   Nsta = (1.0-sqrt((double)(DIVbunbo-DIVbunshi)/DIVbunbo)) * Nmollib;  
   Nend = (1.0-sqrt((double)(DIVbunbo-DIVbunshi-1)/DIVbunbo)) * Nmollib;  
   if (DIVbunshi==0){ Nsta = 0;}
   if (DIVbunshi==(DIVbunbo-1)){ Nend = Nmollib;}

   an = &HeadLibMol;
   while (an->next!=NULL){
     an = an->next;
     printf("%s %s\n",an->name,an->molform);
   } 
   if ((filetype_osim=='S')||(filetype_osim=='D')){ Malloc_FLOAT2DMAP(&Dmap,Nmollib);}

   /** write similarities one pair by one pair */
   printf("#write_similarities_in_list_format()-->'%s' start\n",osimfile);

   if (osimfile[0]=='\0') sprintf(osimfile,"osl");
   if (osimfile[0]!='\0') fpo = fopen(osimfile,"w");

   if (CalMCS=='T'){ 
     printf("#READ ALL THE MOLECULE\n"); 
     fpL = NULL;
     an = &HeadLibMol;
     while (an->next!=NULL){
       an = an->next;
       an->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));
       sprintf(fname,"%s/%s",idirlib,an->name);
       if (Nlibfile==0){ 
          sprintf(fname,"%s/%s",idirlib,an->name);
          ok = Read_MOLECULE(fname,an->mol,molL.filetype,molL.atmhet,molL.BondType);
       }
       else {
         String_from_LINENODE_num(fname,&HeadLibFile,an->num_file);
         if (strcmp(fname,fname0)!=0){
           if (fpL!= NULL) fclose(fpL);
           fpL = fopen(fname,"r"); 
         }
       fseek(fpL,an->offset_libfile,SEEK_SET);
       ok = Read_One_Molecule_from_one_Library_file(fpL,an->mol,molL.filetype,molL.atmhet,'L');
       sprintf(fname0,"%s",fname);
       }
    }
  if (fpL != NULL) fclose(fpL);  
  }
 
   printf("#START CHECK SIMILRITIES\n");
   Write_Similarities_in_List_format_Header(fpo,&HeadLibMol,idirlib);
   fprintf(fpo,"#[numA:1] [nameA:2] [numB:3] [nameB:4] [sim_mcs:5] [sim_pair:6] [sim_pair_norm:7] [sim_one:8]\n");
   n = 0; 
   an = &HeadLibMol;
   while (an->next!=NULL){
     an = an->next;
     if ((Nsta <= an->num) && (an->num <Nend)){
       bn = an;
       while (bn->next!=NULL){
         bn = bn->next;
          
          sim_one = Tanimoto_bwn_oneatom_descriptors(an->oneatom_descriptor,bn->oneatom_descriptor);
          sim_pair = similarity_bwn_atompair_descriptors(an->atompair_descriptor,bn->atompair_descriptor);
          sim_pair_norm = similarity_bwn_atompair_descriptors_oneatom_normalized(an->atompair_descriptor, bn->atompair_descriptor);
          sim_ring = similarity_bwn_ring_descriptors(an->ring_descriptor,bn->ring_descriptor);
         /* printf("%s %s ",an->name,bn->name); */
           sim_oneringA = similarity_bwn_oneatom_ring_descriptors(an->oneatom_descriptor,an->ring_descriptor,bn->oneatom_descriptor,bn->ring_descriptor,'A');
           sim_oneringI = similarity_bwn_oneatom_ring_descriptors(an->oneatom_descriptor,an->ring_descriptor,bn->oneatom_descriptor,bn->ring_descriptor,'I');

             sim_mcs = 0.0;

             if (CalMCS=='F'){ 
               if ((filetype_osim=='D')||(filetype_osim=='S')) Dmap.map[an->num][an->num] = Dmap.map[bn->num][an->num] = sim_pair;
             }
  
             /*** cal MCS **/ 
             if (CalMCS=='T'){
               Match_Two_Molecules(an->mol,bn->mol,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
               if (optMlist.next != NULL){
                 sim_mcs  = Tanimoto_Coefficient(optMlist.next->Npair,an->mol,bn->mol);
               }
               Free_MATCHlist(&optMlist); 
               if ((filetype_osim=='D')||(filetype_osim=='S')) Dmap.map[an->num][bn->num] =  Dmap.map[bn->num][an->num] = sim_mcs;
             }
         
/*
              printf("%s %s NheavyA %d NheavyB %d Ncommon %d sim_mcs %f\n", an->name,bn->name,an->mol->Nheavyatom,bn->mol->Nheavyatom,optMlist.next->Npair,sim_mcs);
 */ 
              Get_Part_Of_Line(buffA,an->name,strlen(an->name)-16,strlen(an->name)); 
              Get_Part_Of_Line(buffB,bn->name,strlen(bn->name)-16,strlen(bn->name)); 
              if ((filetype_osim=='L')&&(osimfile[0]!='\0')){fprintf(fpo,"%d %s %d %s %f %f %f %f\n",
                       an->num,buffA,bn->num,buffB,sim_mcs,sim_pair,sim_pair_norm,sim_one);}
         ++n;
        if ((n%1000000)==0){printf("[%5.1f %%] %d/%.0lf\n",100.0*(float)n/Nmolpair,n,Nmolpair);}
      } /* bn */
    } /* Nsta <= an->num < Nend */ 
  } /* an */
  
  Write_Similarities_in_List_format_Tail(fpo);
  printf("#write_similarities_in_list_format()-->'%s' done.\n",osimfile);
  fclose(fpo);

  if (filetype_osim =='D')
    Write_Distance_FLOAT2DMAP_in_Phylip_format(osimfile,&Dmap,&HeadLibMol,'D',OMatWithProperty);
  if (filetype_osim =='S')
    Write_Distance_FLOAT2DMAP_in_Phylip_format(osimfile,&Dmap,&HeadLibMol,'S',OMatWithProperty);
}


 /**********************************************/
 /*** MODE 'a': all-vs-all search            ***/
 /** (not keep library descriptor in memory) ***/
 /**********************************************/

 if (MODE=='a'){
   HeadLibMol.num  = 0;
   HeadLibMol.next = NULL;
   descriptorA  = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor); 
   descriptorB  = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor); 
   if (ideslib[0]!='\0'){ 
     Read_Headers_of_AtomPair_Descriptor(ideslib,idirlib,&des_filemode,&HeadLibFile); 
     Nmollib = Read_AtomPair_Descriptor_into_LIBMOL(ideslib,&HeadLibMol,des_filemode,'>'); 
   }
   else { printf("#ERROR: the option -ides is not specified.\n"); exit(1);}

   printf("Nmolib %d\n",Nmollib);
   /*
   an = &HeadLibMol;
   while (an->next!=NULL){
     an = an->next;
     printf("%d %s\n",an->num,an->name);
   }   
   */ 

   if (ideslib[0]!='\0') fpl = fopen(ideslib,"r");

   Nmolpair = Nmollib * (Nmollib + 1)/2; 

   Nsta = (1.0-sqrt((double)(DIVbunbo-DIVbunshi)/DIVbunbo)) * Nmollib;  
   Nend = (1.0-sqrt((double)(DIVbunbo-DIVbunshi-1)/DIVbunbo)) * Nmollib;  
   if (DIVbunshi==0) Nsta = 0;
   if (DIVbunshi==(DIVbunbo-1)) Nend = Nmollib;


   printf("#Nmollib %d Npair %.0lf DIVbunshi %d DIVbunbo %d Nsta %d Nend %d\n",Nmollib,Nmolpair,DIVbunshi, DIVbunbo, Nsta, Nend); 

   /** write similarities one pair by one pair */
   if (ideslib[0]=='\0') {printf("#ERROR:option -isl is required for -M a.\n"); exit(1);} 

   fpl = fopen(ideslib,"r");
   
   if (osimfile[0]=='\0') {printf("#ERROR:option -osl is required for -M a.\n"); exit(1);} 
   fpo = fopen(osimfile,"w");
   
   printf("#write_similarities_in_list_format()-->'%s' start\n",osimfile);
   if (DIVbunshi==0) Write_Similarities_in_List_format_Header(fpo,&HeadLibMol,"",idirlib);
   n = 0; 
   an = &HeadLibMol;
   while (an->next!=NULL){
     an = an->next;
     a = 0; 
     if ((Nsta <= an->num) && (an->num <Nend)){
     bn = an;
     while (bn->next!=NULL){
       bn = bn->next;
       if (ThreTanimotoOneAtom>0.0){
         sim_one = Tanimoto_bwn_oneatom_descriptors(an->oneatom_descriptor,bn->oneatom_descriptor);
       }
       if ((ThreTanimotoOneAtom<0.0)||(sim_one>=ThreTanimotoOneAtom)){
         if (a==0){
           fseek(fpl,an->offset_descfile,SEEK_SET); 
           Read_AtomPair_Descriptor_into_LIBMOL_One_by_One(fpl,&HeadLibMol,descriptorA,des_filemode,'D'); 
           a = 1;
         }
         fseek(fpl,bn->offset_descfile,SEEK_SET); 
         Read_AtomPair_Descriptor_into_LIBMOL_One_by_One(fpl,&HeadLibMol,descriptorB,des_filemode,'D'); 
         sim_pair = similarity_bwn_atompair_descriptors(descriptorA,descriptorB);
         if ((ThreTanimotoAtomPair<0.0)||(sim_pair>=ThreTanimotoAtomPair)){
           if ((CalMCS!='T')||(ThreTanimotoMCS<0.0)) fprintf(fpo,"%d %d %f\n",an->num,bn->num,sim_pair);

           /*** cal MCS **/ 
           if ((CalMCS=='T')&&(ThreTanimotoMCS>=0.0)){
             if (an->mol==NULL){
               an->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));
               sprintf(fname,"%s/%s",idirlib,an->name);
               Read_MOLECULE(fname,an->mol,molL.filetype,molL.atmhet,molL.BondType);
              }
             if (bn->mol==NULL){
               bn->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));
               sprintf(fname,"%s/%s",idirlib,bn->name);
               Read_MOLECULE(fname,bn->mol,molL.filetype,molL.atmhet,molL.BondType);
             }
            Match_Two_Molecules(an->mol,bn->mol,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
            sim_mcs = 0.0; 
            if (optMlist.next != NULL){
              sim_mcs  = Tanimoto_Coefficient(optMlist.next->Npair,an->mol,bn->mol);
              if (ThreTanimotoMCS<0.0)
                fprintf(fpo,"%d %d %f %f %f\n",an->num,bn->num,sim_mcs,sim_pair,sim_one);
              else if (sim_mcs >= ThreTanimotoMCS)
                fprintf(fpo,"%d %d %f\n",an->num,bn->num,sim_mcs);
            }
           }
         }
       }
       if ((n%1000000)==0){printf("[%d/(%d-%d)] [%5.1f %%] %d/%.0lf\n",an->num,Nsta,Nend,100.0*(float)n/Nmolpair,n,Nmolpair);}
       ++n;
       if (bn->mol!=NULL) {Free_MOLECULE(bn->mol); free(bn->mol); bn->mol = NULL;}
     } /* bn */
    }
     if (an->mol!=NULL) {Free_MOLECULE(an->mol); free(an->mol); an->mol = NULL;}
   } /* an */
     
    Write_Similarities_in_List_format_Tail(fpo);
    printf("#write_similarities_in_list_format()-->'%s' done.\n",osimfile);
    fclose(fpo);
 }


 /*****************************************************************/ 
 /*** MODE 'M': Library Search for the multiple query molecules ***/
 /*****************************************************************/ 
 if (MODE=='M'){ 
   if (PAR.ConnectGraphType=='S'){
     ThreTanimotoMCS = ThreTanimotoOneAtom = ThreTanimotoAtomPair = 0.0;  
   }



   /** [1] Read only headers of the descriptor file **/
   if (ideslib[0]!='\0') {
     HeadLibFile.next = NULL;
     HeadLibFile.num = -1; 
     Nmollib = Read_Headers_of_AtomPair_Descriptor(ideslib,idirlib,&des_filemode,&HeadLibFile,&(molL.filetype));
     /* Show_LINENODEs(&HeadLibFile); */
     Nlibfile = Number_of_LINENODEs(&HeadLibFile);
   }

   descriptorL = (unsigned char *)malloc(sizeof(unsigned char)*PAR.max_atompair_descriptor); 

   /** [2] Read all the query molecules **/
   if (idirquery[0]!='\0') Nmolquery = Get_List_of_LIBMOLs_Under_Dir(idirquery,&HeadQueryMol,"");
   printf("#READ ALL THE QUERY MOLECULES\n"); 
   fpL = NULL;
   qn = &HeadQueryMol;
   while (qn->next!=NULL){
     qn = qn->next;
     qn->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));
     sprintf(fname,"%s/%s",idirquery,qn->name);
     Read_MOLECULE(fname,qn->mol,molQ.filetype,molQ.atmhet,molL.BondType);
     printf("#[%d] read '%s' Natom %d Nbond %d\n",qn->num,fname,qn->mol->Natom,qn->mol->Nbond); 
     if (PAR.type_atompair=='D')        Set_DualAtomPair_Descriptor(qn->mol);
     else if (PAR.type_atompair == 'N') Set_AtomPair_Descriptor(qn->mol);
   }
  printf("#Nmolquery %d\n",Nmolquery);
  /** [3] Read Library List One by One and compare with molQ by atompair descriptor **/
  printf("## DO ATOMPAIR DESCRIPTOR SEARCH ##\n"); 
  Nmcs = 0;
  Nmol = 0;
  if (ideslib[0]!='\0') fpl = fopen(ideslib,"r");
  if (fpl==NULL) {printf("#ERROR:Can't open description lib file '%s'\n",ideslib); exit(1);}
  HeadLibMol.num  = 0;
  HeadLibMol.next = NULL;
  
  do{
    xn = Read_AtomPair_Descriptor_into_LIBMOL_One_by_One(fpl,descriptorL,des_filemode,'A'); 

    if(xn!=NULL){
      /*  printf(">%d '%s' %d\n",xn->num,xn->name,xn->offset_libfile); */ 
      xn->mul_tanimoto_mcs      =  (float *)malloc(sizeof(float)*Nmolquery);
      xn->mul_tanimoto_atompair =  (float *)malloc(sizeof(float)*Nmolquery);
      xn->mul_tanimoto_oneatom  =  (float *)malloc(sizeof(float)*Nmolquery);
      
      xn->tanimoto_atompair = xn->tanimoto_oneatom = xn->tanimoto_mcs = 0.0;

      qn = &HeadQueryMol;
      while (qn->next != NULL){

        qn = qn->next;

        accept = 1;

        if (PAR.ConnectGraphType=='S'){
           if ((Substructure_condition_oneatom_descriptors(qn->mol->oneatom_descriptor,xn->oneatom_descriptor)==0) || 
              (substructure_condition_atompair_descriptors(qn->mol->atompair_descriptor,descriptorL)==0)){accept = 0;}
         }

         if (accept==1){
             xn->mul_tanimoto_oneatom[qn->num]  = Tanimoto_bwn_oneatom_descriptors(qn->mol->oneatom_descriptor,xn->oneatom_descriptor);
             xn->mul_tanimoto_atompair[qn->num] = similarity_bwn_atompair_descriptors_oneatom_normalized(qn->mol->atompair_descriptor, descriptorL);
         } 
      } /* qn */

      xn->tanimoto_oneatom  = fusion_scores(Nmolquery,xn->mul_tanimoto_oneatom,FusionScDescriptor);
      xn->tanimoto_atompair = fusion_scores(Nmolquery,xn->mul_tanimoto_atompair,FusionScDescriptor);
      xn->score_for_sort = xn->tanimoto_atompair;


      accept = 1;
      if (((ThreTanimotoOneAtom >=0.0)&&(xn->tanimoto_oneatom  < ThreTanimotoOneAtom)) || 
          ((ThreTanimotoAtomPair>=0.0)&&(xn->tanimoto_atompair < ThreTanimotoAtomPair)) ){accept = 0;}
       }

      if (accept==0){
        free(xn->name); free(xn->molform); 
        free(xn->mul_tanimoto_mcs); free(xn->mul_tanimoto_atompair); free(xn->mul_tanimoto_oneatom); 
        HeadLibMol.next = xn->next;
        if (xn->next != NULL) xn->next->prev = &HeadLibMol;
        free(xn);
      }

      if (accept==1){ Nmcs += 1; }

       Nmol += 1; 
       if ((Nmol%5000)==0){ 
          sprintf(buffA,"Step I:descriptor search: %5.1f %% %d/%d",100.0*(float)Nmol/Nmollib,Nmol,Nmollib);
          write_progress_message(oprgfile,buffA);
       }

  } while (xn != NULL); 

  printf("#Number_of_LIBMOL_after_descriptor_search %d\n",Number_of_LIBMOL(&HeadLibMol));

/*
 printf(">INIT\n"); 
 Show_LIST_of_LIBMOL(&HeadLibMol);
 
 reverse_LIBMOL_list(&HeadLibMol);
 
 printf(">REVERSED\n"); 
 Show_LIST_of_LIBMOL(&HeadLibMol);
 */

  Merge_Sort_Double_Linked_List_LIBMOL(&HeadLibMol,'D');

 /*
  xn = &HeadLibMol;
  while (xn->next != NULL){
    xn = xn->next;
    printf("[%d] %f ",xn->num,xn->tanimoto_atompair);
    for (i=0;i<Nmolquery;++i) printf(" %d:%f",i,xn->mul_tanimoto_atompair[i]);
    printf("\n");
  }
  */
 
  /** Calculate MCS-based similarity  **/
  sprintf(buffA,"Step I: descriptor search: %5.1f %% %d/%d",100.0*Nmol/Nmollib,Nmol,Nmollib);
  write_progress_message(oprgfile,buffA);

  printf("## DO MCS COMPARISON for [%d/%d] over tanimoto_atompair %f (NmaxMCS %d)##\n",Nmcs,Nmollib,ThreTanimotoAtomPair,NmaxMCS); 

  if ((NmaxMCS>0)&&(Nmcs>NmaxMCS)) Nmcs = NmaxMCS;  

  Nmol = 0; 
  xn = &HeadLibMol;
  fpL = NULL;
  fname[0] = fname0[0] = '\0';
  Ncalmatch  =  0;
  while (xn->next !=NULL){
    xn = xn->next;
    xn->score_for_sort = 0.0;
    accept = 0;
    if ((NmaxMCS < 0) || (Ncalmatch < NmaxMCS)){

      sprintf(molL.name,"%s",xn->name);
      sprintf(molL.filename,"%s/%s",idirlib,xn->name);
 
      if (Nlibfile==0){ 
        ok = Read_MOLECULE(molL.filename,&molL,molL.filetype,molL.atmhet,molL.BondType);
      }
      else {
       String_from_LINENODE_num(fname,&HeadLibFile,xn->num_file);
       if (strcmp(fname,fname0)!=0){
         if (fpL!= NULL) fclose(fpL);
         fpL = fopen(fname,"r"); 
       }
       fseek(fpL,xn->offset_libfile,SEEK_SET);
       ok = Read_One_Molecule_from_one_Library_file(fpL,&molL,molL.filetype,molL.atmhet,'L');
       if ((xn->name[0]=='\0')&&(molL.name[0]!='\0')) sprintf(xn->name,"%s",molL.name); 
       sprintf(fname0,"%s",fname);
     }


      if (ok>0){ 
        /* printf(">DoMCS '%s' score_for_sort %f\n",xn->name,xn->score_for_sort); */
        xn->Nheavyatom = molL.Nheavyatom; 
        xn->molform = (char *)malloc(sizeof(char)*(strlen(molL.molform)+1));
        sprintf(xn->molform,"%s",molL.molform);
        if (molL.Nheavyatom>0){
           xn->doneMCS = 1;
           Ncalmatch += 1;
           qn = &HeadQueryMol;

           while (qn->next != NULL){
             qn = qn->next;
             Match_Two_Molecules(qn->mol,&molL,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
             xn->mul_tanimoto_mcs[qn->num] = 0.0;
             if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E') && (optMlist.next->Npair>0)){
                if (qn->num==0) Copy_MATCH_Info_to_LIBMOL(xn,optMlist.next,qn->mol,&molL);
                xn->mul_tanimoto_mcs[qn->num] = Tanimoto_Coefficient(optMlist.next->Npair,&molL,qn->mol);
                accept = 1;
             }

           Free_MATCHlist(&optMlist);

          } /* qn */

          xn->tanimoto_mcs = fusion_scores(Nmolquery,xn->mul_tanimoto_mcs,FusionScMCS);
          xn->score_for_sort = xn->tanimoto_mcs;
/*
          printf("[mcs_cal %d] [%d] %s Natom %d Npair %d tanimoto_atompair %f tanimoto_mcs %f select_dis %f\n",
            Ncalmatch,xn->num,xn->name,molQ.Nheavyatom,xn->Npair,xn->tanimoto_atompair,xn->tanimoto_mcs,xn->select_dis);
 */   
          Free_MOLECULE(&molL);
        }
     }

     Nmol += 1;
     if ((Nmol%50)==0){ 
       sprintf(buffA,"Step II: MCS search: %5.1f %% %d/%d",100.0*(float)Nmol/Nmcs,Nmol,Nmcs);
       write_progress_message(oprgfile,buffA);
     }
   }

     if (accept==0){ 
       /* printf("REMOVE %s\n",xn->name); */
       if (xn->prev != NULL) xn->prev->next = xn->next;
       if (xn->next != NULL) xn->next->prev = xn->prev;
     }

  } /* while (xn) */
 
  Merge_Sort_Double_Linked_List_LIBMOL(&HeadLibMol,'D');

  if (oschfile[0] != '\0'){
   Write_Search_Result(oschfile,HeadQueryMol.next->mol,idirlib,ideslib,&HeadLibMol,Nmollib,ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,&HeadLibFile);  
  }

   write_search_result_for_multiple_query_search("woops.out",&HeadQueryMol,&HeadLibMol,
         ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,Nlibfile,&HeadLibFile,Nmax_outmol);  
 
   if (omolsfile[0] != '\0'){
   Write_Similar_Molecules_To_One_File(omolsfile,"",idirlib,&HeadLibMol,
         ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,molL.filetype,Nlibfile,&HeadLibFile,Nmax_outmol);  
  }

  if (omolsdir[0] != '\0'){
   Write_Similar_Molecules_Under_Dir(omolsdir,"",idirlib,&HeadLibMol,
     ThreTanimotoMCS,ThreTanimotoOneAtom,ThreTanimotoAtomPair,molL.filetype,Nlibfile,&HeadLibFile,Nmax_outmol);  
  }


  write_progress_message(oprgfile,"FINISH_CALCULATION");
}

 /*********************************************************************/ 
 /*** MODE 'G': Comparison bwn Two Groups: Query-library vs Library ***/
 /*********************************************************************/ 
 if (MODE=='G'){ 

  /** [1] Read Descriptors for the library */
   if (ideslib[0]!='\0'){ 
     Nmollib = Read_Headers_of_AtomPair_Descriptor(ideslib,idirlib,&des_filemode,&HeadLibFile,&(molL.filetype));
     Nmollib = Read_AtomPair_Descriptor_into_LIBMOL(ideslib,&HeadLibMol,des_filemode,'A');
   }

  /** [2] Read List of Query library molecules **/
  if (idirlibQ[0]!='\0') NmollibQ = Get_List_of_LIBMOLs_Under_Dir(idirlibQ,&HeadLibMolQ,"");
  Nsta = (NmollibQ*DIVbunshi)/DIVbunbo;
  Nend = (NmollibQ*(DIVbunshi+1))/DIVbunbo;
  printf("#NlibmolQ %d div %d/%d Nsta %d Nend %d\n",NmollibQ,DIVbunshi,DIVbunbo,Nsta,Nend);

  /** [3] Compare Query-library vs library */
   if (osimfile[0]=='\0') sprintf(osimfile,"sim.out");
   fpo = fopen(osimfile,"w");

   printf("#write_similarity_list -->'%s'\n",osimfile);
   fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
   fprintf(fpo,"#Nmollib  %d\n",Nmollib);
   fprintf(fpo,"#NmollibQ %d\n",NmollibQ);
   fprintf(fpo,"#div %d/%d Nsta %d Nend %d\n",DIVbunshi,DIVbunbo,Nsta,Nend);
   fprintf(fpo,"#[numQ] [nameQ] [numL] [nameL] [similarity]");
   qn = &HeadLibMolQ;
   while (qn->next != NULL){
     qn = qn->next;
     if ((Nsta<= qn->num) && (qn->num < Nend)){
       /* qn->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE)); */
       sprintf(fname,"%s/%s",idirlibQ,qn->name);
       Read_MOLECULE(fname,&molQ,molQ.filetype,molQ.atmhet,molL.BondType); 

       printf(">query[%d/%d] %s read '%s' Natom %d Nbond %d\n",qn->num,NmollibQ,qn->name,fname,molQ.Natom,molQ.Nbond); 
       if (molQ.Nheavyatom>0){
         if (PAR.type_atompair=='D') Set_DualAtomPair_Descriptor(&molQ);
         else if (PAR.type_atompair == 'N') Set_AtomPair_Descriptor(&molQ);

         xn = &HeadLibMol;
         while (xn->next != NULL){ 
           xn = xn->next;
           xn->tanimoto_oneatom = xn->tanimoto_atompair = xn->tanimoto_mcs = 0.0;
           accept = 1;


           if (PAR.ConnectGraphType=='I'){
             if ((Isomorphic_condition_oneatom_descriptors(molQ.oneatom_descriptor,xn->oneatom_descriptor)==0) || 
               (isomorphic_condition_atompair_descriptors(molQ.atompair_descriptor,xn->atompair_descriptor)==0)){accept = 0;}
           }

          if (PAR.ConnectGraphType=='S'){
             if ((Substructure_condition_oneatom_descriptors(molQ.oneatom_descriptor,xn->oneatom_descriptor)==0) || 
                (substructure_condition_atompair_descriptors(molQ.atompair_descriptor,xn->atompair_descriptor)==0)){accept = 0;}
           }

           if (accept==1){
             xn->tanimoto_oneatom = Tanimoto_bwn_oneatom_descriptors(molQ.oneatom_descriptor,xn->oneatom_descriptor);
             xn->tanimoto_atompair = similarity_bwn_atompair_descriptors_oneatom_normalized(molQ.atompair_descriptor,xn->atompair_descriptor);
             xn->score_for_sort = xn->tanimoto_atompair; 

            if (((ThreTanimotoOneAtom >=0.0)&&(xn->tanimoto_oneatom  < ThreTanimotoOneAtom)) || 
               ((ThreTanimotoAtomPair>=0.0)&&(xn->tanimoto_atompair < ThreTanimotoAtomPair)) ){accept = 0;}

           }

           if (accept==1){
             sprintf(molL.name,"%s",xn->name);
             sprintf(molL.filename,"%s/%s",idirlib,xn->name);
             ok = Read_MOLECULE(molL.filename,&molL,molL.filetype,molL.atmhet,molL.BondType);
             if (molL.Nheavyatom>0){
               xn->doneMCS = 1;
               Match_Two_Molecules(&molQ,&molL,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
               Ncalmatch += 1;
               if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E') && (optMlist.next->Npair>0)){
                  xn->tanimoto_mcs   = Tanimoto_Coefficient(optMlist.next->Npair,&molQ,&molL);
                  xn->score_for_sort = xn->tanimoto_mcs; 
                  accept = 1;
               }
               Free_MATCHlist(&optMlist);
             }
             Free_MOLECULE(&molL);
            if (xn->tanimoto_mcs>=ThreTanimotoMCS){
              fprintf(fpo,"%d %s %d %s %f\n",qn->num,qn->name,xn->num,xn->name,xn->tanimoto_mcs);
             }
           }
          } /* xn */
        } /* if (molQ.NNheavyatom >0 */

        Free_MOLECULE(&molQ);

     } /* (Nsta<= qn->num) && (qn->num < Nend)) */

   } /* qn */

   printf("#write_similarity_list -->'%s'\n",osimfile);
   Set_END_DATE();
   fprintf(fpo,"#DATE_START %s\n",PAR.START_DATE);
   fprintf(fpo,"#DATE_END   %s\n",PAR.END_DATE);
   fprintf(fpo,"#COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);
   fclose(fpo);
 }

  return(1);
} /* end of main() */ 















float fusion_scores(N,scores,type)
  int N;
  float *scores;
  char type;   /* 'A'verage, 'X':maximum */
{
  int i;
  float sc;

  if (type=='A'){
    sc = 0.0;
    for (i=0;i<N;++i){ sc += scores[i]; }
    sc = sc/N;
  }
  else {
    sc = scores[0];
    for (i=1;i<N;++i){ 
      if (scores[i]>sc) sc = scores[i]; 
    }
  }
  return(sc);
}


void write_search_result_for_multiple_query_search(ofname,HeadQueryMol,HeadLibMol,
       thre_tanimoto_mcs,thre_tanimoto_oneatom,thre_tanimoto_atompair,Nlibfile,HeadLibFile,Nmax_outmol)
 char   *ofname;
 struct LIBMOL *HeadQueryMol;
 struct LIBMOL *HeadLibMol;
 float  thre_tanimoto_mcs;
 float  thre_tanimoto_oneatom;
 float  thre_tanimoto_atompair;
 int    Nlibfile;
 struct LINENODE *HeadLibFile;
 int    Nmax_outmol; /* maximum number of output molecules */
{
  FILE *fpo;
  int Nquerymol,q;
  struct LIBMOL *qn,*xn;

  printf("#write_search_result_for_multiple_query_search() --> '%s'\n",ofname);
  fpo = fopen(ofname,"w");
  Nquerymol = Number_of_LIBMOL(HeadQueryMol);
  qn = HeadQueryMol;
  while (qn->next != NULL){
    qn = qn->next;
    fprintf(fpo,"#QUERY %d %s\n",qn->num,qn->name);
  } 
  
  xn = HeadLibMol;
  while (xn->next != NULL){
    xn = xn->next;
    fprintf(fpo,"%-5d %-20s %.3f ",xn->num,xn->name,xn->tanimoto_mcs);
    for (q=0;q<Nquerymol;++q){
      fprintf(fpo," %.3f",xn->mul_tanimoto_mcs[q]);
    }
    fprintf(fpo,"\n");
  }
  fclose(fpo);


} /* end of write_search_result_for_multiple_query_search() */




void write_progress_message(oprgfile,message)
  char *oprgfile;
  char *message;
{
  FILE *fpo;

  if (oprgfile[0] != '\0'){
    fpo = fopen(oprgfile,"w");
    fprintf(fpo,"%s\n",message);
    fclose(fpo);
  }
}
