/*
 
 <lkcombu.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================


  for obtaining maximum common substructure
  by build-up procedures using molecular library
  
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "ioLINE.h"
#include "ioSDF.h"
#include "ioKCF.h"
#include "ioMOL2.h"
#include "ioPDB.h"
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
#include "ringblock.h"
#include "RingDesc.h"
#include "options.h"

struct PARAMETERS PAR;

/** FUNCTIONS (LOCAL) **/
static int Read_Library_MOLECULE_Non_Sequential();

int main(argc,argv)
 int argc;
 char **argv;
{
  struct MOLECULE molQ,molL;
  int L,k,a,n,r,ok;
  char MODE,DetailOptHelp;
  struct MATCH optMlist;

  /** variables for library */
  char ilistlib[MAX_FILENAME],idirlib[MAX_FILENAME],iheadlib[MAX_FILENAME],idirmlib[MAX_FILENAME];
  struct LINENODE HeadLibFile;
  struct LIBMOL  HeadLibMol,HeadLibMolSelected;
  int    Nlibfile;
  struct LINENODE *ln;
  int    Nmollib;

  /** variables for query library */
  char ilistlibQ[MAX_FILENAME],idirlibQ[MAX_FILENAME],iheadlibQ[MAX_FILENAME],idirmlibQ[MAX_FILENAME];
  struct LINENODE HeadLibFileQ;
  struct LIBMOL  HeadLibMolQ;
  int    NlibfileQ;
  int    NmollibQ;

  int    Nmax_outmol; 
  char odirlib[MAX_FILENAME],oschfile[MAX_FILENAME],osimfile[MAX_FILENAME], osmifile[MAX_FILENAME];
  char filetype_osim; 
  int  Nsimilar_pair; 
  char ocommon_rasfile[MAX_FILENAME];
  char omolsfile[MAX_FILENAME],omolsdir[MAX_FILENAME]; 
  char ScanOpeType,ActInactComp,KeepLargestConGraph;
  char ofname[MAX_FILENAME],buff[MAX_FILENAME];
  float ThreTanimotoOneAtom,ThreTanimotoMCS,tanimoto,f;
  struct LIBMOL *an,*bn,*qn,*xn;
  struct FLOAT2DMAP Dmap;
  char  OMatWithProperty,doMCS; 
  int   DIVbunshi, DIVbunbo,Nsta,Nend;
  int   Ncalmatch,Nhitmol;
  struct PERMUTATION *permu; 
  char fname[MAX_FILENAME],fname0[MAX_FILENAME];
  FILE *fpL,*fpo;
  int Nmollib_nonredun;
  char comment_str[512];
  char *string;
  struct LUNIT HeadLunit;
  struct OPTION_VALUE *ov;

 
  /****** SET INITIAL PARAMETERS ***********/
  Set_Default_Global_PAR(); 
  PAR.PROGRAM_TYPE = 'l'; 
  MODE = 'S';
  Initialize_MOLECULE(&molQ);
  Initialize_MOLECULE_string(&molQ);
  Initialize_MOLECULE(&molL);
  Initialize_MOLECULE_string(&molL);
  molQ.filename[0] = molL.filename[0] = '\0';
  molQ.filetype = molL.filetype = '-'; 
  molQ.atmhet = molL.atmhet =  'B';
  molQ.BondType = molL.BondType =  'B';
  molQ.chain = molL.chain = '-';
  molQ.resi[0] = molL.resi[0] = '\0'; 
  molQ.SetStereo3D = molL.SetStereo3D = 'F';

  ocommon_rasfile[0] =  '\0';
  sprintf(oschfile,"sch.out");
  filetype_osim = 'N';
  sprintf(osimfile,"sim.out");

  ilistlib[0] = idirlib[0] = idirmlib[0] = iheadlib[0] = '\0';
  Nlibfile = Nmollib = 0;

  ilistlibQ[0] = idirlibQ[0] = idirmlibQ[0] = iheadlibQ[0] = '\0';
  NlibfileQ = NmollibQ = 0;

  odirlib[0] = osmifile[0] = '\0';
  ThreTanimotoOneAtom = ThreTanimotoMCS = 0.0;
  OMatWithProperty = 'F';
  DIVbunshi = 0; DIVbunbo = 1;
  ScanOpeType = '-';

  ActInactComp = 'F';
  
  HeadLibFile.next = NULL;
  HeadLibFile.num  = -1;
  
  HeadLibMol.next = NULL;
  HeadLibMol.num  = -1;

  HeadLibFileQ.next = NULL;
  HeadLibFileQ.num  = -1;
  
  HeadLibMolQ.next = NULL;
  HeadLibMolQ.num  = -1;

  omolsfile[0] = omolsdir[0] = '\0'; 
  Nmax_outmol = -1;
  Nmollib_nonredun = 10;

  KeepLargestConGraph = 'F';


  DetailOptHelp = '-';
  if((argc>1) && ((strcmp(argv[1],"-h")==0) || (strcmp(argv[1],"-help")==0)) )         DetailOptHelp = 'D';
  if((argc>1) && ((strncmp(argv[1],"-hm",3)==0) || (strncmp(argv[1],"-HM",3)==0)) ) DetailOptHelp = 'M';

  if ((argc<2)&&(DetailOptHelp!=1)){
    printf("lkcombu <options>\n");
    printf("  %s\n",KCOMBU_STANDS_FOR);
    printf("  'lkcombu' is for query-library search and all-vs-all comparison.\n");
    printf("  %s LastModified:%s\n",CODED_BY,LAST_MOD_DATE);
    printf("** Simple Usage for 'lkcombu' **\n");
    printf(">for query-library search\n");  
    printf("  $lkcombu -M S -Q [query_mol] -ifl  [files of lib multi mols] -fL [S|P|2] -osc [output_file]\n");
    printf("  $lkcombu -M S -Q [query_mol] -idl  [dir   of lib mols]       -fL [S|P|2] -osc [output_file]\n");
    printf("  $lkcombu -M S -Q [query_mol] -idml [dir   of lib multi mols] -fL [S|P|2] -osc [output_file]\n");
    printf(">for all-vs-all comparison\n");
    printf("  $lkcombu -M A -ifl  [file of lib. multi mols] -fL [S|P|2] -osm [output_file]\n");
    printf("  $lkcombu -M A -idl  [dir of lib. mols]        -fL [S|P|2] -osm [output_file]\n");
    printf("  $lkcombu -M A -idml [dir of lib. multi mols]  -fL [S|P|2] -osm [output_file]\n");
    printf(" (-fL : file formats.'P'db,'S'df,'K'cf, '2':MOL2)\n");
    printf(">for comparison bwn two groups \n");
    printf("  $lkcombu -M G -idlq [dir of query lib. mols] -idl [dir of lib mols] -fQ [S|P|2] -fL [S|P|2] -osm [output_file]\n");
    printf(">for randomly selecting library\n");
    printf("  $lkcombu -M R -idl [dir of lib. mols] -nlr [num_select_mol] -tam [threMCS] -oms [out_multi_mol] -omsd [out_mol_dir]\n");
    printf(">for convert liblary into SMILES strings\n");
    printf("  $lkcombu -M C -ifl [file of lib. multi mols] -osmi [out_smiles_str_file] -kplcn [T/F]\n");
    printf("  $lkcombu -M C -idl [dir  of lib. mols]       -osmi [out_smiles_str_file] -kplcn [T/F]\n");
    printf(">for showing options in detail\n");
    printf("  $lkcombu -h\n");
    printf(">for showing options for MCS\n");
    printf("  $lkcombu -hmcs\n");
    }
  if (DetailOptHelp=='D'){
    printf("lkcombu <options>\n");
    printf("  for obtaining maximum common substructures by build-up method.\n");
    printf("  %s\n",KCOMBU_STANDS_FOR); 
    printf("   coded by T.Kawabata. LastModified:%s\n",LAST_MOD_DATE);
    printf("<mode option for 'lkcombu'>\n");
    printf(" -M   : MODE 'S'earch_for_query_vs lib 'A'll_vs_all 'G':comparison bwn two groups [%c]\n",MODE);
    printf("      :      'R'andom_select_lib, 'C'onverting molecule types[%c]\n",MODE);
    printf("<options for input library for '-M S','-M A','-M G','-M R'>\n");
    printf(" -ifl : input files for all the library molecules (file1:file2:..) [%s]\n","");
    printf(" -idl : input directory for library molecules   [%s]\n",idirlib);
    printf(" -idml: input directory for multi-mol-library molecular file(multi-mol-one-file) [%s]\n",idirmlib);
    printf(" -ill : input list of library molecules[%s]\n",ilistlib);
    printf(" -ihl : required header string of input file (combined with -idl option) [%s]\n",iheadlib);
    printf(" -o   : output calculation process 'T' or 'F' [%c]\n",PAR.OutCalProcess);
    printf(" -fL  : file formats.'P'db,'S'df,'K'cf, '2':MOL2 [%c]\n",molL.filetype);    
    printf(" -aL  : AtomHetero type. 'A'tom 'H'etatm 'B'oth[%c]\n",molL.atmhet);    
    printf("<options for '-M S','-M A','-M G'>\n"); 
    printf(" -tao : Tanimoto threshold value for one atom descriptor (filtering) [%f]\n",ThreTanimotoOneAtom);
    printf(" -tam : Tanimoto threshold value for MCS (output)                 [%f]\n",ThreTanimotoMCS);
    printf("<options for '-M A' and '-M G'>\n"); 
    printf(" -osm : output file for similarities  [%s]\n",osimfile);
    printf(" -fsm : file formats for osm. 'L':simple pairwise list.'N':list with molecular names\n");
    printf("      :                       'D':distance matrix.  'S'imilarity matrix [%c]\n",filetype_osim);
    printf(" -div : Job division. {0..Ndivision-1}/(Ndivision). [%d/%d]\n",DIVbunshi,DIVbunbo);
    printf("<options for '-M S'>\n"); 
    printf(" -Q   : query molecule Q (*.sdf|*.mol2|*.pdb|*.kcf|*.smi)[%s]\n",molQ.filename);
    printf(" -fQ  : file formats.'P'db,'S'df,'K'cf, '2':MOL2, 'M':SMILES [%c]\n",molQ.filetype);    
    printf(" -aQ  : AtomHetero type. 'A'tom 'H'etatm 'B'oth[%c]\n",molQ.atmhet);    
    printf(" -bQ  : Bond type. 'B':consider bonds, otherwise, not consider bonds[%c]\n",molQ.BondType);
    printf(" -osc : output search result file   [%s]\n",oschfile);
    printf(" -oms : output file for similar molecules     [%s]\n",omolsfile);
    printf(" -omsd: output dir  for files for similar molecules [%s]\n",omolsdir);
    printf(" -nms : max number of output similar molecules. -1:don't care. [%d]\n",Nmax_outmol);
    printf(" -ste3Q : Setup Stereo Parity from 3D for mol Q (T or F)[%c]\n",molQ.SetStereo3D);
    printf(" -ste3L : Setup Stereo Parity from 3D for mol L (T or F)[%c]\n",molL.SetStereo3D);
    printf(" -ocras : output rasmol script file for commmon chemical structuctures [%s]\n",ocommon_rasfile);
    printf("<options for '-M A'>\n"); 
    printf(" -ocras : output rasmol script file for commmon chemical structuctures [%s]\n",ocommon_rasfile);
    printf(" -wpr : ourput dis/sim matrix with class names ('T' or 'F')[%c]\n",OMatWithProperty);
    printf(" -act : active-vs-inactive comparison. avoid 'inactive'-'inactive'.   ('T' or 'F')[%c]\n",ActInactComp);
    printf("<options for '-M G'>\n"); 
    printf(" -iflq  : input files for the query library molecules (file1:file2:..) [%s]\n","");
    printf(" -idlq  : input directory for query library molecules   [%s]\n",idirlibQ);
    printf(" -idmlq : input directory for multi-mol query library molecular file(multi-mol-one-file) [%s]\n",idirmlibQ);
    printf(" -illq  : input list of query library molecules[%s]\n",ilistlibQ);
    printf(" -kplcn : keep only the largest connected graph (T or F)[%c]\n",KeepLargestConGraph);
    printf("<options for '-M C'>\n"); 
    printf(" -sca  : output. 'C': only connected parts(desalting) 'D':with DABRL label 'R':Ring plane and DABL sphere [%c]\n",ScanOpeType);
    printf(" -odl  : output directory for library molecules   [%s]\n",odirlib);
    printf(" -osmi : output file for SMILES string file  [%s]\n",osmifile);
    printf("<options for 'lkcombu -M R'>\n");
    printf(" -nlr : number of selected library molecule for non-redundant library [%d]\n",Nmollib_nonredun); 
  }
 
  if (DetailOptHelp=='M'){
    Show_Option_Help_MCS();
  }

  if ((argc<2)||(DetailOptHelp!='-')) exit(1);

 /****** READ ARGUMENTS ***********/
 PAR.COMMAND[0] = '\0';
 for (k=0;k<argc;++k){ 
   if (k>0) strcat(PAR.COMMAND," ");
   strcat(PAR.COMMAND,argv[k]);
}

/****** READ ARGUMENTS ***********/
 PAR.OptValHead.next = NULL;
 Read_Options_From_Arguments(argc,argv,PAR.COMMAND,&(PAR.OptValHead));

 ov = &(PAR.OptValHead);
 while (ov->next != NULL){
   ov = ov->next;
   if (Set_Global_PAR_from_OPTION_VALUE(ov)==0){
           if (strcmp(ov->opt,"M")==0){ MODE = ov->val[0];}
      else if (strcmp(ov->opt,"Q")==0){ sprintf(molQ.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"ifl")==0){ 
        Nlibfile = Make_LibFileList_from_Colon_Separated_String(&HeadLibFile,ov->val);
        Show_LINENODEs(&HeadLibFile);
      }
      else if (strcmp(ov->opt,"iflq")==0){ 
        NlibfileQ = Make_LibFileList_from_Colon_Separated_String(&HeadLibFileQ,ov->val);
        Show_LINENODEs(&HeadLibFileQ);
      }
      else if (strcmp(ov->opt,"idml")==0){ sprintf(idirmlib,"%s",ov->val);}
      else if (strcmp(ov->opt,"idl")==0) { sprintf(idirlib,"%s",ov->val);}
      else if (strcmp(ov->opt,"ihl")==0) { sprintf(iheadlib,"%s",ov->val);}
      else if (strcmp(ov->opt,"ill")==0) { sprintf(ilistlib,"%s",ov->val);}
      else if (strcmp(ov->opt,"odl")==0) { sprintf(odirlib,"%s",ov->val);}
      else if (strcmp(ov->opt,"osm")==0) { sprintf(osimfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"fsm")==0) { filetype_osim = ov->val[0];}
      else if (strcmp(ov->opt,"wpr")==0) {OMatWithProperty = ov->val[0];}
      else if (strcmp(ov->opt,"osc")==0) { sprintf(oschfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"ocras")==0) { sprintf(ocommon_rasfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"tao")==0) { ThreTanimotoOneAtom = atof(ov->val);}
      else if (strcmp(ov->opt,"tam")==0) { ThreTanimotoMCS   = atof(ov->val);}
      else if (strcmp(ov->opt,"sca")==0) {  ScanOpeType = ov->val[0];}
      else if (strcmp(ov->opt,"act")==0) { ActInactComp = ov->val[0];}
      else if (strcmp(ov->opt,"nlr")==0) {  Nmollib_nonredun = atoi(ov->val);}
      else if (strcmp(ov->opt,"fQ")==0) { molQ.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"fL")==0) { molL.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"aQ")==0) { molQ.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"aL")==0) { molL.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"bQ")==0) { molQ.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"bL")==0) { molL.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"idlq")==0) { sprintf(idirlibQ,"%s",ov->val);}
      else if (strcmp(ov->opt,"idmlq")==0) { sprintf(idirmlibQ,"%s",ov->val);}
      else if (strcmp(ov->opt,"oms")==0)  { sprintf(omolsfile,"%s",ov->val); }
      else if (strcmp(ov->opt,"omsd")==0)  { sprintf(omolsdir,"%s",ov->val); }
      else if (strcmp(ov->opt,"omsi")==0)  { sprintf(osmifile,"%s",ov->val); }
      else if (strcmp(ov->opt,"kplcn")==0) { KeepLargestConGraph = ov->val[0]; }
      else if (strcmp(ov->opt,"nms")==0) { Nmax_outmol = atoi(ov->val);}
      else if (strcmp(ov->opt,"ste3Q")==0) { molQ.SetStereo3D = ov->val[0];}
      else if (strcmp(ov->opt,"ste3L")==0) { molL.SetStereo3D = ov->val[0];}
      else if (strcmp(ov->opt,"div")==0) { 
         a = First_Appear_Pos_in_String(ov->val,'/');
         L = strlen(ov->val);
         if ((L>1)&&(a>0)){ 
           Get_Part_Of_Line(buff,ov->val,a+1,L-1); DIVbunbo  = atoi(buff);
           Get_Part_Of_Line(buff,ov->val,0,a-1);   DIVbunshi = atoi(buff);
         }
       }
     else { 
       printf("#ERROR:Can't understand option %s\n",ov->val); 
       exit(1);
     }
   }
 }

/*
 k = 0;
 while (k<argc){
    if ((argv[k][0]=='-') && ((k+1)<argc)){
     L = strlen(argv[k]);
     if (Set_Global_PAR_from_ARGV(argv[k],argv[k+1])==1){ ++k;}
     else{

           if ((L==2)&&(argv[k][1]=='M')) {++k; MODE = argv[k][0];}
      else if ((L==2)&&(argv[k][1]=='Q')) {++k; sprintf(molQ.filename,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='f')&&(argv[k][3]=='l')) {
        ++k; 
        Nlibfile = Make_LibFileList_from_Colon_Separated_String(&HeadLibFile,argv[k]);
        Show_LINENODEs(&HeadLibFile);
      }
      else if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='f')&&(argv[k][3]=='l')&&(argv[k][4]=='q')) {
        ++k; 
        NlibfileQ = Make_LibFileList_from_Colon_Separated_String(&HeadLibFileQ,argv[k]);
        Show_LINENODEs(&HeadLibFileQ);
      }

      else if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='d')&&(argv[k][3]=='m')&&(argv[k][4]=='l')) {++k; sprintf(idirmlib,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='d')&&(argv[k][3]=='l')) {++k; sprintf(idirlib,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='h')&&(argv[k][3]=='l')) {++k; sprintf(iheadlib,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='l')&&(argv[k][3]=='l')) {++k; sprintf(ilistlib,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='d')&&(argv[k][3]=='l')) {++k; sprintf(odirlib,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='s')&&(argv[k][3]=='m')) {++k; sprintf(osimfile,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='f')&&(argv[k][2]=='s')&&(argv[k][3]=='m')) {++k; filetype_osim = argv[k][0];}
      else if ((L==4)&&(argv[k][1]=='w')&&(argv[k][2]=='p')&&(argv[k][3]=='r')) {++k; OMatWithProperty = argv[k][0];}
      else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='s')&&(argv[k][3]=='c')) {++k; sprintf(oschfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='o')&&(argv[k][2]=='c')&&(argv[k][3]=='r')&&(argv[k][4]=='a')&&(argv[k][5]=='s'))
       {++k; sprintf(ocommon_rasfile,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='t')&&(argv[k][2]=='a')&&(argv[k][3]=='o')) {++k; ThreTanimotoOneAtom = atof(argv[k]);}
      else if ((L==4)&&(argv[k][1]=='t')&&(argv[k][2]=='a')&&(argv[k][3]=='m')) {++k; ThreTanimotoMCS   = atof(argv[k]);}
      else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='c')&&(argv[k][3]=='a')) {++k; ScanOpeType = argv[k][0];}
      else if ((L==4)&&(argv[k][1]=='a')&&(argv[k][2]=='c')&&(argv[k][3]=='t')) {++k; ActInactComp = argv[k][0];}
      else if ((L==4)&&(argv[k][1]=='n')&&(argv[k][2]=='l')&&(argv[k][3]=='r')) {++k; Nmollib_nonredun = atoi(argv[k]);}
      else if ((L==3)&&(argv[k][1]=='f')&&(argv[k][2]=='Q')) {++k; molQ.filetype = MoleculeFileType(argv[k]);}
      else if ((L==3)&&(argv[k][1]=='f')&&(argv[k][2]=='L')) {++k; molL.filetype = MoleculeFileType(argv[k]);}
      else if ((L==3)&&(argv[k][1]=='a')&&(argv[k][2]=='Q')) {++k; molQ.atmhet = argv[k][0];}
      else if ((L==3)&&(argv[k][1]=='a')&&(argv[k][2]=='L')) {++k; molL.atmhet = argv[k][0];}
      else if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='d')&&(argv[k][3]=='l')&&(argv[k][4]=='q')) {++k; sprintf(idirlibQ,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='i')&&(argv[k][2]=='d')&&(argv[k][3]=='m')&&(argv[k][4]=='l')&&(argv[k][5]=='Q')) {++k; sprintf(idirmlibQ,"%s",argv[k]);}

      else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='m')&&(argv[k][3]=='s')) {++k; sprintf(omolsfile,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='m')&&(argv[k][3]=='s')&&(argv[k][4]=='d')) {++k; sprintf(omolsdir,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='s')&&(argv[k][3]=='m')&&(argv[k][4]=='i')) {++k; sprintf(osmifile,"%s",argv[k]); }
   
      else if ((L==6)&&(argv[k][1]=='k')&&(argv[k][2]=='p')&&(argv[k][3]=='l')&&(argv[k][4]=='c')&&(argv[k][5]=='n')) 
         {++k; KeepLargestConGraph = argv[k][0]; }
      else if ((L==4)&&(argv[k][1]=='n')&&(argv[k][2]=='m')&&(argv[k][3]=='s'))
         {++k; Nmax_outmol = atoi(argv[k]); }
      else if ((L==6)&&(argv[k][1]=='s')&&(argv[k][2]=='t')&&(argv[k][3]=='e')&&(argv[k][4]=='3')&&(argv[k][5]=='Q'))
         {++k; molQ.SetStereo3D = argv[k][0];}
      else if ((L==6)&&(argv[k][1]=='s')&&(argv[k][2]=='t')&&(argv[k][3]=='e')&&(argv[k][4]=='3')&&(argv[k][5]=='L'))
          {++k; molL.SetStereo3D = argv[k][0];}
      else if ((L==4)&&(argv[k][1]=='d')&&(argv[k][2]=='i')&&(argv[k][3]=='v')){ 
         ++k; a = First_Appear_Pos_in_String(argv[k],'/');
         L = strlen(argv[k]);
         if ((L>1)&&(a>0)){ 
           Get_Part_Of_Line(buff,argv[k],a+1,L-1); DIVbunbo  = atoi(buff);
           Get_Part_Of_Line(buff,argv[k],0,a-1);   DIVbunshi = atoi(buff);
         }
       }
      else { printf("#ERROR:Can't understand option %s\n",argv[k]); exit(1);}
       }
    }
    ++k;
 } 
*/


printf("#COMMAND '%s'\n",PAR.COMMAND);

 srand(PAR.SeedRand);

 
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



 /********************************/ 
 /*** MODE 'S': Library search ***/
 /********************************/
 if (MODE=='S'){ 
   
  /** [1] Read molecule Q **/
  Read_MOLECULE(molQ.filename,&molQ,molQ.filetype,molQ.atmhet,molQ.BondType);
  printf("#QUERY_MOLECULE: Natom %d '%s'\n",molQ.Natom, molQ.name);
  /** [2] Read Library List  **/
  if (idirmlib[0]!='\0'){
    Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlib,&HeadLibFile,"");
  } 

 if (Nlibfile==0){
    if (ilistlib[0]!='\0')     Nmollib = Read_List_of_LIBMOLs(ilistlib,&HeadLibMol);
    else if (idirlib[0]!='\0') Nmollib = Get_List_of_LIBMOLs_Under_Dir(idirlib,&HeadLibMol,iheadlib);
  }
  else{
    ln = &HeadLibFile;
    while (ln->next != NULL){
      ln = ln->next;
      n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMol,&(molL.filetype),ln->num);
      Nmollib += n;
    }
  }

  printf("#Nmollib %d\n",Nmollib); 
  /** [3] Compare proQ vs each molecule in the Library List  **/

  xn = &HeadLibMol;
  Ncalmatch = 0;
  Nhitmol   = 0;
  fpL = NULL;
  fname[0] = fname0[0] = '\0';

  while (xn->next !=NULL){
    xn = xn->next;
    if ((xn->num%100)==0) printf("[%d/%d]\n",xn->num,Nmollib);
    xn->score_for_sort = xn->select_dis = 0.0; 
    xn->tanimoto_mcs = xn->tanimoto_atompair = xn->tanimoto_oneatom = 0.0; 
    xn->Npair = 0;  
    xn->doneMCS = 0;
    sprintf(molL.name,"%s",xn->name);
    sprintf(molL.filename,"%s/%s",idirlib,xn->name);
    /* printf("[%d] '%s'\n",xn->num,xn->name); fflush(stdout); */

    if (Nlibfile==0){
      ok = Read_MOLECULE(molL.filename,&molL,molL.filetype,molL.atmhet,molL.BondType);
      xn->num_file =  xn->offset_libfile = 0;
    }
    else{
      String_from_LINENODE_num(fname,&HeadLibFile,xn->num_file);
      /* printf("fname '%s' fname0 '%s'\n",fname,fname0); */
      if (strcmp(fname,fname0)!=0){
        if (fpL!= NULL) fclose(fpL);
        fpL = fopen(fname,"r");
       }
      xn->offset_libfile = ftell(fpL);
      if (xn->offset_libfile<0){
        printf("WOOPS %s offset %ld\n",xn->name,xn->offset_libfile); 
        ok = 0;
      }
      else{
        ok = Read_One_Molecule_from_one_Library_file(fpL,&molL,molL.filetype,molL.atmhet,'L');
        sprintf(fname0,"%s",fname);
      }
    }
    
    if (ok>0){ 
     /* printf("#molL.name '%s' Nheavyatom %d\n",molL.name,molL.Nheavyatom); fflush(stdout); */
      if (xn->name[0]=='\0') sprintf(xn->name,"%s",molL.name);
      xn->Nheavyatom = molL.Nheavyatom; 
      xn->molform = (char *)malloc(sizeof(char)*(strlen(molL.molform)+1));
      sprintf(xn->molform,"%s",molL.molform);
      doMCS = 0;
      if ((PAR.ConnectGraphType != 'S') && (PAR.ConnectGraphType != 'I')){
        xn->tanimoto_oneatom =  Max_Tanimoto_from_Natomtype(&molQ,&molL);
        if (xn->tanimoto_oneatom >=ThreTanimotoOneAtom) doMCS = 1;
      }
      else if ((PAR.ConnectGraphType=='S')
              &&(Substructure_condition_oneatom_descriptors(molQ.oneatom_descriptor,molL.oneatom_descriptor)==1)) doMCS = 1;
      else if ((PAR.ConnectGraphType=='I')
              &&(Isomorphic_condition_oneatom_descriptors(molQ.oneatom_descriptor,molL.oneatom_descriptor)==1)) doMCS = 1; 
  
      if ((molL.Nheavyatom>0)&&(doMCS==1 )){
        xn->doneMCS = 1;
        Match_Two_Molecules(&molQ,&molL,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
        Ncalmatch += 1;
        if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E')){
           Copy_MATCH_Info_to_LIBMOL(xn,optMlist.next,&molQ,&molL);
           if (xn->tanimoto_mcs >= ThreTanimotoMCS){
             Nhitmol += 1;
             Increment_Nmatch_for_Matched_Atoms(optMlist.next,&molQ,&molL);
           }
        }
        if (PAR.OutCalProcess=='T')
           printf("[%d] %s Nheavyatom %d-%d Npair %d tanimoto %f select_dis %f\n",
           xn->num,xn->name,molQ.Nheavyatom,xn->Nheavyatom,xn->Npair,xn->tanimoto_mcs,xn->select_dis);
        Free_MATCHlist(&optMlist);
       } /* Nheavyatom > 0 */

     Free_MOLECULE(&molL);
    } /* ok */

  } /* xn */
  
  printf("#Ncalmacth %d\n",Ncalmatch);
  if (fpL != NULL) fclose(fpL);


   /*** [4] Output Results ***/
   Merge_Sort_Double_Linked_List_LIBMOL(&HeadLibMol,'D');
   if (oschfile[0] != '\0'){
    Write_Search_Result(oschfile,&molQ,idirlib,"",&HeadLibMol,Nmollib,Nmollib,
           ThreTanimotoMCS,ThreTanimotoOneAtom,-1.0,&HeadLibFile);
   }
 

   if (omolsfile[0] != '\0'){
      Write_Similar_Molecules_To_One_File(omolsfile,"",idirlib,&HeadLibMol,
         ThreTanimotoMCS,ThreTanimotoOneAtom,-1.0,molL.filetype,Nlibfile,&HeadLibFile,Nmax_outmol);
    }

   if (omolsdir[0] != '\0'){
     Write_Similar_Molecules_Under_Dir(omolsdir,"",idirlib,&HeadLibMol,
        ThreTanimotoMCS,ThreTanimotoOneAtom,-1.0,molL.filetype,Nlibfile,&HeadLibFile,Nmax_outmol);
  }

  if (ocommon_rasfile[0]!='\0'){
    Write_Nmatch_rasmol_Script(ocommon_rasfile,&molQ,Nhitmol,"");
  }
}


 /*********************************************/
 /*** MODE 'A': all-vs-all search           ***/
 /*** (keep all the molecule in the memory) ***/ 
 /*********************************************/
 
 if (MODE=='A'){ 
  comment_str[0] = '\0';

  /** [1] Read List of library molecules **/
  if (idirmlib[0]!='\0'){
    Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlib,&HeadLibFile,"");
  }
 
  if (Nlibfile ==0){
         if (ilistlib[0]!='\0') Nmollib = Read_List_of_LIBMOLs(ilistlib,&HeadLibMol);
    else if (idirlib[0] !='\0') Nmollib = Get_List_of_LIBMOLs_Under_Dir(idirlib,&HeadLibMol,iheadlib);
  }
  else{
    ln = &HeadLibFile;
    while (ln->next != NULL){
      ln = ln->next;
      n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMol,&(molL.filetype),ln->num);
      Nmollib += n;
    }
  }
  printf("#Nmollib %d\n",Nmollib); 

  Cal_Nsta_Nend_from_Bunshi_Bunbo_for_all_vs_all(Nmollib,DIVbunshi,DIVbunbo,&Nsta,&Nend);

  sprintf(comment_str,"division %d/%d Nmollib %d Nsta %d Nend %d",DIVbunshi, DIVbunbo, Nmollib,Nsta, Nend);
  printf("#%s\n",comment_str);

 /*
  an = &HeadLibMol;
  while (an->next !=NULL){
    an = an->next;
    printf("'%s'\n",an->name);
  }
 */

  /** [2] Read all molecules **/
  an = &HeadLibMol;
  fname[0] = fname0[0] = '\0';
  fpL = NULL;
 
  while (an->next !=NULL){
    an = an->next;
    an->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));  
    an->mol->filetype = molL.filetype;
    an->mol->atmhet   = molL.atmhet;
    an->mol->BondType = molL.BondType;
    an->mol->chain    = molL.chain;
    an->mol->resi[0]  = '\0';
    an->mol->name[0]  = '\0';
    an->mark = 0;
    if (an->name[0]!='\0') sprintf(an->mol->name,"%s",an->name);
    sprintf(an->mol->filename,"%s/%s",idirlib,an->name);
    if (Nlibfile==0){
       ok = Read_MOLECULE(an->mol->filename,an->mol,an->mol->filetype,an->mol->atmhet,an->mol->BondType);
       an->num_file =  an->offset_libfile = 0;
    } 
    else{
      String_from_LINENODE_num(fname,&HeadLibFile,an->num_file);
      if (strcmp(fname,fname0)!=0){
           if (fpL!= NULL) fclose(fpL);
           fpL = fopen(fname,"r"); 
         }
       an->offset_libfile = ftell(fpL);
       ok = Read_One_Molecule_from_one_Library_file(fpL,an->mol,molL.filetype,molL.atmhet,'L');
       sprintf(fname0,"%s",fname);
    }
    if (an->mol->name[0]!='\0'){ 
      free(an->name); an->name = (char *)malloc(sizeof(char)*(strlen(an->mol->name)+1));
      sprintf(an->name,"%s",an->mol->name);
    }
     if ((Nmollib<100)||((an->num%100)==0)) 
       printf(">read [%d/%d] '%s' '%s' Natom %d\n",an->num,Nmollib,an->mol->filename,an->mol->name,an->mol->Natom); 
   }
  
  if (fpL != NULL) fclose(fpL);

  /** [3] all-vs-all  comparison **/
  if (((filetype_osim=='D')||(filetype_osim=='S'))&&(osimfile[0]!='\0')) Malloc_FLOAT2DMAP(&Dmap,Nmollib); 
  fpo = NULL; 
  if (((filetype_osim=='L')||(filetype_osim=='N'))&&(osimfile[0]!='\0')){ 
    printf("#write_similarity_list() -->'%s'\n",osimfile);
    fpo = fopen(osimfile,"w");
    Write_Similarities_in_List_format_Header(fpo,&HeadLibMol,idirlib);
  }
 
  an = &HeadLibMol;
  bn = NULL;
  while (an->next !=NULL){
    an = an->next;  
    f = (float)an->num/(float)Nmollib;
    if ((Nsta <= an->num) && (an->num <Nend)){
      if (((filetype_osim=='D')||(filetype_osim=='S'))&&(osimfile[0]!='\0')) Dmap.map[an->num][an->num] = 1.0;
      printf("[%3d/%d](%5.1f%%) %s\n",an->num,Nmollib,f*(2.0-f)*100.0,an->mol->core_molname);
      an->mark = 1;
      bn = an;
      while (bn->next !=NULL){
         bn = bn->next;
         if (Max_Tanimoto_from_Natomtype(an->mol,bn->mol) < ThreTanimotoOneAtom) { ok = 0; tanimoto = 0.0;}
         if ((ActInactComp=='T') && (strcmp(an->class,"inactive")==0) && (strcmp(bn->class,"inactive")==0)){
          ok = 0; tanimoto = -1.0;
          /* printf("#omit calculation %s %s vs %s %s\n",an->name,an->class,bn->name,bn->class); */
          } 

         doMCS = 0;
         tanimoto = 0.0;

         if ((PAR.ConnectGraphType != 'S') && (PAR.ConnectGraphType != 'I')){
           an->tanimoto_oneatom =  Max_Tanimoto_from_Natomtype(an->mol,bn->mol);
          if (an->tanimoto_oneatom >=ThreTanimotoOneAtom) doMCS = 1;
         }
         else if ((PAR.ConnectGraphType=='S')
                  &&(Substructure_condition_oneatom_descriptors(an->mol->oneatom_descriptor,bn->mol->oneatom_descriptor)==1)) doMCS = 1;
         else if ((PAR.ConnectGraphType=='I')
                  &&(Isomorphic_condition_oneatom_descriptors(an->mol->oneatom_descriptor,bn->mol->oneatom_descriptor)==1)) doMCS = 1;

         if (doMCS==1){ 
           Match_Two_Molecules(an->mol,bn->mol,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
           Increment_Nmatch_for_Matched_Atoms(optMlist.next,an->mol,bn->mol);
           tanimoto = Tanimoto_Coefficient(optMlist.next->Npair,an->mol,bn->mol);
           if ((tanimoto >= ThreTanimotoMCS)&&(osimfile[0]!='\0')){ 
             if ((filetype_osim=='D')||(filetype_osim=='S')) Dmap.map[an->num][bn->num] = Dmap.map[bn->num][an->num] = tanimoto;
             if (filetype_osim=='L') fprintf(fpo,"%d %d %f\n",an->num,bn->num,tanimoto);
             if (filetype_osim=='N') fprintf(fpo,"%d %s %d %s %f\n",an->num,an->mol->core_molname,bn->num,bn->mol->core_molname,tanimoto);
           } 
         }

         Free_MATCHlist(&optMlist);
       } /* bn */ 
     } 
   } /* an */

   if (osimfile[0]!='\0'){
     if (filetype_osim=='D')
       Write_Distance_FLOAT2DMAP_in_Phylip_format(osimfile,&Dmap,&HeadLibMol,'D',OMatWithProperty);
     if (filetype_osim=='S')
       Write_Distance_FLOAT2DMAP_in_Phylip_format(osimfile,&Dmap,&HeadLibMol,'S',OMatWithProperty);
     if ((filetype_osim=='L')||(filetype_osim=='N')){
       printf("#write_similarity_list() -->'%s'\n",osimfile);
       fclose(fpo);
     }
   }
 
   /* output the most common structure */ 
   if (ocommon_rasfile[0]!='\0'){
     an = &HeadLibMol;
     n = 0;
     while (an->next !=NULL){
      an = an->next;
      if (an->mol->Nmatch>n) {bn = an; n = an->mol->Nmatch;}
      printf("'%s' '%-20s' Nmatch %4d\n",an->mol->filename,an->mol->name,an->mol->Nmatch); 
     }
     sprintf(buff,"#most common molecule '%s' '%s' Nmatch %d",bn->mol->filename,bn->mol->name,bn->mol->Nmatch); 
     printf("%s\n",buff);
     Write_Nmatch_rasmol_Script(ocommon_rasfile,bn->mol,Nmollib,buff);
   }
 
  }


 /*************************************************/
 /*** MODE 'G': comparison between two groups    **/
 /**  (comparison bwn library and query library) **/
 /*************************************************/
 if (MODE=='G'){ 
  
  /** [1] Read List of library molecules **/
  if (idirmlib[0]!='\0'){
    Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlib,&HeadLibFile,"");
  }
 
  if (Nlibfile ==0){
         if (ilistlib[0]!='\0')  Nmollib = Read_List_of_LIBMOLs(ilistlib,&HeadLibMol);
    else if (idirlib[0] !='\0')  Nmollib = Get_List_of_LIBMOLs_Under_Dir(idirlib,&HeadLibMol,iheadlib);
  }
  else{
    ln = &HeadLibFile;
    while (ln->next != NULL){
      ln = ln->next;
      n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMol,&(molL.filetype),ln->num);
      Nmollib += n;
    }
  }
  printf("#Nmollib %d\n",Nmollib); 


  /** [2] Read List of Query library molecules **/
  printf("NlibfileQ %d\n",NlibfileQ);
  if (idirmlibQ[0]!='\0'){
    Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlibQ,&HeadLibFileQ,"");
  }
 
  if (NlibfileQ ==0){
         if (ilistlibQ[0]!='\0')  NmollibQ = Read_List_of_LIBMOLs(ilistlibQ,&HeadLibMolQ);
    else if (idirlibQ[0] !='\0')  NmollibQ = Get_List_of_LIBMOLs_Under_Dir(idirlibQ,&HeadLibMolQ,iheadlibQ);
  }
  else{
    printf("read file\n");
    ln = &HeadLibFileQ;
    while (ln->next != NULL){
      ln = ln->next;
      n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMolQ,&(molQ.filetype),ln->num);
      NmollibQ += n;
    }
  }
  printf("#NmollibQ %d\n",NmollibQ); 
 
  /*
  qn  = &HeadLibMolQ;
  while (qn->next != NULL){
    qn = qn->next;
    printf("[qn] num %d name '%s'\n",qn->num,qn->name);
  }  
  */
 
  /** [3] Read all molecules in the library **/
  an = &HeadLibMol;
  fname[0] = fname0[0] = '\0';
  fpL = NULL;
 
  while (an->next !=NULL){
    an = an->next;
    an->mol = (struct MOLECULE*)malloc(sizeof(struct MOLECULE));  
    if (an->name[0]!='\0') sprintf(an->mol->name,"%s",an->name);
    sprintf(an->mol->filename,"%s/%s",idirlib,an->name);
    an->mol->filetype = molL.filetype;
    an->mol->atmhet   = molL.atmhet;
    an->mol->BondType = molL.BondType;
    an->mol->chain    = molL.chain;
    an->mol->resi[0]  = '\0';
    an->mark = 0;
    if (Nlibfile==0){
       ok = Read_MOLECULE(an->mol->filename,an->mol,an->mol->filetype,an->mol->atmhet,an->mol->BondType);
       an->num_file =  an->offset_libfile = 0;
    } 
    else{
      String_from_LINENODE_num(fname,&HeadLibFile,an->num_file);
      if (strcmp(fname,fname0)!=0){
           if (fpL!= NULL) fclose(fpL);
           fpL = fopen(fname,"r"); 
         }
       an->offset_libfile = ftell(fpL);
       ok = Read_One_Molecule_from_one_Library_file(fpL,an->mol,molL.filetype,molL.atmhet,'L');
       sprintf(fname0,"%s",fname);
    }
    if (an->mol->name[0]!='\0'){ 
      free(an->name); an->name = (char *)malloc(sizeof(char)*(strlen(an->mol->name)+1));
      sprintf(an->name,"%s",an->mol->name);
    }
    if ((an->num%100)==0)
       printf(">read [%d/%d] '%s' '%s' Natom %d\n",an->num,Nmollib,an->mol->filename,an->mol->name,an->mol->Natom); 
   }
  
  if (fpL != NULL) fclose(fpL);

  /** [4] Read query molecule one by one, and compare all the library molecule **/
  Nsta = (NmollibQ*DIVbunshi)/DIVbunbo;
  Nend = (NmollibQ*(DIVbunshi+1))/DIVbunbo;
  printf("#NlibmolQ %d div %d/%d Nsta %d Nend %d\n",NmollibQ,DIVbunshi,DIVbunbo,Nsta,Nend);

  fname[0] = fname0[0] = '\0';
  fpL = fpo = NULL;
  Nsimilar_pair = 0; 
  if (((filetype_osim=='L')||(filetype_osim=='N'))&&(osimfile[0]!='\0')){ 
    printf("#write_similarity_list() -->'%s'\n",osimfile);
    fpo = fopen(osimfile,"w");
    Write_Two_Groups_Similarities_in_List_format_Header(fpo,&HeadLibMolQ,&HeadLibMol);
  }

 
  qn = &HeadLibMolQ;
  while (qn->next !=NULL){
    qn = qn->next;
    molQ.Natom = 0;
    sprintf(molQ.filename,"%s/%s",idirlibQ,qn->name);
    if (Nlibfile==0){
       ok = Read_MOLECULE(molQ.filename,&molQ,molQ.filetype,molQ.atmhet,molQ.BondType);
       qn->num_file =  qn->offset_libfile = 0;
    } 
    else{
      String_from_LINENODE_num(fname,&HeadLibFileQ,qn->num_file);
      if (strcmp(fname,fname0)!=0){
           if (fpL!= NULL) fclose(fpL);
           fpL = fopen(fname,"r"); 
         }
       if (fpL!=NULL){
         qn->offset_libfile = ftell(fpL);
         ok = Read_One_Molecule_from_one_Library_file(fpL,&molQ,molQ.filetype,molQ.atmhet,'L');
         sprintf(fname0,"%s",fname);
       }
    }

    if ((molQ.Natom > 0) && (Nsta<= qn->num) && (qn->num < Nend)){
      if (molQ.name[0] != '\0'){
        free(qn->name);
        qn->name = (char *)malloc(sizeof(char)*(strlen(molQ.name)+1));
        sprintf(qn->name,"%s",molQ.name);
      }
      printf(">query [%d/%d] name '%s' filename '%s' Natom %d\n",qn->num,NmollibQ,molQ.name,molQ.filename,molQ.Natom); 

      an = &HeadLibMol;
      while (an->next != NULL){
        an = an->next;

        doMCS = 0;
        if ((PAR.ConnectGraphType != 'S') && (PAR.ConnectGraphType != 'I')){
          an->tanimoto_oneatom =  Max_Tanimoto_from_Natomtype(&molQ,an->mol);
          if (an->tanimoto_oneatom >=ThreTanimotoOneAtom) doMCS = 1;
        }
        else if ((PAR.ConnectGraphType=='S')
                 &&(Substructure_condition_oneatom_descriptors(molQ.oneatom_descriptor,an->mol->oneatom_descriptor)==1)) doMCS = 1;
        else if ((PAR.ConnectGraphType=='I')
                 &&(Isomorphic_condition_oneatom_descriptors(molQ.oneatom_descriptor,an->mol->oneatom_descriptor)==1)) doMCS = 1;

        if (doMCS==1) {
          Match_Two_Molecules(&molQ,an->mol,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
          an->tanimoto_mcs = Tanimoto_Coefficient(optMlist.next->Npair,&molQ,an->mol);
          /* fprintf(fpo,"%s %s %f\n",qn->name,an->name,tanimoto); */
          if ((ThreTanimotoMCS<=0.0)||(an->tanimoto_mcs>=ThreTanimotoMCS)){
            if (osimfile[0]!='\0'){ 
              if (filetype_osim=='L') fprintf(fpo,"%d %d %f\n",qn->num,an->num,an->tanimoto_mcs); 
              if (filetype_osim=='N') fprintf(fpo,"%d %s %d %s %f\n",qn->num,qn->name,an->num,an->name,an->tanimoto_mcs); 
              Nsimilar_pair += 1;
            }
          }
          Free_MATCHlist(&optMlist);
        }
      }
      Free_MOLECULE(&molQ);
    }
  }
  printf("#Nsimilar_pair %d\n",Nsimilar_pair); 
  if (fpL != NULL) fclose(fpL);
  if (osimfile[0]!='\0'){ 
     printf("#write_similarity_list() -->'%s'\n",osimfile);
     fclose(fpo); 
   }
   
 }



 /*************************************************************/
 /*** MODE 'R': randomly selecting non-redundant molecules  ***/
 /*************************************************************/
 
 if (MODE=='R'){ 


  /** [1] Read Library List  **/
   if (idirmlib[0]!='\0'){
     Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlib,&HeadLibFile,"");
   }

 
   if (Nlibfile ==0){
          if (ilistlib[0]!='\0')  Nmollib = Read_List_of_LIBMOLs(ilistlib,&HeadLibMol);
     else if (idirlib[0] !='\0')  Nmollib = Get_List_of_LIBMOLs_Under_Dir(idirlib,&HeadLibMol,iheadlib);
   }
  else{
    ln = &HeadLibFile;
    while (ln->next != NULL){
      ln = ln->next;
      n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMol,&(molL.filetype),ln->num);
      Nmollib += n;
    }
  }

  printf("#Nmollib %d\n",Nmollib); 


   HeadLibMolSelected.next = NULL;
   HeadLibMolSelected.num  = -1;
   xn = &HeadLibMol;
   while (xn->next != NULL){
     xn = xn->next;
     xn->subnum = 0;
     if ((xn->num%100)==0) printf("[%d] %s\n",xn->num,xn->name);
   }
  printf("#Nmollib %d\n",Nmollib); 
   
   /*** [2] Select LIBMOL randomly ***/

   fpL = NULL;
   fname[0] = '\0';

   for (n=0;n<Nmollib_nonredun;++n){
     ok = 0;  
     while (ok==0){
       r = rand() % Nmollib;
       xn = Get_LIBMOL_by_number(&HeadLibMol,r);
       if (xn->subnum==0){
         ok = 1;
         sprintf(molL.filename,"%s/%s",idirlib,xn->name);
         Read_Library_MOLECULE_Non_Sequential(xn,&molL,Nlibfile,&HeadLibFile,fpL,fname,fname0);
         if (ThreTanimotoMCS>0.0){
           an = &HeadLibMolSelected;
           while ((an->next!=NULL)&&(ok==1)){
             an = an->next;
             Match_Two_Molecules(&molL,an->mol,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,"");
             tanimoto = 0.0;
             if (optMlist.next != NULL)
               tanimoto = Tanimoto_Coefficient(optMlist.next->Npair,&molL,an->mol);
             if (tanimoto>ThreTanimotoMCS){ 
               printf("#%s vs %s tanimoto %f\n",xn->name,an->name,tanimoto);
               ok = 0;
               xn->subnum = 2; /* already chosen, but similar to previously selected molecule */ 
             }
             Free_MATCHlist(&optMlist);
           }
         }
         if (ok==1){
           xn->subnum = 1;
           bn = Add_LIBMOL_to_LIBMOLlist(xn,&HeadLibMolSelected);
           bn->mol = (struct MOLECULE *)malloc(sizeof(struct MOLECULE)); 
           Malloc_MOLECULE(bn->mol,molL.Natom,molL.Nbond);
           Copy_MOLECULE(bn->mol,&molL); 
           Free_MOLECULE(&molL); 
           printf("[%d/%d/%d] select %d %s\n",n,Nmollib_nonredun,Nmollib,xn->num,xn->name);
         }
       }
     }
   }

   /*** [3] Output ***/
   if (omolsfile[0] != '\0'){
      Write_Similar_Molecules_To_One_File(omolsfile,"",idirlib,&HeadLibMolSelected,
         -1.0,-1.0,-1.0,molL.filetype,Nlibfile,&HeadLibFile,-1);
    }

   if (omolsdir[0] != '\0'){
     Write_Similar_Molecules_Under_Dir(omolsdir,"",idirlib,&HeadLibMolSelected,
        -1.0,-1.0,-1.0,molL.filetype,Nlibfile,&HeadLibFile,-1);
  }

 }




 /****************************************************************/
 /*** MODE 'C': do converting operation  for molecular library ***/
 /****************************************************************/
 if (MODE=='C'){ 

   /**  Read Library List  **/
  if (idirmlib[0]!='\0'){
    Nlibfile =  Get_List_of_LibFiles_Under_Dir(idirmlib,&HeadLibFile,"");
  }

  if (Nlibfile ==0){
         if (ilistlib[0]!='\0') Nmollib = Read_List_of_LIBMOLs(ilistlib,&HeadLibMol);
    else if (idirlib[0]!='\0')  Nmollib = Get_List_of_LIBMOLs_Under_Dir(idirlib,&HeadLibMol,iheadlib);
  }
  else{
    ln = &HeadLibFile;
    while (ln->next != NULL){
      ln = ln->next;
      n = Get_List_of_LIBMOLs_from_one_library_file(ln->line,&HeadLibMol,&(molL.filetype),ln->num);
      Nmollib += n;
    }
  }

  printf("#Nmollib %d\n",Nmollib); 
  fpo = NULL;

  if (osmifile[0] != '\0'){
     if (strcmp(osmifile,"stdout")==0) fpo = stdout; 
     else{ 
       fpo = fopen(osmifile,"w");
       /* fprintf(fpo,"#COMMAND '%s'\n",PAR.COMMAND);*/
     }
  }

  /** Read all molecules **/
   an = &HeadLibMol;
   fpL = NULL;
   fname[0] = fname0[0] = '\0';

   while (an->next !=NULL){
     an = an->next;  
     /* printf(">%s\n",molL.name); */
     sprintf(molL.filename,"%s/%s",idirlib,an->name);
     ok = 0;

     if (Nlibfile==0){
       ok = Read_MOLECULE(molL.filename,&molL,molL.filetype,molL.atmhet,molL.BondType);
       if (molL.name[0]=='\0') sprintf(molL.name,"%s",an->name);

     } 
     else{
      String_from_LINENODE_num(fname,&HeadLibFile,an->num_file);
      if (strcmp(fname,fname0)!=0){
           if (fpL!= NULL) fclose(fpL);
           fpL = fopen(fname,"r"); 
         }

       an->offset_libfile = ftell(fpL);
       ok = Read_One_Molecule_from_one_Library_file(fpL,&molL,molL.filetype,molL.atmhet,'L');
       sprintf(fname0,"%s",fname);
     }
     if (ok>0){
       /* printf(">'%s' Natom %d\n",molL.name,molL.Natom); fflush(stdout); */
       if (osmifile[0] != '\0'){
         make_LUNITs_from_Molecule_by_DFS(&HeadLunit,&molL);
         string = (char *)malloc(sizeof(char)*(length_string_from_LUNITs(&molL,&HeadLunit)));
         make_string_from_LUNITs(&molL,string,&HeadLunit,KeepLargestConGraph); 
         if ((an->num %100)==0){
           printf("#[%d/%d] '%s' '%s'\n",an->num,Nmollib,string,molL.name);
         }
         fprintf(fpo,"%s %s\n",string,molL.name);
         Free_LUNITs(&HeadLunit);
         free(string);
       } 

       if (odirlib[0]!='\0'){
         n = Delete_Disconnected_Parts_Of_Molecule(&molL);
         if (n>0) printf("#DELETE ISOLATED %d ATOMS FOR '%s'\n",n,an->name);
         sprintf(ofname,"%s/%s",odirlib,an->name);
         Write_SDF_Molecule(ofname,&molL,"");
         Read_and_Append_SDF_Annotations(molL.filename,ofname);
       }

       if (PAR.GenEquivAtomPermu=='T'){
         PAR.N = Number_of_EquivMap_for_Heavyatoms(&molL);
         Make_All_Permutation_of_Equivalent_Heavyatoms(&molL); 
         printf("%-30s Nheavyatom %3d Nequivpair %3d %6.2f%% Npermu %d\n",
           molL.filename,molL.Nheavyatom,PAR.N,100.0*2.0*PAR.N/(molL.Nheavyatom)/(molL.Nheavyatom-1),molL.Npermu);
         permu = &(molL.permuhead);
         while (permu->next != NULL){
           permu = permu->next;
           show_atom_permutation(&molL,permu->index_heavy);
         }
       }


       if (ScanOpeType=='D'){ 
         if (odirlib[0]!='\0'){
           sprintf(ofname,"%s/%s",odirlib,an->name);
           Write_PDB_Molecule(ofname,&molL,'w','A',1,-1,'D','-',""); 
          }
       }

       if (ScanOpeType=='R'){ 
         if (odirlib[0]!='\0'){
           sprintf(ofname,"%s/%s",odirlib,an->name);
           Write_DAMLatom_and_RING_in_PDB(ofname,&molL,'w','A',1,""); 
         }
       }

     Free_MOLECULE(&molL);
    } /* ok > 0 */
   }

    if ((osmifile[0]!= '\0') && (fpo != stdout)){
      printf("write_SMILES() --> '%s'\n",osmifile);
      fclose(fpo);
    }
 }
 
 return(1);
} /* end of main() */ 









/**********

  FUNCTIONS

***********/



int Read_Library_MOLECULE_Non_Sequential(xn,mol,Nlibfile,HeadLibFile,fpL,fname,fname0)
  struct LIBMOL *xn;
  struct MOLECULE *mol;
  int Nlibfile;
  struct LINENODE *HeadLibFile;
  FILE *fpL;
  char *fname,*fname0;
{
  int ok;
  /* 
  printf("#Read_Library_MOLECULE_Non_Sequential(xn:%d offset %d mol,Nlibfile:%d HeadLibFile,fpL,fname:'%s' fname0:'%s')\n",xn->num,xn->offset_libfile,Nlibfile,fname,fname0);
  */ 
  ok = 0;
  if (Nlibfile==0){
     ok = Read_MOLECULE(mol->filename,mol,mol->filetype,mol->atmhet,mol->BondType);
     xn->num_file =  xn->offset_libfile = 0;
  }
  else{
     String_from_LINENODE_num(fname,HeadLibFile,xn->num_file);
     if ((fpL==NULL)||(strcmp(fname,fname0)!=0)){
        if (fpL!= NULL) fclose(fpL);
        fpL = fopen(fname,"r");
        if (fpL == NULL){ printf("#ERROR:Can't open '%s'\n",fname); exit(1);}
     }
    if (fpL != NULL){
      fseek(fpL,xn->offset_libfile,SEEK_SET);
      ok = Read_One_Molecule_from_one_Library_file(fpL,mol,mol->filetype,mol->atmhet,'L');
      xn->name = (char *)malloc(sizeof(char)*(strlen(mol->name)+1));
      sprintf(xn->name,"%s",mol->name);
    }
    sprintf(fname0,"%s",fname);
  }
  return(ok);
} 
