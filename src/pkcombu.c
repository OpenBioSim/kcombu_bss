/*
 
 <pkcombu.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

  for obtaining maximum common substructure
  by build-up procedures only for pairwise comparison.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "ioSDF.h"
#include "ioPDB.h"
#include "ioMOL2.h"
#include "ioSMILES.h"
#include "gen2D.h"
#include "ioLib.h"
#include "match.h"
#include "io_match.h"
#include "buildup.h"
#include "clique.h"
#include "c_clique.h"
#include "PCAfit.h"
#include "MrgSrtMATCH.h"
#include "qRMS.h"
#include "molprop.h"
#include "molring.h"
#include "molpermu.h"
#include "moltopodis.h"
#include "ringblock.h"
#include "RingDesc.h"
#include "outPS.h"
#include "options.h"
#include "ioLINE.h"
#include "stereo_check.h"
#include "energies.h"

struct PARAMETERS PAR;

/** FUNCTIONS (LOCAL) **/

int main(argc,argv)
 int argc;
 char **argv;
{
  struct MOLECULE molA,molB;
  int k;
  char DetailOptHelp;
  struct MATCH optMlist,*om;
  char orasfile[MAX_FILENAME],orasfileA[MAX_FILENAME],orasfileB[MAX_FILENAME];
  char opdbfileAB[MAX_FILENAME];
  char omcs_sdffile[MAX_FILENAME],oamchfile[MAX_FILENAME],oAmchfile[MAX_FILENAME], iamchfile[MAX_FILENAME];
  char opsfile[MAX_FILENAME],osmilesfileA[MAX_FILENAME],osmilesfileB[MAX_FILENAME]; 
  char Superimp_molA, Superimp_molB;
  char Paste2D_molA,Paste3D_molA, Gen2D_molA, Gen2D_molB, AnnotationSmiles_molA;
  char KeepLargestConGraph,KeyValueOut,RerankMatchByRMSD,CompEleCharge;
  char Rescale3D_molB; 
  char ofname[MAX_FILENAME],newname_molA[MAX_FILENAME],*string;
  char buff_string[MAX_FILENAME];
  float tanimoto;
  int    Natmmatch;
  double gA[3],gB[3],Rmat[3][3],raxisang[4],tvec[3];
  float RMSD;
  int okA, okB;
  int RankMatchOutput; 
  struct LINENODE HeadComment;
  struct LUNIT HeadLunit;
  char CheckStereoParity;
  char CheckSelfClash_A, CheckSelfClash_B;
  struct OPTION_VALUE *ov;

  
  HeadComment.next = NULL;
  /****** SET INITIAL PARAMETERS ***********/
  Set_Default_Global_PAR(); 
  PAR.PROGRAM_TYPE = 'p'; 
  Superimp_molA =  Superimp_molB = 'F'; 
  Paste2D_molA = Paste3D_molA = 'F'; 
  Initialize_MOLECULE(&molA);
  Initialize_MOLECULE_string(&molA);
  Initialize_MOLECULE(&molB);
  Initialize_MOLECULE_string(&molB);
  molA.filename[0] = molB.filename[0] = '\0';
  molA.filetype = molB.filetype = '-'; 
  molA.atmhet = molB.atmhet = 'B';
  molA.BondType = molB.BondType = 'B'; 
  molA.opdbfile[0]  = molB.opdbfile[0] = '\0';  
  molA.osdffile[0]  = molB.osdffile[0] = '\0';  
  molA.omol2file[0] = molB.omol2file[0] = '\0';  
  molA.chain = molB.chain = '-'; 
  molA.resi[0] = molB.resi[0] = '\0'; 
  osmilesfileA[0] = osmilesfileB[0] = '\0'; 
  omcs_sdffile[0] =  oamchfile[0] = iamchfile[0] = oAmchfile[0] =  '\0';
  opdbfileAB[0]= orasfile[0] = orasfileA[0] = orasfileB[0] = opsfile[0] = '\0'; 
   
  PAR.max_ring_size  = 30;
  PAR.max_block_size = 40;
  PAR.max_ring_descriptor  = PAR.max_ring_size + PAR.max_block_size  - 3;
  RankMatchOutput = 1;

  PAR.ps_label   = 'M';
  PAR.ps_hydro   = 'F';
  PAR.ps_circle  = 'M';
  PAR.ps_planeB  = 'T';
  PAR.ps_supAonB = 'T';
  PAR.ps_color   = 'C';
  PAR.ps_circ_outline   = 'F';
  PAR.ps_care_bondtype  = 'T';

  DetailOptHelp = '-';

  Gen2D_molA = Gen2D_molB = 'F';
  newname_molA[0] = '\0';
  AnnotationSmiles_molA = 'F';
  KeepLargestConGraph = 'F';
  Rescale3D_molB = 'F';

  molA.SetStereo3D = 'F'; molB.SetStereo3D = 'F';

  KeyValueOut = 'F';
  RerankMatchByRMSD = 'F';
  CompEleCharge = 'F';

  CheckStereoParity = 'F';
  CheckSelfClash_A = CheckSelfClash_B = 'F';

  om = NULL; 


  if((argc>1) && ((strcmp(argv[1],"-h")==0) || (strcmp(argv[1],"-help")==0)) ){DetailOptHelp = 'D';}
  if((argc>1) && ((strncmp(argv[1],"-hm",3)==0) || (strncmp(argv[1],"-HM",3)==0)) ){ DetailOptHelp = 'M';}

  if (argc<2){
    printf("pkcombu <options>\n");
    printf("  %s\n",KCOMBU_STANDS_FOR);
    printf(" 'pkcombu' is for pairwise comparison.\n"); 
    printf("  %s LastModified:%s \n",CODED_BY,LAST_MOD_DATE);
    printf("** Simple Usage for 'pkcombu' **\n");
    printf(">for default pairwise comparison.\n");
    printf("  $pkcombu -A [molecule_fileA] -B [molecule_fileB] -oam [atom-matching file]\n");
    printf(">for Topologically constrained D-MCS (TD-MCS) with max_diff_of_dis = 1\n");
    printf("  $pkcombu  -A [molecule_fileA] -B [molecule_fileB] -con T -mtd 1\n");
    printf(">for exact C-MCS\n");
    printf("  $pkcombu -A [molecule_fileA] -B [molecule_fileB] -con C -alg X\n");
    printf(">for checking chirality matches bwn 3D molA and 2D molB(*.sdf) \n");
    printf("  $pkcombu -A [3D-molecule_fileA] -B [2D-molecule_fileB.sdf] -ste3A T -chkste C\n");
    printf(">for showing options in detail\n");
    printf("  $pkcombu -h\n");
    printf(">for showing options for MCS\n");
    printf("  $pkcombu -hmcs\n");
   }
  if (DetailOptHelp=='D'){
    printf("pkcombu <options>\n");
    printf("  for obtaining maximum common substructures by build-up method.\n");
    printf("  %s\n",KCOMBU_STANDS_FOR); 
    printf("   coded by T.Kawabata. LastModified:%s\n",LAST_MOD_DATE);
    printf("<input options for 'pkcombu'>\n"); 
    printf(" -A   : molecule A (molA)(*.sdf|*.mol2|*.pdb|*.kcf|*.smi|*.cif)[%s]\n",molA.filename);
    printf(" -B   : molecule B (molB)(*.sdf|*.mol2|*.pdb|*.kcf|*.smi|*.cif)[%s]\n",molB.filename);
    printf(" -fA,-fB : file formats.'P'db,'S'df,'K'cf, '2':MOL2, 'C':CCD CIF 'M':SMILES, 'm':SMILES_STDIN [%c %c]\n",molA.filetype,molB.filetype);    
    printf(" -aA,-aB : AtomHetero type. 'A'tom 'H'etatm 'B'oth[%c%c]\n",molA.atmhet,molB.atmhet);    
    printf(" -bA,-bB : Bond type. 'B':consider bonds, otherwise, not consider bonds[%c%c]\n",molA.BondType,molB.BondType);
    printf("<output options for 'pkcombu'>\n"); 
    printf(" -oam   : output file for atom matchings [%s]\n",oamchfile);
    printf(" -oAm   : output file for atom matchings for all the candidate matches[%s]\n",oAmchfile);
    printf(" -opdbA : output PDB  file for molecule A [%s]\n",molA.opdbfile);
    printf(" -osdfA : output SDF  file for molecule A [%s]\n",molA.osdffile);
    printf(" -omol2A: output MOL2 file for molecule A [%s]\n",molA.omol2file);
    printf(" -osmiA : output SMILES  file for molecule A [%s]\n",osmilesfileA);
    printf(" -orasA : output rasmol script for atom type for moleculeA [%s]\n",orasfileA);
    printf(" -opdbB : output PDB  file for molecule B [%s]\n",molB.opdbfile);
    printf(" -osdfB : output SDF  file for molecule B [%s]\n",molB.osdffile);
    printf(" -omol2B: output MOL2 file for molecule B [%s]\n",molB.omol2file);
    printf(" -osmiB : output SMILES  file for molecule B [%s]\n",osmilesfileB);
    printf(" -opdbAB: output PDB file for rotated molA and fixed molB[%s]\n",opdbfileAB);
    printf(" -orasB : output rasmol script for atom type for moleculeB [%s]\n",orasfileB);
    printf(" -oras  : output rasmol script for the best match [%s]\n",orasfile);
    printf(" -omcs  : output maximum common substructue in SDF [%s]\n",omcs_sdffile);
    printf(" -o     : output calculation process 'T' or 'F' [%c]\n",PAR.OutCalProcess);
    printf(" -rk    : rank of output atom MATCH by '-oam' option [%d]\n",RankMatchOutput);
    printf(" -nameA : new name for molecule A [%s]\n",newname_molA);
    printf(" -kplcn : keep only the largest connected subgraph (T or F)[%c]\n",KeepLargestConGraph);
    printf(" -anosmiA : add SMILES string as property annotation for molecule A (T or F)[%c]\n",AnnotationSmiles_molA);
    printf(" -iam   : input file for atom matchings [%s]\n",iamchfile);
    printf("<output options for generating XYZ coodinartes>\n"); 
    printf(" -sup3A  : Do    3D-superimpose  molA onto fixed molB (T or F)[%c]\n",Superimp_molA);
    printf(" -sup3B  : Do    3D-superimpose  molB onto fixed molA (T or F)[%c]\n",Superimp_molB);
    printf(" -pas2   : Paste 2D coordinates to molA from fixed molB and building hydrogen XY (T or F)[%c]\n",Paste2D_molA);
    printf(" -pas3   : Paste 3D coordinates to molA from fixed molB (T or F)[%c]\n",Paste3D_molA);
    printf(" -gen2dA : generate 2D coordinates for molecule A (T or F)[%c]\n",Gen2D_molA);
    printf(" -gen2dB : generate 2D coordinates for molecule B (T or F)[%c]\n",Gen2D_molB);
    printf(" -rscB   : rescaling coordinates of moleculeB (T or F)[%c]\n",Rescale3D_molB);
    printf("<output options for PostScript>\n"); 
    printf(" -ops    : output PostScript file atom matching [%s]\n",opsfile);
    printf(" -pslab  : Label on atoms. 'M'atched_atom_number, 'T':atomtype, 'N':atom number,'A':molA_num,'B':molB_num\n");
    printf("         :                 'a':atom name,'E':element,'F':no label. [%c]\n",PAR.ps_label);
    printf(" -pscir  : Draw filled circle of atoms. 'A'll 'M'arked_atom, 'H'eavy_chain_full, hydrogen_half [%c]\n",PAR.ps_circle);
    printf(" -pshyd  : Hydrogen atoms. 'T':show, 'F':not show [%c]\n",PAR.ps_hydro); 
    printf(" -psplB  : rotate molB on plane [%c]\n",PAR.ps_planeB); 
    printf(" -pssupA : superimpose molA on molB [%c]\n",PAR.ps_supAonB); 
    printf(" -pscol  : color scheme. 'M'onochrome,'C'olor [%c]\n",PAR.ps_color); 
    printf(" -pscline: draw circle outline (T or F) [%c]\n",PAR.ps_circ_outline); 
    printf(" -psbtype: care bond type (single,double) for drawing bonds (T or F) [%c]\n",PAR.ps_care_bondtype);
    printf("<other options>\n"); 
    printf(" -ste3A : Setup Stereo Parity from 3D for mol A (T or F)[%c]\n",molA.SetStereo3D); 
    printf(" -ste3B : Setup Stereo Parity from 3D for mol B (T or F)[%c]\n",molB.SetStereo3D); 
    printf(" -KV    : key-value-style stdout. (T or F)[%c]\n",KeyValueOut);
    printf(" -rkrm  : reranking MATCHes by RMSD (recommended to used with '-per a'). (T or F)[%c]\n",RerankMatchByRMSD);
    printf(" -elch  : compare electric charges (T or F)[%c]\n",CompEleCharge);
    printf(" -chkste:Check stereo parity. 'C':just checking, 'D':delete unmacth chiral atom pairs.[%c]\n",CheckStereoParity);
    printf(" -clshA : Check atomic clashes within the molecule A (T or F)[%c]\n",CheckSelfClash_A);
    printf(" -clshB : Check atomic clashes within the molecule A (T or F)[%c]\n",CheckSelfClash_B);

  }
  if (DetailOptHelp=='M'){ 
    Show_Option_Help_MCS();
  }

  if ((argc<2)||(DetailOptHelp!='-')) exit(1);

 /****** READ ARGUMENTS ***********/
 PAR.OptValHead.next = NULL;
 Read_Options_From_Arguments(argc,argv,PAR.COMMAND,&(PAR.OptValHead));

 ov = &(PAR.OptValHead);
 while (ov->next != NULL){
   ov = ov->next;
   if (Set_Global_PAR_from_OPTION_VALUE(ov)==0){
           if (strcmp(ov->opt,"A")==0)  { sprintf(molA.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"B")==0)  { sprintf(molB.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"T")==0)  { sprintf(molA.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"R")==0)  { sprintf(molB.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"oras")==0)  { sprintf(orasfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"orasA")==0)  { sprintf(orasfileA,"%s",ov->val);}
      else if (strcmp(ov->opt,"orasB")==0)  { sprintf(orasfileB,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdbAB")==0)  { sprintf(opdbfileAB,"%s",ov->val);}

      else if (strcmp(ov->opt,"opdbA")==0)  { sprintf(molA.opdbfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"osdfA")==0)  { sprintf(molA.osdffile,"%s",ov->val);}
      else if (strcmp(ov->opt,"omol2A")==0)  { sprintf(molA.omol2file,"%s",ov->val);}
      else if (strcmp(ov->opt,"osmiA")==0)  { sprintf(osmilesfileA,"%s",ov->val);}

      else if (strcmp(ov->opt,"opdbB")==0)  { sprintf(molB.opdbfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"osdfB")==0)  { sprintf(molB.osdffile,"%s",ov->val);}
      else if (strcmp(ov->opt,"omol2B")==0) { sprintf(molB.omol2file,"%s",ov->val);}
      else if (strcmp(ov->opt,"osmiB")==0)  { sprintf(osmilesfileB,"%s",ov->val);}

      else if (strcmp(ov->opt,"ops")==0)  { sprintf(opsfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"omcs")==0) { sprintf(omcs_sdffile,"%s",ov->val);}
      else if (strcmp(ov->opt,"rkrm")==0) { RerankMatchByRMSD = ov->val[0];}

      else if (strcmp(ov->opt,"oam")==0) { sprintf(oamchfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"iam")==0) { sprintf(iamchfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"oAm")==0) { sprintf(oAmchfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"sup3A")==0) { Superimp_molA = ov->val[0];}
      else if (strcmp(ov->opt,"sup3B")==0) { Superimp_molB = ov->val[0];}
      else if (strcmp(ov->opt,"pas2")==0) {  Paste2D_molA = ov->val[0];}
      else if (strcmp(ov->opt,"pas3")==0) {  Paste3D_molA = ov->val[0];}
      else if (strcmp(ov->opt,"fA")==0) { molA.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"fB")==0) { molB.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"aA")==0) { molA.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"aB")==0) { molB.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"bA")==0) { molA.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"bB")==0) { molB.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"rk")==0) { RankMatchOutput = atoi(ov->val);}
      else if (strcmp(ov->opt,"gen2dA")==0) { Gen2D_molA = ov->val[0];}
      else if (strcmp(ov->opt,"gen2dB")==0) { Gen2D_molB = ov->val[0];}
      else if (strcmp(ov->opt,"elch")==0) { CompEleCharge = ov->val[0];}
      else if (strcmp(ov->opt,"pslab")==0) { PAR.ps_label = ov->val[0];}
      else if (strcmp(ov->opt,"pshyd")==0) { PAR.ps_hydro = ov->val[0];}
      else if (strcmp(ov->opt,"pscir")==0) { PAR.ps_circle = ov->val[0];}
      else if (strcmp(ov->opt,"psplB")==0) { PAR.ps_planeB = ov->val[0];}
      else if (strcmp(ov->opt,"pssupA")==0) { PAR.ps_supAonB = ov->val[0];}
      else if (strcmp(ov->opt,"pscol")==0) { PAR.ps_color = ov->val[0];}
      else if (strcmp(ov->opt,"pscline")==0) { PAR.ps_circ_outline = ov->val[0];}
      else if (strcmp(ov->opt,"psbtype")==0) { PAR.ps_care_bondtype = ov->val[0];}
      else if (strcmp(ov->opt,"anosmiA")==0) { AnnotationSmiles_molA = ov->val[0];} 
      else if (strcmp(ov->opt,"kplcn")==0) { KeepLargestConGraph = ov->val[0];}
      else if (strcmp(ov->opt,"nameA")==0) { sprintf(newname_molA,"%s",ov->val);}
      else if (strcmp(ov->opt,"rscB")==0) {  Rescale3D_molB = ov->val[0];}
      else if (strcmp(ov->opt,"ste3A")==0) { molA.SetStereo3D = ov->val[0];}
      else if (strcmp(ov->opt,"ste3B")==0) { molB.SetStereo3D = ov->val[0];}
      else if (strcmp(ov->opt,"KV")==0)  { KeyValueOut = ov->val[0];}
      else if (strcmp(ov->opt,"chkste")==0)  { CheckStereoParity = ov->val[0];}
      else if (strcmp(ov->opt,"clshA")==0)  { CheckSelfClash_A = ov->val[0];}
      else if (strcmp(ov->opt,"clshB")==0)  { CheckSelfClash_B = ov->val[0];}
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

  /** [1] Read molecule A **/
  okA = Read_MOLECULE(molA.filename,&molA,molA.filetype,molA.atmhet,molA.BondType);
  if (okA==0) {
     printf("#ERROR: file for molA '%s' (type '%c') does not contain any atoms.\n",molA.filename,molA.filetype);
     exit(1);
   }

   if (KeepLargestConGraph == 'T'){
     Set_Connected_Structure(&molA);
     Keep_Only_Largest_Connected_Structure(&molA);
   }


  if (molA.osdffile[0] !='\0'){ Write_SDF_Molecule(molA.osdffile,&molA,""); }
  if (osmilesfileA[0]!='\0'){ Write_SMILES_Molecule(osmilesfileA,&molA,'-'); }

  if (newname_molA[0] != '\0'){ sprintf(molA.name,"%s",newname_molA);}

   if (AnnotationSmiles_molA == 'T'){
     make_LUNITs_from_Molecule_by_DFS(&HeadLunit,&molA);
     string = (char *)malloc(sizeof(char)*(length_string_from_LUNITs(&molA,&HeadLunit)));
     make_string_from_LUNITs(&molA,string,&HeadLunit,KeepLargestConGraph); printf("#string '%s'\n",string);
     Add_key_value_string_to_KEYVALUEs("SMILES",string,&(molA.HeadProperty));
     Free_LUNITs(&HeadLunit);
     free(string);
  }

  if (orasfileA[0]!='\0'){ Write_atomtype_in_rasmol_script(orasfileA,&molA);}


  /** [2] Read molecule B **/
  okB = Read_MOLECULE(molB.filename,&molB,molB.filetype,molB.atmhet,molA.BondType);
  if (okB==0) {
     printf("#ERROR: file for molB '%s' (type '%c') does not contain any atoms.\n",molB.filename,molB.filetype);
     exit(1);
   }

   if (KeepLargestConGraph == 'T'){
     Set_Connected_Structure(&molB);
     Keep_Only_Largest_Connected_Structure(&molB);
   }

  /*
  Show_Ring_Structure(&molB,&(molB.HeadRing));
  */

  if ((PAR.GenEquivAtomPermu=='A')||(PAR.GenEquivAtomPermu=='N')||(PAR.GenEquivAtomPermu=='B')){
    Make_All_Permutation_of_Equivalent_Heavyatoms(&molA); 
    Make_All_Permutation_of_Equivalent_Heavyatoms(&molB); 
    printf("#PERMUTAION molA %d molB %d\n",molA.Npermu,molB.Npermu);
  }
  if (PAR.GenEquivAtomPermu=='a'){
    Make_All_Permutation_of_Equivalent_Heavyatoms(&molA);
    printf("#PERMUTAION molA %d\n",molA.Npermu);
  }

  if (orasfileB[0]!='\0'){ Write_atomtype_in_rasmol_script(orasfileB,&molB);}

   /** [3] Matching Two molecules (molA and molB) **/
  Match_Two_Molecules(&molA,&molB,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,iamchfile);
  printf("#CheckStereoParity '%c'\n", CheckStereoParity );

  if (CheckStereoParity != 'F'){
    check_stereo_parity_for_MATCH(optMlist.next,&molA,&molB,CheckStereoParity,"");
  }
 
  if (PAR.GenEquivAtomPermu=='A'){ Make_Permutated_MATCHlist_from_Top_MATCH(&optMlist,&molA,&molB);}
  if (PAR.GenEquivAtomPermu=='a'){ Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation(&optMlist,&molA,&molB);}
  if (PAR.GenEquivAtomPermu=='N'){ Make_NonRedundant_MATCHlist_by_Permutation(&optMlist,&molA,&molB);}

  if (RerankMatchByRMSD=='T') {Rerank_MATCHlist_by_transformed_RMSD(&optMlist,&molA,&molB);}


   /** [4] Output Matching  **/

   if (Gen2D_molA=='T'){ Generate_2D_Coordinates(&molA);}
   if (Gen2D_molB=='T'){ Generate_2D_Coordinates(&molB);}
   if (Rescale3D_molB == 'T'){ Rescaling_Corrdinates_of_Molecule(&molB); }

   /* Print_MATCHlist(&optMlist,"optMlist(final)"); */

   Write_MATCHlist_in_Horizontal_Element("-",&optMlist,&molA,&molB);
   /*
   Write_MATCHlist_in_Horizontal_NumStr("-",&optMlist,&molA,&molB);
   Write_MATCHlist_in_Horizontal_Number("-",&optMlist,&molA,&molB);
   */

  RMSD = tanimoto = 0.0; 
 
  if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E')){
    Natmmatch = optMlist.next->Npair;

    /** decide 'om' for output MATCH */
    om = optMlist.next;
    k = 1;
    while ((om->next!=NULL)&&(k<RankMatchOutput)){
      om = om->next;
      k += 1;
    }

 /*
    tanimoto = Tanimoto_Coefficient(Natmmatch,molA.Nheavyatom,molB.Nheavyatom);
  */   
    tanimoto = Tanimoto_Coefficient(Natmmatch,&molA,&molB);
    RMSD = Calculate_CRMS_MATCH_Quaternion(om,&molA,&molB,gA,gB,Rmat,"AonB");

    Assign_Aligned_AtomNumber_to_tFactor(om,&molA,'A');
    Assign_Aligned_AtomNumber_to_tFactor(om,&molB,'B');
    printf("#RMSD no_superimpose %.2f superimpose %.2f\n",CRMS_MATCHed_Molecules_without_Superimpose(om,&molA,&molB),RMSD);

    if (Superimp_molA=='T'){
       RMSD = Calculate_CRMS_MATCH_Quaternion(om,&molA,&molB,gA,gB,Rmat,"AonB");
       Rotate_Molecule(&molA,gA,gB,Rmat);
       make_rot_axis_angle_from_rot_matrix(raxisang,Rmat);
       make_tvec_from_rmat_gorig_gtra(tvec, Rmat, gA, gB);
       sprintf(buff_string,"TRANS_A_ON_B rot %lf %lf %lf %lf tvec %lf %lf %lf",raxisang[0],raxisang[1],raxisang[2],raxisang[3],tvec[0],tvec[1],tvec[2]); 
       printf("#%s\n",buff_string); 
       Add_string_to_LINENODEs(buff_string,&HeadComment);
    }

    if (Superimp_molB=='T'){
       RMSD = Calculate_CRMS_MATCH_Quaternion(om,&molA,&molB,gA,gB,Rmat,"BonA");
       Rotate_Molecule(&molB,gB,gA,Rmat);
       make_rot_axis_angle_from_rot_matrix(raxisang,Rmat);
       make_tvec_from_rmat_gorig_gtra(tvec, Rmat, gB, gA);
       sprintf(buff_string,"TRANS_B_ON_A rot %lf %lf %lf %lf tvec %lf %lf %lf",raxisang[0],raxisang[1],raxisang[2],raxisang[3],tvec[0],tvec[1],tvec[2]); 
       printf("#%s\n",buff_string); 
       Add_string_to_LINENODEs(buff_string,&HeadComment);
    }

  

    if (Paste2D_molA=='T'){ Paste_Coordinates_into_matched_molA_atoms(om,&molA,&molB,'2');}
    if (Paste3D_molA=='T'){ Paste_Coordinates_into_matched_molA_atoms(om,&molA,&molB,'3');}


    if (molA.osdffile[0] !='\0'){ Write_SDF_Molecule(molA.osdffile,&molA,""); }
    if (molA.opdbfile[0] !='\0'){ Write_PDB_Molecule(molA.opdbfile,&molA,'w','A',1,-1,'-','-',""); }
    if (molA.omol2file[0]!='\0'){ Write_MOL2_Molecule(molA.omol2file,&molA,'w',""); }
    if (osmilesfileA[0]!='\0'){ Write_SMILES_Molecule(osmilesfileA,&molA,'-'); }

    if (molB.osdffile[0] !='\0'){ Write_SDF_Molecule(molB.osdffile,&molB,""); }
    if (molB.opdbfile[0] !='\0'){ Write_PDB_Molecule(molB.opdbfile,&molB,'w','B',1,-1,'-','-',""); }
    if (molB.omol2file[0]!='\0'){ Write_MOL2_Molecule(molB.omol2file,&molB,'w',""); }
    if (osmilesfileB[0]!='\0'){ Write_SMILES_Molecule(osmilesfileB,&molB,'-'); }

    if (opdbfileAB[0]!='\0'){ 
      Write_PDB_Molecule(opdbfileAB,&molA,'w','A',-1,1,'-','-',""); 
      Write_PDB_Molecule(opdbfileAB,&molB,'a','B',-1,molA.Natom+1,'-','-',""); 
    }
 

    if (opsfile[0] != '\0'){Write_MATCH_in_PostScript(opsfile,om,&molA,&molB);}
 
    if (orasfile[0] != '\0'){
      sprintf(ofname,"%s.ras",orasfile);
      Write_Superimposed_MATCH_in_rasmol_script(ofname,om,&molA,&molB);
      sprintf(ofname,"%s-A.ras",orasfile);
      Write_MATCH_in_rasmol_script(ofname,om,&molA,&molB,'A','M','L');
      sprintf(ofname,"%s-B.ras",orasfile);
      Write_MATCH_in_rasmol_script(ofname,om,&molA,&molB,'B','M','L');
    }

    if (oamchfile[0] != '\0'){ Write_MATCHlist(oamchfile,&optMlist,&molA,&molB,'O',RankMatchOutput,&HeadComment);}
  }

  else {Natmmatch = 0; tanimoto = 0.0;}

  if (oAmchfile[0] != '\0'){ Write_MATCHlist(oAmchfile,&optMlist,&molA,&molB,'A',0,&HeadComment);}

  if ((omcs_sdffile[0] != '\0') && (optMlist.next != NULL)){
     Write_MATCH_in_SDF(omcs_sdffile,om,&molA,&molB); 
  } 
  printf("#COMMENT A %s B %s NhvatmA %d NhvatmB %d Natmmatch %d tanimoto %f\n",molA.filename,molB.filename,molA.Nheavyatom,molB.Nheavyatom,Natmmatch,tanimoto);
  
  /* printf("M %d N %d O %d P %d\n",PAR.M,PAR.N,PAR.O,PAR.P); */

  if (CompEleCharge=='T'){Compare_Electric_Charges(om,&molA,&molB);}

  if (KeyValueOut=='T'){
    printf("NHEAVYATOM_A %d\n",molA.Nheavyatom);
    printf("NHEAVYATOM_B %d\n",molB.Nheavyatom);
    printf("NATOM_MATCH  %d\n",Natmmatch);
    printf("TANIMOTO     %f\n",tanimoto);
    printf("RMSD         %f\n",RMSD);
    printf("RMSD_NO_SUPERIMPOSE %f\n",CRMS_MATCHed_Molecules_without_Superimpose(om,&molA,&molB));
  }

  if (CheckSelfClash_A =='T'){
    printf("#Nclash_atom_pair_A %d\n",count_selfclash(&molA));
    printf("#Eselfclash_A %f\n",energy_selfclash(&molA));
    if (count_selfclash(&molA)>0){ show_selfclash(&molA,"_A"); }
  }

  if (CheckSelfClash_B =='T'){
    printf("#Nclash_atom_pair_B %d\n",count_selfclash(&molB));
    printf("#Eselfclash_B %f\n",energy_selfclash(&molB));
    if (count_selfclash(&molB)>0){ show_selfclash(&molB,"_B"); }
  }



  return(1);

} /* end of main() */ 
