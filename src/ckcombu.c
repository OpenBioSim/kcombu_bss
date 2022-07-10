/*
 
 <ckcombu.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

  for obtaining maximum common substructure
  by build-up procedures for converting molecular format 

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
#include "ioPDB.h"
#include "ioMOL2.h"
#include "ioSMILES.h"
#include "kcombu_io_CIF.h"
#include "MOLECULE_from_CCD_CIF.h"
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
#include "energies.h"
#include "stereo_check.h"

struct PARAMETERS PAR;

/** FUNCTIONS (LOCAL) **/

int main(argc,argv)
 int argc;
 char **argv;
{
  struct MOLECULE molA;
  struct OPTION_VALUE *ov;
  int okA;
  char DetailOptHelp,KeyValueOut;
  char *string; 
  char orasfile[MAX_FILENAME], opsfile[MAX_FILENAME], osmilesfileA[MAX_FILENAME];
  char opdbfileAdamrl[MAX_FILENAME];
  char  Gen2D_molA, RandAtomOrder, AnnotationSmiles,CheckSelfClash;
  char KeepLargestConGraph,CheckStereoParity3D;
  char newname_molA[128];
  struct LUNIT HeadLunit;

  /****** SET INITIAL PARAMETERS ***********/
  Set_Default_Global_PAR(); 
  PAR.PROGRAM_TYPE = 'p'; 
  Initialize_MOLECULE(&molA);
  Initialize_MOLECULE_string(&molA);
  molA.filetype =  '-'; 
  molA.atmhet =  'B';
  molA.chain =  '-'; 
  molA.BondType =  'B'; 
  molA.resi[0] =  '\0'; 
  molA.maxNatom = -1;

  osmilesfileA[0] = '\0'; 
  orasfile[0] = opsfile[0] = '\0'; 
  opdbfileAdamrl[0] = '\0';
 
  PAR.max_ring_size  = 30;
  PAR.max_block_size = 40;
  PAR.max_ring_descriptor  = PAR.max_ring_size + PAR.max_block_size  - 3;

  PAR.ps_label  = 'F';
  PAR.ps_circle = 'A';
  PAR.ps_hydro  = 'F';
  PAR.ps_color   = 'C';
  PAR.ps_circ_outline   = 'F';
  PAR.ps_care_bondtype   = 'T';


  DetailOptHelp = '-';

  Gen2D_molA = 'F';
  RandAtomOrder = 'F';
  AnnotationSmiles = 'F';
  CheckStereoParity3D = 'F';
  newname_molA[0] = '\0';

  KeepLargestConGraph = 'F';

  CheckSelfClash = 'F';
  KeyValueOut = 'F';


  if((argc>1) && ((strcmp(argv[1],"-h")==0) || (strcmp(argv[1],"-help")==0)) ){         DetailOptHelp = 'D';}
  if((argc>1) && ((strncmp(argv[1],"-hm",3)==0) || (strncmp(argv[1],"-HM",3)==0)) ){ DetailOptHelp = 'M';}

  if (argc<2){
    printf("ckcombu <options>\n");
    printf("  %s\n",KCOMBU_STANDS_FOR);
    printf(" 'ckcombu' is for converting molecular file format.\n"); 
    printf("  %s LastModified:%s \n",CODED_BY,LAST_MOD_DATE);
    printf("** Simple Usage for 'ckcombu' **\n");
    printf(">for converting into SDF format.\n");
    printf("  $ckcombu -A [molecule_file] -osdf [output SDF file]\n");
    printf(">for converting into SMILES format.\n");
    printf("  $ckcombu -A [molecule_file] -osmi [output SMILES file]\n");
    printf(">for making 2D-coordinates and writing in SDF format.\n");
    printf("  $ckcombu -A [molecule_file] -gen2d T -osdf [output SDF file]\n");
    printf(">for showing options in detail\n");
    printf("  $ckcombu -h\n");
   }
  if (DetailOptHelp=='D'){
    printf("ckcombu <options>\n");
    printf("  for obtaining maximum common substructures by build-up method.\n");
    printf("  %s\n",KCOMBU_STANDS_FOR); 
    printf("  'ckcombu' is for converting molecular file format.\n"); 
    printf("   coded by T.Kawabata. LastModified:%s\n",LAST_MOD_DATE);
    printf("<input options for 'ckcombu'>\n"); 
    printf(" -A   : molecule A (molA)(*.sdf|*.mol2|*.pdb|*.kcf|*.smi|*.cif)[%s]\n",molA.filename);
    printf(" -fA  : file formats.'P'db,'S'df,'K'cf, '2':MOL2, 'C':CCD CIF, 'M':SMILES 'm':SMILES_STDIN.[%c]\n",molA.filetype);
    printf(" -aA  : AtomHetero type. 'A'tom 'H'etatm 'B'oth[%c]\n",molA.atmhet);
    printf(" -bA  : 'B':consider bonds, otherwise, not consider bonds[%c]\n",molA.BondType);
    printf(" -maxNatomA: Number of maximum atoms for molA. [%d]\n",molA.maxNatom);
    printf("<output options for 'ckcombu'>\n"); 
    printf(" -opdb  : output PDB  file for molecule A [%s]\n",molA.opdbfile);
    printf(" -osdf  : output SDF  file for molecule A [%s]\n",molA.osdffile);
    printf(" -omol2 : output MOL2 file for molecule A [%s]\n",molA.omol2file);
    printf(" -osmi  : output SMILES  file for molecule A [%s]\n",osmilesfileA);
    printf(" -opdbdamrl: output PDB file for molecule A with chain_ids of DAMRL classification.[%s]\n",opdbfileAdamrl);
    printf(" -oras  : output rasmol script for the best match [%s]\n",orasfile);
    printf(" -ops   : output PostScript file for the molecule [%s]\n",opsfile);
    printf(" -o     : output calculation process 'T' or 'F' [%c]\n",PAR.OutCalProcess);
    printf(" -gen2d : generate 2D coordinates for molecule A (T or F)[%c]\n",Gen2D_molA);
    printf("<output options for PostScript>\n"); 
    printf(" -pslab : Label on atoms 'M'atched_atom_number, 'T':atomtype, 'N':atom number, 'a':atom name,'E':element [%c]\n",PAR.ps_label);
    printf(" -pscir  : Draw filled circle of atoms. 'A'll 'M'arked_atom, 'H'eavy_chain_full, hydrogen_half [%c]\n",PAR.ps_circle);
    printf(" -pshyd : Hydrogen atoms. 'T':show, 'F':not show [%c]\n",PAR.ps_hydro);
    printf(" -pscol  : color scheme. 'M'onochrome,'C'olor [%c]\n",PAR.ps_color);
    printf(" -pscline: draw circle outline (T or F) [%c]\n",PAR.ps_circ_outline);
    printf(" -psbtype: care bond type (single,double) for drawing bonds (T or F) [%c]\n",PAR.ps_care_bondtype);
    printf("<output options for 2D depiction>\n"); 
    printf(" -corv : Do collision resolving for 2D (T or F)[%c]\n",PAR.CollisionResolve2D);
    printf("<other output options>\n"); 
    printf(" -clsh  : Check atomic clashes within the moleculeA (T or F)[%c]\n",CheckSelfClash);
    printf(" -ratm  : randomize order of ATOMs (T or F)[%c]\n",RandAtomOrder);
    printf(" -anosmi: add SMILES string as property annotation (T or F)[%c]\n",AnnotationSmiles);
    printf(" -nameA : new name for molecule A [%s]\n",newname_molA);
    printf(" -sp2ub : change space to underbar for key string (T or F)[%c]\n",PAR.CHANGE_KEY_SPACE_TO_UNDERBAR);
    printf(" -no2su : change 'NOTE' into 'SUPPLIERS' (T or F)[%c]\n",PAR.CHANGE_NOTE_TO_SUPPLIERS);
    printf(" -orksm : output rank-draw RasMol Script for SMILES (T or F)[%c]\n",PAR.OutRankForSMILES);
    printf(" -kplcn : keep only the largest connected subgraph (T or F)[%c]\n",KeepLargestConGraph);
    printf(" -chkste3D: check stereo parity of 3D molecule with assigned parity (T or F)[%c]\n",CheckStereoParity3D);
    printf(" -KV    : key-value-style stdout. (T or F)[%c]\n",KeyValueOut);
  }
  if (DetailOptHelp=='M'){ 
    Show_Option_Help_MCS();
  }



  if ((argc<2)||(DetailOptHelp!='-')){ exit(1);}

 /****** READ ARGUMENTS ***********/

 PAR.OptValHead.next = NULL;
 Read_Options_From_Arguments(argc,argv,PAR.COMMAND,&(PAR.OptValHead));
 
 ov = &(PAR.OptValHead);
 while (ov->next != NULL){
   ov = ov->next;
   if (Set_Global_PAR_from_OPTION_VALUE(ov)==0){
           if (strcmp(ov->opt,"A")==0){sprintf(molA.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"o")==0){PAR.OutCalProcess = ov->val[0];}
      else if (strcmp(ov->opt,"KV")==0) { KeyValueOut = ov->val[0];}
      else if (strcmp(ov->opt,"oras")==0)  { sprintf(orasfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdb")==0)  { sprintf(molA.opdbfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"osdf")==0)  { sprintf(molA.osdffile,"%s",ov->val);}
      else if (strcmp(ov->opt,"omol2")==0)  { sprintf(molA.omol2file,"%s",ov->val);}
      else if (strcmp(ov->opt,"osmi")==0)  { sprintf(osmilesfileA,"%s",ov->val);}
      else if (strcmp(ov->opt,"ops")==0)  {  sprintf(opsfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"fA")==0)  { molA.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"aA")==0)  { molA.atmhet   = ov->val[0];}
      else if (strcmp(ov->opt,"bA")==0)  { molA.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"gen2d")==0) { Gen2D_molA = ov->val[0];}
      else if (strcmp(ov->opt,"pslab")==0) { PAR.ps_label = ov->val[0];}
      else if (strcmp(ov->opt,"pscir")==0) { PAR.ps_circle = ov->val[0];}
      else if (strcmp(ov->opt,"pshyd")==0) { PAR.ps_hydro = ov->val[0];}
      else if (strcmp(ov->opt,"pscol")==0) { PAR.ps_color = ov->val[0];}
      else if (strcmp(ov->opt,"pscline")==0) {PAR.ps_circ_outline = ov->val[0];}
      else if (strcmp(ov->opt,"psbtype")==0) {PAR.ps_care_bondtype = ov->val[0];}
      else if (strcmp(ov->opt,"corv")==0) { PAR.CollisionResolve2D = ov->val[0];}
      else if (strcmp(ov->opt,"ratm")==0) { RandAtomOrder = ov->val[0];}
      else if (strcmp(ov->opt,"anosmi")==0) {AnnotationSmiles = ov->val[0];}
      else if (strcmp(ov->opt,"opdbdamrl")==0)  { sprintf(opdbfileAdamrl,"%s",ov->val);}
      else if (strcmp(ov->opt,"nameA")==0) {sprintf(newname_molA,"%s",ov->val);}
      else if (strcmp(ov->opt,"sp2ub")==0) {PAR.CHANGE_KEY_SPACE_TO_UNDERBAR = ov->val[0];}
      else if (strcmp(ov->opt,"no2su")==0) {PAR.CHANGE_NOTE_TO_SUPPLIERS = ov->val[0];}
      else if (strcmp(ov->opt,"orksm")==0) { PAR.OutRankForSMILES = ov->val[0];}
      else if (strcmp(ov->opt,"chkste3D")==0) { CheckStereoParity3D = ov->val[0];}
      else if (strcmp(ov->opt,"kplcn")==0) { KeepLargestConGraph = ov->val[0];}
      else if (strcmp(ov->opt,"clsh")==0) {  CheckSelfClash = ov->val[0];}
      else if (strcmp(ov->opt,"maxNatomA")==0) {  molA.maxNatom = atoi(ov->val);}
      else { printf("#ERROR:Can't understand option %s\n",ov->opt); exit(1);}
   }
 }
printf("#COMMAND '%s'\n",PAR.COMMAND);

 
 if ((PAR.ConnectGraphType=='T') || (PAR.ConnectGraphType=='t')){
   if (PAR.maxDIFtopodis <0) PAR.maxDIFtopodis = 0;
 }
 
  Set_max_atomtype_by_atomtype_class();
  srand(PAR.SeedRand);
  /** [1] Read molecule A **/
  okA = Read_MOLECULE(molA.filename,&molA,molA.filetype,molA.atmhet,molA.BondType);
  if (okA==0) {
    printf("#ERROR: file for molA '%s' (type '%c') does not contain any atoms.\n",molA.filename,molA.filetype);
    exit(1);
  }

  if ((KeepLargestConGraph == 'T')&&(molA.BondType=='B')){
    /* added by T.Kawabata (2018/09/07) */
    Set_Connected_Structure(&molA);
    Keep_Only_Largest_Connected_Structure(&molA); 
    Set_MOLECULE(&molA,'B');
  }
   

/*
  Show_KEYVALUEs(&(molA.HeadProperty));
 */
  if (CheckSelfClash=='T'){
    printf("#Nclash_atom_pair_A %d\n",count_selfclash(&molA));
    printf("#Eselfclash_A %f\n",energy_selfclash(&molA));
    if (count_selfclash(&molA)>0){ show_selfclash(&molA,"_A"); }
  }



  if (RandAtomOrder=='T'){
    randomize_atom_order_of_MOLECULE(&molA);
  }

  if (Gen2D_molA=='T'){ 
    Generate_2D_Coordinates(&molA);
  }
  
  if (newname_molA[0] != '\0') sprintf(molA.name,"%s",newname_molA);

  if (AnnotationSmiles == 'T'){
     make_LUNITs_from_Molecule_by_DFS(&HeadLunit,&molA);
     string = (char *)malloc(sizeof(char)*(length_string_from_LUNITs(&molA,&HeadLunit)));
     make_string_from_LUNITs(&molA,string,&HeadLunit,KeepLargestConGraph); printf("#string '%s'\n",string);
     Add_key_value_string_to_KEYVALUEs("SMILES",string,&(molA.HeadProperty));
     Free_LUNITs(&HeadLunit);
     free(string);
  }

  if (CheckStereoParity3D == 'T'){
    check_3D_stereo_parity_vs_assigned_parity(&molA);
   /*
    Set_Stereo_Parity_of_Molecule_from_3D(&molA,'C'); 
    Set_Stereo_Parity_of_Molecule_from_2D_BondStereo(&molA,'C');
    */
  }
  if (molA.osdffile[0] !='\0'){ Write_SDF_Molecule(molA.osdffile,&molA,""); }
  if (osmilesfileA[0]!='\0'){ Write_SMILES_Molecule(osmilesfileA,&molA); }
  if (molA.opdbfile[0] !='\0'){ Write_PDB_Molecule(molA.opdbfile,&molA,'w','A',-1,-1,'-','-',""); }

  if (opdbfileAdamrl[0] !='\0'){ 
    Set_atomtype_of_Atoms(&molA,'D');
    Write_PDB_Molecule(opdbfileAdamrl,&molA,'w','A',-1,-1,'D','-',"");
  }

  if (molA.omol2file[0]!='\0'){ Write_MOL2_Molecule(molA.omol2file,&molA,'w',""); }
  if (opsfile[0] != '\0'){ Write_Molecule_in_PostScript(opsfile,&molA);}

  if (KeyValueOut=='T'){
    printf("MOLFORM    %s\n",molA.molform);
    printf("NATOM      %d\n",molA.Natom);
    printf("NHEAVYATOM %d\n",molA.Nheavyatom);
    printf("NRING      %d\n",molA.Nring);
    printf("NRINGBLOCK %d\n",molA.Nringblock);
  }
  return(1);

} /* end of main() */ 
