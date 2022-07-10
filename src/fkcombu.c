/*
 
 <fkcombu.c> 

==============================================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under the three-clause BSD License, see LICENSE.txt.
==============================================================================

  'fkcombu' is for pairwise comparison and flexible molecular alignment.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include "globalvar.h"
#include "2DMAP.h"
#include "molecule.h"
#include "ioLINE.h"
#include "ioSDF.h"
#include "ioPDB.h"
#include "ioMOL2.h"
#include "ioLib.h"
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
#include "transform.h"
#include "energies.h"
#include "steep_descen.h"
#include "tpdisrest.h"
#include "stamp_transf.h"
#include "stereo_adjust.h"
#include "AtmPairDesc.h"
#include "ringblock.h"
#include "RingDesc.h"
#include "DualAtmPair.h"
#include "options.h"
#include "outPS.h"
#include "ioSMILES.h"

struct PARAMETERS PAR;


struct MCS_TYPE{
  char   string[64];       /* String from command. such as 'C:0:K', 'T:1:E' */
  char   ConnectGraphType; /* Connection of MCS. 'C'onnected, 'D'isconnected, 'T'D-MCS,'t'C-MCS */
  int    maxDIFtopodis;    /* accepted maximum difference of topological distance (number of bonds in the shortest path) */
  char   atomtype_class;   /* atomtype classification */
};



int main(argc,argv)
 int argc;
 char **argv;
{
  struct OPTION_VALUE *ov;
  struct MOLECULE molT,molR,molP,molTo;
  struct TRANSFORM_VARIABLE traT;
  struct MOLECULE *molTtrial;            /* array of transformed molT [NinitialGradFit] (malloc later) */ 
  struct TRANSFORM_VARIABLE *traTtrial;  /* array of transformed variable for molT [NinitialGradFit] (malloc later) */ 

  int t,k,n,nn,ns,ok,r,okT,okR,okP,Nam;
  char DetailOptHelp;
  struct MATCH optMlist,optMlist0;
  struct MATCH *m_model;
  char orasfile[MAX_FILENAME], opdbfileTR[MAX_FILENAME],opdbfileTP[MAX_FILENAME],opdbfileTRP[MAX_FILENAME],SeparateOutT;
  char omcs_sdffile[MAX_FILENAME],oamchfile[MAX_FILENAME],iamchfile[MAX_FILENAME],oAmchfile[MAX_FILENAME];
  char ofname[MAX_FILENAME],buff[MAX_FILENAME],otransfile[MAX_FILENAME];
  float tanimoto,RMSD;
  int   Natom_match;
  double g1[3],g2[3],Rmat[3][3];
  char KeyValueOut,EnergyStrategy,SearchStrategy,RigidFlexibleType;
  char RgdMatchFitType,CnfStampingType,SteepestDescentType,ChiralSp3FitType,CnfStampRingType,CheckRingConfType,RgdPcaFitType;
  char InitChiralChange,CnfEselfclashMin,DoMCS, CnfRandInitType,RgdRandInitType;
  int   NiterateGradFit, NinitialGradFit, NoutconfGradFit;
  float maxRotBondAngleInit,rmsd;
  int   n_min,RankAtomMatchFor3D, MaxAtomMatchFor3D, Ntarget_conf;
  char SortValueForTargets,wmode,PlaneDecision,FixPlaneAngle,FixStampAngle;
  char  comment[1024];
  int  *sindex; 
  char PermuRmsdMin,AdjustStereo2D,CheckStereo2D,CheckBondAngle,CheckStereo3D,AdjustMatchStereo3D;
  struct LINENODE HeadComment;
  char opsfileT[MAX_FILENAME],opsfile[MAX_FILENAME];
  char substring[MAX_FILENAME];
  int  Wsta[100],Wend[100],Nword;
  float ProbFlipRingInit;
  int  minNselfclash;
  char PyramidInversion,AllReflectChiral,FoldHexaRing; 
  int  Ndiff_chiral_sp3,Ndiff_chiral_sp3_init;
  int  Ndiff_ring_conf,Ndiff_ring_conf_init;
  char newname_molT[128];
  char new_resi_molT[4];  /* new residue name     */
  char new_rnum_molT[6];  /* new residue number   */
  char new_chain_molT;    /* new chain identifier */
  
  int    Nmcs_types;
  struct MCS_TYPE mcs_types[100];
  char KeepLargestConGraph;

 /****** SET INITIAL PARAMETERS ***********/
  Set_Default_Global_PAR(); 
  substring[0] = '\0';
  KeyValueOut = 'T';
  PAR.PROGRAM_TYPE = 'p'; 


  Initialize_MOLECULE(&molT);
  Initialize_MOLECULE_string(&molT);

  Initialize_MOLECULE(&molR);
  Initialize_MOLECULE_string(&molR);

  Initialize_MOLECULE(&molP);
  Initialize_MOLECULE_string(&molP);

  Initialize_MOLECULE(&molTo);
  Initialize_MOLECULE_string(&molTo);

  molT.filename[0] = molR.filename[0] = '\0';
  molT.filetype = molR.filetype = '-'; 
  molT.atmhet = molR.atmhet = 'B';
  molT.opdbfile[0]  = molR.osdffile[0] = '\0';  
  molT.chain = molR.chain = '-';
  molT.resi[0] = molR.resi[0] = '\0';
  molT.BondType = molR.BondType = 'B';

  molT.SetStereo3D = molR.SetStereo3D = 'T';

  KeepLargestConGraph = 'F';  
 
  molP.filename[0] =  '\0';
  molP.filetype = '-'; 
  molP.atmhet = 'A';
  molP.Natom = molP.Nheavyatom = molP.Nbond = 0;
  molP.chain = '-';
  molP.resi[0] = '\0';
  molP.BondType = 'F';
 
  omcs_sdffile[0] =  oamchfile[0] = iamchfile[0] = oAmchfile[0] =  '\0';
  opdbfileTR[0]= opdbfileTP[0] = opdbfileTRP[0] = otransfile[0] = '\0'; 
  orasfile[0] = '\0';
  opsfileT[0] = opsfile[0] = '\0'; 
  PAR.maxTRA = 0.1;
  PAR.maxROT = 10.0; /* degree */
  PAR.maxRB  = 10.0; /* degree */
 
  EnergyStrategy     = 'A'; 
  SearchStrategy  = 'F'; 
  RigidFlexibleType = 'F'; 
  DoMCS = 'T';
  RgdPcaFitType = 'F';
  RgdMatchFitType = 'T';
  CnfStampingType = 'T';
  SteepestDescentType = 'T';
  ChiralSp3FitType = 'F';
  CnfStampRingType = 'T';
  CheckRingConfType = 'F';
  CnfRandInitType     = 'T';
  RgdRandInitType     = 'T';
  CnfEselfclashMin = 'F';
  NinitialGradFit = 10;
  NoutconfGradFit = 1;
  NiterateGradFit = 100;
  maxRotBondAngleInit  = 180.0; /* degree */
  InitChiralChange = 'F'; 
  RankAtomMatchFor3D = -1;
  MaxAtomMatchFor3D  =  1;

  PAR.max_separation = 10;
  PAR.simtype_descriptor='T';
  PAR.type_atompair = 'F';
  
  PAR.max_ring_size  = 30;
  PAR.max_block_size = 40;

  PAR.max_ring_descriptor  = PAR.max_ring_size + PAR.max_block_size  - 3;

  PAR.max_separation_cyclic  = 10;
  PAR.max_separation_acyclic = 10;


  PAR.WEatommatch =  1.0;
  PAR.WEselfclash =  1.0;
  PAR.WEprotclash =  1.0;
  PAR.WEprotatrct =  0.0;
  PAR.WEvolmovlap =  0.0;
  PAR.WEtpdisrest =  0.0;

  PAR.p_volmovlap     = 2.82842712474619009760; /* p = 2*sqrt(2) */
  PAR.atomtype_volmovlap = 'I';
  PlaneDecision = 'X';
  FixStampAngle  = 'R';
  FixPlaneAngle = 'I';
  PermuRmsdMin  = 'F';
  AdjustStereo2D = 'F';
  CheckStereo2D = 'F';
  CheckStereo3D = 'F';
  CheckBondAngle = 'F';
  AdjustMatchStereo3D = 'T';
  PAR.HydrogenOptimize = 'T';
  PyramidInversion = 'T';
  AllReflectChiral = 'T';
  FoldHexaRing = 'T'; 
  ProbFlipRingInit = 0.0;

  SortValueForTargets = 'E';
  PAR.string_tpdisrest[0][0] = PAR.string_tpdisrest[1][0] =  PAR.string_tpdisrest[2][0] = '\0';
  SeparateOutT = 'F';

  HeadComment.next = NULL;
   
  newname_molT[0] = '\0';
  new_resi_molT[0] = '\0'; 
  new_rnum_molT[0] = '\0';
  new_chain_molT   = ' ';

  Ndiff_chiral_sp3 = Ndiff_chiral_sp3_init = 0;
  Ndiff_ring_conf  = Ndiff_ring_conf_init = 0;
  molTtrial = NULL;
  traTtrial = NULL;
  sindex = NULL;
  nn = 0;
  n_min = 0;

  DetailOptHelp = '-';
  if((argc>1) && ((strcmp(argv[1],"-h")==0) || (strcmp(argv[1],"-help")==0)) )         DetailOptHelp = 'D';
  if((argc>1) && ((strncmp(argv[1],"-hm",3)==0) || (strncmp(argv[1],"-HM",3)==0)) ) DetailOptHelp = 'M';


  mcs_types[0].ConnectGraphType =  PAR.ConnectGraphType;;
  mcs_types[0].maxDIFtopodis    =  PAR.maxDIFtopodis;   
  mcs_types[0].atomtype_class   =  PAR.atomtype_class;

  Nmcs_types = 0;

  if (argc<2){
    printf("fkcombu <options>\n");
    printf("  %s\n",KCOMBU_STANDS_FOR);
    printf(" 'fkcombu' is for pairwise comparison and flexible molecular alignment\n"); 
    printf("  %s LastModified:%s\n",CODED_BY,LAST_MOD_DATE);
    printf("** Simple Usage for 'fkcombu' **\n");
    printf(">Transforming of the target molecule onto the reference.\n");
    printf("  $fkcombu -T [target molecule] -R [reference molecule] -opdbT [output target in PDB]\n");
    printf(">Transforming of the target molecule onto the reference avoiding the receptor.\n");
    printf("  $fkcombu -T [target molecule] -R [reference molecule] -P [receptor protein molecule] -opdbT [output target in PDB]\n");
    printf(">Assigning higher number of '-nini' (number of initial conformations) leads to a lower energy with a larger computation.\n");
    printf("  $fkcombu -nini 20 -T [target molecule] -R [reference molecule] -P [receptor molecule] -opdbT [output target in PDB]\n");
    printf(">Transforming of the target molecule onto the reference with volume-overlap by rigid fitting \n");
    printf("  $fkcombu -T [target molecule] -R [reference molecule] -opdbT [output target in PDB] -E V -S R\n");
    printf(">Transforming of the target molecule onto the reference with volume-overlap by flexible fitting \n");
    printf("  $fkcombu -T [target molecule] -R [reference molecule] -opdbT [output target in PDB] -E V -S F\n");
    printf(">Transforming of the target molecule stereo chilarity onto the reference with 2D stereo chilarity\n");
    printf("  $fkcombu -T [target molecule] -R [reference molecule] -opdbT [output target in PDB] -adjste2D T\n");
    printf(">for showing options in detail\n");
    printf("  $fkcombu -h\n");
    printf(">for showing options for MCS\n");
    printf("  $fkcombu -hmcs\n");
  }
  if (DetailOptHelp=='D'){
    printf(">for showing all the options\n");
    printf("  $fkcombu -h\n");
    printf("<options for input molecules>\n");
    printf(" -T   :target    molecule T (molT) (*.sdf|*.mol2|*.pdb|*.kcf)[%s]\n",molT.filename);
    printf(" -R   :reference molecule R (molR) (*.sdf|*.mol2|*.pdb|*.kcf)[%s]\n",molR.filename);
    printf(" -P   :protein   receptor P (molP) (*.sdf|*.mol2|*.pdb|*.kcf)[%s]\n",molP.filename);
    printf(" -fT, -fR, -fP : file formats.'P'db,'S'df,'K'cf, '2':MOL2, 'C':CCD CIF [%c%c%c]\n",molT.filetype,molR.filetype,molP.filetype);
    printf(" -aT, -aR, -aP : AtomHetero types. 'A'tom 'H'etatm 'B'oth for PDB[%c%c%c]\n",molT.atmhet,molR.atmhet,molP.atmhet);
    printf(" -bT, -bR, -bP : Bond type. 'B':consider bonds, otherwise, not consider bonds[%c%c%c]\n",molT.BondType,molR.BondType,molP.BondType);
    printf(" -chT,-chR,-chP: ChainIDs for PDB [%c%c%c]\n",molT.chain,molR.chain,molP.chain);
    printf(" -resT, -resR     : residue name for molecule T,R[%s%s]\n",molT.resi,molR.resi);
    printf("<options for general stargety of transformation>\n");
    printf(" -E  : Energy. 'A'tom-match, 'V'olume-overlap,         'X':by detailed options.[%c]\n",EnergyStrategy);
    printf(" -S  : Search. 'F'lexible,'R'igid,      'N':do nothing 'X':by detailed options [%c]\n",SearchStrategy);
    printf("<options for transformation procedures>\n");
    printf(" -mcs     : Do MCS. ('T' or 'F'). If 'F', do not use any atom match [%c]\n",DoMCS); 
    printf(" -iam     : input file for atom matching [%s]\n",iamchfile);
    printf(" -SD      : Gradient-based Steepest Descent fitting (T or F) [%c]\n",SteepestDescentType); 
    printf(" -rgdflx  : Rigid or Flexible.  R'igid-body(do not set up rotational bond),'F':lexible (set up rotational bond) [%c]\n",RigidFlexibleType);
    printf(" -cfstp   : Conf Change by stamping rotatable bond into molT from molR (T or F) [%c]\n",CnfStampingType); 
    printf(" -cfesc   : Conf Change by minimizing only Eselfclash (only for -niini >=2,recommend for volume-flexible fitting) [%c]\n",CnfEselfclashMin);
    printf(" -cfstprng: Stamp nonplanar 5- or 6-ring conformations (T or F) [%c]\n",CnfStampRingType); 
    printf(" -rgmch   : rigid-body rmsd-min rigid-body fitting of matched atoms bwn molT and molR (T or F) [%c]\n",RgdMatchFitType); 
    printf(" -rgpca   : rigid-body PCA-Evolmovlap fitting (T or F) [%c]\n",RgdPcaFitType); 
    printf("<output options for transformed molecular file>\n"); 
    printf(" -KV      : key-value-style stdout. (T or F)[%c]\n",KeyValueOut);
    printf(" -opdbT   : output PDB file  for molT [%s]\n",molT.opdbfile);
    printf(" -osdfT   : output SDF file  for molT [%s]\n",molT.osdffile);
    printf(" -omol2T  : output MOL2 file for molT [%s]\n",molT.omol2file);
    printf(" -opdbTR  : output PDB file  for rotated molT and fixed molR[%s]\n",opdbfileTR);
    printf(" -opdbTP  : output PDB file  for rotated molT and receptor protein[%s]\n",opdbfileTP);
    printf(" -opdbTRP : output PDB file  for rotated molT and fixed molR and receptor protein[%s]\n",opdbfileTRP);
    printf(" -otra    : output TRANSFORM variable [%s]\n",otransfile); 
    printf(" -nameT   : new name for molecule T [%s]\n",newname_molT);
    printf(" -newresT : new residue name for molT ([resi_name]:[chain]:[resi_number]) for PDB output[%s:%c:%s]\n",new_resi_molT,new_chain_molT,new_rnum_molT);
    printf(" -sepoutT : output multiple conformations with separated files (*_1.*, *_2,*,...)  (T or F) [%c]\n",SeparateOutT);
    printf("<options for adjust 3D stereo into 2D stereo reference>\n");
    printf(" -adjste2D : adjusting 3D stereo into 2D stereo reference (T or F)[%c]\n",AdjustStereo2D); 
    printf(" -pyramid  : Do pyramid inversion for stereo-parity agreement (-ste) (T or F). [%c]\n",PyramidInversion); 
    printf(" -reflect  : Do all-atom Reflection after exchange and pyramid (T or F) [%c]\n",AllReflectChiral);  
    printf(" -foldhex  : fold hexa-ring with 2-disagree atoms after exchange and pyramid(T or F)[%c]\n",FoldHexaRing);   
    printf("<options for transformations>\n");
    printf(" -chsp3  : Chilarity fit for sp3  atoms (T or F) [%c]\n",ChiralSp3FitType); 
    printf(" -hyop   : Hydrogen-position optimize (T or F)[%c]\n",PAR.HydrogenOptimize); 
    printf("<options about MCS atom matches >\n");
    printf(" -rkam  : rank of employing atom match for the 3D-modelling (-1:not assigned)  [%d]\n",RankAtomMatchFor3D); 
    printf(" -ops   : output PostScript file atom matching [%s]\n",opsfile);
    printf(" -opsT  : output PostScript file for molT [%s]\n",opsfileT);
    printf("<options about multiple types of MCS atom matches >\n");
    printf(" -mcs[x] : for [x]-th type for multiple MCS matches (x=0,1,2,..). [con]:[mtd]:[at]:\n");
    printf("           ex) -mcs0 C:-1:K -mcs1 T:1:E -mcs2 T:1:K \n");
    printf(" -mxam   : maximum number of atom matches for the 3D-modelling [%d]\n",MaxAtomMatchFor3D); 
    printf(" -sort : sort value for multiple target 3D models. 'E':Etotal, 'V':tanimoto_volume [%c]\n",SortValueForTargets); 
    printf("<options for chirarity check and MCS adjustment.>\n");
    printf(" -chkrng : Check nonplanar 5- or 6-ring conformations (T or F) [%c]\n",CheckRingConfType); 
    printf(" -chste2D: checking chirarity agreement stereo parities of molR 2D SDF file(T or F)[%c]\n",CheckStereo2D); 
    printf(" -chste3D: checking chirarity agreement between two 3D conformations (T or F)[%c]\n",CheckStereo3D); 
    printf(" -chbang : checking bond angle agreement between two 3D conformations (T or F)[%c]\n",CheckBondAngle); 
    printf(" -amste3D: adjust atom match to agree 3D stereo chirality by permutations (T or F)[%c]\n",AdjustMatchStereo3D); 
    printf(" -permsd : Permutations for getting minimum RMSD. 'T':permute only target,'B'oth_permute, 'F'alse[%c]\n",PermuRmsdMin); 
    printf("<detailed options for steepest-descent optimization>\n");
    printf(" -nini  : num. of random initial conformatinon. (nini==1, only keep the input conf). [%d]\n",NinitialGradFit);
    printf(" -nout  : num. of output conformations [%d]\n",NoutconfGradFit);
    printf(" -sort  : sort value for multiple target 3D models. 'E':Etotal, 'V':tanimoto_volume [%c]\n",SortValueForTargets); 
    printf(" -cfrnd : do randomize conformation    before steepest descent ('T' or 'F')[%c]\n",CnfRandInitType);
    printf(" -rgrnd : do randomize rigid-body pose before steepest descent ('T' or 'F')[%c]\n",RgdRandInitType);
    printf(" -xroini: initial maximum step for rotational bond (degree) [%f]\n",maxRotBondAngleInit);
    printf(" -ich   : initial chirarity change ('T' or 'F') [%c]\n",InitChiralChange);
    printf(" -pflp  : probabality of flip ring fragment for initial conf (0.0..1.0)[%f]\n",ProbFlipRingInit); 
    printf(" -xtr   : maximum step for translation      for optimization (angstrom)[%f]\n",PAR.maxTRA); 
    printf(" -xro   : maximum step for rotation         for optimization (degree)  [%f]\n",PAR.maxROT); 
    printf(" -xrb   : maximum step for rotational bond  for optimization (degree)  [%f]\n",PAR.maxRB); 
    printf(" -nig   : Number of iteration for gradient-based fitting [%d]\n",NiterateGradFit); 
    printf(" -pln   : plane decision for rotational bond. '3':by 3D structure of molT. '2':by 2D strucutre of molT with hydrogens. 'X':if hydrogens are included then 2D, otherwise 3D. [%c]\n",PlaneDecision);
    printf(" -fixstp: fixing rule for stamped bond for steepest descent initials. 'F':fix. 'R'otatable. [%c]\n",FixStampAngle);
    printf(" -fixpln: fixing rule for plane   bond for steepest descent initials. 'F':fix. 'I'nversion(allow 180 degree rotation). 'R'otatable.[%c]\n",FixPlaneAngle);
    printf("<options for potential energy>\n");
    printf(" -weam  : Weight for Eatommatch [%f]\n",PAR.WEatommatch); 
    printf(" -wesc  : Weight for Eselfclash [%f]\n",PAR.WEselfclash); 
    printf(" -wepc  : Weight for Eprotclash [%f]\n",PAR.WEprotclash); 
    printf(" -wepa  : Weight for Eprotatrct [%f]\n",PAR.WEprotatrct); 
    printf(" -wevo  : Weight for Evolmovlap [%f]\n",PAR.WEvolmovlap); 
    printf(" -wetp  : Weight for Etpdisrest [%f]\n",PAR.WEtpdisrest);
    printf(" -tolc  : tolerant distance for clash [%f]\n",PAR.tole_dis_clash);
    printf(" -dxc   : maximum distance for clash energy [%f]\n",PAR.Dmax_clash);
    printf(" -pvo   : parameter p for Evolmovlap [%f]\n",PAR.p_volmovlap);
    printf(" -avo   : atomtype score for Evolmovlap. 'I'gnore,'C'are [%c]\n",PAR.atomtype_volmovlap);
    printf(" -tppair[x]: [x]-th target-protein atom pairs for Etpdisrest. (Tatom_num):(Patom_num):(Dupper). [x] is 0 or 1 or 2. [%s %s %s]\n",
       PAR.string_tpdisrest[0], PAR.string_tpdisrest[1], PAR.string_tpdisrest[2]);
    printf("<options for stereo parity adjusting '-adjste3D T'>\n"); 
    printf(" -miflp: minimum angle for flip of fragment for -pflp>0.0. [%f]\n",PAR.minAngleFlipOfFragment);  
    printf("<other options>\n"); 
    printf(" -oam  : output file for atom matchings [%s]\n",oamchfile);
    printf(" -oAm  : output file for atom matchings for all the candidate matches[%s]\n",oAmchfile);
    printf(" -oras : output rasmol script for the best match [%s]\n",orasfile);
    printf(" -omcs : output maximum common substructue in SDF [%s]\n",omcs_sdffile);
    printf(" -o    : output calculation process 'T' or 'F' [%c]\n",PAR.OutCalProcess);
    printf(" -kplcn : keep only the largest connected subgraph (T or F)[%c]\n",KeepLargestConGraph);
     
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
           if (strcmp(ov->opt,"T")==0) { sprintf(molT.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"A")==0) { sprintf(molT.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"R")==0) { sprintf(molR.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"B")==0) { sprintf(molR.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"P")==0) { sprintf(molP.filename,"%s",ov->val);}
      else if (strcmp(ov->opt,"E")==0) { EnergyStrategy  = ov->val[0];}
      else if (strcmp(ov->opt,"S")==0) { SearchStrategy  = ov->val[0];}

      else if (strcmp(ov->opt,"rgdflx")==0)  { RigidFlexibleType = ov->val[0];}
      else if (strcmp(ov->opt,"K")==0)   { KeyValueOut = ov->val[0];}
      else if (strcmp(ov->opt,"ap")==0)  { PAR.type_atompair = ov->val[0];}
      else if (strcmp(ov->opt,"rkam")==0)  { RankAtomMatchFor3D = atoi(ov->val);}
      else if (strcmp(ov->opt,"mxam")==0)   { MaxAtomMatchFor3D = atoi(ov->val);}
      else if (strcmp(ov->opt,"sort")==0)   { SortValueForTargets = ov->val[0];}
      else if (strcmp(ov->opt,"oras")==0)   { sprintf(orasfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdbTR")==0)  { sprintf(opdbfileTR,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdbTP")==0)  { sprintf(opdbfileTP,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdbTRP")==0)  { sprintf(opdbfileTRP,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdbT")==0)   { sprintf(molT.opdbfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"opdbR")==0)   { sprintf(molR.opdbfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"osdfT")==0)   { sprintf(molT.osdffile,"%s",ov->val);}
      else if (strcmp(ov->opt,"omol2T")==0)   { sprintf(molT.omol2file,"%s",ov->val);}
      else if (strcmp(ov->opt,"omcs")==0)   { sprintf(omcs_sdffile,"%s",ov->val);}
      else if (strcmp(ov->opt,"oam")==0)  { sprintf(oamchfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"iam")==0)  { sprintf(iamchfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"oAm")==0) { sprintf(oAmchfile,"%s",ov->val);}
      else if (strcmp(ov->opt,"xtr")==0)  { PAR.maxTRA = atof(ov->val);}
      else if (strcmp(ov->opt,"xro")==0)  { PAR.maxROT = atof(ov->val);}
      else if (strcmp(ov->opt,"xrb")==0)  { PAR.maxRB  = atof(ov->val);}
      else if (strcmp(ov->opt,"mcs")==0)  { DoMCS  = ov->val[0];}
      else if (strcmp(ov->opt,"rgmch")==0) { RgdMatchFitType  = ov->val[0];}
      else if (strcmp(ov->opt,"cfstp")==0)  { CnfStampingType  = ov->val[0];}
      else if (strcmp(ov->opt,"chsp3")==0)   { ChiralSp3FitType  = ov->val[0];}
      else if (strcmp(ov->opt,"cfstprng")==0)  { CnfStampRingType = ov->val[0];}
      else if (strcmp(ov->opt,"chkrng")==0) { CheckRingConfType = ov->val[0];}
      else if (strcmp(ov->opt,"SD")==0)  { SteepestDescentType  = ov->val[0];}
      else if (strcmp(ov->opt,"rgpca")==0)   { RgdPcaFitType  = ov->val[0];}
      else if (strcmp(ov->opt,"fixstp")==0) { FixStampAngle = ov->val[0];}
      else if (strcmp(ov->opt,"fixpln")==0)  { FixPlaneAngle = ov->val[0];}
      else if (strcmp(ov->opt,"pln")==0)  {  PlaneDecision = ov->val[0];}
      else if (strcmp(ov->opt,"adjste2D")==0)   { AdjustStereo2D  = ov->val[0];  molR.SetStereo3D = 'F';}
      else if (strcmp(ov->opt,"nini")==0)  { NinitialGradFit  = atoi(ov->val);}
      else if (strcmp(ov->opt,"nig")==0)  {  NiterateGradFit  = atoi(ov->val);}
      else if (strcmp(ov->opt,"nout")==0)  { NoutconfGradFit = atoi(ov->val);} 
      else if (strcmp(ov->opt,"xroini")==0)  {  maxRotBondAngleInit   = atof(ov->val);}
      else if (strcmp(ov->opt,"ich")==0)  {  InitChiralChange = ov->val[0];}
      else if (strcmp(ov->opt,"fT")==0)  {molT.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"fR")==0)  {molR.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"fP")==0)  {molP.filetype = MoleculeFileType(ov->val);}
      else if (strcmp(ov->opt,"aT")==0)  { molT.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"aR")==0)  { molR.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"aP")==0)  { molP.atmhet = ov->val[0];}
      else if (strcmp(ov->opt,"bT")==0)  { molT.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"bR")==0)  { molR.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"bP")==0)  { molP.BondType = ov->val[0];}
      else if (strcmp(ov->opt,"sep")==0)  { PAR.max_separation = atoi(ov->val);} 
      else if (strcmp(ov->opt,"pvo")==0)  { PAR.p_volmovlap = atof(ov->val);} 
      else if (strcmp(ov->opt,"avo")==0)   { PAR.atomtype_volmovlap = ov->val[0];} 
      else if (strcmp(ov->opt,"sepc")==0)   { PAR.max_separation_cyclic  = atoi(ov->val);} 
      else if (strcmp(ov->opt,"sepa")==0)   { PAR.max_separation_acyclic = atoi(ov->val);} 
      else if (strcmp(ov->opt,"weam")==0)   { PAR.WEatommatch = atof(ov->val);} 
      else if (strcmp(ov->opt,"wesc")==0)   { PAR.WEselfclash = atof(ov->val);} 
      else if (strcmp(ov->opt,"wepc")==0)   { PAR.WEprotclash = atof(ov->val);} 
      else if (strcmp(ov->opt,"wepa")==0)   { PAR.WEprotatrct = atof(ov->val);} 
      else if (strcmp(ov->opt,"wevo")==0)   { PAR.WEvolmovlap = atof(ov->val);} 
      else if (strcmp(ov->opt,"wetp")==0)   { PAR.WEtpdisrest = atof(ov->val);} 
      else if (strcmp(ov->opt,"tolc")==0)   { PAR.tole_dis_clash  = atof(ov->val);} 
      else if (strcmp(ov->opt,"cfcnd")==0)   { CnfRandInitType = ov->val[0];} 
      else if (strcmp(ov->opt,"rgrnd")==0)   { RgdRandInitType = ov->val[0];} 
      else if (strcmp(ov->opt,"cfesc")==0)   { CnfEselfclashMin = ov->val[0];} 
      else if (strcmp(ov->opt,"dxc")==0)  { PAR.Dmax_clash  = atof(ov->val);} 
      else if (strcmp(ov->opt,"chT")==0)  { molT.chain = ov->val[0];} 
      else if (strcmp(ov->opt,"chR")==0)  { molR.chain = ov->val[0];} 
      else if (strcmp(ov->opt,"chP")==0)  { molP.chain = ov->val[0];} 
      else if (strcmp(ov->opt,"resT")==0)  { sprintf(molT.resi,"%s",ov->val);} 
      else if (strcmp(ov->opt,"resR")==0)  { sprintf(molR.resi,"%s",ov->val);} 
      else if (strcmp(ov->opt,"ops")==0)  {  sprintf(opsfile,"%s",ov->val);} 
      else if (strcmp(ov->opt,"opsT")==0)  {  sprintf(opsfileT,"%s",ov->val);} 
      else if (strcmp(ov->opt,"otra")==0)  {  sprintf(otransfile,"%s",ov->val);} 
      else if (strcmp(ov->opt,"hyop")==0)  {  PAR.HydrogenOptimize = ov->val[0];} 
      else if (strcmp(ov->opt,"sub")==0)  {  sprintf(substring,"%s",ov->val);} 
      else if (strcmp(ov->opt,"pflp")==0)  {  ProbFlipRingInit = atof(ov->val);} 
      else if (strcmp(ov->opt,"pyramid")==0)  {  PyramidInversion = ov->val[0];} 
      else if (strcmp(ov->opt,"reflect")==0) { AllReflectChiral = ov->val[0];} 
      else if (strcmp(ov->opt,"foldhex")==0) {  FoldHexaRing = ov->val[0];} 
      else if (strcmp(ov->opt,"miflp")==0)  { PAR.minAngleFlipOfFragment = atoi(ov->val);} 
      else if (strcmp(ov->opt,"tppair0")==0)  { 
        sprintf(PAR.string_tpdisrest[0],"%s",ov->val); printf("#%s\n",PAR.string_tpdisrest[0]);
        if (PAR.WEtpdisrest==0.0){ PAR.WEtpdisrest = 1.0;}
        /* printf("#PAR.string_tpdisrest[0] '%s'\n",PAR.string_tpdisrest[0]); */
      } 
      else if (strcmp(ov->opt,"tppair1")==0)  { 
        sprintf(PAR.string_tpdisrest[1],"%s",ov->val); printf("#%s\n",PAR.string_tpdisrest[1]);
        if (PAR.WEtpdisrest==0.0){ PAR.WEtpdisrest = 1.0;}
      } 
      else if (strcmp(ov->opt,"tppair2")==0)  { 
        sprintf(PAR.string_tpdisrest[2],"%s",ov->val); printf("#%s\n",PAR.string_tpdisrest[2]);
      } 
      else if ((strncmp(ov->opt,"mcs",3)==0) &&
         ( (ov->opt[3]=='0')|| (ov->opt[3]=='1')|| (ov->opt[3]=='2')|| (ov->opt[3]=='3')|| (ov->opt[3]=='4')||
           (ov->opt[3]=='5')|| (ov->opt[3]=='6')|| (ov->opt[3]=='7')|| (ov->opt[3]=='8')|| (ov->opt[3]=='9')   ) ) {
         buff[0] = ov->opt[3]; buff[1] = '\0';
         n = atoi(buff);
         if (Nmcs_types != n) { printf("#ERROR:IMPROPER ORDER -mcs[x] ('%s %s').\n",ov->opt,ov->val); exit(1);}
         if ((n+1)>Nmcs_types){ Nmcs_types = n + 1;}
         Split_to_Words(ov->val,':',&Nword,Wsta,Wend,100); 
         Get_Part_Of_Line(buff,ov->val,Wsta[0],Wend[0]);  mcs_types[n].ConnectGraphType = buff[0];
         Get_Part_Of_Line(buff,ov->val,Wsta[1],Wend[1]);  mcs_types[n].maxDIFtopodis    = atoi(buff);
         Get_Part_Of_Line(buff,ov->val,Wsta[2],Wend[2]);  mcs_types[n].atomtype_class   = buff[0];
      }
      else if (strcmp(ov->opt,"chste2D")==0) { CheckStereo2D = ov->val[0];} 
      else if (strcmp(ov->opt,"chste3D")==0) { CheckStereo3D = ov->val[0];} 
      else if (strcmp(ov->opt,"chbang")==0) { CheckBondAngle = ov->val[0];} 
      else if (strcmp(ov->opt,"amste3D")==0) { AdjustMatchStereo3D = ov->val[0];} 
      else if (strcmp(ov->opt,"permsd")==0) { PermuRmsdMin = ov->val[0];} 
      else if (strcmp(ov->opt,"nameT")==0) {  sprintf(newname_molT,"%s",ov->val);}
      else if (strcmp(ov->opt,"newresT")==0) {
        Split_to_Words(ov->val,':',&Nword,Wsta,Wend,100); 
        Get_Part_Of_Line(new_resi_molT,ov->val,Wsta[0],Wend[0]);
        Get_Part_Of_Line(buff,ov->val,Wsta[1],Wend[1]); new_chain_molT = buff[0];
        Get_Part_Of_Line(new_rnum_molT,ov->val,Wsta[2],Wend[2]);
      }
      else if (strcmp(ov->opt,"sepoutT")==0) { SeparateOutT = ov->val[0]; }
      else if (strcmp(ov->opt,"kplcn")==0) { KeepLargestConGraph = ov->val[0]; }
      else { printf("#ERROR:Can't understand option '%s'.\n",ov->opt); exit(1);}
   }
 }



 srand(PAR.SeedRand);


 if ((PAR.ConnectGraphType=='T') || (PAR.ConnectGraphType=='t')){
   if (PAR.maxDIFtopodis <0) PAR.maxDIFtopodis = 0;
 }
  
 
/* Default Weighting for MCS */

 if ((PAR.WeightScheme=='D') && ((PAR.AlgoType=='B')||(PAR.AlgoType=='b'))){ 
     PAR.Wneiatm  = 1.0; 
     PAR.Wextcon  = 1.0; 
     PAR.Wdistan  = 0.0; 
     PAR.Wncompo  = 0.0;
     PAR.Wtopodis = 0.0; 
     if ((PAR.ConnectGraphType=='T')||(PAR.ConnectGraphType=='t')||(PAR.ConnectGraphType=='D')){ PAR.Wtopodis = 1.0;} 
 }


/* Setup options by EnergyStrategy(-E) and SearchStrategy (-T) */
  if ((EnergyStrategy=='A') && (SearchStrategy == 'F')){
    RigidFlexibleType   = 'F';
    DoMCS               = 'T'; 
    CnfStampingType     = 'T';
    CnfEselfclashMin    = 'F';
    RgdMatchFitType     = 'T';
    RgdPcaFitType       = 'F';
    SteepestDescentType = 'T'; 
    CnfRandInitType     = 'T'; 
    RgdRandInitType     = 'F'; 
    PAR.WEatommatch     = 1.0; 
    PAR.WEvolmovlap     = 0.0; 
  }
  else if ((EnergyStrategy=='A') && (SearchStrategy == 'R')){
    RigidFlexibleType   = 'R';
    DoMCS               = 'T'; 
    CnfStampingType     = 'F';
    CnfEselfclashMin    = 'F';
    RgdMatchFitType     = 'T';
    RgdPcaFitType       = 'F';
    SteepestDescentType = 'T'; 
    CnfRandInitType     = 'F'; 
    RgdRandInitType     = 'F'; 
    PAR.WEatommatch     = 1.0; 
    PAR.WEvolmovlap     = 0.0; 
  }
  else if ((EnergyStrategy=='V') && (SearchStrategy == 'F')){
    RigidFlexibleType   = 'F';
    DoMCS               = 'F'; 
    CnfStampingType     = 'F';
    CnfEselfclashMin    = 'T';
    RgdMatchFitType     = 'F';
    RgdPcaFitType       = 'T';
    SteepestDescentType = 'T'; 
    CnfRandInitType     = 'T'; 
    RgdRandInitType     = 'T'; 
    PAR.WEatommatch     = 0.0; 
    PAR.WEvolmovlap     = 1.0; 
  }
  else if ((EnergyStrategy=='V') && (SearchStrategy == 'R')){
    RigidFlexibleType   = 'R';
    DoMCS               = 'F'; 
    CnfStampingType     = 'F';
    CnfEselfclashMin    = 'F';
    RgdMatchFitType     = 'F';
    RgdPcaFitType       = 'T';
    SteepestDescentType = 'T'; 
    CnfRandInitType     = 'F'; 
    RgdRandInitType     = 'T'; 
    PAR.WEatommatch     = 0.0; 
    PAR.WEvolmovlap     = 1.0; 
  }


  else if (SearchStrategy=='N'){
    RigidFlexibleType   = 'R';
    /* DoMCS               = 'F'; */
    DoMCS               = 'T'; 
    CnfStampingType     = 'F';
    CnfEselfclashMin    = 'F';
    RgdMatchFitType     = 'F';
    RgdPcaFitType       = 'F';
    SteepestDescentType = 'F';
    CnfRandInitType     = 'F'; 
    RgdRandInitType     = 'F';
    NiterateGradFit = 0;
    NinitialGradFit = 1;
    NoutconfGradFit = 1;
    /* printf("#SearchStrategy %c\n",SearchStrategy); */
  }

  /* Read Options Again */
  ov = &(PAR.OptValHead);
  while (ov->next != NULL){
    ov = ov->next;
    if (Set_Global_PAR_from_OPTION_VALUE(ov)==0){
           if (strcmp(ov->opt,"mcs")==0)  { DoMCS  = ov->val[0];}
      else if (strcmp(ov->opt,"cfstp")==0)   { CnfStampingType  = ov->val[0];}
      else if (strcmp(ov->opt,"cfesc")==0)   { CnfEselfclashMin = ov->val[0];} 
      else if (strcmp(ov->opt,"rgmch")==0)   { RgdMatchFitType  = ov->val[0];}
      else if (strcmp(ov->opt,"rgpca")==0)   { RgdPcaFitType  = ov->val[0];}
      else if (strcmp(ov->opt,"SD")==0)      { SteepestDescentType  = ov->val[0];}
      else if (strcmp(ov->opt,"cfcnd")==0)   { CnfRandInitType = ov->val[0];} 
      else if (strcmp(ov->opt,"rgrnd")==0)   { RgdRandInitType = ov->val[0];} 
    }
  }

  printf("#CnfEselfclashMin %c CnfRandInitType %c\n",CnfEselfclashMin,CnfRandInitType);

  if (AdjustStereo2D=='T'){
    RgdMatchFitType = CnfStampRingType = ChiralSp3FitType = CnfStampingType = 'F'; SteepestDescentType = 'T';
    PAR.WEselfclash = 1.0; 
    PAR.WEatommatch =  PAR.WEprotclash =  PAR.WEprotatrct = PAR.WEvolmovlap = 0.0;
  }
  if (CheckStereo2D=='T'){
    RgdMatchFitType = CnfStampRingType = ChiralSp3FitType = CnfStampingType =  SteepestDescentType = 'F';
  }

  if (DoMCS=='F'){ RgdMatchFitType =  CnfStampingType =  CnfStampRingType = ChiralSp3FitType = 'F'; }

  Set_max_atomtype_by_atomtype_class();

  PAR.cal_topodis = 'T';

  m_model     = NULL;
  tanimoto    = 0.0;
  Natom_match = 0;

  if (NoutconfGradFit>NinitialGradFit){ NoutconfGradFit = NinitialGradFit;}
 
  /********************************/ 
  /** [1] Read molecule T(arget) **/
  /********************************/ 
  okT = Read_MOLECULE(molT.filename,&molT,molT.filetype,molT.atmhet,molT.BondType);
  if (okT==0) {
    printf("#ERROR: file for molT '%s' (type '%c') does not contain any atoms.\n",molT.filename,molT.filetype); 
    exit(1);}

  if (KeepLargestConGraph == 'T'){
    Set_Connected_Structure(&molT);
    Keep_Only_Largest_Connected_Structure(&molT);
  }


  if ((new_resi_molT[0]!='\0') && (new_rnum_molT[0]!='\0')){Set_Resi_Chain_Rnum_to_Molecule(&molT,new_resi_molT,new_chain_molT,new_rnum_molT);}

  printf("#molT Nheavyatom %d\n",molT.Nheavyatom);

  if (newname_molT[0] != '\0') sprintf(molT.name,"%s",newname_molT);

  Set_Topological_Distance(&molT,'H');  

  printf("#type atom_pair %c\n",PAR.type_atompair);

  if (opsfileT[0]!='\0') {Write_Molecule_in_PostScript(opsfileT,&molT); exit(1);}

  /***********************************/
  /** [2] Read molecule R(eference) **/
  /***********************************/
  okR = 0;
  if (molR.filename[0] != '\0'){
    okR = Read_MOLECULE(molR.filename,&molR,molR.filetype,molR.atmhet,molR.BondType);
  }
  else{
    printf("#ERROR: the reference molecule (-R) is necessary. If you don't have it, prepare the probe sphere file in the putative active site.\n");
    exit(1);
  }
  if (okR==0) {
    printf("#ERROR: file for molR '%s' (type '%c') does not contain any atoms.\n",molR.filename,molR.filetype); 
    exit(1);}

  if (KeepLargestConGraph == 'T'){
    Set_Connected_Structure(&molR);
    Keep_Only_Largest_Connected_Structure(&molR);
  }


  if ((PAR.GenEquivAtomPermu=='A')||(PAR.GenEquivAtomPermu=='N')||(PAR.GenEquivAtomPermu=='B')){
    Make_All_Permutation_of_Equivalent_Heavyatoms(&molT); 
    Make_All_Permutation_of_Equivalent_Heavyatoms(&molR); 
    printf("#PERMUTAION molT %d molR %d\n",molT.Npermu,molR.Npermu);
  }
  if (PAR.GenEquivAtomPermu=='a'){
    Make_All_Permutation_of_Equivalent_Heavyatoms(&molT); 
    printf("#PERMUTAION molT %d\n",molT.Npermu);
  }

  /****************************************************************/
  /** [3] Matching Two molecules (molT and molR) (if DoMCS=='T') **/
  /****************************************************************/
  printf("#DoMCS %c Nmcs_types %d\n",DoMCS, Nmcs_types);
 
  optMlist.next = NULL;

  if (DoMCS == 'T'){

   /* Do Multiple Type of MCSs */
    if (Nmcs_types>0){
      for (n=0;n<Nmcs_types;++n){
        printf("#MCS_TYPES[%d] con '%c' mtd %d at %c\n",n,mcs_types[n].ConnectGraphType,mcs_types[n].maxDIFtopodis,mcs_types[n].atomtype_class);
        optMlist0.next = NULL;
        PAR.ConnectGraphType = mcs_types[n].ConnectGraphType;
        PAR.maxDIFtopodis    = mcs_types[n].maxDIFtopodis;
        PAR.atomtype_class   = mcs_types[n].atomtype_class;
 
        PAR.Wtopodis = 0.0; 
        if ((PAR.ConnectGraphType=='T')||(PAR.ConnectGraphType=='t')||(PAR.ConnectGraphType=='D')){ PAR.Wtopodis = 1.0; }

        Set_max_atomtype_by_atomtype_class();
        Set_atomtype_of_Atoms(&molT,mcs_types[n].atomtype_class);
        Set_atomtype_of_Atoms(&molR,mcs_types[n].atomtype_class);
        Set_Nneighbor_atomtype(&molT);
        Set_Nneighbor_atomtype(&molR);
        Match_Two_Molecules(&molT,&molR,&optMlist0,PAR.AlgoType,mcs_types[n].ConnectGraphType,mcs_types[n].maxDIFtopodis,PAR.Nkeep,iamchfile);
        Write_MATCHlist_in_Horizontal_Element("-",&optMlist0,&molT,&molR);
        m_model = &optMlist0;
        while ((m_model->next != NULL)&&(m_model->next->nodetype!='E')){
          m_model  = m_model->next;
          m_model->atomtype_class = mcs_types[n].atomtype_class;
          Add_MATCH_to_MATCHlist(m_model,&optMlist);
        }
      }
       Sort_MATCHlist_by_Npair_and_select_dis(&optMlist); 
    }


   /* Do Single Type of MCS (default) */
    else if (Nmcs_types==0){
      Match_Two_Molecules(&molT,&molR,&optMlist,PAR.AlgoType,PAR.ConnectGraphType,PAR.maxDIFtopodis,PAR.Nkeep,iamchfile);
      Write_MATCHlist_in_Horizontal_Element("-",&optMlist,&molT,&molR);
    }
 
    if (PAR.GenEquivAtomPermu=='A'){Make_Permutated_MATCHlist_from_Top_MATCH(&optMlist,&molT,&molR);}
    if (PAR.GenEquivAtomPermu=='a'){Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation(&optMlist,&molT,&molR);}
    if (PAR.GenEquivAtomPermu=='N'){Make_NonRedundant_MATCHlist_by_Permutation(&optMlist,&molT,&molR);}
 
   if (CheckBondAngle=='T') {Check_Bond_Angle_Difference(&molT, &molR,optMlist.next);}

   if ((CheckStereo3D=='T')||(AdjustMatchStereo3D=='T')){
     Ndiff_chiral_sp3_init = Check_3D_sp3_Chirality(&molT, &molR, optMlist.next,'T');
   }

   if (AdjustMatchStereo3D=='T'){
     if (Ndiff_chiral_sp3_init>0){
       Make_All_Permutation_of_Equivalent_Heavyatoms(&molT);
       Make_All_Permutation_of_Equivalent_Heavyatoms(&molR);
       Ndiff_chiral_sp3_init = Modify_MATCH_For_3D_sp3_Chirality_with_Permutation(&molT, &molR, optMlist.next);
     }
   }

     Write_MATCHlist_in_Horizontal_Element("-",&optMlist,&molT,&molR);
  
    if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E')){
      Natom_match = optMlist.next->Npair;
      tanimoto = Tanimoto_Coefficient(Natom_match,&molT,&molR);
      Assign_Aligned_AtomNumber_to_tFactor(optMlist.next,&molT,'A');
      Assign_Aligned_AtomNumber_to_tFactor(optMlist.next,&molR,'B');
  
  
      if (orasfile[0] != '\0'){
       sprintf(ofname,"%s.ras",orasfile);
       Write_Superimposed_MATCH_in_rasmol_script(ofname,optMlist.next,&molT,&molR);
       sprintf(ofname,"%s-A.ras",orasfile);
       Write_MATCH_in_rasmol_script(ofname,optMlist.next,&molT,&molR,'A','M','L');
       sprintf(ofname,"%s-B.ras",orasfile);
       Write_MATCH_in_rasmol_script(ofname,optMlist.next,&molT,&molR,'B','M','L');
      }
  
     if (oamchfile[0] != '\0'){
       Write_MATCHlist(oamchfile,&optMlist,&molT,&molR,'O',RankAtomMatchFor3D,&HeadComment);
     }
   
  
    }
    else {Natom_match = 0; tanimoto = 0.0;}
  
    if (oAmchfile[0] != '\0') Write_MATCHlist(oAmchfile,&optMlist,&molT,&molR,'A',0,&HeadComment);
    if (opsfile[0] != '\0')   Write_MATCH_in_PostScript(opsfile,optMlist.next,&molT,&molR);
  
    if ((omcs_sdffile[0] != '\0') && (optMlist.next != NULL))
       Write_MATCH_in_SDF(omcs_sdffile,optMlist.next,&molT,&molR); 
  
    if ((optMlist.next != NULL) && (optMlist.next->nodetype != 'E')){ 
     /*
     printf("#initial_rmsd %f\n",rmsd_atommatch(&molT,&molR,optMlist.next));
     printf("#length of optMlist %d\n",Length_of_MATCHlist(&optMlist));
     */
     if ((RankAtomMatchFor3D>1) && (RankAtomMatchFor3D <= Length_of_MATCHlist(&optMlist))){
        m_model = optMlist.next;
        n = 1;
        while ((n<RankAtomMatchFor3D)&&(m_model->next != NULL)){
          m_model = m_model->next;
          n += 1;
        } 
        optMlist.next = NULL;
        Add_MATCH_to_MATCHlist(m_model,&optMlist);
     } 
     else{ m_model = optMlist.next;}
    }
 
   Nam = 0;
   m_model = &optMlist;
   while ((m_model->next != NULL)&&(m_model->next->nodetype!='E')){
     m_model = m_model->next;
     m_model->num = Nam;
     Nam += 1;   
  }

 
  } /* DoMCS=='T' */


   /* If DoMCS=='F' or no atom match is found, just superimpose Gcenters */ 
   if ((DoMCS=='F')||(optMlist.next == NULL) || (optMlist.next->nodetype == 'E')){ 
    Superimpose_Just_Two_Gcenters(&molT,&molR);
    rmsd = 0.0;
    m_model = (struct MATCH *)malloc(sizeof(struct MATCH));
    m_model->Npair = 0;
    m_model->num  = 0;
    optMlist.next = NULL;
    m_model->ConnectGraphType = PAR.ConnectGraphType;
    m_model->maxDIFtopodis    = PAR.maxDIFtopodis;
    m_model->atomtype_class   = PAR.atomtype_class;
    Add_MATCH_to_MATCHlist(m_model,&optMlist);
  }

  /************************************/ 
  /*** [5] Read Molecule P(rotein)  ***/
  /************************************/ 

  if (molP.filename[0] != '\0'){
    okP = Read_MOLECULE(molP.filename,&molP,molP.filetype,molP.atmhet,molP.BondType);
    if (molP.Natom==0){
      printf("#ERROR: file for moleculeR '%s' does not contain any atoms.\n",molP.filename);
      exit(1);
    } 
   printf("#read_receptor_file('%s') Natom %d\n",molP.filename,molP.Natom); 

   setup_variables_from_string_tpdisrest(&molT,&molP);
 }

   if ((PAR.string_tpdisrest[0][0] != '\0') ||(PAR.string_tpdisrest[1][0] != '\0')||(PAR.string_tpdisrest[2][0] != '\0')){
      if (molP.filename[0] == '\0'){
         printf("#ERROR: target-protein atom pairs option -tppair[x] requires the receptor protein ('-P' option)\n");
         exit(1);
      }
   }

   if (RigidFlexibleType=='F'){
     Set_Rotatable_Bond(&molT);  
     Set_Rotatable_Bond(&molR);  
     Set_TRANSFORM_VARIABLE(&traT,&molT,'R');
     if (PlaneDecision=='X'){
       if (molT.Natom > molT.Nheavyatom) {Set_Type_of_Rotatable_Bond(&traT,&molT,m_model,'2');}
                                    else {Set_Type_of_Rotatable_Bond(&traT,&molT,m_model,'3');}
     } 
     else{
       Set_Type_of_Rotatable_Bond(&traT,&molT,m_model,PlaneDecision);
     }
   }
   else{
     Set_TRANSFORM_VARIABLE(&traT,&molT,'-');
     traT.Nrbond = 0;
   } 

  if (SteepestDescentType=='T'){  
    molTtrial   = (struct MOLECULE*)malloc(sizeof(struct MOLECULE)*NinitialGradFit*MaxAtomMatchFor3D);
    traTtrial   = (struct TRANSFORM_VARIABLE*)malloc(sizeof(struct TRANSFORM_VARIABLE)*NinitialGradFit*MaxAtomMatchFor3D);
    sindex      = (int *)malloc(sizeof(int)*NinitialGradFit*MaxAtomMatchFor3D);
    for (n=0;n<NinitialGradFit*MaxAtomMatchFor3D;++n){
      molTtrial[n].atoms = NULL;
      molTtrial[n].bonds = NULL;
      molTtrial[n].Natom = 0;
      molTtrial[n].Nbond = 0;
      traTtrial[n].Natom = 0;
      sindex[n] = n;
    } 
 }

  Malloc_MOLECULE(&molTo,molT.Natom,molT.Nbond);
  Copy_MOLECULE(&molTo,&molT);

  /**************************************************/
  /*** [6] TRANSFORM MOLECULE BY EACH ATOM MATCH  ***/
  /**************************************************/
   m_model = &optMlist;
   Nam = 0; 
   Ntarget_conf = 0; 

   while ((m_model->next != NULL) && (m_model->next->nodetype != 'E') && (Nam < MaxAtomMatchFor3D)){
     m_model = m_model->next;
     Copy_MOLECULE(&molT,&molTo);
     if (m_model != NULL){ printf("#>ATOM_MATCH [%d]\n",m_model->num); } 
    if ((PermuRmsdMin=='T')||(PermuRmsdMin=='B')){
      if (PermuRmsdMin=='T'){
        Make_All_Permutation_of_Equivalent_Heavyatoms(&molT); 
        Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation(&optMlist,&molT,&molR); 
      }
      if (PermuRmsdMin=='B'){
        Make_All_Permutation_of_Equivalent_Heavyatoms(&molT); 
        Make_All_Permutation_of_Equivalent_Heavyatoms(&molR); 
       /* Make_Permutated_MATCHlist_from_Top_MATCH(&optMlist,&molT,&molR); */
      }

      m_model = chose_MATCH_with_min_rmsd(&molT,&molR,&optMlist);
    }

     /*** (6A) TRANSFORM MOLECULE BY NON-GRADIENT METHODS ***/

     if (RgdPcaFitType=='T'){
       Rotate_TargetMolecule_by_minimum_energy_PCA(&molT,&molR,&molP,m_model);
     }
  
     if (DoMCS=='T'){
        

       if (CheckRingConfType=='T'){Ndiff_ring_conf_init =  Check_Five_or_Six_Atom_Ring_Geometry(m_model,&molT,&molR); }
       if (CnfStampRingType=='T'){Stamp_Five_or_Six_Atom_Ring_Geometry(m_model,&molT,&molR); }
  
       if (ChiralSp3FitType=='T'){ Ndiff_chiral_sp3 = Stamp_sp3_Chirality(&molTo,&molR,m_model); }
  
       if (CnfStampingType=='T'){  
         traT.Nrbond_stamp = 0; 
         for (r=0;r<traT.Nrbond;++r){
            Stamp_RotAngle_for_Reference_Dihedral_Angle(r,&traT,m_model,&molT,&molR);
            if (traT.rbAngle[r]!=0.0){
              Transform_MOLECULE_by_RotBond(r,traT.rbAngle[r],&traT,&molT);
              traT.Nrbond_stamp += 1;
            }
            /* printf("rmsd_match %f\n",rmsd_atommatch(&molT,&molR,m_model)); */
         }
         printf("#Nrbond_stamp %d Nrbond_no_stamp %d\n",traT.Nrbond_stamp, traT.Nrbond - traT.Nrbond_stamp);
       }
  
        if ((RgdMatchFitType=='T')||(CnfStampRingType=='T')||(ChiralSp3FitType=='T')||(CnfStampingType=='T')){
         RMSD = Calculate_CRMS_MATCH_Quaternion(m_model,&molT,&molR,g1,g2,Rmat,"AonB");
         Rotate_Molecule(&molT,g1,g2,Rmat);
        }
     }
  
    if (molP.filename[0] != '\0'){
      n = mark_ligand_neighbor_receptor_atoms_for_Eprotclash(&molP,&molT,PAR.Dmax_clash*1.2);
      printf("#molP->Natom %d Natom_nei_molT %d\n",molP.Natom,n);
    }
  
     if (CheckStereo2D=='T'){
       if (PAR.AtomStereoFromBondStereo != 'F'){
         Set_Atomic_Stereo_Parity_from_2D_BondStereo(&molR,PAR.AtomStereoFromBondStereo);
       }
       /* Ndisagree = Check_sp3_Chirality_for_Stereo_Parity(&molT, &molR, m_model); */
       printf("CRASH_ATOM_PAIR %d\n",count_selfclash(&molT));
       if (count_selfclash(&molT)>0){ show_selfclash(&molT,"_T");}
       if (molT.osdffile[0] !='\0'){Write_SDF_Molecule(molT.osdffile,&molT,comment);}
       exit(1);
     }
  
     if (AdjustStereo2D=='T'){
       if (PAR.AtomStereoFromBondStereo!='F'){ 
         Set_Atomic_Stereo_Parity_from_2D_BondStereo(&molR,PAR.AtomStereoFromBondStereo);
       }
       Adjust_3D_Stereo_for_2D_Stereo_Reference(&molT, &molR, m_model, PyramidInversion,AllReflectChiral,FoldHexaRing);
       if (molT.osdffile[0] !='\0'){Write_SDF_Molecule(molT.osdffile,&molT,comment);}
       if (molT.opdbfile[0] !='\0'){Write_PDB_Molecule(molT.opdbfile,&molT,comment);}
       exit(1);

    } 
  
    traT.Etotal     = energy_total(&molT,&molR,&molP,m_model);
    traT.Eatommatch = energy_atommatch(&molT,&molR,m_model);
    traT.Eselfclash = energy_selfclash(&molT);
    traT.Eprotclash = energy_protclash(&molT,&molP);
    traT.Eprotatrct = energy_protatrct(&molT,&molP);
    traT.Evolmovlap = energy_volmovlap(&molT,&molR);
    traT.rmsd_match = fixed_rmsd_atommatch(&molT,&molR,m_model);
    traT.tanimoto_volume =  -energy_volmovlap(&molT,&molR)/(-energy_volmovlap(&molT,&molT)-energy_volmovlap(&molR,&molR) + energy_volmovlap(&molT,&molR));
    
  
    printf("#Etotal %f Eatommatch %f Eselfclash %f Eprotclash %f Eprotatrct %f rmsd_match %f\n",
     traT.Etotal, traT.Eatommatch, traT.Eselfclash,traT.Eprotclash, traT.Eprotatrct,traT.rmsd_match);
  
   
    /*** (6B) GRADIENT-BASED OPTIMIZATION ***/
    if (SteepestDescentType=='T'){  
       printf("#DO GRADIENT-BASED OPTIMIZATION\n");
  
       minNselfclash = -1;
  
       for (n=0;n<NinitialGradFit;++n){

         nn = Nam*NinitialGradFit + n;
         Malloc_MOLECULE(&(molTtrial[nn]),molT.Natom,molT.Nbond);
         Copy_MOLECULE(&(molTtrial[nn]),&molT);
         Malloc_TRANSFORM_VARIABLE(&(traTtrial[nn]),traT.Nrbond,traT.Natom);
         Copy_TRANSFORM_VARIABLE(&(traTtrial[nn]),&traT);
         traTtrial[nn].Ninit = nn;
         traTtrial[nn].num_atom_match = Nam;
         Malloc_MATCH(&(traTtrial[nn].atom_match),m_model->Npair); 
         Copy_MATCH(&(traTtrial[nn].atom_match),m_model); 
         Ntarget_conf += 1; 

         /** NOTE: if n==0, random transformation is not applied (the same as the input conformations) !! **/
         if (n>0){
           if (CnfRandInitType=='T'){
             Random_Transform_MOLECULE_by_RotBond(&traTtrial[nn],&(molTtrial[nn]),maxRotBondAngleInit,FixStampAngle,FixPlaneAngle);
             if (ProbFlipRingInit>0.0){Random_Change_Flip_Of_Fragments(&(molTtrial[nn]),ProbFlipRingInit);} 
             if (InitChiralChange=='T'){Random_Change_sp3_Chirality(&(molTtrial[nn]),0.5);}
           
             if (CnfEselfclashMin=='T'){ 
               steepest_descent_search_using_only_Eselfclash(NiterateGradFit,&traTtrial[nn],&(molTtrial[nn]),&molR,&molP,m_model,FixStampAngle,FixPlaneAngle);
             }

           }
  
           if ((RgdMatchFitType!='T')&&(RgdRandInitType=='T')){
             if (PAR.string_tpdisrest[0][0] != '\0'){ 
               Transform_MOLECULE_by_Random_3PointPairs_with_TPDISREST_restraint(&(traTtrial[nn]),&(molTtrial[nn]),&molR,&molP);
             }
             else{
               if ((1<=n)&&(n<=4)){
                 Rotate_TargetMolecule_by_PCA(&(molTtrial[nn]),&molR,&molP,n-1,m_model);
               }
               else{
                 Transform_MOLECULE_by_Random_3PointPairs_with_TPDISREST_restraint(&(traTtrial[nn]),&(molTtrial[nn]),&molR,&molP);
               }
             } 
           } 
         }
  
         mark_ligand_neighbor_receptor_atoms_for_Eprotclash(&molP,&molT,PAR.Dmax_clash*1.2);
         cal_energy_for_transform(&(traTtrial[nn]),&(molTtrial[nn]),&molR,&molP,m_model);
         traTtrial[nn].Etotal_init = traTtrial[nn].Etotal;
  
         
         for (k=0;k<NiterateGradFit;++k){
           if (RgdMatchFitType=='T'){
             RMSD = Calculate_CRMS_MATCH_Quaternion(m_model,&(molTtrial[nn]),&molR,g1,g2,Rmat,"AonB");
             Rotate_Molecule(&(molTtrial[nn]),g1,g2,Rmat);
           }
           cal_Force_and_Torque(&traTtrial[nn], &(molTtrial[nn]), &molR, &molP,m_model); 
           cal_dEd_rbAngle(&traTtrial[nn],&(molTtrial[nn]),&molR, &molP, m_model);
           mark_ligand_neighbor_receptor_atoms_for_Eprotclash(&molP,&molTtrial[nn],PAR.Dmax_clash*1.2);
           ok = steepest_descent_golden_linear_search(&traTtrial[nn], &(molTtrial[nn]), &molR, &molP, m_model, FixStampAngle,FixPlaneAngle);
           if (ok==0) k = 1000000;
           if (RgdMatchFitType=='T'){
             RMSD = Calculate_CRMS_MATCH_Quaternion(m_model,&(molTtrial[nn]),&molR,g1,g2,Rmat,"AonB");
             Rotate_Molecule(&(molTtrial[nn]),g1,g2,Rmat); 
           }
           /* printf("#STEP [%d] ene %10.5f\n",k,traTtrial[nn].Etotal); */
         } /* k */
 
         mark_ligand_neighbor_receptor_atoms_for_Eprotclash(&molP,&(molTtrial[nn]),PAR.Dmax_clash*1.2);
         cal_energy_for_transform(&(traTtrial[nn]),&(molTtrial[nn]),&molR,&molP,&(traTtrial[nn].atom_match));
         if ((minNselfclash<0)||(traTtrial[nn].Nselfclash<minNselfclash)){
            minNselfclash = traTtrial[nn].Nselfclash;
         }
         printf("#>trial[%3d] init_ene %10.5f --> fin_ene %10.5f\n",nn,traTtrial[nn].Etotal_init,traTtrial[nn].Etotal);
  
       } /* n */

     } /* end of SteepestDescentType=='T' */
  
     Nam += 1; 

  } /* m_model */


 /*** [7] Procedures for ending  (Sort molTtrial[] and choose the best one as molT)***/

  if (SteepestDescentType=='T'){
    printf("#Ntarget_conf %d\n",Ntarget_conf);
    Bubble_Sort_Array_of_TRANSFORM_VARIABLE(Ntarget_conf, sindex, traTtrial, SortValueForTargets);
    for (ns=0;ns<Ntarget_conf;++ns){
      n = sindex[ns];
      if (molTtrial[n].Natom>0){
      printf("#[rank %2d] num_atom_match %2d ninit %2d MCS %c%d%c Etotal %+7.2f Eam %+7.2f Esc %+7.2f Epc %+7.2f Epa %7.2f Evo %+6.2f Etp %+6.2f tani_vol %5.3f rmsd_ref %5.3f rmsd_rank1 %5.3f",
       ns+1,traTtrial[n].num_atom_match,n, 
       traTtrial[n].atom_match.ConnectGraphType, traTtrial[n].atom_match.maxDIFtopodis, traTtrial[n].atom_match.atomtype_class,
       traTtrial[n].Etotal, traTtrial[n].Eatommatch, traTtrial[n].Eselfclash, traTtrial[n].Eprotclash,traTtrial[n].Eprotatrct,
       traTtrial[n].Evolmovlap, traTtrial[n].Etpdisrest, traTtrial[n].tanimoto_volume,
       traTtrial[n].rmsd_match,CRMS_Same_Molecules_without_Superimpose(&molTtrial[sindex[0]],&molTtrial[n]));
      }
      for (t=0;t<3;++t){
        mk_short_string_status_tpdisrest(buff,t,&(molTtrial[n]),&molP);
        printf(" %s",buff);
      } 
      printf("\n");
    }

    /** special care for 'AdjustStereoParity' (if no clash, accept candidate 0 (non-randomized initial) **/
    if ((AdjustStereo2D=='T')&&(count_selfclash(&(molTtrial[0]))==0)){
      printf("#conf 0 is zero clash : %d Etotal %f\n",count_selfclash(&(molTtrial[0])),traTtrial[0].Etotal);
      n_min = 0;
    }
 
    n_min = sindex[0];
    printf("# The best conformation [%d] ene_min %f Nselfclash %d\n",n_min,traTtrial[n_min].Etotal,traTtrial[n_min].Nselfclash); 
 
    Copy_MOLECULE_XYZ(&molT,&(molTtrial[n_min]));
    Copy_TRANSFORM_VARIABLE(&traT,&(traTtrial[n_min]));
    mark_ligand_neighbor_receptor_atoms_for_Eprotclash(&molP,&molT,PAR.Dmax_clash*1.2);
    show_energy_total(&molT,&molR,&molP,&(traTtrial[n_min].atom_match));
 }



  sprintf(buff,"Etotal     %f",traT.Etotal);     Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"Eatommatch %f",traT.Eatommatch); Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"Eselfclash %f",traT.Eselfclash); Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"Eprotclash %f",traT.Eprotclash); Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"Eprotatrct %f",traT.Eprotatrct); Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"Evolmovlap %f",traT.Evolmovlap); Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"rmsd_match %f",traT.rmsd_match); Add_string_to_LINENODEs(buff,&HeadComment);
  sprintf(buff,"%d",count_selfclash(&molT)); Add_key_value_string_to_KEYVALUEs("CRASH_ATOM_PAIR",buff,&(molT.HeadProperty));


  /*** [8] WRITE TRANSFORMED MOLT ***/
  n = mark_ligand_neighbor_receptor_atoms_for_Eprotclash(&molP,&molT,PAR.Dmax_clash*1.2);
  show_protclash(&molT,&molP);
  m_model = optMlist.next;
  show_energy_total(&molT,&molR,&molP,m_model);

  mk_string_status_tpdisrest(buff,&molT,&molP); printf("%s\n",buff);


 if ((DoMCS=='T')&&(SteepestDescentType=='T')){
   sprintf(comment,"[rank 1] ninit %d mcs %c%d%c Nmatch %2d tani %.4f Etotal %+5.2f Eam %+5.2f Esc %+5.2f Epc %+5.2f Epa %+5.2f Evo %+5.2f rmsd_match %5.3f tani_vol %.4f",
     n_min, traTtrial[n_min].atom_match.ConnectGraphType, traTtrial[n_min].atom_match.maxDIFtopodis, 
     traTtrial[n_min].atom_match.atomtype_class, m_model->Npair, Tanimoto_Coefficient(m_model->Npair,&molT,&molR),
     traTtrial[n_min].Etotal, traTtrial[n_min].Eatommatch, traTtrial[n_min].Eselfclash, traTtrial[n_min].Eprotclash, traTtrial[n_min].Eprotatrct,traTtrial[n_min].Evolmovlap,
     traTtrial[n_min].rmsd_match,traTtrial[n_min].tanimoto_volume
   );
 }

 printf("#COMMENT %s\n",comment); 

 if ((CheckStereo3D=='T')||(AdjustMatchStereo3D=='T')){
   Ndiff_chiral_sp3 = Check_3D_sp3_Chirality(&molT, &molR,m_model,'T');
 }
 if ((CheckRingConfType=='T')||(CnfStampRingType=='T')){
   Ndiff_ring_conf =  Check_Five_or_Six_Atom_Ring_Geometry(m_model,&molT,&molR);
 }

 if (KeyValueOut=='T'){
   printf("NHEAVYATOM_T %d\n",molT.Nheavyatom);
   printf("NHEAVYATOM_R %d\n",molR.Nheavyatom);
   printf("NROTBOND_T %d\n",traT.Nrbond);
   printf("NATOM_MATCH  %d\n",Natom_match);
   printf("TANIMOTO     %f\n",tanimoto);
   printf("Etotal     %f\n",traT.Etotal);
   printf("Eatommatch %f\n",traT.Eatommatch);
   printf("Eselfclash %f\n",traT.Eselfclash);
   printf("Eprotclash %f\n",traT.Eprotclash);
   printf("Eprotatrct %f\n",traT.Eprotatrct);
   printf("Evolmovlap %f\n",traT.Evolmovlap);
   printf("RMSD_MATCH %f\n",traT.rmsd_match);
   printf("TANIMOTO_VOL %f\n",traT.tanimoto_volume);
   if ((CheckStereo3D=='T')||(AdjustMatchStereo3D=='T')) {
      printf("NDIFF_CHIRAL_SP3_INIT %d\n",Ndiff_chiral_sp3_init);
      printf("NDIFF_CHIRAL_SP3      %d\n",Ndiff_chiral_sp3);
   }
   if ((CheckRingConfType=='T')||(CnfStampRingType=='T')) {
      printf("NDIFF_RING_CONF_INIT %d\n",Ndiff_ring_conf_init);
      printf("NDIFF_RING_CONF      %d\n",Ndiff_ring_conf);
   }

 }

 if (molT.osdffile[0] !='\0'){
   Write_SDF_Molecule(molT.osdffile,&molT,comment);
 }

 if ((molT.opdbfile[0] !='\0')||(opdbfileTP[0]!='\0')){ 
   if (SteepestDescentType=='T'){
     for (ns=0;ns<NoutconfGradFit;++ns){
       n = sindex[ns];
       sprintf(comment,"[rank %d] ninit %d mcs %c%d%c Nmatch %2d tani %.4f Etotal %+5.2f Eam %+5.2f Esc %+5.2f Epc %+5.2f Epa %+5.2f Evo %+5.2f Etp %+5.2f rmsd_match %5.3f tani_vol %.4f rmsd_rank1 %5.3f",
         ns+1,n, traTtrial[n].atom_match.ConnectGraphType, traTtrial[n].atom_match.maxDIFtopodis, traTtrial[n].atom_match.atomtype_class,
         m_model->Npair, Tanimoto_Coefficient(m_model->Npair,&molT,&molR),
         traTtrial[n].Etotal, traTtrial[n].Eatommatch, traTtrial[n].Eselfclash, traTtrial[n].Eprotclash,traTtrial[n].Eprotatrct,
         traTtrial[n].Evolmovlap, traTtrial[n].Etpdisrest,
         traTtrial[n].rmsd_match, traTtrial[n].tanimoto_volume, CRMS_Same_Molecules_without_Superimpose(&molTtrial[sindex[0]],&molTtrial[n]));
       if (PAR.string_tpdisrest[0][0]!='\0'){
         mk_string_status_tpdisrest(buff,&(molTtrial[n]),&molP); strcat(comment," "); strcat(comment,buff);
       }
       if ((SeparateOutT=='T')&&(molT.opdbfile[0] != '\0')) {
         make_multiple_conf_fname_from_single_conf_fname(ofname,molT.opdbfile,ns+1,NoutconfGradFit);
         Write_PDB_Molecule(ofname,&molTtrial[n],'w','A',ns+1,-1,'-','-',comment); 
       }
       else{         
         if (ns==0) {wmode = 'w';} else {wmode = 'a';}
         if (molT.opdbfile[0] != '\0') {Write_PDB_Molecule(molT.opdbfile,&molTtrial[n],wmode,'A',ns+1,-1,'-','-',comment); }
         if (opdbfileTP[0] != '\0') {Write_PDB_Molecule(opdbfileTP,&molTtrial[n],wmode,'A',ns+1,-1,'-','-',comment); }
       }
     }
   }
 }

 if (molT.omol2file[0]!='\0'){
    for (ns=0;ns<NoutconfGradFit;++ns){
        n = sindex[ns];
        if (SeparateOutT=='T'){
          make_multiple_conf_fname_from_single_conf_fname(ofname,molT.omol2file,ns+1,NoutconfGradFit);
          Write_MOL2_Molecule(ofname,&molTtrial[n],'w');
        }
        else{
          if (ns==0) wmode = 'w'; else wmode = 'a';
          Write_MOL2_Molecule(molT.omol2file,&molTtrial[n],wmode);
        }
    }
 } 

 
 if (molR.opdbfile[0] !='\0'){Write_PDB_Molecule(molR.opdbfile,&molR,'w','-',-1,-1,'-','-',comment); }

 if (opdbfileTR[0]!='\0'){ 
   Write_PDB_Molecule(opdbfileTR,&molT,'w','T',1,1,           '-','-',""); 
   Write_PDB_Molecule(opdbfileTR,&molR,'a','R',2,molT.Natom+1,'-','-',""); 
 }
 if (opdbfileTP[0]!='\0'){ 
   /* Write_PDB_Molecule(opdbfileTP,&molT,'w','T',-1,1,           '-','-',"");  */
   if (molP.filename[0]!='\0') Write_PDB_Molecule(opdbfileTP,&molP,'a','P',3,molT.Natom + molR.Natom+1,'-','-',""); 
 }
 if (opdbfileTRP[0]!='\0'){ 
   Write_PDB_Molecule(opdbfileTRP,&molT,'w','T',1,1,'-','-',""); 
   Write_PDB_Molecule(opdbfileTRP,&molR,'a','R',2,molT.Natom+1,'-','-',""); 
   if (molP.filename[0]!='\0') Write_PDB_Molecule(opdbfileTRP,&molP,'a','P',3,molT.Natom + molR.Natom+1,'-','-',""); 
 }

   if (oamchfile[0] != '\0'){
     Write_MATCHlist(oamchfile,&optMlist,&molT,&molR,'O',RankAtomMatchFor3D,&HeadComment); 
   }

   if (oAmchfile[0] != '\0'){
    Write_MATCHlist(oAmchfile,&optMlist,&molT,&molR,'A',0,&HeadComment);
   }

  if (otransfile[0] != '\0'){Write_TRANSFORM_VARIABLE(otransfile,&traT,&molT);}

  if (count_selfclash(&molT)>0){ show_selfclash(&molT,"_T");}

  
  Free_MOLECULE(&molTo);
  
  return(1);

} /* end of main() */ 





