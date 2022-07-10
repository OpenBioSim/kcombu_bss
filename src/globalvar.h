/*
 
 <globalvar.h>
 
  Definition of global parameters 
 
*/

#define LAST_MOD_DATE     "2022/06/10"
#define CODED_BY          "coded by Takeshi Kawabata."
#define KCOMBU_STANDS_FOR "KCOMBU stands for 'K'emical structure 'COM'parison using 'B'uild-'U'p algorithm."
#define MAX_PROPERTY  30 
#define MAX_LENGTH_LINE 1024 
#define MAX_FILENAME     512 

#include "ReadOpt.h"

/** PARAMETERS: for globally defined important parameters **/
struct PARAMETERS{
  /** BASIC PARAMETERS **/ 
  char   COMMAND[MAX_LENGTH_LINE];
  struct OPTION_VALUE OptValHead;
  char   COMMENT[MAX_LENGTH_LINE];
  char   PROGRAM_TYPE;     /* 'p':pkcombu, 'l':lkcombu, 'd':dkcombu, 'f':fkcombu, 'c':ckcombu */ 
 
  /** Date and Computation time **/ 
  char   START_DATE[MAX_FILENAME];  /* Start Date */
  char   END_DATE[MAX_FILENAME];    /* End   Date */
  double START_TIME_SEC;   /* Date on starting time ( in second) */
  double COMP_TIME_SEC_PROGRESS;  /* Computation time in progress (in second); */
  double COMP_TIME_SEC;           /* Final Computation time (in second); */
  double MAX_COMP_TIME_SEC;/* Maximum Computational time (in second) if negative, don't care */

  /*** IO of molecular property ***/
  char   READ_MOLECULAR_PROPERTY; /* 'T':read molecular property stored in SDF or MOL2 files */ 
  int    Nproperty;               /* Number of focused properties */
  char   PROP_NAME[MAX_PROPERTY][64]; /* name of i-th property */ 
  int    PROP_NUM[MAX_PROPERTY];      /* number of i-th property */
  float  PROP_MIN[MAX_PROPERTY];      /* minimum value  for i-th property */
  float  PROP_MAX[MAX_PROPERTY];      /* maximum value  for i-th property */
  char   PROP_STR[MAX_PROPERTY][64];  /* string pattern for i-th property */
  char   CHANGE_KEY_SPACE_TO_UNDERBAR;  /* change ' ' to '_' in the key string */
  char   CHANGE_NOTE_TO_SUPPLIERS;     /* change 'NOTES' into 'SUPPLIERS' */

  /** for MCS style description **/
  char   AlgoType;         /* Algorithm type. 'B'uild-up, 'X':exact */
  char   atomtype_class;   /* atomtype classification */ 
  char   bondtype_class;   /* bondtype classfication. */
  char   DAMRLassign;      /* 'H'ydrogen(2d), '3'D, 'X':pdb->3D, others->hydrogen */ 
  int    max_atomtype;     /* maximum number of atomtype  (< MAX_ATOMTYPE defined in 'molecule.h') */ 
  char   ConnectGraphType; /* Connection of MCS. 'C'onnected, 'D'isconnected, 'T'D-MCS,'t'C-MCS */
  int    maxDIFtopodis;    /* accepted maximum difference of topological distance (number of bonds in the shortest path) */ 
  int    maxNcomponent;    /* accepted maximum number of connected component. -1: not assign maximum number */
  int    maxNatom_in_ring; /* maximum number of atom in one ring (default is 20). */
  int    maxNatom_in_aromatic_ring;  /* maximum number of atom in one 'aromatic' ring (default is 8). */
  char   includeHinMCS;    /* include Hydrogen atom matching for MCS (T or F) */ 

  /** for Build-up calculation **/
  int    Nkeep;            /* maximum number of keeping atom correspondences */
  char   WeightScheme;     /* weighting scheme. 'D'efault weighting. 'S'pecified by user */
  float  Wneiatm;          /* Weight for diff. of neighboring atoms    (for 2D) */
  float  Wextcon;          /* Weight for diff .of extended connectivity (for 2D) */
  float  Wtopodis;         /* Weight for diff .of topological distance  (for 2D) */
  float  Wrank;            /* Weight for rank by unique extended connectivity (for 2D) */
  float  Wncompo;          /* Weight for number of connected component  (for 2D) */
  float  Wdistan;          /* Weight for diff .of distance                        (for 3D)*/
  char   levelEC_Dextcon;  /* Level of EC for scoring. '0' EC0, '1':EC1, '2':EC2 for Dextcon */ 
  int    maxD_Dtopodis;    /* maximum topological distance for Dtopodis */ 
  char   NonRedSelectDis;  /* Non Redandancy Check by Selection Distance Score. ('T' or 'F') */
 
 /** for molecule descriptor */ 
  int    max_separation;          /* maximum separation for atompair_descriptor[] */
  int    max_atompair_descriptor; /* maximum length of atompair_descriptor[] */
  char   type_atompair;           /* 'T':ringblock-dependent, otherwise:not dependent */ 
  int    ring_atomtype;           /* number of ring atom types */ 
  int    max_ring_size;           /* maximum Natom of ring for descriptors [] */
  int    max_block_size;          /* maximum Natom of block for descriptors [] */
  int    max_ring_descriptor;     /* maximum length of ring_descriptor[] */
 
  int    max_separation_cyclic;   /* maximum separation for atompair_descriptor[] only for cyclic bond */
  int    max_separation_acyclic;  /* maximum separation for atompair_descriptor[] only for acyclic bond */
  char   simtype_descriptor;      /* 'T'animoto 'C'ahart */
 
 /** for obtaining MCS for various ways*/ 
  char   GenEquivAtomPermu;     /* Equivalent Atom Permutation. 'F':false, 'A':make all, 'N'on_redundant */
  int    BKmaxNpair;            /* maximum number of atom pairs for Bron-Kebosch algorithm */
  char   calc_finish;           /* 'N'ot start caculation, 'B'reak before finished, 'F'inished. */ 

 /** for flexible transforming molecule */
  float tole_dis_clash;         /* tolarant distance for clash */ 
  float tole_dis_atrct;         /* tolarant distance for attraction */ 
  float alpha_sigmoid_clash;    /* Alpha parameter for sigmoidal clash energy function  */ 
  float WEatommatch;            /* Weight for Eatommatch */
  float WEselfclash;            /* Weight for Eselfclash */
  float WEprotclash;            /* Weight for Eprotclash (clash     with receptor protein) */
  float WEprotatrct;            /* Weight for Eprotatrct (atraction with receptor protein) */
  float WEvolmovlap;            /* Weight for Evolmovlap (volume overlap) */
  float WEtpdisrest;            /* Weight for Etpdisrest (target-protein distance restraint) */
  float Dmax_clash;             /* maximum distance considering for Eselfclash and Eprotclash */
  float maxTRA;                 /* maximum step for translational vector dT  (angstrom) */
  float maxROT;                 /* maximum step for rotational vector    dW (degree)    */
  float maxRB;                  /* maximum step for angle for rotational bond (degree) */
  float p_volmovlap;            /* parameter p for volmovlap */ 
  char  atomtype_volmovlap;     /* atom type for Evolmovlap. 'I'gnore,'C'are  */ 
  char  HydrogenOptimize;       /* Hydrogen-position Optimize. 'T' or 'F'  */  
  float minAngleFlipOfFragment; /* minimum angle for CA and DB for FlipOfFragment */
  char  string_tpdisrest[3][100]; /* [Tarom_num_tpdisrest]:[Patom_num_tpdisrest]:[Dupper_tpdisrest] */ 
  int   Tatom_num_tpdisrest[3];   /* Target  Atom number (0,1,2,... ) for Etpdistrest */
  int   Patom_num_tpdisrest[3];  /* Protein Atom number (0,1,2,... ) for Etpdistrest */
  float Dupper_tpdisrest[3];     /* Dupper for Etpdistrest */


  /** other parameters */ 
  int    SeedRand;         /* Seed number for rand() */ 
  char   cal_topodis;      /* if 'T', calculate topological distance      */
  float  Tole_distan;      /* Torelance for distance difference (for 3D)  */
  int    TotalNatompair;   /* Total Number of possible atomic pair        */
  char   OutCalProcess;    /* if 'T', print calculation process to stdout */
  char   molnamekey[64];   /* key for molecular name, such as 'NAMIKI_ID' or 'LIGANDBOX_ID'.   */

  /** for permutation **/
  int   maxPermutations;    /* maximum number of calculating permutations at the mol->permuhead. */

  /* for PostScript output */
  char ps_label;         /* 'M'atched_atom_number, 'T':atomtype, 'N':atom number, 'A':atom name,'E':element */
  char ps_circle;        /* 'A'll, 'M'arked atoms */
  char ps_hydro;         /* Hydrogen atoms. 'T':show, 'F':not show */
  char ps_planeB;        /* rotate molA on PCA plane ('T' or 'F') */ 
  char ps_supAonB;       /* superimpose molA on molB ('T' or 'F') */ 
  char ps_color;         /* Color scheme ('M'onochrome, 'C'pk )   */
  char ps_circ_outline;  /* Draw Circle Outline ('T' or 'F')      */
  char ps_care_bondtype; /* Care bond type ('T' or 'F') (single,double,triple) */
  /* for 2D depiction */
  float LengthOfBond2D;     /* Length of bond for 2D depiction */
  char  CollisionResolve2D; /* Do CollisionResolving 'T' or 'F' */
  /** for SMILES **/
  char OutRankForSMILES;  /* if 'T', output rasmol script for SMILES rank */
  
  /** for chiralities **/
  char AtomStereoFromBondStereo;


 /** parameters for many purposes */ 
  int    subint;           /* sub integer for many purposes */
  char   subchar;          /* sub character for many purposes */
  int    M;            /* count for many purpose */
  int    N;            /* count for many purpose */
  int    O;            /* count for many purpose */
  int    P;            /* count for many purpose */

};

/** GLOBAL VARIABLES **/

extern struct PARAMETERS PAR;

