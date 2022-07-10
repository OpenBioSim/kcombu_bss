/*

 <molecule.h> 

 definition of molecules.
 
*/ 


#define MAX_ATOMTYPE         68
#define MAX_NELEMENT         87 
#define MAX_MOL_FILENAME    512


struct ATOM{
  int   num;          /* Atom number (0,1,2,... ) */ 
  int   num_in_file;  /* Atom number specified in file. (Normally >=1. not always starts with 1.) */ 
  int   num_heavy;    /* Atom number (0,1,2,... ) only for heavy atoms (-1 for hydrogen) */ 
  float Pos[3];       /* X,Y,Z coordinates of ATOM */
  char  element[3];   /* chemical element, such as 'C','O','Zn','Fe' */
  char  one_char_ele;  /* one character element. If element has two char (such as Zn, Fe, ..) one_char_ele == 'X' */  
  char  atomtype[6];   /* atomtype for matching. hydrogen has '' atomtype. */ 
  char  atomname[6];   /* atomname specified in the file (length=4) (only for PDB) */ 
  char  resi[4];       /* residue name     (only for PDB) */
  char  rnum[6];       /* residue number   (only for PDB) */
  char  chain;         /* chain identifier (only for PDB) */ 
  char  atmhet;        /* 'A'tom,'H'etatm  (only for PDB) */
  char  ring;          /* 'r'ing, 'p'lane-ring, '-':not-ring */
  char  ringblock;     /* char to show block. 'a','b','c',...  */
  char  aromatic;      /* 'A':in the 3,4,5-member ring, 'a': in the ring with not more than 8 atoms, '-':otherwise */
  int   constr_num;    /* connected structure number (0,1,2,3,...). */
  char  DAMRL;         /* 'D'oner, 'A'cceptor, 'M'ixed donor and acceptor, a'R'omatic, a'L'yphatic */ 
  char  orbital;       /* '2':sp2 geometry (trigonal) '3':sp3 geometry (tetrahedral) */ 
  float tFactor;       /* temperature factor (only for PDB) */
  /** for SDF file **/
  int  mass_diff;     /* (-3, -2, -1, 0, 1, 2, 3, 4) Difference from mass in periodic table. */
  char char_charge;    /* '0' = uncharged or value other than these, '1' = +3, '2' = +2, '3' = +1, 
                          '4' = doublet radical, '5' = -1, '6' = -2, '7' = -3 */
  char  stereo_parity; /* '0':not stereo, '1':odd, '2':even, '3':either or unmarked stereo center */
  /** for MOL2 file **/
  float charge;           /* value for partial charge */
  char  mol2atomtype[6];  /* atom type in mol2 format such as "C2.2","C2.ar","O.co2"*/ 
  float radius;        /* van der Waals radius */  
 /**  variables for connection **/
  int   Nneighbor;         /* Number of connected atoms       */ 
  int   Nnei_heavy;        /* Number of connected heavy atoms */ 
  int   *neighbors;        /* array of Atom Number(num) of neighbor_atoms (malloc later)*/ 
  /**  variables for select_distance **/
  int   NneiC,NneiN,NneiO,NneiS,NneiP,NneiH;  /* Number of neighbor for C,N,O,S,P,H atoms */
  int   Nnei_atomtype[MAX_ATOMTYPE];          /* Number of neighbor for each atomtype */ 
  int   EC0,EC1,EC2,EC3,EC4;   /* Number of extended connectivity (only for heavy atoms) */ 
  int   rank;                  /* atom rank using unique Extended Connectivity [0,1,2,...] */ 
  /**  other variables **/
  int    Nmatch;        /* Number of matched other molecule for this atom */
  int    subnum;         /* sub number  for various purpose */
  float  subfloat;       /* sub float number  for various purpose */
  int    mark;           /* mark for various purpose */
};


struct BOND{
  int   num;        /* Number of bonds (0,1,2,... ) */ 
  int  anum1;       /* connected atom number 1  (0,1,2,..) */ 
  int  anum2;       /* connected atom number 2  (0,1,2,..) */ 
  char bondtype;    /*  '1':single,'2':double '3':triple 'r':aromatic */
  char bondstereo;  /*'0' not stereo (default), '1':Up, '4':Either, '6'Down */
};



struct RING{
/* this structure is also used for the "block" of rings */
  int   num;                   /* Number of rings (0,1,2,... ) */ 
  int   Natom;                 /* Number of member ATOMS */
  int   *num_atoms;            /* [PAR.maxNatom_in_ring] :array of numbers of member ATOMS */
                               /* for SSSR-ring, num_atoms[] is 'ring'-order, but for ring-block, not 'ring'-order.*/ 
  char  planarity;             /* 'p':almost plane, 'n':ot plane */
  char  aromatic;              /* 'A'romatic */
  float cpos[3];               /* center position of the plane */ 
  float nvec[3];               /*  normal vector  for the plane */ 
  float rg;                   /* radius of gyration */
  struct RING *prev, *next;    /* variables for linked list of RINGs */
};



struct CHAR2DMAP{
 int  N;              /* Number of grids */
 char **map; /* map[N][N] (malloc later) */
 char  malloced;      /* 0:not_malloced, 1:already_malloced */
 char  calculated;    /* 0:not_calculated, 1:already_calculated */
};


struct UNSIGNED_CHAR2DMAP{
 int  N;              /* Number of grids */
 unsigned char **map; /* map[N][N] (malloc later) */
 char  malloced;      /* 0:not_malloced, 1:already_malloced */
 char  calculated;    /* 0:not_calculated, 1:already_calculated */
};


struct FLOAT2DMAP{
 int  N;          /* Number of grids */
 float     **map; /* map[N][N] (malloc later) */
 char  malloced;  /* 0:not_malloced, 1:already_malloced */
 char  calculated;    /* 0:not_calculated, 1:already_calculated */
};


struct INT2DMAP{
 int  N;          /* Number of grids */
 int   **map; /* map[N][N] (malloc later) */
 char  malloced;  /* 0:not_malloced, 1:already_malloced */
 char  calculated;    /* 0:not_calculated, 1:already_calculated */
};


struct PERMUTATION{
  int   N;              /* length of index[] (usuary, mol->Nheavyatom) */
  int   *index_heavy;   /* index_heavy[0...N-1] : heavy atom number for permutated atom */
  struct PERMUTATION *next, *prev;
};


struct KEYVALUE{
  int  num;    /* number of KEYVALUE (0,1,2,....) */
  char *key;   /* string of 'key'   (malloc later) */
  char *value; /* string of 'value' (malloc later) */
  struct KEYVALUE *next, *prev;
};


struct MOLECULE{
  char   name[MAX_MOL_FILENAME];           /* name of molecule */ 
  char   core_molname[MAX_MOL_FILENAME]; 
    /* 
     If this molecule is 'one-molecule-from-one-file'.   
       then core_molname[] is the last part of filename. ex) filename:'SDF/ATP.sdf' --> core_filename:'ATP.sdf'.
     Otherwise,(if this molecule is 'many-molecule-from-one-file')  
       then core_molname[] is the same as name[].
   */ 
  char   filetype;                         /* 'P'db, 'S'df, '2':mol2, 'K'cf */ 
  char   filename[MAX_MOL_FILENAME];       /* filename */ 
  char   osdffile[MAX_MOL_FILENAME];       /* output file in sdf */ 
  char   opdbfile[MAX_MOL_FILENAME];       /* output file in pdb */ 
  char   omol2file[MAX_MOL_FILENAME];      /* output file in mol2 */ 
  char   comment[MAX_MOL_FILENAME];        /* comment of molecule (2nd line of SDF file)*/
  struct ATOM      *atoms;        /* array of ATOMs (0,1,2,...,Natom-1) (malloc later) */
  struct BOND      *bonds;        /* array of BONDs (0,1,2,,,,.Nbond-1) (malloc later) */
  struct RING      HeadRing;      /* Head of linked list of RINGs (SSSR) */
  struct RING      HeadBlock;     /* Head of linked list of the block of RINGs. One Block may include multiple SSSR rings. */
  int    Natom;                   /* Number of atoms */
  int    maxNatom;                /* predefined maximum Number of atoms. If Natom > maxNatom, the program stops. */
  int    Nbond;                   /* Number of bonds */
  char   chain;                   /* chain identifier (only for PDB)  */ 
  char   resi[4];                 /* residue name     (only for PDB)  */
  int    Nheavyatom;              /* Number of heavy atoms */
  int    Nheavybond;              /* Number of bonds between heavy atoms */
  int    *num_frm_heavy;          /* num    of ATOM = num_frm_heavy[num_heavy] (malloc later) */
  int    *heavy_frm_num;          /* num_heavy of ATOM = heavy_frm_num[num]    (malloc later) */
  struct CHAR2DMAP  conmap;       /* connection map. map value is character.
                        DEFAULT SCHEME   '0':not bonded, '1':single-bond, '2':double-bond. 'r':aromatic bond 
                        ROTATABLE SCHEME '0':not bonded, 'R':rotatable bond, 'N':not rotatable bond  */
  struct FLOAT2DMAP dismap;       /* distance map     */
  char   molform[MAX_MOL_FILENAME];            /* Molecular Formula */
  char   chiral_flag;             /* 0=not chiral, 1=chiral (only for SDF) */ 
  int    Nelement[MAX_NELEMENT];         /* Number of elements */
  unsigned char   oneatom_descriptor[MAX_ATOMTYPE]; /* Natomtype. Number of atomtype class (excluding Hydrogens) */ 
  int    Nring;                   /* Number of ring */
  int    Nringblock;              /* Number of ring block */
  int    Nconstr;                 /* Number of connected structure */
  int    Nmatch;                  /* Number of matched other molecule for this molecule */
  float  Gcen[3];                 /* Position of Gravity of Center */
  float  eigval[3];               /* Eigen Values  for Covariance Matrix (eigval[0]>eival[1]>eigval[2])*/
  float  eigvec[3][3];            /* Eigen Vectors for Covariance Matrix */
  char   atmhet;                  /* read only 'A'tom, only 'H'etam, 'B'oth */
  char   BondType;                /* 'B': guess bonds by distances and cal topo-distances. Otherwise, not consider bonds  */
  struct INT2DMAP  topodismap;  /* topological distance map (number of bonds on the shortest path) */
  struct INT2DMAP  pathmap;     /* path for topological distance p[i][j]:the last node just before j on the path i->j. */
  int    Npermu;                  /* Number of permutaion of equivalent atom match */
  struct PERMUTATION permuhead;   /* head of permutation */
  struct KEYVALUE HeadProperty;   /* head of KEYVALUE-linked list for molecular property */  
  unsigned char *atompair_descriptor; /* descriptor of atom pairs (malloc later)*/ 
  unsigned char *ring_descriptor;     /* descriptor of rings and blocks  (malloc later)*/ 
  char    SetStereo3D;                /* Setup Stereo Parity from 3D conformation (T or F) */ 
};



struct MATCH{
/*
   MATCH: structure of matching atoms.
    example)
     Npair = 6 
     anumA = [0,1,2,3,4,6]
     anumB = [0,1,3,5,7,8]
 */
  int num;                      /* number */
  int Npair;                    /* Number of atom pairs */
  int Npair_malloc;             /* Number of atom pairs malloced */
  /*
  unsigned char *anumA;          atom_numbers for molecule A[Npair] (malloc later) [0..Npair-1] 
  unsigned char *anumB;          atom_numbers for molecule B[Npair] (malloc later) [0..Npair-1] 
  */
  int *anumA;         /* atom_numbers for molecule A[Npair] (malloc later) [0..Npair-1] */
  int *anumB;         /* atom_numbers for molecule B[Npair] (malloc later) [0..Npair-1] */
  float  select_dis;            /* heuristic selection distance (dissimilarity) */
  float  Dneiatm;               /* Dis for diff. of neighboring atoms     (for 2D) */
  float  Dextcon;               /* Dis for diff .of extended connectivity (for 2D) */
  float  Dtopodis;              /* Dis for diff .of topological distance  (for 2D) */
  unsigned char  Natomtype[MAX_ATOMTYPE]; /* Number of atom type class */
  struct MATCH *next, *prev;    /* pointer to 'next' and 'previous' MATCH nodes */
  char   mark;                  /* mark for various purpose */
  int    Ncomponent;            /* Number of component (connected graph) */
  char   nodetype;              /* node type 'I'nitial,'N'ormal,'E'nd. (this 'E' is for use a part of malloced linked list)*/
  float  subfloat;              /* sub float number  for various purpose */
  char   ConnectGraphType;      /* Connection of MCS. 'C'onnected, 'D'isconnected, 'T'D-MCS,'t'C-MCS */
  int    maxDIFtopodis;         /* accepted maximum difference of topological distance (number of bonds in the shortest path) */
  char   atomtype_class;        /* atomtype classification */
};



/*
   
   MATCHlist : linear linked list of "struct MATCH"
 
    >> Normal Linked List << 

       ['I'] -> ['N'] -> ['N'] -> ['N'] -> ['N'] -> NULL
        head     normal   normal   normal   normal
         dummy
            
     >> Partial Linked List using 'type' variables << 
           
     ['I'] -> ['N'] -> ['N'] -> ['E'] -> ['N'] -> NULL 
      head    normal   normal   end
      dummy                     dummy
 */



struct LIBMOL{
  /* A node in library */
  int  num;                   /* number of LIBMOL (0,1,2,3, ....)    */
  char *name;                 /* name of LIBMOL (malloc later)       */
  char *class;                /* class of LIBMOL (malloc later) */
  char *molform;              /* molecular formula (malloc later)  */
  int   Nheavyatom;           /* Number of heavy atoms */
  float select_dis;           /* heuristic seclection distance (dissimilarity) */
  float score_for_sort;       /* Score for sorting */
  float alph_score;           /* Alphabet score converted from name[] (0>1>2..>A>B>C..>a>b>c) */
  char   doneMCS;             /* 1: done MCS comparison, 0: do not */
  struct MOLECULE *mol;       /* pointer for MOLECULE (only for mode 'A') (malloc later)  */ 
  int    Npair;               /* Number of matched atom pairs */
  unsigned char *anumA;       /* atom_numbers for molecule A[Npair] (malloc later) [0..Npair-1] */
  unsigned char *anumB;       /* atom_numbers for molecule B[Npair] (malloc later) [0..Npair-1] */
  int    *num_in_fileA;       /* atom_num_in_file for molecule A[Npair] (malloc later) [0..Npair-1] */
  int    *num_in_fileB;       /* atom_num_in_file for molecule B[Npair] (malloc later) [0..Npair-1] */
  struct LIBMOL *next, *prev;

  unsigned char  *oneatom_descriptor;   /* descriptor of one atom types (Number of atom type class) (malloc later) */
  unsigned char  *atompair_descriptor;  /* descriptor of atom pairs (malloc later) */ 
  unsigned char  *ring_descriptor;      /* descriptor of rings  and blocks (malloc later)*/ 
  char  *property0;                     /* property number 0 */ 
  int    num_file;                      /* file number of this molecule */ 
  long   offset_libfile;                /* file offset for library file */
  long   offset_descfile;               /* file offset for descriptor file */ 
  float  tanimoto_mcs;                  /* Tanimoto Coefficient using mcs */
  float  tanimoto_oneatom;              /* Tanimoto Coefficient using one_atom descriptor*/
  float  tanimoto_atompair;             /* Tanimoto Coefficient using atom_pair descriptor*/
  int    subnum;                        /*  number for various purposes */
  float  *mul_tanimoto_mcs;             /* only for multiple queries (malloc later) */
  float  *mul_tanimoto_atompair;        /* only for multiple queries (malloc later) */
  float  *mul_tanimoto_oneatom;         /* only for multiple queries (malloc later) */
 
  float  RMSD;                          /* RMSD for matched atoms between lib and query */
  double Raxisang[4];                   /* rot_axis[0][1][2] rot_angle[3] for superposition on query */ 
  double Tvec[4];                       /* translation vector [0][1][2]   for superposition on query */ 
  int   mark;                           /* mark for various purpose */
};



struct LINENODE{
  int  num;                         /* number of LINENODE (0,1,2,....) */
  char *line;
  struct LINENODE *next,*prev;
};


