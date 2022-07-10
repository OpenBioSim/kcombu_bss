/*

< transform.h >

*/

struct TRANSFORM_VARIABLE{
  int   Natom;            /* Number of atom */
  float dT[3];            /* translational vector around the center atom.*/
  float RotAxis[3];       /* Rotational Axis */
  float RotAngle;         /* Rotational Angle */ 
  float Force[3];         /* total force */ 
  float Torque[3];        /* Torque */ 
  float **force_atom;     /* force_atom[Natom][3] (malloc later) (force[i] = - dE_dPos[i])*/

  /*** Variables for rotational bond **/
  int   Canum;            /* atom number for Center atom */
  int   Nrbond;           /* Number of rotational bonds */
  int   Nrbond_stamp;     /* Number of stamped rotational bonds */
  char  *rbType;          /* [Nrbond](malloc later) Type of rotational bond.  'M'atched, 'F'ree */ 
  char  *rbStamp;         /* [Nrbond](malloc later) Stamped bond or not 'S':stamped.  */
  char  *rbPlane;         /* [Nrbond](malloc later) Plane   bond or not 'P': plane. rbAngle should be 0 or 180 degree.*/
  float *rbAngle;         /* Rotational angle(radian). rbAngle[Nrbond] (malloc later) */
  int   *rbOanum;         /* atom numbers for Origin atoms(0,1,..,Natom-1). rbOanum[Nrbond] (malloc later) */ 
  int   *rbZanum;         /* atom numbers for Z-axis atoms(0,1,,,,Natom-1). rbZanum[Nrbond] (malloc later) */ 
  unsigned char **rbManum;/* rbManum[Nrbond][mol->Natom-1] (malloc later). */
                          /* if rbManum[k][a]=1, atom a is rotated by the bond k. */
  float *dEd_rbAngle;     /* 1st derivative of energy by rotational Angle*/

  /**  Energies for transformations **/ 
  float Etotal;           /* total energy */
  float Etotal_init;      /* initial total energy */
  float Eatommatch;       /* Eatommatch */
  float Eselfclash;       /* Eselfclash */
  float Eprotclash;       /* Eprotclash (clash with receptor protein) */
  float Eprotatrct;       /* Eprotatrct (attraction with receptor protein) */
  float Evolmovlap;       /* Evolmovlap (non-specific atom matching )*/
  float Etpdisrest;       /* Etpdisrest (target-protein distance restraint)*/
  float rmsd_match;       /* RMSD for matched atoms */
  float tanimoto_volume;  /* Tanimoto Coefficient of Volume */ 
  int   Nselfclash;       /* Number of self-clashing atom pairs */

  int  confnumT, confnumR;  /* Conformer number for Target and Reference (only for 'rkcombu') */
  int  Ninit;               /* Initial conformer number */
  int  num_atom_match;      /* Atom Match Number */

  struct MATCH atom_match;  
};

extern int  Set_TRANSFORM_VARIABLE();
extern void Malloc_TRANSFORM_VARIABLE();
extern void Free_TRANSFORM_VARIABLE();
extern void Copy_TRANSFORM_VARIABLE();
extern void Rotate_TargetMolecule_by_PCA();
extern void Rotate_TargetMolecule_by_minimum_energy_PCA();
extern void Write_TRANSFORM_VARIABLE();
extern void Bubble_Sort_Array_of_TRANSFORM_VARIABLE();
extern void Set_Type_of_Rotatable_Bond();


