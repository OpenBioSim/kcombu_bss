/*
 <gen2D.h>

*/
extern int  Generate_2D_Coordinates();
extern int  set_hydrogen_XY();
extern void Write_SDF_UNIT();
extern void Free_UNITs();
extern void Rescaling_Corrdinates_of_Molecule();


/** STRUCTURE for UNIT OF ATOMS FOR 2D DEPICTION **/
struct UNIT{
  int    num;            /* Number of rings (0,1,2,... ) */
  int    Natom;          /* Number of member ATOMS */
  int    *num_atoms;     /* [Natom] :array of atom numbers */
  int    *Nneighbor;     /* [Natom] :number of neighbors within the unit */
  int    **neighbors;    /* [Natom][Natom] : neighbor atom index number */
  float  **pos;          /* [Natom][2] : XY-coordinates of points */
  char  type;            /* 'R'ing, 'C'hain */
  char  set_localXY;     /* set local XY coordinate. 'T' or 'F' */
  char  set_globalXY;    /* set global XY coordinate. 'T' or 'F' */
  struct UNIT *next,*prev;
};

