/*

 <ioSMILES.h>

*/
/** FUNCTIONS (GLOBAL) **/
extern int  Translate_Lines_into_SMILES_Molecule();
extern int  make_LUNITs_from_Molecule_by_DFS();
extern void make_string_from_LUNITs();
extern int  length_string_from_LUNITs();
extern void Write_SMILES_Molecule();
extern void randomize_atom_order_of_MOLECULE();
extern void Set_Unique_Extended_Connectivity();
extern void Set_Unique_Order_of_Atoms();
extern int  Number_of_LUNITs();
extern void Free_LUNITs();
extern void Set_Connected_Structure();
extern void Keep_Only_Largest_Connected_Structure();
extern void Set_Aromatic_Bond_Atoms();


/** STRUCTURE for Linear Unit **/
struct LUNIT{
  int  num;
  char type;       /* 'A'tom, 'N'umber, 'B'ond,'P'arenthesis, '.':period */
  char str[8];
  int  atom_num;
  struct LUNIT *next,*prev;
};

