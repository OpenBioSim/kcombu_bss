/*
 
 <ioPDB.h> 

*/ 

extern int  Translate_Lines_into_PDB_Molecule();
extern void Write_PDB_Molecule();
extern void Guess_Bonds_from_Atom_XYZ();
extern void Guess_Orbital_from_Atom_XYZ();
extern void Set_Aromatic_Ring_by_3D_Planarity();
extern void Set_DAMRL_from_3D_structure();
extern void Set_DAMRL_from_3D_RING();
extern void Set_Resi_Chain_Rnum_to_Molecule();
extern void Set_Unique_atomname_to_Molecule();
