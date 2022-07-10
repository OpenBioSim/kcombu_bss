/*

 <molring.h>
 
*/

extern void Set_Ring_Structure_with_maxNatom_in_ring();
extern void Set_Property_using_RING();
extern void Show_Ring_Structure();
extern void Set_DAML_ring_ATOMs();
extern void Check_Planarity_of_RING_from_Atom_XYZ();
extern int  Delete_Disconnected_Parts_Of_Molecule();
extern void Write_DAMLatom_and_RING_in_PDB();
extern void Read_DAMLring_Molecule_in_PDB();
extern int Number_of_RING();
extern void Free_RING();
extern void Free_BLOCK();
extern void Malloc_Copy_RINGs();
extern void Set_Rotatable_Bond();
