/*

 <molprop.h>
 
*/

extern char MoleculeFileType();
extern int  Read_MOLECULE();
extern int  Set_MOLECULE();
extern void Initialize_MOLECULE();
extern void Initialize_MOLECULE_string();
extern void Malloc_MOLECULE();
extern void Free_MOLECULE();
extern void Copy_MOLECULE();
extern void Copy_MOLECULE_XYZ();
extern void Copy_ATOM();
extern void Copy_BOND();
extern void Set_Center_of_Gravity();
extern void Set_one_char_ele();
extern void Set_Nheavyatom_Nheavybond();
extern void Set_radius();
extern void Set_Nneighbor_CHNOSP();
extern void Set_Nneighbor_atomtype();
extern void Cal_Distance_Map();
extern void Cal_Distance_Map_for_Heavy_atoms();
extern void Set_Molform();
extern void Set_atomtype_of_Atoms();
extern void Count_atomtype_in_Natomtype();
extern float Max_Tanimoto_from_Nelement();
extern float Max_Tanimoto_from_Natomtype();
extern float Tanimoto_bwn_oneatom_descriptors();
extern int  Substructure_condition_oneatom_descriptors();
extern int  Isomorphic_condition_oneatom_descriptors();
extern int  Set_Extended_Connectivity();
extern void Set_DAMRL_from_2D_structure_with_hydrogen();
extern void Set_Ring_Structure();
extern int  Delete_Disconnected_Parts_Of_Molecule();
extern void Show_conmap();
extern void Set_max_atomtype_by_atomtype_class();
extern int Number_of_atomtype();
extern float Score_bwn_atomtypes();
extern char *atomtype_string_from_number();
extern char molecular_file_type_from_tail_of_file();
extern void Get_Part_Of_Line();
extern void Split_to_Words();
extern void Split_to_Words_With_Double_Splsym();
extern int First_Appear_Pos_in_String();
extern void Remove_Symbol_from_String();
extern void Remove_HeadTail_Space_from_String();
extern void Change_Space_To_Underbar();
extern void Change_String_ToUpper();
extern int Match_Tail_String();
extern int Number_of_Element();
extern int Check_Element_String_Or_Not();
extern void Write_atomtype_in_rasmol_script();
extern char* string_getenv();
