/*

<molpermu.h>
 
*/

/** FUNCTIONS (GLOBAL) **/
extern int  Make_All_Permutation_of_Equivalent_Heavyatoms();
extern int  Make_EquivMap_for_Heavyatoms();
extern int  Number_of_EquivMap_for_Heavyatoms();
extern void show_atom_permutation();
extern void show_all_the_atom_permutations();
extern int  Number_of_PERMUTATION();
extern int Make_Permutated_MATCHlist_from_Top_MATCH();
extern int Make_Permutated_MATCHlist_from_Top_MATCH_by_only_molA_permutation();
extern int Make_NonRedundant_MATCHlist_by_Permutation();
extern int Check_Equivalence_bwn_Two_MATCHes_by_Permutation();
extern void Free_PERMUTATION();
extern void Rerank_MATCHlist_by_transformed_RMSD();
