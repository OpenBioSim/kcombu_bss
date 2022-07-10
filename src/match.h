/*
 
 <match.h>

*/


/** FUNCTIONS (GLOBAL) **/
extern int  Match_Two_Molecules();
extern int  Make_All_Atom_Pair_MATCH();
extern int  Make_All_Atom_Pair_MATCH_for_Substructure();
extern float cal_select_dis_for_MATCH();
extern void cal_select_dis_for_MATCHlist();
extern void Malloc_MATCH();
extern void Free_MATCH();
extern void Malloc_MATCHlist();
extern void Add_MATCH_to_MATCHlist();
extern int Add_MATCH_to_MATCHlist_Keep_Only_Best();
extern void Free_MATCHlist();
extern void Copy_MATCH();
extern void Copy_MATCHlist();
extern int Length_of_MATCHlist();
extern int Nmark_in_MATCHlist();
extern float Tanimoto_Coefficient();
extern int Smaller_Int();
extern int Larger_Int();
extern void Keep_Only_Nonredundant_MATCH();
extern void add_hydrogen_pairs_to_MATCH();
extern float Compare_Electric_Charges();
extern void Sort_MATCHlist_by_Npair_and_select_dis();
