/*
 
 <ioList.h>

*/


/*** FUNCTIONS (GLOBAL) ***/
extern int Make_LibFileList_from_Colon_Separated_String();
extern int Get_List_of_LibFiles_Under_Dir();
extern int Get_List_of_LIBMOLs_from_one_library_file();
extern int Read_One_Molecule_from_one_Library_file();
extern int Read_One_Entry_as_LINENODEs();
extern int Read_List_of_LIBMOLs();
extern int Get_List_of_LIBMOLs_Under_Dir();
extern void Copy_MATCH_Info_to_LIBMOL();
extern void Write_Search_Result();
extern int  Write_Similar_Molecules_To_One_File();
extern int  Write_Similar_Molecules_Under_Dir();
extern void Write_Clustering_Res_List_Of_NODEs();
extern void Cal_Nsta_Nend_from_Bunshi_Bunbo_for_all_vs_all();
extern void Write_Similarities_in_List_format_Header();
extern void Write_Two_Groups_Similarities_in_List_format_Header();
extern void Write_Similarities_in_List_format_Tail();
extern int  Number_of_LIBMOL();
extern struct LIBMOL *Get_LIBMOL_by_number();
extern struct LIBMOL *Add_LIBMOL_to_LIBMOLlist();
extern void Increment_Nmatch_for_Matched_Atoms();
extern void Write_Nmatch_rasmol_Script();
extern void insert_LIBMOL_into_score_decreasing_proper_position();
extern void delete_the_first_LIBMOL_node();
extern void reverse_LIBMOL_list();
extern void Show_LIST_of_LIBMOL();
