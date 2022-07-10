/*
 
 <ioLINE.h>

*/
extern int Read_LINENODEs();
extern void Add_string_to_LINENODEs();
extern void Free_LINENODEs();
extern void Show_LINENODEs();
extern int  Number_of_LINENODEs();
extern void String_from_LINENODE_num();

extern void Add_key_value_string_to_KEYVALUEs();
extern void Show_KEYVALUEs();
extern void Free_KEYVALUEs();
extern int  Get_Value_from_KEYVALUE();
extern void Write_Property_for_One_Molecule();
extern void Write_Headers_of_Property_File();
extern int  wildcard_match();
extern int  substring_match();
extern int match_properties_and_patterns();
