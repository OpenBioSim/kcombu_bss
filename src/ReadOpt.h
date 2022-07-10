/*

 <ReadOpt.h>

*/


struct OPTION_VALUE{
 char *opt;   /* key string   of option */
 char *val;   /* value string of option */
 int  Lopt;   /* Length of opt[] string */
 int  num;
 char accept; /* accept = 1, this is already translated */
 struct OPTION_VALUE *next,*prev;
};


/*** FUNCTIONS (GLOBAL) ***/
extern void Show_Options();
extern void Read_Options_From_Arguments();
extern void Read_Options_From_File();
extern void Make_COMMAND_String_From_OptVal_List();
