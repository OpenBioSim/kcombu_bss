/*
 < outPS.h>

*/

extern void Write_Molecule_in_PostScript();
extern void Write_MATCH_in_PostScript();

struct PS_PARAMETER{
  char text;  /* 'M'atched_atom_number, 'T':atomtype, 'N':atom number, 'A':atom name,'E':element */
};
