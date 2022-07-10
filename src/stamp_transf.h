/*

< stamp_transf.h >

functions to transform molecular conformation.

*/

extern int  Stamp_RotAngle_for_Reference_Dihedral_Angle();
extern int  Stamp_sp3_Chirality();
extern int  Check_3D_sp3_Chirality();
extern int  Modify_MATCH_For_3D_sp3_Chirality_with_Permutation();
extern int  Stamp_Five_or_Six_Atom_Ring_Geometry();
extern int  Check_Five_or_Six_Atom_Ring_Geometry();
extern void Superimpose_Just_Two_Gcenters();
extern int  Random_Change_sp3_Chirality();
extern int  Random_Change_Flip_Of_Fragments();
extern int  Change_Flip_Of_Fragment();
extern int  Check_Bond_Angle_Difference();
