/*

 <AtmPairDesc.h>
 
*/

extern void Set_AtomPair_Descriptor();
extern void Set_Ring_Descriptor();
extern void Write_Headers_of_AtomPair_Descriptor();
extern int  Read_Headers_of_AtomPair_Descriptor();
extern void Write_Footers_of_AtomPair_Descriptor();
extern void Write_AtomPair_Descriptor();
extern void write_atompair_descriptor_for_one_molecule();
extern float similarity_bwn_atompair_descriptors();
extern float similarity_bwn_atompair_descriptors_oneatom_normalized();
extern int substructure_condition_atompair_descriptors();
extern int isomorphic_condition_atompair_descriptors();
extern int Read_AtomPair_Descriptor_into_LIBMOL();
extern struct LIBMOL* Read_AtomPair_Descriptor_into_LIBMOL_One_by_One();
