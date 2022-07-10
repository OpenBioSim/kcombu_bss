/*

< energies.h >

functions to transform molecular conformation.

*/

extern int  Transform_MOLECULE_by_RotBond();
extern void Transform_MOLECULE_by_Translation_Rotation();

extern float energy_total();
extern void  show_energy_total();
extern void  cal_energy_for_transform();
extern float energy_atommatch();
extern float energy_selfclash();
extern float energy_protclash();
extern int   mark_ligand_neighbor_receptor_atoms_for_Eprotclash();
extern float energy_protatrct();
extern float energy_volmovlap();
extern float energy_tpdisrest();
extern float sphere_volume_overlap();
extern void  show_selfclash();
extern int   count_selfclash();
extern void  show_protclash();
extern float fixed_rmsd_atommatch();
extern struct MATCH  *chose_MATCH_with_min_rmsd();
extern void  cal_Force_and_Torque();
extern void  cal_dEd_rbAngle();
extern void  Random_Transform_MOLECULE_by_RotBond();
