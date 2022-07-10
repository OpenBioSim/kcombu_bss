/*
 <qRMS.h>

*/


struct POS_ARRAY{
  int N;
  int Nmalloc;
  double **pos; /* [0..N-1][3] (malloc later) */
};



extern float Calculate_CRMS_MATCH_Quaternion();
extern float Calculate_CRMS_bwn_POS_ARRAYs();
extern float CRMS_MATCHed_Molecules_without_Superimpose();
extern float CRMS_Same_Molecules_without_Superimpose();
extern float Rotate_Molecule();
extern float RMSD_bwn_Rotated_and_Init();
extern void calculate_gcenter();
extern void Malloc_POS_ARRAY();
extern void Free_POS_ARRAY();
extern int  make_rot_axis_angle_from_rot_matrix();
extern void make_rot_matrix_from_rot_axis_angle();
extern void make_tvec_from_rmat_gorig_gtra();
