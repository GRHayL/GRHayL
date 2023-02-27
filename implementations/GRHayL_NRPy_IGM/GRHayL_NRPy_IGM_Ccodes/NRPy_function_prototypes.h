void MoL_malloc_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_free_memory_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_malloc_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_free_memory_non_y_n_gfs(const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_step_forward_in_time(griddata_struct *restrict griddata, const REAL dt);
void bcstruct_set_up(const paramstruct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct);
void apply_bcs_outerradiation_and_inner(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                        const REAL custom_wavespeed[NUM_EVOL_GFS],
                                        const REAL custom_f_infinity[NUM_EVOL_GFS],
                                        REAL *restrict gfs, REAL *restrict rhs_gfs);
void apply_bcs_inner_only(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs);
void apply_bcs_outerextrap_and_inner(const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs);
REAL find_timestep(const paramstruct *restrict params, REAL *restrict xx[3], const REAL CFL_FACTOR);
void xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]);
void set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],paramstruct *restrict params, REAL *restrict xx[3]);
void Cart_to_xx_and_nearest_i0i1i2(const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]);
void Cart_to_xx_and_nearest_i0i1i2_global_grid_center(const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]);
void prims_to_cons(const paramstruct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, REAL *restrict evol_gfs);
void Balsara1(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Balsara2(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Balsara3(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Balsara4(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Balsara5(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Cylindrical_Explosion(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Magnetic_Rotor(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void Loop_Advection(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict evol_gfs, REAL *restrict auxevol_gfs);
void initial_data(const char *initial_data_option, const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict evol_gfs,REAL *restrict auxevol_gfs);
void calculate_metric_gfs(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict auxevol_gfs);
void interpolate_metric_gfs_to_cell_faces(const paramstruct *params,REAL *auxevol_gfs,const int flux_dirn);
void compute_B_and_Bstagger_from_A(const paramstruct *restrict params, const REAL *evol_gfs, REAL *auxevol_gfs);
int main(int argc, const char *argv[]);
void set_Cparameters_to_default(paramstruct *restrict params);
