#ifndef LAGRANGE_MIN_H
#define LAGRANGE_MIN_H
void lagrange_min(PARAMS *p, double * coords);
double K(PARAMS * p, double * coords);
double U(PARAMS * p, double * coords);
void f_pot(PARAMS * p, double * coords,int part_num, int time_slice, double * f_out);
void f_kin(PARAMS * p, double * coords,int part_num, int time_slice, double * f_out);
double f_kin_x(PARAMS * p, double * coords,int  dim, int part_num, int time_slice);
void  df(PARAMS * p, double * coords, int part_num, int time_slice, double * df_out);
int gsl_lagrange_min(PARAMS * p, double * coords);
int gsl_lagrange_min_conj(PARAMS * p, double * coords);
void fsq_min(PARAMS *p, double * coords);
void relax(PARAMS * p, double * coords, int local);
void test_df(PARAMS * p, double * coords);
double u_tot(PARAMS * p, double * coords, int time_slice);
void init_relax(PARAMS * p, double * coords, double radius);
//void init_relax(PARAMS * p, double * coords, double radius, double (* u_f)(PARAMS *, double *, int, int));
double u_tot_verbose(PARAMS * p, double * coords, int time_slice,int verbose);
double total_path_error(PARAMS *p, double * coords, int verbose);
void set_functions(PARAMS * p);
double u_tot_verbose(PARAMS * p, double * coords, int time_slice,int verbose);
int monte_relax_ends(PARAMS * p, double * coords, double (* u_f)(PARAMS *, double *, int, int), int verbose);
#endif// LAGRANGE_MIN_H
