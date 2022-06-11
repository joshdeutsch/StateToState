#ifndef LAGRANGE_MIN_H
#define LAGRANGE_MIN_H
void lagrange_min(PARAMS *p, double * coords);
double K(PARAMS * p, double * coords);
double U(PARAMS * p, double * coords);
void f_pot(PARAMS * p, double * coords,int part_num, int time_slice, double * f_out);
void f_kin(PARAMS * p, double * coords,int part_num, int time_slice, double * f_out);
double f_kin_x(PARAMS * p, double * coords,int  dim, int part_num, int time_slice);
int gsl_lagrange_min(PARAMS * p, double * coords);
int gsl_lagrange_min_conj(PARAMS * p, double * coords);
#endif// LAGRANGE_MIN_H
