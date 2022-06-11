#ifndef PARAMS_H
#define PARAMS_H
typedef struct params_{
   int num_its, num_time_slices,num_parts,num_dim, num_follow, num_runs,num_runs_in,monte_its,relax_its,local;
   double f_ampl,eps,L,invL,error,dx,dtsq,beta_init,step_size,tol,f_err;
   double min_out;
   unsigned int seed, seed2;
   double * tmp_coords;
   double (* err_function)(struct params_ *p, double * coords, int t, int verbose);
   double (* total_err_function)(struct params_ *p, double * coords, int verbose);
   void (* relax_once)(struct params_ *p, double * coords, int verbose);
} PARAMS;
#endif// PARAMS_H
