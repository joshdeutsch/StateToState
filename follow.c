#include <math.h>
#include <float.h> 
#include <assert.h> 
#include "defs.h"
#include "params.h"
#include "follow.h"
#include "json.h"

//#define DEBUG
//#define PERIODIC_BC

typedef struct  monte_ {
   int num_its_per_slice, num_accepted;
   int ti, tf;
   int verbose;
   int it_max;
   long long count_sweeps;
   int sample_size;
   int num_in_sample;
   int num_since_last;
   double  frac_accepted;
   double u_max;
   double energy;
   double beta;
   double step_size;
   double ave_e_sweep;
   double ave_e, ave_e_old;
   void (* monte_relax_once)(PARAMS * p, double * coords, struct monte_ * monte, double (* u_f)(PARAMS *, double *, int, int), int verbose);
} MONTE;


#define LOCAL

#ifdef LOCAL
#define ERR_FUNCTION two_sided_local_error
#define TOTAL_ERR_FUNCTION total_path_error_local
#else
#define ERR_FUNCTION two_sided_path_error
#define TOTAL_ERR_FUNCTION total_path_error
#endif


2d_arr(3,5)++;
arr[3+10*5]++

#define 2d_arr(arr,i,j) arr[i+N*j]
#define unpacked(arr,dim,part_num,time_slice) (arr)[(dim) + p->num_dim*((part_num) + p->num_parts*(time_slice))]
#define unpacked_vec(arr,dim,part_num_i,part_num_j) (arr)[(dim) + p->num_dim*(p->num_dim*((part_num_i)+p->num_parts*(part_num_j)))]
#define unpacked_mat(arr,dim_comp,dim_deriv,part_num_i,part_num_j) (arr)[(dim_comp) + p->num_dim*((dim_deriv)+p->num_dim*((part_num_i)+p->num_parts*(part_num_j)))]

void f_spring(PARAMS * p, double * rs, double * f_out)
{
   int dim;
   for(dim=0;dim < p->num_dim;dim++)
   {
      //f_out[dim] =  0;
      f_out[dim] = p->f_ampl*(rs[dim]);
   }
   //f_out[0] = p->f_ampl*(rs[0]);
}

double u_power(PARAMS * p, double * rs)
{
   int dim;
   double rsq = 0.0;
   for(dim=0;dim < p->num_dim;dim++)
   {
      double dx = rs[dim];
#ifdef PERIODIC_BC
      double whole = dx*p->invL;
      double rem = whole - floor(whole);
      if (rem > 0.5) rem = rem-1.0;
      rem *= p->L;
      rsq += SQR(rem);
#else
      rsq += SQR(dx);
#endif
      
   }
   double ins = 1 + rsq;
   //double ins = rsq;
   ins = QUART(ins);
   //ins = (ins);

   double u = 1.0/ins;
   return -p->f_ampl *  u;
}

double u(PARAMS * p, double * rs)
{
    return u_power(p, rs);
    //return u_osc(p, rs);
    //return u_spring(p, rs);
}

double u_tot(PARAMS * p, double * coords, int time_slice)
{
   int part_num, other_part_num,dim;
   double tot_u = 0.0;
   double * x_tmp;
   alloca_(x_tmp,p->num_dim);
   for(part_num=0;part_num < p->num_parts;part_num++)
   {
      for(other_part_num=0;other_part_num < part_num;other_part_num++)
      {
	 for(dim=0;dim < p->num_dim;dim++)
	 {
	    x_tmp[dim] = unpacked(coords,dim,other_part_num,time_slice)-unpacked(coords,dim,part_num,time_slice);
	 }

	 tot_u += u(p,x_tmp);
      }
   }
   return tot_u;
}

double u_tot_verbose(PARAMS * p, double * coords, int time_slice,int verbose)
{
   return u_tot(p, coords, time_slice);
}

//double monte_relax_slice(PARAMS * p, double * coords, int time_slice,double radius, int * num_accepted, double beta,int num_its, double (* u_f)(PARAMS *, double *, int,int))
double monte_relax_slice(PARAMS * p, double * coords, int time_slice, MONTE * monte, double (* u_f)(PARAMS *, double *, int,int))
{
   int part_num;
   int it=0;
   double *x_tmp;
   alloca_(x_tmp,p->num_dim);
   double u_0 = u_f(p,coords,time_slice,0);
   double u;
   if (monte->verbose) printf("in monte_relax_slice:time_slice = %d, u_0 = %lf\n",time_slice, u_0);
   //monte->num_accepted = 0;
   while(it < monte->num_its_per_slice)
   {
      int dim;
      part_num = p->num_parts*ran();
      for(dim=0;dim < p->num_dim;dim++)
      {
	 x_tmp[dim] = unpacked(coords,dim,part_num,time_slice); 
	 unpacked(coords,dim,part_num,time_slice) = x_tmp[dim] + monte->step_size*(2*ran()-1);
      }
      u = u_f(p, coords, time_slice,0);
      if (ran() > exp(-monte->beta*(u-u_0)))//rejected
      {
	 for(dim=0;dim < p->num_dim;dim++)
	 {
	    unpacked(coords,dim,part_num,time_slice) = x_tmp[dim];
	 }
      } 
      else 
      {
	 u_0 = u;
	 monte->num_accepted++;
      }
      if (u < monte->u_max)
      {
	 u_0=u;
	 printf("iterations = %d, u for time_slice %d = %lf is < u_max = %lf\n",it,time_slice, u, monte->u_max);
	 break;
      } 
      monte->num_since_last++;
      assert (it < monte->it_max);
      it++;
      if (it%1000 == 0 || monte->verbose) printf("in monte_relax_slice:it = %d, u = %lf\n",it, u);
   }
   monte->energy = u_0;
   return u_0;
}

//void monte_relax(PARAMS * p, double * coords, double radius, double beta, int num_its,int ti, int tf,  double (* u_f)(PARAMS *, double *, int, int), int verbose)
void monte_relax_once_local(PARAMS * p, double * coords, MONTE * monte, double (* u_f)(PARAMS *, double *, int, int), int verbose)
{
   int time_slice;
   double ave_u = 0;
   for(time_slice = monte->ti;time_slice < monte->tf; time_slice++)
   {
     double u;
     //u = monte_relax_slice(p, coords, time_slice, radius, num_its, beta, &frac_accepted, u_f);
     u = monte_relax_slice(p, coords, time_slice,  monte,  u_f);
     ave_u += u;
     monte->count_sweeps++;
     if (verbose) printf("monte_relax: time_slice = %d, u=%lf\n",time_slice, u);
   }
   monte->ave_e_sweep = ave_u/(monte->tf-monte->ti);
}

void monte_relax_once_follow(PARAMS * p, double * coords, MONTE * monte, double (* u_f)(PARAMS *, double *, int, int), int verbose)
{
   int time_slice;
   double ave_u = 0;
   int num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
   for(time_slice=0; time_slice < num_time_slices; time_slice += p->num_follow)
   {
     double u;
     u = monte_relax_slice(p, coords, time_slice,  monte,  u_f);
     ave_u += u;
     monte->count_sweeps++;
     if (verbose) printf("monte_relax: time_slice = %d, u=%lf\n",time_slice, u);
   }
   int num_slices = 2*num_time_slices/p->num_follow -1;
   monte->ave_e_sweep = ave_u/num_slices;
}


int monte_relax_ends(PARAMS * p, double * coords, double (* u_f)(PARAMS *, double *, int, int), int verbose)
{
	MONTE monte;
	monte.step_size = 0.2;
	monte.num_its_per_slice = 1e9;
	monte.ti = 1;
	monte.tf = p->num_time_slices - 1;
	monte.beta = 100.0;
	monte.verbose = 0;
	monte.it_max = 1e7;
	monte.u_max = 1.0;

	int time_slice;
	int i;
	int num_its = 0;
	int * mark;
	int finished = 0;
	double U;
	alloca_(mark,2);
	memset_(mark,0,2);
	while (finished  == 0 && num_its < monte.it_max)
	{
		finished = 1;
		for(i=0;i < 2;i++)
		{
			time_slice = i*(p->num_time_slices-1);
			if (mark[i] != 1) 
			{
				U = monte_relax_slice(p, coords, time_slice,  &monte,  u_f);
				finished = 0;
			}
			if (U < monte.u_max) 
			{
				mark[i] = 1;
			}
		}
		num_its++;
	}
	if (num_its >= monte.it_max)
	{
		printf("couldn't relax ends after %d total iterations\n",monte.it_max);
		return 0;
	}
	return 1;
}


void init_relax(PARAMS * p, double * coords, double radius) 
{
   double beta = 100.0; 
   int num_its = 10000000;
   MONTE monte;
   monte.step_size = radius;
   monte.num_its_per_slice = 1e9;
   monte.ti = 1;
   monte.tf = p->num_time_slices - 1;
   monte.beta = beta;
   monte.verbose = 0;
   monte.it_max = 1e5;
   monte.u_max = 1.0;
   monte.monte_relax_once = monte_relax_once_local;

   monte.monte_relax_once(p, coords, &monte, u_tot_verbose, 1);
}

double local_error(PARAMS *p, double * coords, int t)
{
   double dt = 1.0/p->num_time_slices;
   double dtsq = SQR(dt);
   int dim,part_num;
   double err=0.0, diff=0.0;
   double * f_out;
   alloca_(f_out,p->num_dim);

   for(part_num=0;part_num < p->num_parts;part_num++)
   {
      f_pot(p, coords,part_num, t, f_out);
      for(dim=0; dim < p->num_dim;dim++)
      {
            diff = unpacked(coords,dim,part_num,t+1) - 2*unpacked(coords,dim,part_num,t)
		                                   +unpacked(coords,dim,part_num,t-1) - dtsq*f_out[dim]; 
            err += SQR(diff);
      }
   }
   return err;
}


double total_path_error_local(PARAMS *p, double * coords, int verbose)
{
   int time_slice,t,dim;
   double error = 0;
   double err_max = 0;

   int indx=-1;
   for(time_slice=1; time_slice < p->num_time_slices-1; time_slice++)
   {
      double err =  local_error(p, coords, time_slice);
      if (err > err_max)
      {
	 err_max = err;
	 indx = time_slice;
      } 

      error += err;
   }

   if (verbose)
      printf("In total_path_error_local: err_max = %lf, indx = %d\n",err_max, indx);
   return error;
}
void adjust_monte_params(PARAMS * p, MONTE * monte, double * coords, int verbose)
{
   monte->ave_e += monte->ave_e_sweep;
   if (monte->num_in_sample >= monte->sample_size)
   {
      monte->frac_accepted = monte->num_accepted*1.0/monte->num_since_last;
      if (monte->frac_accepted < 0.3)
      {
	 monte->step_size *= 0.9;
	 if (verbose) 
	 {
	    printf("monte->frac_accepted = %lf < 0.3 step_size = %lf, total_path_error_local = %lf\n",monte->frac_accepted,monte->step_size
             ,total_path_error_local(p, coords, 1));
	 }
      }
      else if (monte->frac_accepted > 0.7)
      {
	 monte->step_size *= 1.1;
	 if (verbose) printf("monte->frac_accepted = %lf > 0.7 step_size = %lf\n",monte->frac_accepted,monte->step_size);
      }
      monte->ave_e /= monte->sample_size;

      if (verbose) printf("monte->beta = %lf monte->ave_e = %lf\n", monte->beta, monte->ave_e);
      monte->ave_e_old = monte->ave_e;
      monte->ave_e_old = monte->ave_e;
      monte->ave_e = 0;
      monte->num_in_sample = 0;
      monte->num_since_last = 0;
      monte->num_accepted = 0;
   }
   monte->num_in_sample++;
}


void f_soft(PARAMS * p, double * rs, double * f_out)
{
	/*
	 * U = -f_ampl (1-r^2)^4 r<1, 0 otherwise
	 * f = f_ampl 8(1-r^2)^3 r r<1, 0 otherwise, refers to force on other particle
	 */
	int dim;
	double rsq = 0.0;
	for(dim=0;dim < p->num_dim;dim++)
	{
#ifdef PERIODIC_BC
		double dx = rs[dim];
		double whole = dx*p->invL;
		double rem = whole - floor(whole);
		if (rem > 0.5) rem = rem-1.0;
		f_out[dim] = rem*p->L;
#else
		f_out[dim] = rs[dim];
#endif

		rsq += SQR(f_out[dim]);
	}
	if (rsq >= 1.0)
	{
		memset_(f_out,0,p->num_dim);
		return;
	}
	double ins = 1 - rsq;
	ins *= SQR(ins);
	for(dim=0;dim < p->num_dim;dim++)
	{
		f_out[dim] *= p->f_ampl * ins;
		//f_out[dim] *= p->f_ampl * 8.0/ins;
	}
}

void f_power(PARAMS * p, double * rs, double * f_out)
{
   int dim;
   double rsq = 0.0;
   for(dim=0;dim < p->num_dim;dim++)
   {
#ifdef PERIODIC_BC
      double dx = rs[dim];
      double whole = dx*p->invL;
      double rem = whole - floor(whole);
      if (rem > 0.5) rem = rem-1.0;
      f_out[dim] = rem*p->L;
#else
      f_out[dim] = rs[dim];
#endif

      rsq += SQR(f_out[dim]);
   }
   //double ins = 1 + rsq;
   double ins = rsq;
   //ins *= QUART(ins);
   ins *= (ins);
  // ins *= sqrt(ins);
   if (ins < 1e-16)
   {
	   printf("in f_power rsq = %lf, dim=%d,rs = [%lf,%lf]\n",rsq, dim, rs[0],rs[1]);
	   ins = 1e-16;
   }
   double inv8ins = 8.0/ins;
   for(dim=0;dim < p->num_dim;dim++)
   {
      f_out[dim] *= p->f_ampl * inv8ins;
      //f_out[dim] *= p->f_ampl * 8.0/ins;
   }
}

void f(PARAMS * p, double * rs, double * f_out)
{
   //return f_soft(p,rs,f_out);
   return f_power(p,rs,f_out);
   //return f_osc(p,rs,f_out);
   //return f_spring(p,rs,f_out);
}

void f_ext(PARAMS * p, double * coords,int part_num, int time_slice, double * f_out)
{
   int dim;
   double * x_tmp;
   alloca_(x_tmp,p->num_dim);
   for(dim=0;dim < p->num_dim;dim++)
   {
      x_tmp[dim] = unpacked(coords,dim,part_num,time_slice);
   }
   f_spring(p,x_tmp,f_out);
   //f_power(p,x_tmp,f_out);
}


void f_pot(PARAMS * p, double * coords,int part_num, int time_slice, double * f_out)
{
   int dim;
   double * f_tmp;
   double * x_tmp;
   double * x_self;
   alloca_(f_tmp,p->num_dim);
   alloca_(x_self,p->num_dim);
   alloca_(x_tmp,p->num_dim);
   memset_(f_out,0,p->num_dim);
   memset_(f_tmp,0,p->num_dim);

   for(dim=0;dim < p->num_dim;dim++)
   {
      x_self[dim] = unpacked(coords,dim,part_num,time_slice);
   }

   int other_part_num;
#if 1
   for (other_part_num =0; other_part_num < p->num_parts;other_part_num++)
   {
      if (other_part_num == part_num) continue;
      for(dim=0;dim < p->num_dim;dim++)
      {
         x_tmp[dim] = unpacked(coords,dim,other_part_num,time_slice)-x_self[dim];//unpacked(coords,dim,part_num,time_slice);
      }
#ifdef DEBUG
      double diff = fabs(x_tmp[0]) + fabs(x_tmp[1]);
      if (diff < 1e-10) printf("other = %d, part_num = %d, [x_tmp[0]=%lf, x_tmp[1]=%lf]\n",other_part_num,part_num,x_tmp[0],x_tmp[1]);
      assert(diff > 1e-10);
#endif
      f(p,x_tmp,f_tmp);
      for(dim=0;dim < p->num_dim;dim++)
      {
         f_out[dim] += f_tmp[dim];
      }
   }
#else
   f_ext(p, coords, part_num, time_slice, f_tmp);
   for(dim=0;dim < p->num_dim;dim++)
   {
      f_out[dim] += f_tmp[dim];
   }
#endif
}

double partial_path_error(PARAMS *p, double * coords, int time_slice, int verbose)
{
   double dt = 1.0/p->num_time_slices;
   double dtsq = SQR(dt);
   int t,dim,part_num;
   double err=0.0;
   double * f_out;
   double * coords_tmp = p->tmp_coords;
   alloca_(f_out,p->num_dim);
   /*
   0 is fixed        9 is fixed
     1 var
       2 slaved
         3 slaved
           4,5 var
               6 slaved
                 7 slaved
   | |     | |     | | 
   0 1 2 3 4 5 6 7 8 9
     ->2     ->6
       ->3     ->7
         ->4     ->8
             
   */
   if (time_slice%p->num_follow == 0)
   {
      int shift = p->num_dim*p->num_parts*time_slice;
      memcpy_(coords_tmp+shift,coords+shift,2*p->num_dim*p->num_parts);
      for(t = time_slice; t < time_slice+p->num_follow; t++)
      {
	 for(part_num=0; part_num < p->num_parts; part_num++)
	 {
	    f_pot(p, coords_tmp,part_num, t+1, f_out);
	    for(dim=0;dim < p->num_dim;dim++)
	    {
	       unpacked(coords_tmp,dim,part_num,t+2) = 2*unpacked(coords_tmp,dim,part_num,t+1)
		                                         -unpacked(coords_tmp,dim,part_num,t) + dtsq*f_out[dim]; 
	    }
	 }
      }

      int t_init = time_slice+ p->num_follow;
      int t_final = time_slice+ p->num_follow+2;
      if (t_final == p->num_time_slices) t_init++;
      for(t = t_init ; t < t_final; t++)
      {
	 for(part_num=0; part_num < p->num_parts; part_num++)
	 {
            //f_pot(p, coords,part_num, t-1, f_out);
	    for(dim=0;dim < p->num_dim;dim++)
	    {
	       double diff = unpacked(coords,dim,part_num,t) - unpacked(coords_tmp,dim,part_num,t); 
	       err += SQR(diff);
	       if (verbose && fabs(diff) > 0) printf("in partial_path_error: t=%d diff=%lf\n",t,diff);
	    }
	 }
      }
      t_init = time_slice+ p->num_follow;
      t_init = time_slice+ p->num_follow;
      if (0 && t_final < p->num_time_slices)
      {
	 for(part_num=0; part_num < p->num_parts; part_num++)
	 {
	    for(dim=0;dim < p->num_dim;dim++)
	    {
	       double diff = unpacked(coords,dim,part_num,t_init+1) - unpacked(coords_tmp,dim,part_num,t_init+1)
	                     -(unpacked(coords,dim,part_num,t_init) - unpacked(coords_tmp,dim,part_num,t_init)); 
	       //err += SQR(diff*(p->num_time_slices));
	       //err += SQR(diff*(10));
	       err += SQR(diff*(2));//worked well for rep.5part with (r^2)^4 potential
	       if (verbose && fabs(diff) > 0) printf("in partial_path_error: t=%d diff=%lf\n",t,diff);
	    }
	 }
      }

   }
   else
   {
      printf("timeslice value %d seems to be off\n",time_slice);
   }
   return err;
}

double two_sided_path_error(PARAMS *p, double * coords, int time_slice,int verbose)
{
    double error_m=0;
    if (time_slice - p->num_follow >= 0)
        error_m = partial_path_error(p, coords, time_slice-p->num_follow, verbose);
    double error_p = partial_path_error(p, coords, time_slice, verbose);
    return error_m+error_p;
}

double two_sided_local_error(PARAMS *p, double * coords, int t, int verbose)
{
   if (t == 1)
      return  local_error(p,coords,t) + local_error(p,coords,t+1);
   else if (t == p->num_time_slices -2)
      return  local_error(p,coords,t) + local_error(p,coords,t-1);
   else
      return  local_error(p,coords,t) + local_error(p,coords,t-1) + local_error(p,coords,t+1);
}


double total_path_error(PARAMS *p, double * coords, int verbose)
{
   int time_slice,t,dim;
   double dt = 1.0/p->num_time_slices;
   double dtsq = SQR(dt);
   double error=0.0;
   int num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
   /*
   0 is fixed        9 is fixed
   | |     | |     | | 
   0 1 2 3 4 5 6 7 8 9
     ->2     ->6
       ->3     ->7
         ->4     ->8
   */
   for(time_slice=0; time_slice < num_time_slices; time_slice += p->num_follow)
   {
      double err = partial_path_error(p, coords, time_slice, verbose);
      if (verbose) printf("in total_path_error time_slice=%d, err=%lf\n",time_slice,err);
      error +=  err;
   }
   if (verbose)
   {
       printf("total_path_error = %lf\n",error);
   }
   
   return error;
}


/* 
 * dx/dt = - dV/dx
 */

void relax_once_local(PARAMS * p, double * coords, int verbose)
{
   int time_slice,ts0,i,dim,part_num;
   double * error;
   double dx = p->dx;
   int num_time_slices = p->num_time_slices;
   newarr_(error,2*p->num_dim*p->num_parts*(num_time_slices-1));

   for(i=1; i < num_time_slices-1; i++)
   {
      time_slice = i;
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
         {
            double tmp = unpacked(coords,dim,part_num,time_slice);
            unpacked(coords,dim,part_num,time_slice) = tmp + dx;
            double error_pos =  p->err_function(p, coords, time_slice,0);
            unpacked(coords,dim,part_num,time_slice) = tmp - dx;
            double error_neg =  p->err_function(p, coords, time_slice,0);
            unpacked(error,dim,part_num,i-1) =  0.5*(error_pos - error_neg);
            unpacked(coords,dim,part_num,time_slice) = tmp;
         }
      }
   }
   double norm = p->eps/dx;
   for(i=1; i < num_time_slices-1; i++)
   {
      time_slice = i;
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
         {
            unpacked(coords,dim,part_num,time_slice) -= norm * unpacked(error,dim,part_num,i-1);
            //unpacked(coords,dim,part_num,time_slice) += 0.001*(2*ran()-1);
         }
      }
   }
   freearr_(error);
}

void reverse_coords(PARAMS * p, double * coords, int verbose)
{
   int part_num, dim, t;

   for(t=0;t < p->num_time_slices/2;t++)
   {
      int t_rev = p->num_time_slices-1-t;
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
	 {
	    double a = unpacked(coords,dim,part_num,t);
	    double b = unpacked(coords,dim,part_num,t_rev);
	    unpacked(coords,dim,part_num,t_rev) = a;
	    unpacked(coords,dim,part_num,t) = b;
	 }
      }
   }
}

void relax_once_follow(PARAMS * p, double * coords, int verbose)
{
   int time_slice,ts0,i,dim,part_num;
   double * error;
   double dx = p->dx;
   int num_follow = p->num_time_slices/p->num_follow;
   int num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
   int num = 2*(num_time_slices/p->num_follow);
   newarr_(error,2*p->num_dim*p->num_parts*(num-1));

   for(i=1; i < num; i++)
   {
      int j = i/2;
      ts0 = j*p->num_follow;
      time_slice = ts0 + i%2;
      double error0 =  p->err_function(p, coords, ts0, 0);
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
         {
            double tmp = unpacked(coords,dim,part_num,time_slice);

            unpacked(coords,dim,part_num,time_slice) = tmp + dx;
            double error_pos =  p->err_function(p, coords, ts0, 0);
            unpacked(coords,dim,part_num,time_slice) = tmp - dx;
            double error_neg =  p->err_function(p, coords, ts0, 0);
#ifdef FIRST_ORDER_DERIV
            unpacked(error,dim,part_num,i-1) =  0.5*(error_pos - error_neg);
#else

            unpacked(coords,dim,part_num,time_slice) = tmp + 2*dx;
            double error_pos2 =  p->err_function(p, coords, ts0, 0);
            unpacked(coords,dim,part_num,time_slice) = tmp - 2*dx;
            double error_neg2 =  p->err_function(p, coords, ts0, 0);

            unpacked(error,dim,part_num,i-1) = (-error_pos2+error_neg2 + 8*(error_pos - error_neg))/12.0;
#endif

            unpacked(coords,dim,part_num,time_slice) = tmp;
         }
      }
   }
   double norm = p->eps/dx;
   for(i=1; i < num; i++)
   {
      int j = i/2;
      time_slice = j*p->num_follow + i%2;
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
         {
            unpacked(coords,dim,part_num,time_slice) -= norm * unpacked(error,dim,part_num,i-1);
            //unpacked(coords,dim,part_num,time_slice) += 0.001*(2*ran()-1);
         }
      }
   }
   freearr_(error);
}

void set_functions(PARAMS * p)
{

	if (p->local)
	{
		p->err_function  = two_sided_local_error;
		p->total_err_function =  total_path_error_local;
		p->relax_once = relax_once_local;
	}
	else
	{
		p->err_function  = two_sided_path_error;
		p->total_err_function =  total_path_error;
		p->relax_once = relax_once_follow;

	}

}

void relax(PARAMS * p, double * coords, int local)
{
    int i;
    int num_time_slices = (p->num_time_slices/p->num_follow)*p->num_follow;
    double dt = 1.0/num_time_slices;
    double err;
    double err_old = 1e0;
    double * coords_tmp;
    memcpyarr_(coords_tmp,coords);
    MONTE * monte;
    newarr_(monte,1);

    set_functions(p);
    if (p->monte_its > 0)
    {
       double radius = 0.1;
       printf("entering monte relax ...\n");
       monte->step_size = radius;
       monte->num_its_per_slice = 1;
       monte->ti = 1;
       monte->beta = p->beta_init;
       monte->tf = p->num_time_slices - 1;
       monte->it_max = p->monte_its;
       monte->verbose = 0;
       monte->u_max = -1e9;
       monte->sample_size = 100;
       if (local)
       {
	       monte->monte_relax_once = monte_relax_once_local;
       }
       else
       {
	       monte->monte_relax_once = monte_relax_once_follow;
       }
       for(i=0; i < monte->it_max;i++)
       {
	  //monte_relax(p, coords, radius, beta, 2, 1, p->num_time_slices-1, p->err_function, 0);
          monte->beta = p->beta_init*pow((i+10000)/10000.0,1.0);
	  //printf("monte->beta = %lf\n",monte->beta);
	  monte->monte_relax_once(p, coords, monte, p->err_function,0);
	  adjust_monte_params(p,monte, coords, 1);
       }
    }
    
    {
       printf("entering relax ...\n");
       for(i=0; i < p->relax_its;i++)
       {
	  //printf("in relax: i=%d\n",i);
	  memcpy_(coords_tmp,coords,p->num_dim*p->num_parts*p->num_time_slices);
	  p->relax_once(p, coords,1);
	  int verbose;
	  verbose = (i%100000) == 0;
	  err = p->total_err_function(p, coords,verbose);
	  if (0)
	  {
	     if (err/err_old > 1.2)
	     {
		memcpy_(coords,coords_tmp,p->num_dim*p->num_parts*p->num_time_slices);
		p->eps *= 0.9;
		printf ("err = %lf > err_old = %lf :changing p->eps to %.*e\n",err, err_old,12,p->eps);
	     }
	     else 
	     {
		p->eps *= 1.01;
	     }
	  }
	  err_old = err;

	  if (i%10000 == 0) 
	  {
	    printf("err = %.*e, p->eps = %.*e\n",10, err, 10, p->eps);
	  }
	  if (err < p->error)
	  {
	     printf("err = %lf p->error = %lf,  breaking\n",err,p->error);
	     break;
	  }
	  if (p->eps < 1e-8)
	  {
	     printf("p->eps = %lf breaking\n",p->eps);
	     break;
	  }
	  //reverse_coords(p, coords, 0);

       }
    }
    printf("err = %lf\n",err);
    if (!local) memcpy_(coords,p->tmp_coords,size_(coords)-1);
    printf("coords[%d] = %lf\n",p->num_time_slices-1,coords[p->num_time_slices-1]);
}

void expand_vars1(PARAMS * p, double * min_vars, double * coords)
{

   int i,time_slice, part_num, dim;
   int num_active_slices = p->num_time_slices/p->num_follow;

   for(time_slice=0,i=0;time_slice < p->num_time_slices; time_slice += p->num_follow, i++)
   {
      for(part_num=0; part_num < p->num_parts;part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
         {
            int indx = dim + p->num_dim*(part_num+p->num_parts*i);
            unpacked(coords,dim,part_num,time_slice) = min_vars[indx];
         }
      }
   }
}

