
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include "defs.h"
#include "params.h"
#include "follow.h"
#include "lagrange_min.h"

typedef struct {
   PARAMS * p;
   double * coords;
} GSL_PARAMS;


#define unpacked(arr,dim,part_num,time_slice) arr[(dim) + p->num_dim*((part_num) + p->num_parts*(time_slice))]

void expand_vars(const gsl_vector *v, void *void_params)
{
	GSL_PARAMS * gsl_params = (GSL_PARAMS *) void_params;
	double * coords = gsl_params->coords;
	PARAMS * p = gsl_params->p;

	if (p->local)
	{
		int t;
		for(t=1;t < p->num_time_slices-1;t++)
		{
			int part_num, dim;
			for(part_num=0; part_num < p->num_parts; part_num++)
			{
				for(dim=0;dim < p->num_dim;dim++)
				{
					int indx = dim + p->num_dim*(part_num+p->num_parts*(t-1));
					unpacked(coords,dim,part_num,t) = gsl_vector_get(v, indx);
				}
			}
		}

	}
	else
	{
		int num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
		int num = 2*(num_time_slices/p->num_follow);
		int i;
		for(i=1; i < num; i++)
		{
			int j = i/2;
			int time_slice, part_num, dim;
			int ts0 = j*p->num_follow;
			time_slice = ts0 + i%2;
			for(part_num=0; part_num < p->num_parts; part_num++)
			{
				for(dim=0;dim < p->num_dim;dim++)
				{
					int indx = dim + p->num_dim*(part_num+p->num_parts*(i-1));
					unpacked(coords,dim,part_num,time_slice) = gsl_vector_get(v, indx);
				}
			}

		}
	}
}

double my_f (const gsl_vector *v, void *void_params)
{
   GSL_PARAMS * gsl_params = (GSL_PARAMS *) void_params;
   double * coords = gsl_params->coords;
   PARAMS * p = gsl_params->p;
   int time_index, part_num,dim;
   int indx=0;
   expand_vars(v, void_params);
   int verbose=0;
   //double error =  total_path_error(p, coords, verbose);
   double error =  p->total_err_function(p, coords, verbose);

   return error;
}

void my_df (const gsl_vector *v, void *void_params, gsl_vector *df)
{
   GSL_PARAMS * gsl_params = (GSL_PARAMS *) void_params;
   double * coords = gsl_params->coords;
   PARAMS * p = gsl_params->p;
   int part_num,dim;
   int i, indx=0;
   int num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
   int num;
   double dx = p->dx;
   double invdx = 1.0/dx;

   if (p->local)
   {
	   num = p->num_time_slices-1;
   }
   else
   {
	   num = 2*(num_time_slices/p->num_follow);
   }
   expand_vars(v, void_params);

   int j, ts0, time_slice;
   for(i=1; i < num; i++)
   {
      if (p->local)
      {
	      ts0 = time_slice = i;
      }
      else
      {
	      j = i/2;
	      ts0 = j*p->num_follow;
	      time_slice = ts0 + i%2;
      }
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
         for(dim=0;dim < p->num_dim;dim++)
         {
            double tmp = unpacked(coords,dim,part_num,time_slice);
	    unpacked(coords,dim,part_num,time_slice) = tmp + dx;
	    double error_pos =  p->err_function(p, coords, ts0, 0);
	    //double error_pos =  p->total_err_function(p, coords, 0);
            unpacked(coords,dim,part_num,time_slice) = tmp - dx;
            double error_neg =  p->err_function(p, coords, ts0, 0);
            //double error_neg =  p->total_err_function(p, coords, 0);
	    //double diff = error_pos-error_neg;
	    //double diff_t = error_pos_t-error_neg_t;

#ifdef FIRST_ORDER_DERIV
            double de_dx =  0.5*(error_pos - error_neg)*invdx;
#else
            unpacked(coords,dim,part_num,time_slice) = tmp + 2*dx;
            double error_pos2 =  p->err_function(p, coords, ts0, 0);
            //double error_pos2 =  p->total_err_function(p, coords, 0);
            unpacked(coords,dim,part_num,time_slice) = tmp - 2*dx;
            double error_neg2 =  p->err_function(p, coords, ts0, 0);
            //double error_neg2 =  p->total_err_function(p, coords, 0);

            double de_dx =  (-error_pos2+error_neg2 + 8*(error_pos - error_neg))/(12.0*dx);
#endif

            int indx = dim + p->num_dim*(part_num+p->num_parts*(i-1));
            gsl_vector_set(df, indx, de_dx);

            unpacked(coords,dim,part_num,time_slice) = tmp;
         }
      }
   }

}

void my_fdf (const gsl_vector *x, void *gsl_params, double *f, gsl_vector *df)
{
  *f = my_f(x, gsl_params);
  my_df(x, gsl_params, df);
}

int gsl_lagrange_min_conj(PARAMS * p, double * coords)
{
   size_t iter = 0;
   int status=-1234;
   GSL_PARAMS  gsl_params;

   const gsl_multimin_fdfminimizer_type *T;
   //const gsl_multimin_fminimizer_type *T;
   gsl_multimin_fdfminimizer *s;

   gsl_params.p= p;
   gsl_params.coords = coords;

   gsl_vector *x;
   gsl_multimin_function_fdf my_func;
   //int n = (p->num_time_slices-2)*p->num_dim*p->num_parts;
   int num_time_slices=0;
   int num;
   int n;

   if (p->local)
   {
	   num = p->num_time_slices-1;
	   n = p->num_dim*p->num_parts*(num-1);
   }
   else
   {
	   num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
	   num = 2*(num_time_slices/p->num_follow);
	   n = p->num_dim*p->num_parts*(num-1);
   }

   my_func.n = n;
   my_func.f = my_f;
   my_func.df = my_df;
   my_func.fdf = my_fdf;
   my_func.params = &gsl_params;

   x = gsl_vector_alloc(n);
   /* Starting point, x = (2.90,3.498) */

   int time_index, part_num,dim;
   int i,j, time_slice, ts0;

   for(i=1; i < num; i++)
   {
      if (p->local)
      {
	      ts0 = time_slice = i;
      }
      else
      {
	      j = i/2;
	      ts0 = j*p->num_follow;
	      time_slice = ts0 + i%2;
      }
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
	 for(dim=0;dim < p->num_dim;dim++)
	 {
	    int indx = dim + p->num_dim*(part_num+p->num_parts*(i-1));
            gsl_vector_set(x, indx, unpacked(coords,dim,part_num,time_slice));
	 }
      }
   }


   T = gsl_multimin_fdfminimizer_conjugate_fr;//seems to be the fastest
   //T = gsl_multimin_fdfminimizer_vector_bfgs2;
   //T = gsl_multimin_fdfminimizer_vector_bfgs;
   //T = gsl_multimin_fdfminimizer_conjugate_pr;
   //T = gsl_multimin_fdfminimizer_steepest_descent;
   s = gsl_multimin_fdfminimizer_alloc (T, n);

   gsl_multimin_fdfminimizer_set (s, &my_func, x, p->step_size, p->tol);

   do
   {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      if (iter%100 == 0) printf("iter = %ld, status = %d\n",iter,status);

      if (status)
      {
         printf("iter = %ld, status changed %d\n",iter,status);
         break;
      }

      //status = gsl_multimin_test_gradient (s->gradient, p->error);
      double nrm = gsl_blas_dnrm2(s->gradient);

     // if (nrm < p->f_err)
      if (s->f < p->f_err)
      {
         status =GSL_SUCCESS;
         //printf("iter = %ld, nrm = %f < f_err = %lf,total_path_error = %lf\n",iter,nrm,p->f_err,p->total_err_function(p,coords,0));
         printf("iter = %ld, nrm = %f < error_in = %lf,total_path_error = %lf,s->f=%.*e\n",iter,nrm,p->error,p->total_err_function(p,coords,0),12,s->f);
      }
      else
      {
         status = GSL_CONTINUE;
      }
      if (iter%100 == 0)
      { 
         printf("iter = %ld, nrm = %f, s->f=%.*e,  error_in = %lf\n",iter,nrm,12,s->f,p->error);
      }

      if (status == GSL_SUCCESS)
      {
         printf ("Minimum found\n");
         printf("iter = %ld, nrm = %f, s->f=%.*e,  error_in = %lf\n",iter,nrm,12,s->f,p->error);
      }

      //printf ("%5d %.5f %.5f %10.5f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->f);

   }
   while (status == GSL_CONTINUE && iter < p->num_its);

   p->min_out = s->f;

   if (p->local)
   {

	   for(time_index=1;time_index < p->num_time_slices-1;time_index++)
	   {
		   for(part_num=0;part_num < p->num_parts;part_num++)
		   {
			   for(dim=0;dim < p->num_dim;dim++)
			   {
				   int indx = dim + p->num_dim*(part_num+p->num_parts*(time_index-1));
				   unpacked(coords,dim,part_num,time_index) = gsl_vector_get(s->x, indx);
				   indx++;
			   }
		   }
	   }
   }
   else
   {
        memcpy_(coords,p->tmp_coords,size_(coords)-1);
   }


   gsl_multimin_fdfminimizer_free (s);
   gsl_vector_free (x);

   return status == GSL_SUCCESS;
}


int gsl_lagrange_min(PARAMS * p, double * coords)
{
   size_t iter = 0;
   int status;
   GSL_PARAMS  gsl_params;

   const gsl_multimin_fminimizer_type *T;
   gsl_multimin_fminimizer *s=NULL;

   gsl_params.p= p;
   gsl_params.coords = coords;

   gsl_vector *x, *ss;
   gsl_multimin_function my_func;
   int num_time_slices = p->num_follow * (p->num_time_slices/p->num_follow);
   int num = 2*(num_time_slices/p->num_follow);
   int n = p->num_dim*p->num_parts*(num-1);

   my_func.n = n;
   my_func.f = my_f;
   my_func.params = &gsl_params;

   x = gsl_vector_alloc(n);
   /* Starting point, x = (5,7) */

   int time_index, part_num,dim;
   int i=0;

   for(i=1; i < num; i++)
   {
      int j = i/2;
      int time_slice, part_num, dim;
      int ts0 = j*p->num_follow;
      time_slice = ts0 + i%2;
      for(part_num=0; part_num < p->num_parts; part_num++)
      {
	 for(dim=0;dim < p->num_dim;dim++)
	 {
	    int indx = dim + p->num_dim*(part_num+p->num_parts*(i-1));
            gsl_vector_set(x, indx, unpacked(coords,dim,part_num,time_slice));
	 }
      }
   }



   ss = gsl_vector_alloc (n);
   gsl_vector_set_all (ss, 1.0);
   T = gsl_multimin_fminimizer_nmsimplex2;
   s = gsl_multimin_fminimizer_alloc (T, n);

   gsl_multimin_fminimizer_set (s, &my_func, x, ss);

   do
   {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (iter%1000 == 0) printf("iter = %ld, status = %d\n",iter,status);

      if (status)
      {
         printf("iter = %ld, status changed %d\n",iter,status);
         break;
      }

      double size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, p->f_err);

      double nrm=0;
      if (iter%1000 == 0)
      { 
         printf("iter = %ld, num_its = %d,  size = %f,s->f=%.*e, error_in = %lf\n",iter,p->num_its,size,12,s->fval,p->error);
      }

      if (status == GSL_SUCCESS)
      {
         printf ("Minimum found\n");
         printf("iter = %ld, nrm = %f, s->f=%.*e,  error_in = %lf\n",iter,nrm,12,s->fval,p->error);
      }

      //printf ("%5d %.5f %.5f %10.5f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->fval);

   }
   while (status == GSL_CONTINUE && iter < p->num_its);

   memcpy_(coords,p->tmp_coords,size_(coords)-1);
   if (0)
   {
	   int indx = 0;
	   for(time_index=1;time_index < p->num_time_slices-1;time_index++)
	   {
		   for(part_num=0;part_num < p->num_parts;part_num++)
		   {
			   for(dim=0;dim < p->num_dim;dim++)
			   {
				   unpacked(coords,dim,part_num,time_index) = gsl_vector_get(s->x, indx);
				   indx++;
			   }
		   }
	   }
   }


   gsl_multimin_fminimizer_free (s);
   gsl_vector_free(ss);
   gsl_vector_free (x);

   return status;
}

#undef unpacked
