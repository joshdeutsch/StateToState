#include <stdio.h>
#include <math.h> 
#include <float.h> 
#include <time.h> 
#include "defs.h"
#include "params.h"
#include "follow.h"
#include "json.h"
#include "analyze.h"

char **  main_args(int * argc)
{
   FILE * args_in = fopen("params","r");
   char ** argv = 0;
   int i=1, j=0;
   int last_break=0;
   newarr_(argv,2);
   newarr_(argv[i],1);


   while (1)
   {
      char c = getc(args_in);
      if (c == EOF) break;
      if ((c == ' ' || c == '\t' || c == '\n'))
      {
         if (!last_break)
         {
            incrarr_(argv[i]);
            argv[i][j] = '\0';
            i++;
            incrarr_(argv);
            j=0;
//            printf("i=%d\n",i);
            last_break=1;
         }
      }
      else
      {
         last_break=0;
         incrarr_(argv[i]);
         argv[i][j]=c;
         j++;
//         printf("j=%d\n",j);
      }
   }
   argc[0] = i;
   return argv;
}

#define unpacked(arr,dim,part_num,time_slice) arr[dim + p->num_dim*(part_num + p->num_parts*(time_slice))]
double * init_coords_slice(PARAMS * p)
{
   int part_indx, i=0;
   int dim_indx;
   double * x;
   newarr_(x,(p->num_dim*p->num_parts));

   for(part_indx=0;part_indx < p->num_parts;part_indx++) 
   {
      for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
      {
         int time_slice=0;
         x[i] = p->L*ran();
	 i++;
      }
   }
   return x;
}

double * init_coords(PARAMS * p,int ran_offset, double * x_i, double * x_f)
{
   double * coords;
   int arr_size = p->num_time_slices*p->num_dim*p->num_parts;
   printf("in init_coords: p->num_time_slices = %d, arr_size = %d\n",p->num_time_slices,arr_size);
   newarr_(coords,arr_size);
   int part_indx, i=0;
   int created_i=0;
   int created_f=0;
   int dim_indx;

   if (x_i == 0) 
   {
	   x_i = init_coords_slice(p);
	   created_i = 1;
   }
   if (x_f == 0) 
   {
	   x_f = init_coords_slice(p);
	   created_f = 1;
   }

   srandom(p->seed2+ran_offset);
   i=0;
   int time_slice;
   for(part_indx=0;part_indx < p->num_parts;part_indx++) 
   {
      for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
      {

         //srandom(p->seed2);
	 if (0)
	 {
		 int n =  p->num_follow * (p->num_time_slices/p->num_follow);
		 double dt = 1.0/n;
		 double omega = sqrt(-p->f_ampl);
		 //2(1-cos(wdt))=k dt^2 ->  cos(wdt) = 1- k dt^2/2, w^2 dt^2 = k dt^2
		 omega = (1.0/dt)*acos(1+0.5*p->f_ampl*dt*dt);
		 double ampl = x_f[i]/sin(dt*(p->num_time_slices-1)*omega);
		 for(time_slice=0;time_slice < p->num_time_slices;time_slice++) 
		 {
			 unpacked(coords,dim_indx,part_indx,time_slice) = ampl*sin(dt*omega * time_slice);
		 }
	 }
	 else
	 {
	    for(time_slice=0;time_slice < p->num_time_slices;time_slice++) 
	    {
	       double x  = time_slice*(x_f[i]-x_i[i])/(p->num_time_slices-1) + x_i[i];
	       x  = p->L*ran();
	       //printf("x = %lf\n",x);
	       unpacked(coords,dim_indx,part_indx,time_slice) = x;
	    }
	 }
	 unpacked(coords,dim_indx,part_indx,0) = x_i[i];
	 //printf("dim_indx=%d,part_indx=%d,p->num_time_slices=%d\n",
         unpacked(coords,dim_indx,part_indx,p->num_time_slices-1) = x_f[i];
         i++;
      }
   }
   p->invL = 1.0/p->L;
   if (created_i) freearr_(x_i);
   if (created_f) freearr_(x_f);
   return coords;
}

void rand_in_circle(int num_dim, double * center, double radius, double * r_out)
{
   int dim_indx;
   while (1)
   {
      double sum=0;
      for(dim_indx=0;dim_indx < num_dim; dim_indx++)
      {
         r_out[dim_indx] = 2*ran()-1.0;
         sum +=SQR(r_out[dim_indx]);
      }
      if (sum < 1)
         break;
   }
      for(dim_indx=0;dim_indx < num_dim; dim_indx++)
      {
         r_out[dim_indx] = radius*r_out[dim_indx] + center[dim_indx];
      }
}

void make_time_slice(PARAMS * p, int n_begin, int n_end, double * center, double radius, double ** coord_t_slice)
{
   int part_num; 
   for(part_num = n_begin; part_num < n_end;part_num++)
   {
      rand_in_circle(p->num_dim, center, radius,  coord_t_slice[part_num]);
   }
}

void make_init_final_coords(PARAMS * p, double ** c_i, double ** c_f)
{
   int n= p->num_parts;
   int n2 = p->num_parts/2;
   int n4 = n2/2;
   int n34 = 3*n4;
   double * center_1;
   double * center_2;
   double r_small = 0.1*p->L;
   double r_big = 0.5*p->L;


   alloca_(center_1,p->num_dim);
   alloca_(center_2,p->num_dim);

   double c_x = p->L*0.25;
   center_1[0] = c_x;
   center_1[1] = c_x;

   c_x = p->L*0.75;
   center_2[0] = c_x;
   center_2[1] = c_x;

   make_time_slice(p, 0, n2, center_1, r_small, c_i);
   make_time_slice(p, n2, p->num_parts, center_2, r_big, c_i);
   make_time_slice(p, 0, n2, center_1, r_big, c_f);
   //make_time_slice(p, 0, n4, center_1, r_small, c_f);
   //make_time_slice(p, n4, n2, center_1, r_big, c_f);
   make_time_slice(p, n2, p->num_parts, center_2, r_small, c_f);

}

double * init_coords_2_circs(PARAMS * p)
{
   double ** c_i, ** c_f;
   int arr_size = p->num_time_slices*p->num_dim*p->num_parts;
   double * coords;
   newarr_(coords,arr_size);
   alloca_2d(c_i,p->num_parts,p->num_dim);
   alloca_2d(c_f,p->num_parts,p->num_dim);
   make_init_final_coords(p, c_i, c_f);
   int dim_indx, part_indx;

   int n2 = p->num_parts/2;
#if 0
   for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
   {
      SWAP(c_i[1][dim_indx],c_i[n2+1][dim_indx]);
      SWAP(c_f[0][dim_indx],c_f[n2][dim_indx]);
   }
#endif

         

   for(part_indx=0;part_indx < p->num_parts;part_indx++)
   {
      for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
      {
         int time_slice=0;
         double x_i = c_i[part_indx][dim_indx];
         double x_f = c_f[part_indx][dim_indx];
         for(time_slice=0;time_slice < p->num_time_slices;time_slice++) 
         {
            double x  = time_slice*(x_f-x_i)/(p->num_time_slices-1) + x_i;
            unpacked(coords,dim_indx,part_indx,time_slice) = x;
         }
      }
   }

   p->invL = 1.0/p->L;
   return coords;
}

double * init_coords_2_circs_1(PARAMS * p)
{
   double * coords;
   double * center_i;
   double * center_f;
   double * pos_i;
   double * pos_f;
   double r_init,r_fin, c_init_x,c_fin_x;
   int arr_size = p->num_time_slices*p->num_dim*p->num_parts;
   newarr_(coords,arr_size);
   int part_indx;
   alloca_(center_i,p->num_dim);
   alloca_(center_f,p->num_dim);
   alloca_(pos_i,p->num_dim);
   alloca_(pos_f,p->num_dim);
   for(part_indx=0;part_indx < p->num_parts;part_indx++) 
   {
      int dim_indx;
      if (part_indx < p->num_parts/2)
      {
         r_init = p->L * 0.5;
         c_init_x = p->L*0.25;
         r_fin= p->L * 0.1;
         c_fin_x = p->L * 0.10;
      }
      else
      {
         r_init = p->L * 0.1;
         c_init_x = p->L*0.75;
         r_fin= p->L * 0.5;
         c_fin_x = p->L * 0.75;
      }


      for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
      {
         center_i[dim_indx] = c_init_x;
         center_f[dim_indx] = c_fin_x;
      }
      rand_in_circle(p->num_dim, center_i, r_init,  pos_i);
      rand_in_circle(p->num_dim, center_f, r_fin,  pos_f);

      for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
      {
         int time_slice=0;
         double x_i = pos_i[dim_indx];
         double x_f = pos_f[dim_indx];
         for(time_slice=0;time_slice < p->num_time_slices;time_slice++) 
         {
            double x  = time_slice*(x_f-x_i)/(p->num_time_slices-1) + x_i;
            unpacked(coords,dim_indx,part_indx,time_slice) = x;
         }
      }
   }
   p->invL = 1.0/p->L;
   return coords;
}

#undef unpacked

int main(int argc1, char * argv1[])
{
   int i,argc;
   int all_succeed = 0;
   int local = 0;
   int num_follow = 0;
   PARAMS * p;
   newarr_(p,1);
   char ** argv =  main_args(&argc);

   clock_t begin = clock();
   for (i = 1; i < argc; )
   {
      if (argv[i][0] == '#')
      {
         i++;
      }
      else if (strstr (argv[i], "-num"))
      {
         if (strstr (argv[i], "part"))
         {
            if (sscanf(argv[i+1],"%d", &p->num_parts) != 1) 
            {
               printf("num_parts messed up\n");
               exit(2); 
            }
         }
         else if (strstr (argv[i], "tim") && strstr (argv[i], "sli"))
         {
            if (sscanf(argv[i+1],"%d", &p->num_time_slices) != 1) 
            {
               printf("num_time_slices messed up\n");
               exit(2); 
            }
         }
         else if (strstr (argv[i], "fol") && strstr (argv[i], "ow"))
         {
            if (sscanf(argv[i+1],"%d", &p->num_follow) != 1) 
            {
               printf("num_follow messed up\n");
               exit(2); 
            }
            num_follow = p->num_follow;
         }
	 if (strstr (argv[i], "run"))
         {
            if (sscanf(argv[i+1],"%d", &p->num_runs) != 1)
            {
               printf("num_runs messed up\n");
               exit(2);
            }
         }

         if (strstr (argv[i], "dim"))
         {
            if (sscanf(argv[i+1],"%d", &p->num_dim) != 1) 
            {
               printf("num_dim messed up\n");
               exit(2); 
            }
         }
         if (strstr (argv[i], "its"))
         {
            if (sscanf(argv[i+1],"%d", &p->num_its) != 1) 
            {
               printf("num_its messed up\n");
               exit(2); 
            }
         }
         i += 2;
      }
      else if (strstr (argv[i], "-all") && strstr(argv[i],"suc"))
      {
         all_succeed=1;
         i += 1;
      }
      else if (!strcmp (argv[i], "-local"))
      {
         local = p->local=1;
         i += 1;
      }


      else if (strstr (argv[i], "-mont") && strstr(argv[i],"it"))
      {
         if (sscanf(argv[i+1],"%d", &p->monte_its) != 1) 
         {
            printf("monte_its messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (strstr (argv[i], "-rela") && strstr(argv[i],"it"))
      {
         if (sscanf(argv[i+1],"%d", &p->relax_its) != 1) 
         {
            printf("relax_its messed up\n");
            exit(2); 
         }
         i += 2;
      }


      else if (!strcmp (argv[i], "-seed"))
      {
         if (sscanf(argv[i+1],"%d", &p->seed) != 1) 
         {
            printf("seed messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (!strcmp (argv[i], "-seed2"))
      {
         if (sscanf(argv[i+1],"%d", &p->seed2) != 1) 
         {
            printf("seed2 messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (strstr (argv[i], "-f") && strstr (argv[i], "amp"))
      {
         if (sscanf(argv[i+1],"%lf", &p->f_ampl) != 1) 
         {
            printf("f_ampl  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (strstr (argv[i], "-step") && strstr (argv[i], "size"))
      {
         if (sscanf(argv[i+1],"%lf", &p->step_size) != 1) 
         {
            printf("step_size  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (!strcmp (argv[i], "-tol"))
      {
         if (sscanf(argv[i+1],"%lf", &p->tol) != 1) 
         {
            printf("tol  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (strstr (argv[i], "-f") && strstr(argv[i],"err"))
      {
         if (sscanf(argv[i+1],"%lf", &p->f_err) != 1) 
         {
            printf("f_err  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (!strcmp (argv[i], "-error"))
      {
         if (sscanf(argv[i+1],"%lf", &p->error) != 1) 
         {
            printf("eps  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (strstr (argv[i], "-beta") && strstr(argv[i], "ini"))
      {
         if (sscanf(argv[i+1],"%lf", &p->beta_init) != 1) 
         {
            printf("beta_init  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (!strcmp (argv[i], "-dx"))
      {
         if (sscanf(argv[i+1],"%lf", &p->dx) != 1) 
         {
            printf("dx messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (!strcmp (argv[i], "-eps"))
      {
         if (sscanf(argv[i+1],"%lf", &p->eps) != 1) 
         {
            printf("eps  messed up\n");
            exit(2); 
         }
         i += 2;
      }
      else if (!strcmp (argv[i], "-L"))
      {
         if (sscanf(argv[i+1],"%lf", &p->L) != 1) 
         {
            printf("L messed up\n");
            exit(2); 
         }
         i += 2;
      }


      else
      {
         printf("don't understand command line: %s\n",argv[i]);
         exit(3);
      }
   }
   p->num_runs_in = p->num_runs;
   double eps = p->eps;
   int num_time_slices = p->num_time_slices;
   if (1)
   {
	   num_time_slices = 2 +  p->num_follow * (p->num_time_slices/p->num_follow);
   }
   int num_successes=0;
   if (num_time_slices != p->num_time_slices)
   {
      printf("changing num_time_slices from %d to %d\n",p->num_time_slices, num_time_slices);
      p->num_time_slices = num_time_slices;
   }
   set_functions(p);

     double ** coords_arr=0;
     //newarr_(coords_arr,p->num_runs);
     srandom(p->seed);
     double * coords_beg = init_coords_slice(p);
     double * coords_end = init_coords_slice(p);
     double * coords = init_coords(p,0, coords_beg, coords_beg);
     monte_relax_ends(p, coords, u_tot_verbose,0);

     char * filename = "coords_out.json";
     char * filestatsname = "stats_out.d";
     FILE * statsfile= fopen(filestatsname,"w");

     for(i=0; i < p->num_runs_in; i++)
     {
        p->local = local;
        p->num_follow = num_follow;
	double * coords = init_coords(p,1000*(i+1), coords_beg, coords_beg);
	//double * coords = init_coords_2_circs(p);
	init_relax(p, coords, 0.2);
	double rgsq_init = rgsq(p,coords,0);
	memcpyarr_(p->tmp_coords,coords);
	printf("size of coords = %ld\n",(size_(coords)));

	relax(p, coords, 0);
	int success =  gsl_lagrange_min_conj  (p, coords);
	if (1)
	{

		if (p->local)
		{
		    printf("changing from local algorithm to follow algorithm\n");
		    p->local = 0;
		    p->num_follow = 2;
                    set_functions(p);
		}
	}
	while (0 && (!p->local && !success && p->num_follow < p->num_time_slices/2))
	{
		printf("***num_follow changed from %d to %d running:\n",p->num_follow, p->num_follow*2);
		p->num_follow *= 2;
	        success =  gsl_lagrange_min_conj  (p, coords);
	}
        //gsl_lagrange_min(p, coords);//slower

	//   for(i=0;i < size_(coords);i++) printf("%lf\n",coords[i]);
	//json_out(filename,p,coords);
	if (success || p->num_runs_in == 1 || all_succeed)
	{
		appendarr_(coords_arr,coords);
		//coords_arr[num_successes] = coords;

		if (success) num_successes++;
	}
	double min_rg = min_rgsq(p,coords);
//	double min_rgsq(PARAMS * p, double * coords);
	printf("i = %d, rgsq for t=0 %lf, min_f = %lf, min_rg = %lf\n",i, rgsq_init,p->min_out,min_rg);
	fprintf(statsfile, "i = %d, rgsq for t=0 %lf, min_f = %.*e, min_rg = %lf\n",i, rgsq_init,14,p->min_out,min_rg);
     }

     p->num_runs = num_successes;
     p->local = local;
     p->num_follow = num_follow;
     //resizearr_(coords_arr,p->num_runs);
     unique_paths(p,coords_arr);
     printf("number of successful runs = %d\n",p->num_runs);
     json_arr_out(filename, p, coords_arr);

     free2darr_(coords_arr);

     clock_t end = clock();
     double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
     printf("total time = %lf sec\n",time_spent);
     //fclose(filename);
     fclose(statsfile);

}

