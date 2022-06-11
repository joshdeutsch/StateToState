#include <stdio.h>
#include <float.h> 
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include "defs.h"
#include "params.h"
//#include "f_min.h"

#include "json.h"

#define unpacked(arr,dim,part_num,time_slice) arr[dim + p->num_dim*(part_num + p->num_parts*time_slice)]

void json_params_out(FILE * fout,PARAMS * p)
{
   fprintf(fout,"\{");
   fprintf(fout,"\"num_its\":%d,",p->num_its);
   fprintf(fout,"\"num_time_slices\":%d,",p->num_time_slices);
   fprintf(fout,"\"num_parts\":%d,",p->num_parts);
   fprintf(fout,"\"num_dim\":%d,",p->num_dim);
   fprintf(fout,"\"local\":%d,",p->local);
   fprintf(fout,"\"num_follow\":%d,",p->num_follow);
   fprintf(fout,"\"f_ampl\":%lf,",p->f_ampl);
   fprintf(fout,"\"eps\":%lf,",p->eps);
   fprintf(fout,"\"dx\":%lf,",p->dx);
   fprintf(fout,"\"tol\":%lf,",p->tol);
   fprintf(fout,"\"step_size\":%lf,",p->step_size);
   fprintf(fout,"\"error\":%lf,",p->error);
   fprintf(fout,"\"seed\":%d,",p->seed);
   fprintf(fout,"\"seed2\":%d,",p->seed2);
   fprintf(fout,"\"L\":%lf",p->L);
   fprintf(fout,"}");
}



char * read_to_buf(FILE * f_in, int  length)
{
	//int fd = open("filename", O_RDONLY);
	int f_d = fileno(f_in);
        length = lseek(f_d, 0, SEEK_END);
        char * mem_to_string = mmap(0, length, PROT_READ, MAP_PRIVATE, f_d, 0);
	char * whole_string=0;
	cpy_pntr_to_arr_(whole_string,mem_to_string,length);
 
	return whole_string;
}

int shift_til_char(char * in, char c)
{
	int i;
	int pos = -1;
	for(i=0;i < size_(in);i++)
	{
		if (in[i] == c)
		{
			pos = -1;
			break;
		}
	}
	return i;
}


double ** json_coords_in(FILE * f_in)
{
	double * coords=0;
	double ** coords_arr=0;
	int num_time_slices=0;
	int num_parts=0;
	int num_dim=0;
	int num_its=0;
	int i;
	int num_read=982;
        num_read = fscanf(f_in,"%*[^\"]\"num_its\":%d,\"num_time_slices\":%d,\"num_parts\":%d,\"num_dim\":%d,",&num_its,&num_time_slices,&num_parts,&num_dim);
	printf("num_read = %d, num_its = %d, num_time_slices = %d, num_parts = %d, num_dim = %d\n", num_read, num_its, num_time_slices, num_parts, num_dim);

	//if (fscanf(f_in,"%d %d [^ ]%d[^ ]%d",&num_its, &num_time_slices,&num_parts,&num_dim) != 4)
	if (num_read != 4)
	{
		printf("couldn't read in num_time_slices, num_parts, and num_dim\n)");
		exit(1);
	}

	int coord_len = num_time_slices*num_parts*num_dim;
	int result= EOF;
	while ((result = fgetc(f_in)) != EOF )
	{
		if (result == '}')
			break;
	}

	double double_in=-1e10;
	while (fscanf(f_in,"%*c%lf",&double_in) != EOF)

	{
		appendarr_(coords,double_in);
		if (size_(coords) == coord_len)
		{
			appendarr_(coords_arr,coords);
			coords = 0;
		}
	}
	return coords_arr;
}

void json_coords_out(FILE * fout, PARAMS * p, double * coords)
{
   int time_slice;
   fprintf(fout,"[");
   for(time_slice=0;time_slice < p->num_time_slices;time_slice++) 
   {
      int part_indx;
      fprintf(fout,"[");
      for(part_indx=0;part_indx < p->num_parts;part_indx++) 
      {
         int dim_indx;
         fprintf(fout,"[");
         for(dim_indx=0;dim_indx < p->num_dim;dim_indx++) 
         {
            fprintf(fout,"%.*f",DBL_DIG, unpacked(coords,dim_indx,part_indx,time_slice));
            if (dim_indx < p->num_dim-1) fprintf(fout,",");
         }
         fprintf(fout,"]");
         if (part_indx < p->num_parts-1) fprintf(fout,",");
      }
     fprintf(fout,"]");
     if (time_slice < p->num_time_slices-1) fprintf(fout,",");
   }
   fprintf(fout,"]");
}

void json_out(char * filename, PARAMS * p, double * coords)
{
   FILE * fout = fopen(filename,"w");
   fprintf(fout,"[");
   json_params_out(fout,p);
   fprintf(fout,",");
   json_coords_out(fout,p,coords);
   fprintf(fout,"]");
   fclose(fout);
}

void json_arr_out(char * filename, PARAMS * p, double ** coords_arr)
{
   FILE * fout = fopen(filename,"w");
   fprintf(fout,"[");
   json_params_out(fout,p);
   int i;
   for(i=0;i < size_(coords_arr);i++)
   {
        fprintf(fout,",");
        json_coords_out(fout,p,coords_arr[i]);
   }
   fprintf(fout,"]");
}
#undef unpacked

#if 0
int main(int argz, char * argv[])
{
	FILE * f_in = fopen("coords_out.json","r");
        double ** coords_arr =  json_coords_in(f_in);
}
#endif
