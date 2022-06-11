#include <stdio.h>
#include <math.h> 
#include <float.h> 
#include <time.h> 
#include "defs.h"
#include "params.h"
#include "follow.h"
#include "json.h"


#define unpacked(arr,dim,part_num,time_slice) arr[dim + p->num_dim*(part_num + p->num_parts*(time_slice))]

double rgsq(PARAMS * p, double * coords, int t)
{
	double rsq = 0.0;
	double rave = 0.0;
	int n = p->num_parts;
	int N = p->num_dim*n;
	int i;

	for(i=0;i < n;i++)
	{
		int dim;
		for(dim=0;dim < p->num_dim;dim++)
		{
		   double x = unpacked(coords,dim,i,t);
		   rsq += SQR(x);
		   rave += x;
		}
	}
	rave /= N;
	rsq /= N;
	double rg2 = rsq - SQR(rave);
	return rg2;
}

double min_rgsq(PARAMS * p, double * coords)
{
	int t;
	double min_r2 = rgsq(p,coords,0);
	for(t=1;t < p->num_time_slices;t++)
	{
		double rsq = rgsq(p,coords,t);
		min_r2 = MIN(min_r2,rsq);
	}
	return min_r2;
}

int compare(PARAMS * p, double * coords1, double * coords2)
{
	int i;
	int t = p->num_time_slices/2;
	double diff = 0.0;
	int part_indx,dim;
	double err = (1.0)*p->num_time_slices * p->f_err/(2*p->num_follow);

	for(part_indx=0;part_indx < p->num_parts;part_indx++) 
	{
		for(dim=0;dim < p->num_dim;dim++) 
		{
			double delta = unpacked(coords2,dim, part_indx,t) - unpacked(coords1,dim, part_indx,t);
			diff += SQR(delta);
		}
	}

	if (diff > 10*err/p->num_time_slices)
	{
		return 0;
	}

	diff =0.0;

	for(t=0;t < p->num_time_slices;t++) 
	{
		for(part_indx=0;part_indx < p->num_parts;part_indx++) 
		{
			for(dim=0;dim < p->num_dim;dim++) 
			{
				double delta = unpacked(coords2,dim, part_indx,t) - unpacked(coords1,dim, part_indx,t);
				diff += SQR(delta);
			}
		}
	}

	if (diff > 10*err)
	{
		return 0;
	}
	return 1;
}

int in_array(int * hay, int needle)
{
	int i;
        if (!hay)
	  return 0;
	for(i=0; i < size_(hay);i++)
	{
		if (hay[i] == needle)
			return 1;
	}
	return 0;
}



int ** unique_paths(PARAMS * p, double ** coords_arr)
{
	int i,j;
	int * unique_index_map;
	//newarr_(unique_index_map,p->num_runs);
	newarr_(unique_index_map,size_(coords_arr));
	setarr_(unique_index_map,-1);

	//for(i=0;i < p->num_runs;i++)
	for(i=0;i < size_(coords_arr);i++)
	{
		for(j=0;j < i;j++)
		{
			if (unique_index_map[j] > -1)
				continue;
			if (compare(p, coords_arr[i],coords_arr[j]))
			{
				unique_index_map[j] = i;
			}
		}
	}

/*
 * 0 1 2 3 4
 * A C B A B
 *
 * i=1
 * i=2 
 * i=3 cmp(3,0) u[0]=3
 * i=4 cmp(4,2) u[2]=4
 *
 * i=0:
 * num_c=0, pl=0
 * incr(clumps) -> num_c=1
 * u[0]=3 > -1 clumps[0]=[3,0]
 * i=1:
 * num_c=pl=1
 * j=0:if in_array(clumps[0],u[1]=-1) false pl=1 incr(clumps) num_c=2 if (u[1]=-1 > -1) False clumps[1]=[1]
 */

	int ** clumps=0;
	//for(i=0;i < p->num_runs;i++)
	for(i=0;i < size_(coords_arr);i++)
	{
		//in existing clump?
		int j;
		int num_clumps;
		if (clumps) num_clumps = size_(clumps);
		else num_clumps = 0;
		int placement = num_clumps;
		for(j=0; j < num_clumps;j++)
		{
			if (in_array(clumps[j],i))
			{
				placement = j;
				break;
			}
		}
		if (placement == num_clumps) incrarr_(clumps);
		if (unique_index_map[i] > -1) appendarr_(clumps[placement],unique_index_map[i]);	
                if (!in_array(clumps[placement],i))
		{
		   appendarr_(clumps[placement],i);	
		}
	}
	if (!clumps)
		return 0;
	if (size_(clumps) < 1)
		return 0;
	printf("num of unique minima = %ld\n",size_(clumps));
	for(i=0;i<size_(clumps);i++)
	{
	    int path_num = clumps[i][0];
	    double min_r2 = min_rgsq(p,coords_arr[path_num]);
	    printf("i %d, clump size %ld,rgsq %lf\n",i, size_(clumps[i]),min_r2);
	}
	return clumps;
}


#undef unpacked
