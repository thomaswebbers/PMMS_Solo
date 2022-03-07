#include "fail.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include <stdio.h>
#include <output.h>
#include <float.h>
#include <time.h>

static double *top_halo = 0;
static double *bottom_halo = 0;
clock_t start_time;
clock_t end_time;

/*
double *extract_top_halo(struct parameters *p, size_t num_colls)
{
	if (!(top_halo = calloc(p->N, sizeof(double)))) die("compute calloc");

	for(size_t i = 0; i < num_colls; i++)
	{
		top_halo[i] = p->tinit[i];	
	}	

	return top_halo;
}

double *extract_bottom_halo(struct parameters *p, size_t num_colls, size_t num_rows)
{
	
	if (!(bottom_halo = calloc(p->N, sizeof(double)))) die("compute calloc");

	for(size_t i = ((num_colls * num_rows) - num_colls); i < ((num_colls * num_rows) - 1); i++)
	{
		bottom_halo[i] = p->tinit[i];	
	}	

	return bottom_halo;
}
*/

void extract_halos_and_tinit(const struct parameters *p, double *current_matrix, double *next_matrix, size_t num_colls, size_t num_rows)
{	
	//top halo grid points for both current and next matrices
	for(size_t i = 0; i < num_colls; i++)
	{
		current_matrix[i] = p->tinit[i];	
		next_matrix[i] = p->tinit[i];	
	}

	//todo add values tinit
	for(size_t i = num_colls; i < ((num_colls * (num_rows + 2)) - num_colls); i++)
	{
		current_matrix[i] = p->tinit[i];	
	}

	//bottom halo grid points for both current and next matrices
	for(size_t i = ((num_colls * (num_rows + 2)) - num_colls); i < ((num_colls * (num_rows + 2)) - 1); i++)
	{
		current_matrix[i] = p->tinit[i];	
		next_matrix[i] = p->tinit[i];	
	}	
}

double *allocate_matrix(const struct parameters *p)
{
	static double *matrix = 0;
	if (!(matrix = calloc((p->N + 2) * p->M, sizeof(double)))) die("compute calloc");
	return matrix;
}

double compute_left_gridpoint(const struct parameters *p, double *matrix, size_t num_colls, size_t index)
{
	size_t top_left_index = index - 1;
	size_t top_middle_index = index - num_colls;
	size_t top_right_index = index - num_colls + 1;
	size_t middle_left_index = index + num_colls - 1;
	size_t middle_right_index = index + 1;
	size_t bottom_left_index = index + 2 * num_colls - 1;
	size_t bottom_middle_index = index - num_colls;
	size_t bottom_right_index = index + 1 + num_colls;

	
	double retained_temp = 0;
	double conducted_temp = 0;
	double new_temp = 0;


	retained_temp = (p->conductivity[index] * matrix[index]);
	conducted_temp = (matrix[top_middle_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[middle_left_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[middle_right_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[bottom_middle_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));

	conducted_temp = conducted_temp + (matrix[top_left_index] * (1/4) * (1/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[top_right_index] * (1/4) * (1/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[bottom_left_index] * (1/4) * (1/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[bottom_right_index] * (1/4) * (1/(sqrt(2)+1)));

	conducted_temp = ((double)(1 - p->conductivity[index])) * conducted_temp;
	new_temp = retained_temp + conducted_temp;
	return new_temp;
}

double compute_gridpoint(const struct parameters *p, double *matrix, size_t num_colls, size_t col_count, size_t index)
{
	size_t top_left_index = index - num_colls - 1;
	size_t top_middle_index = index - num_colls;
	size_t top_right_index = ((col_count + 1) % num_colls) != 0 ? index - num_colls + 1 : index - 2 * num_colls;
	size_t middle_left_index = index - 1;
	size_t middle_right_index = ((col_count + 1) % num_colls) != 0 ? index + 1 : index - num_colls + 1;
	size_t bottom_left_index = index - 1 + num_colls;
	size_t bottom_middle_index = index - num_colls;
	size_t bottom_right_index = ((col_count + 1) % num_colls) != 0 ? index + 1 + num_colls : index + 1;

	
	double retained_temp = 0;
	double conducted_temp = 0;
	double new_temp = 0;


	retained_temp = (p->conductivity[index] * matrix[index]);
	conducted_temp = (matrix[top_middle_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[middle_left_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[middle_right_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[bottom_middle_index] * (1/4) * (sqrt(2)/(sqrt(2)+1)));

	conducted_temp = conducted_temp + (matrix[top_left_index] * (1/4) * (1/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[top_right_index] * (1/4) * (1/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[bottom_left_index] * (1/4) * (1/(sqrt(2)+1)));
	conducted_temp = conducted_temp + (matrix[bottom_right_index] * (1/4) * (1/(sqrt(2)+1)));

	conducted_temp = ((double)(1 - p->conductivity[index])) * conducted_temp;
	new_temp = retained_temp + conducted_temp;
	return new_temp;
}

void print_tinit_and_conductivity(struct parameters *p, size_t num_colls, size_t num_rows)
{
	
	for (size_t i = 0; i < num_colls * num_rows; ++i)
    {
		printf("Index %li contains temperature %lf\n", i, p->tinit[i]);
	}

	for (size_t i = 0; i < num_colls * num_rows; ++i)
    {
		printf("Index %li contains conductivity value %lf\n", i, p->conductivity[i]);
	}
}

void do_compute(const struct parameters* p, struct results *r)
{
	start_time = clock();

	size_t num_rows = p->N;
	size_t num_colls = p->M;
	//double *top_halo = extract_top_halo(p, num_colls); /*halo values for top half of the cylinder */
	//double *bottom_halo = extract_bottom_halo(p, num_colls, num_rows); /* halo values for bottom half of the cylinder */

	//memory allocated only once and then re-used to reduce overhead
	double *current_matrix = allocate_matrix(p); /* temperature values for current matrix, in row-major order */
	double *next_matrix = allocate_matrix(p); /* temperature values for next matrix, in row-major order */
	//needed to swap next with current in order to save on memory allocation
	double *matrix_swap_pointer = 0;	
	
	extract_halos_and_tinit(p, current_matrix, next_matrix, num_colls, num_rows);


	//todo add iterations
	//todo add case where if this is the first iteration it fills current rather than new
	//TODO add it so that the halo rows top and bottom are only extracted when it is the first iteration

	r->tmin = DBL_MAX;
	r->tmax = DBL_MIN;
	r->maxdiff = 0;
	r->tavg = 0;
	double average_temp = 0;

	for(size_t t = 0; t < p->maxiter; t++)
	{
		size_t col_count = 0;
		double new_point_temperature = 0;
		//i starts at num_colls to skip the top halo row
		//num_rows + 1 to skip the bottom halo row
		for (size_t i = num_colls; i < num_colls * (num_rows + 1); ++i)
    	{
			if((i % num_colls) == 0)
			{
				new_point_temperature = compute_left_gridpoint(p, current_matrix, num_colls, i);
			} else {
				new_point_temperature = compute_gridpoint(p, current_matrix, num_colls, col_count, i);			
			}
			
			if (new_point_temperature < r->tmin) r->tmin = new_point_temperature;	
			if (new_point_temperature > r->tmax) r->tmax = new_point_temperature;	
			if(abs(new_point_temperature - current_matrix[i]) > r->maxdiff) 
			{
				r->maxdiff = (abs(new_point_temperature) - current_matrix[i]);
			}
			average_temp = average_temp + new_point_temperature;

			next_matrix[i] = new_point_temperature;

			col_count = (col_count + 1) % num_colls;
		}
		
		//basic pointer swap in order to save to prevent re-allocation for every iteration
		matrix_swap_pointer = current_matrix;
		current_matrix = next_matrix;
		next_matrix = matrix_swap_pointer;

		r->tavg = average_temp / (p->N * p->M);


		//every k iterations
		if(t % p->period == 0 && t > 0)
		{
			r->niter = t;
			report_results(p, r);
		}
	}

	end_time = clock();
	clock_t time_taken = end_time - start_time;
	double time_in_secs = ((double)time_taken)/CLOCKS_PER_SEC;
	r->time = time_in_secs;

	//todo add timing to the final result 
}


