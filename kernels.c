/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* 
 * Please fill in the following team struct 
 */
team_t team = {
    "hhh",              /* Team name */

    "Donkey-Hao",     /* First member full name */
    "19307130267@fudan.edu.cn",  /* First member email address */

    "",                   /* Second member full name (leave blank if none) */
    ""                    /* Second member email addr (leave blank if none) */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: 4_Loop";
void rotate(int dim, pixel *src, pixel *dst) 
{

	int i, j;
	int tlimi = dim - 3;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < tlimi; i+=4){
			dst[i+k] = src[i*dim+j];
			dst[ i+1+k] = src[(i+1)*dim+j];		
			dst[ i+2+k] = src[(i+2)*dim+j];
			dst[ i+3+k] = src[(i+3)*dim+j];
		}
		for(;i<dim;++i){
			dst[i+k] = src[i*dim+j];  
		}
	}
}
char rotate_NoRIDX_descr[] = "rotate: no ridx & just change i&j";
void rotate_NoRIDX(int dim, pixel *src, pixel *dst) 
{

	int i,j;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < dim; i++){
			 dst[i+k] = src[i*dim+j];  }
	}
}
char rotate_6_descr[] = "rotate: 6_Loop";
void rotate_6(int dim, pixel *src, pixel *dst) 
{

	int i, j;
	int tlimi = dim - 5;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < tlimi; i+=6){
			dst[i+k] = src[i*dim+j];
			dst[ i+1+k] = src[(i+1)*dim+j];		
			dst[ i+2+k] = src[(i+2)*dim+j];
			dst[ i+3+k] = src[(i+3)*dim+j];
			dst[ i+4+k] = src[(i+4)*dim+j];
			dst[ i+5+k] = src[(i+5)*dim+j];
		//	dst[ i+6+k] = src[(i+6)*dim+j];
		//	dst[ i+7+k] = src[(i+7)*dim+j];
		//	dst[ i+8+k] = src[(i+8)*dim+j];
		//	dst[ i+9+k] = src[(i+9)*dim+j];
		//	dst[ i+10+k] = src[(i+10)*dim+j];
		//	dst[ i+11+k] = src[(i+11)*dim+j];
		//	dst[ i+12+k] = src[(i+12)*dim)+j];
			
		}
		for(;i<dim;++i){
			dst[i+k] = src[i*dim+j];  
		}
	}
}
char rotate_8_descr[] = "rotate: 8_Loop";
void rotate_8(int dim, pixel *src, pixel *dst) 
{

	int i, j;
	int tlimi = dim - 7;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < tlimi; i+=8){
			dst[i+k] = src[i*dim+j];
			dst[ i+1+k] = src[(i+1)*dim+j];		
			dst[ i+2+k] = src[(i+2)*dim+j];
			dst[ i+3+k] = src[(i+3)*dim+j];
			dst[ i+4+k] = src[(i+4)*dim+j];
			dst[ i+5+k] = src[(i+5)*dim+j];
			dst[ i+6+k] = src[(i+6)*dim+j];
			dst[ i+7+k] = src[(i+7)*dim+j];
		//	dst[ i+8+k] = src[(i+8)*dim+j];
		//	dst[ i+9+k] = src[(i+9)*dim+j];
		//	dst[ i+10+k] = src[(i+10)*dim+j];
		//	dst[ i+11+k] = src[(i+11)*dim+j];
		//	dst[ i+12+k] = src[(i+12)*dim)+j];
			
		}
		for(;i<dim;++i){
			dst[i+k] = src[i*dim+j];  
		}
	}
}
char rotate_10_descr[] = "rotate: 10_Loop";
void rotate_10(int dim, pixel *src, pixel *dst) 
{

	int i, j;
	int tlimi = dim - 9;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < tlimi; i+=10){
			dst[i+k] = src[i*dim+j];
			dst[ i+1+k] = src[(i+1)*dim+j];		
			dst[ i+2+k] = src[(i+2)*dim+j];
			dst[ i+3+k] = src[(i+3)*dim+j];
			dst[ i+4+k] = src[(i+4)*dim+j];
			dst[ i+5+k] = src[(i+5)*dim+j];
			dst[ i+6+k] = src[(i+6)*dim+j];
			dst[ i+7+k] = src[(i+7)*dim+j];
			dst[ i+8+k] = src[(i+8)*dim+j];
			dst[ i+9+k] = src[(i+9)*dim+j];
		//	dst[ i+10+k] = src[(i+10)*dim+j];
		//	dst[ i+11+k] = src[(i+11)*dim+j];
		//	dst[ i+12+k] = src[(i+12)*dim)+j];
			
		}
		for(;i<dim;++i){
			dst[i+k] = src[i*dim+j];  
		}
	}
}
char rotate_12_descr[] = "rotate: 12_Loop";
void rotate_12(int dim, pixel *src, pixel *dst) 
{

	int i, j;
	int tlimi = dim - 11;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < tlimi; i+=12){
			dst[i+k] = src[i*dim+j];
			dst[ i+1+k] = src[(i+1)*dim+j];		
			dst[ i+2+k] = src[(i+2)*dim+j];
			dst[ i+3+k] = src[(i+3)*dim+j];
			dst[ i+4+k] = src[(i+4)*dim+j];
			dst[ i+5+k] = src[(i+5)*dim+j];
			dst[ i+6+k] = src[(i+6)*dim+j];
			dst[ i+7+k] = src[(i+7)*dim+j];
			dst[ i+8+k] = src[(i+8)*dim+j];
			dst[ i+9+k] = src[(i+9)*dim+j];
			dst[ i+10+k] = src[(i+10)*dim+j];
			dst[ i+11+k] = src[(i+11)*dim+j];
		//	dst[ i+12+k] = src[(i+12)*dim)+j];
			
		}
		for(;i<dim;++i){
			dst[i+k] = src[i*dim+j];  
		}
	}
}
char rotate_16_descr[] = "rotate: 16_Loop";
void rotate_16(int dim, pixel *src, pixel *dst) 
{

	int i, j;
	int tlimi = dim - 15;
	for (j = 0; j < dim; j++){
		int k = dim - 1 - j;
		k = k*dim;
		for (i = 0; i < tlimi; i+=16){
			dst[i+k] = src[i*dim+j];
			dst[ i+1+k] = src[(i+1)*dim+j];		
			dst[ i+2+k] = src[(i+2)*dim+j];
			dst[ i+3+k] = src[(i+3)*dim+j];
			dst[ i+4+k] = src[(i+4)*dim+j];
			dst[ i+5+k] = src[(i+5)*dim+j];
			dst[ i+6+k] = src[(i+6)*dim+j];
			dst[ i+7+k] = src[(i+7)*dim+j];
			dst[ i+8+k] = src[(i+8)*dim+j];
			dst[ i+9+k] = src[(i+9)*dim+j];
			dst[ i+10+k] = src[(i+10)*dim+j];
			dst[ i+11+k] = src[(i+11)*dim+j];
			dst[ i+12+k] = src[(i+12)*dim+j];
			dst[ i+13+k] = src[(i+13)*dim+j];
			dst[ i+14+k] = src[(i+14)*dim+j];
			dst[ i+15+k] = src[(i+15)*dim+j];

			
		}
		for(;i<dim;++i){
			dst[i+k] = src[i*dim+j];  
		}
	}
}
char rotate_16_block_descr[] = "rotate: 16-block";
void rotate_16_block(int dim, pixel *src, pixel *dst)
{
    int i, j, a, b;
    int sdim = dim - 1;
    for (i = 0; i < dim; i += 16)
    {
        for (j = 0; j < dim; j += 16)
        {	
			int ki = i+16;
            for (a = i; a < ki; a++)
            {
				int kj = j+16;
				int ak = a * dim;
                for (b = j; b < kj; b++)
                {
					int k = sdim - b;
					k = k * dim;
                    dst[ k + a ] = src[ak+b];
                }
            }
        }
    }
}
/*********************************************************************
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_rotate_functions() 
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);   
   	add_rotate_function(&rotate, rotate_descr);   
	add_rotate_function(&rotate_NoRIDX, rotate_NoRIDX_descr);   
	add_rotate_function(&rotate_6, rotate_6_descr); 
	add_rotate_function(&rotate_8, rotate_8_descr); 
	add_rotate_function(&rotate_10, rotate_10_descr); 
	add_rotate_function(&rotate_12, rotate_12_descr); 
	add_rotate_function(&rotate_16, rotate_16_descr); 
	add_rotate_function(&rotate_16_block, rotate_16_block_descr); 
   /*  ... Register additional test functions here */
}


/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int min(int a, int b) { return (a < b ? a : b); }
static int max(int a, int b) { return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
    sum->red += (int) p.red;
    sum->green += (int) p.green;
    sum->blue += (int) p.blue;
    sum->num++;
    return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
    current_pixel->red = (unsigned short) (sum.red/sum.num);
    current_pixel->green = (unsigned short) (sum.green/sum.num);
    current_pixel->blue = (unsigned short) (sum.blue/sum.num);
    return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;

    initialize_pixel_sum(&sum);
    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) 
	for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) 
	    accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

/*
 * smooth - Your current working version of smooth. 
 * IMPORTANT: This is the version you will be graded on
 */
char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst) 
{
    naive_smooth(dim, src, dst);
}


/********************************************************************* 
 * register_smooth_functions - Register all of your different versions
 *     of the smooth kernel with the driver by calling the
 *     add_smooth_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_smooth_functions() {
    add_smooth_function(&smooth, smooth_descr);
    add_smooth_function(&naive_smooth, naive_smooth_descr);
    /* ... Register additional test functions here */
}

