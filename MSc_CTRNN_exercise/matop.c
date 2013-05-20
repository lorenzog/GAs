#ifndef _MATOP_C
#define _MATOP_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SFMT.h"

#include "matop.h"

#ifdef _WITH_DEBUG
#define DBG(x);		x
#else
#define DBG(x);
#endif

/// vector operations

inline void zero_vector (int size, double vector[size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		vector[i] = 0.0;
}

/**
 * randomize vectors [-1,1]
 */
inline void random_vector (int size, double vector[size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		//vector[i] = drand48();
		//vector[i] = drand48()-0.5;
		// [-1,1]
		vector[i] = genrand_real1() + genrand_real1() - 1;
}

inline void flatten_vector(double MAX, double MIN, int size, 
		double vector[size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		if ( vector[i] > MAX )
			vector[i] = MAX;
		else
			if ( vector[i] < MIN )
				vector[i] = MIN;
}

inline void normalize_vector(double ZERO, int size, double vector[size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		if ( fabs(vector[size]) < ZERO )
			vector[size] = 0.0;
}

inline void copy_vector(int size, double src[size], double dest[size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		dest[i] = src[i];
}

/**
 * used to sort the vectors; receives pointers to pointers to objects
 */
static int cmpdbl(const void *a, const void *b) {
	//printf("a3: %f\n", *(double*)a);
	// cast and dereferencing
	if ( *(double*)a < *(double*)b )
		return -1;
	else
		if ( *(double*)a > *(double*)b )
			return 1;
		else
			return 0;
}

inline void sort_vector(int size, double vector[size]) {
	int i;
	qsort(vector, size, sizeof(double), cmpdbl);
}

inline void print_vector (int size, double vector[size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		printf("%f ", vector[i]);
	printf("\n");
}

/// matrix operations

void zero_matrix (int size, double matrix[size][size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		zero_vector(size, matrix[i]);
}

void random_matrix (int size, double matrix[size][size])  {
	int i;
	for ( i = 0 ; i < size ; i++ )
		random_vector(size, matrix[i]);
}

void flatten_matrix(double MAX, double MIN, int size, 
		double matrix[size][size]) {
	int i;
	for ( i = 0 ; i < size; i++ )
		flatten_vector(MAX, MIN, size, matrix[i]);
}

void normalize_matrix(double ZERO, int size, double matrix[size][size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		normalize_vector(ZERO, size, matrix[i]);
}

void copy_matrix(int size, double src[size][size], double dest[size][size]) {
	int i;
	for ( i = 0 ; i < size ; i++ )
		copy_vector(size, src[i], dest[i]);
}
void print_matrix(int size, double matrix[size][size]) {
	int i, j;
	for ( i = 0 ; i < size ; i++ )
		print_vector(size, matrix[i]);
}

#endif
