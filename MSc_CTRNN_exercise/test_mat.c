#include <stdio.h>
#include "matop.h"

#include "SFMT.c"

int main ( int argc, char **argv ) {
	if ( argc < 2 )
		return -1;

	int seed = atoi(argv[1]);
	init_gen_rand(seed);
	/*
	double A[3][3];
	random_matrix(3, A);
	printf("A is:\n");
	print_matrix(3, A);
	double B[3][3];
	copy_matrix(3, A, B);
	printf("B is:\n");
	print_matrix(3, B);
	*/

	double X[10];
	random_vector(10, X);
	print_vector(10, X);
	printf("sorted:\n");
	sort_vector(10, X);
	print_vector(10, X);

	printf("\n");
	return 0;
}
