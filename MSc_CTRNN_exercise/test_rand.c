#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "SFMT.c"


int main ( int argc, char **argv ) {
	if ( argc < 2 )
		return -1;

	int seed = atoi(argv[1]);
	printf("seed: %d\n", seed);

	//printf("---\n");

	printf("drand48():\n");

	srand48(seed);
	double rand_1 = drand48();

	printf("normal random: %f\n", rand_1);
	printf("random between -0.5 and 0.5: %f\n", rand_1 - 0.5);
	printf("random between 0.0 and 0.1: %f\n", rand_1 / 10);

	printf("---\n");

	printf("rand():\n");

	srand(seed);
	int rand_2 = rand();

	printf("normal random: %d\n", rand_2);
	printf("random between 1 and 10: %d\n", 
			1 + (int) (10.0 * ( rand_2 / (RAND_MAX + 1.0))));
	printf("random between 0 and 1: %lf\n", 
			(double)( rand_2 / (RAND_MAX + 1.0)));

	printf("---\n");

	printf("SFMT:\n");

	init_gen_rand(seed);

	double rand_3, rand_3_1, rand_3_2;
		
	do {
		rand_3 = genrand_real3();
		printf("normal random: %lf\n", rand_3);
		rand_3_1 = genrand_real1()/10;
		printf("random [0.1,0.01]: %lf\n", rand_3_1);
		assert ( rand_3_1 >= 0.0 && rand_3_1 <= 0.1 );
		rand_3_2 = genrand_real1()/50-0.01;
		printf("random [-0.01, 0,01]: %lf\n", rand_3_2);
		assert ( rand_3_2 >= -0.01 && rand_3_2 <= 0.01 );
	} while ( 1 );

	return 0;
}

