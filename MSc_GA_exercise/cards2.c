/*
 * genetic algorithm exercise
 * Oct 13, 2007
 *
 * 10 cards num. 1 -> 10
 * find a way to
 * divide into 2 piles
 * pile_0 sum near to 36
 * pile_1 multiply near to 360
 *
 * lorenzo grespan
 * lg80@susx.ac.uk
 *
 * compile with: gcc -std=c99 -lm -g -o cards2 cards2.c
 * use (in BASH+GNU/Linux) with: ./cards2 $RANDOM 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define GENOTYPE_LENGTH 10
#define MAX_MUTATIONS	1024	// 2^GENOTYPE_LENGTH
#define TARGET_PILE0	36
#define TARGET_PILE1	360.0

// define as x to see debugging output
#define DBG(x)		
// a little bit more nice output
// (commented in source code to test different order in the fitness functions)
#define INFO(x)		

/**
 * evaluates the individual for the SUM
 */
inline int sum(int *individual) {
		int sum = 0, i;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				// if true return the value + 1 (1 .. 10)
				// note: if not true, add zero
				sum += ( individual[i] == 1 ? (i+1) : 0 );
		return sum;
}

/**
 * evaluates an individual for the PROD
 */
inline int prod(int *individual) {
		int prod = 1, i;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				// if true return the value + 1 (1 .. 10)
				// note: if not true, multiply by 1 (not zero!)
				prod *= ( individual[i] == 1 ? (i+1) : 1 );
		// watcha!
		if ( prod == 1 && individual[0] == 0 )
				prod = 0;

		return prod;
}

/**
 * generates the bit-flip version of a genotype
 */
inline void bitflip(int *original, int *bitflipped) {
		int i;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				bitflipped[i] = ( original[i] == 0 ? 1 : 0 );
}

/**
 * randomly initializes a genotype
 */
inline void rand_init(int *individual) {
		int i;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				individual[i] =  ( rand() % 2 );
}

/**
 * prints the phenotype
 */
inline void print_individual(int *ind) {
		int i;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				if ( ind[i] == 1 )
						printf("%d ", i+1);
				else
						printf("- ");
}

/**
 * basic source-to-dest copy
 */
inline void copy_individual (int *source, int *dest) {
		int i;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				dest[i] = source[i];
}

/*
 * flip one specific bit
 */
inline void bitflip1(int *individual, int index) {
		individual[index] = ( individual[index] == 0 ? 1 : 0 );
}

/**
 * calculate the hash of the genotype as the value (from binary to dec)
 */
inline int hash(int *individual) {
		int i, hashvalue = 0;
		for ( i = 0 ; i < GENOTYPE_LENGTH ; i++ )
				if ( individual[i] == 1 )
						hashvalue += exp2(i);
		return hashvalue;
}

/**
 * mutate the individual with a new one taken from the population
 * but just 1 mutation-step away
 * 
 * @return -1 if no new individual can be found which is 1 gene away
 */
short mutate_1gene(int *individual, int *population) {
		int test_mutation[GENOTYPE_LENGTH];

		// make a backup copy
		copy_individual(individual, test_mutation);

		// breed the possible population
		int max_mutations = GENOTYPE_LENGTH-1;
		int mutations[max_mutations][GENOTYPE_LENGTH];
		int possible_mutations = 0 ;

		DBG(printf("original:\t");print_individual(test_mutation);printf("\n");)

		int i = 0, j;

		// for every bit of the original individual
		for ( j = 0 ; j < GENOTYPE_LENGTH ; j++ )
		{
				// mutate a specific gene
				bitflip1(test_mutation, j);
				// check if such individual has been tested already for fitness
				if ( population[hash(test_mutation)] == 0 )
				{
						// if it's not yet tested put into the pool of possible mutations
						copy_individual(test_mutation, mutations[i]);
						DBG(printf("mutations[%d]:\t", i);print_individual(test_mutation);printf("\n");)
						i++;
				}
				// get the original individual back
				copy_individual(individual, test_mutation);
		}

		DBG(printf("the mutation pool has size of: %d\n", i);)

		if ( i == 0 )
				// no new individual can be found
				return -1;

		// select one random "new" individual
		int new_individual_index = rand() % i;

		// that's our man
		copy_individual(mutations[new_individual_index], individual);

		DBG(printf("chosen mut:\t");print_individual(individual);printf("\n");)

		return 0;
}


/**
 * the actual fitness function
 *
 */
float eval_fit(int* individual_pile0, int* individual_pile1) {
		int sum_val = sum(individual_pile0);
		int prod_val = prod(individual_pile1);

		// partial fitnesses
		float f0, f1;

		// v3
		f0 = log(1+abs(sum_val - TARGET_PILE0));
		f1 = log(1+abs(prod_val - TARGET_PILE1));
		INFO(printf("sum: %d fit: %f 1-fit: %f\n", sum_val, f0, 1-f0);)
		INFO(printf("prod: %d fit: %f 1-fit: %f\n", prod_val, f1, 1-f1);)


		return ( f0 + f1 );
		//return ((1-f0)+(1-f1));
}

/**
 * GA with a single fitness function
 */
int evolve_both(int* first, int* second, int* winner, int* loser, int* population) {
		int round = 0;
		float fitness_first, fitness_second;
		int first_bitflipped[GENOTYPE_LENGTH], second_bitflipped[GENOTYPE_LENGTH];
		do {
				round++;
				// generate the complementary individuals
				bitflip(first, first_bitflipped);
				bitflip(second, second_bitflipped);

				// mark them as tested in the population
				population[hash(first)] = 1;
				population[hash(second)] = 1;

				// evaluate first and its complementary
				fitness_first = eval_fit(first, first_bitflipped);

				DBG(print_individual(first);)
				INFO(printf("fit: %f\n", fitness_first);)

				// evaluate second 
				fitness_second = eval_fit(second, second_bitflipped);

				DBG(print_individual(second);)
				INFO(printf("fit: %f\n", fitness_second);)
								
				// optimal solution?
				if ( fitness_first == 0.0 )
				{
						copy_individual(first, winner);
						copy_individual(second, loser);
						break;
				}
				if ( fitness_second == 0.0 )
				{
						copy_individual(second, winner);
						copy_individual(first, loser);
						break;
				}

				if ( fitness_first < fitness_second )
				{
						// first wins
						INFO(printf ("winner: first\n");)
						// mutate loser
						if ( mutate_1gene(second, population) < 0 )
						{
								// if no more mutations are possible
								copy_individual(first, winner);
								copy_individual(second, loser);
								break;
						}
				}
				else
				{
						// second wins
						INFO(printf ("winner: second\n");)
						// mutate loser
						if ( mutate_1gene(first, population) < 0 )
						{
								copy_individual(second, winner);
								copy_individual(first, loser);
								break;
						}
				}
				INFO(printf("\n");)

		} while ( 1 );

		return round;
}

int main ( int argc, char **argv ) {
		// check if param is given and used to seed
		if ( argc <= 1 )
		{
				printf("Usage: %s <seed>\n", argv[0]);
				return -1;
		}

		// get the seed from the user
		int seed = atoi(argv[1]);
		// initialize the random generator
		srand(seed);

		// the most-used variable around :)
		int i;

		/////

		// initialize the population
		int population[MAX_MUTATIONS];
		for ( i = 0 ; i < MAX_MUTATIONS ; i++ )
				population[i] = 0;

		// generate the first random individual
		int a[GENOTYPE_LENGTH];
		rand_init(a);
		// mark it as known in the population
		population[hash(a)] = 1;

		// generate the second random individual (different)
		int b[GENOTYPE_LENGTH];
		do { 
				rand_init(b);
		} while ( population[hash(b)] != 0 );
		population[hash(b)] = 1;

		int winner[GENOTYPE_LENGTH], loser[GENOTYPE_LENGTH];
		int winner_bitflipped[GENOTYPE_LENGTH], loser_bitflipped[GENOTYPE_LENGTH];

		/////

		// the actual GA
		int round_tot = evolve_both(a, b, winner, loser, population);

		/////
	
		// regenerate the complementary individuals 
		bitflip(winner, winner_bitflipped);
		bitflip(loser, loser_bitflipped);

		INFO(printf("after %d rounds\n", round_tot);)
		INFO(printf("best match: %d, %d\n", sum(winner), prod(winner_bitflipped));)
		INFO(printf("second-best match: %d, %d\n", sum(loser), prod(loser_bitflipped));)

		printf("%d (%d %d)\t", round_tot, sum(winner), prod(winner_bitflipped));
		print_individual(winner);
		printf("\t");
		print_individual(winner_bitflipped);
		printf("\n");


		return 0;
}
