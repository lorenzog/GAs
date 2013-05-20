/**
 * getting started with ctrnn
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <signal.h>

#include "matop.h"
// using .c for the inline functions
#include "SFMT.c"

#ifdef _WITH_DEBUG
#define DBG(x);		x
#else
#define DBG(x);
#endif

#ifdef _WITH_INFO
#define INFO(x);	x
#else
#define INFO(x);	
#endif

#ifdef _WITH_CURSES
#define CUR(x);		x
#else
#define CUR(X);
#endif

#ifdef _WITH_CURSES
#include <ncurses.h>
#endif

// some constants

//const double TIMESTEP = 0.0001;
const double MUTATION_PROB = 0.5;	// gets divided over size
const double RECOMBINATION_PROB = 0.5;
const double THRESHOLD = 0.01;
const double MAX_B = 1.0;
const double MIN_B = -1.0;
const double MAX_W = 1.0;
const double MIN_W = -1.0;
const double ZERO_W = 0.01;
const double ZERO_B = 0.01;

int scr_lastrow, scr_lastcol;

// used to get out of main loop and print data
int user_interaction = 0;


/**
 * signal handling: when ctrl-c, print data
 */
void stats(int signal) {
#ifndef _WITH_CURSES
	user_interaction = 1;
#endif
}

/// randomness

/**
 * returns a random number (0,1)
 */
inline double myrand_01() {
	// returns real [0,1]
	//return genrand_real1();
	// returns real [0,1)
	//return genrand_real2();
	// returns real (0,1)
	return genrand_real3();

	// using linux rand()
	//return (double)( rand_2 / (RAND_MAX + 1.0)))
	
	// using drand48()
	// return drand48();
}

/**
 * returns a random number [-0.01,0.01]
 */
inline double myrand_min001_001() {
	return genrand_real1()/50-0.01;
}

/// math functions for neurons

inline double sigma (double x) {
	return ( 1/(1+exp(-x)) );
}

double big_sigma(int size, double W[size][size], 
		double Y[size], double B[size], int neuron) {
	int j;
	double result = 0.0;
	for ( j = 0 ; j < size ; j++ )
		result += ( W[j][neuron] * sigma(Y[j] + B[j]) );
	//result = result + ( W[j][neuron] * sigma(Y[j] + B[j]) );

	return result;
}

/**
 * returns the distance target-desired fitness
 */
inline double eval_fit(double test, double desired) {
	return fabs(test-desired);
}


/**
 * cross and breed all genes at the same time (creeping rate)
 */
void cross_and_breed_all(int size, 
		double W_winner[size][size], 
		double B_winner[size], 
		double tau_winner,
		double W_loser[size][size], 
		double B_loser[size], 
		double tau_loser) {

	// some constants (needs to be here because of the size constraint
	double mutation_prob = MUTATION_PROB / size;

	int i, j;
	
	// "basic" microbial GA
	for ( i = 0 ; i < size ; i++ )
	{
		// the weights
		for ( j = 0 ; j < size ; j++ )
			if ( myrand_01() < mutation_prob )
				// mutate, not cross
				W_loser[i][j] += myrand_min001_001();
			else
				if ( myrand_01() < RECOMBINATION_PROB )
					// cross, not mutate
					W_loser[i][j] = W_winner[i][j];

		// the biases
		if ( myrand_01() < mutation_prob )
			// mutate, not cross
			B_loser[i] += myrand_min001_001();
		else
			if ( myrand_01() < RECOMBINATION_PROB )
				// cross, not mutate
				B_loser[i] = B_winner[i];	
	}
	
	// the learning rate
	if ( myrand_01() < mutation_prob )
		//tau_loser += myrand_01() - 0.5;
		// XXX
		tau_loser += tau_loser * (myrand_01()-0.5) * 3;
		//tau_loser += myrand_01() * 2.0;
		//tau_loser += myrand_min001_001();
	else
		if ( myrand_01() < RECOMBINATION_PROB )
			tau_loser = tau_winner;

}

short eval_and_mutate(int size, int samples,
		double Y1[size], double B1[size], double W1[size][size],
		double fitness_1[samples], double tau1,
		double Y2[size], double B2[size], double W2[size][size],
		double fitness_2[samples], double tau2
		) {

	int i, j;

	// the total fitness
	double fit_1, fit_2;

	///
	// calculate the fitness as:

	/*
	// the worst result given by a network
	fit_1 = THRESHOLD;
	fit_2 = THRESHOLD;
	for ( i = 0 ; i < samples ; i++ )
	{
		fit_1 = ( fit_1 < fabs(fitness_1[i]) ? fabs(fitness_1[i]) : fit_1 );
		fit_2 = ( fit_2 < fabs(fitness_2[i]) ? fabs(fitness_2[i]) : fit_2 );
	}
	INFO(
			printf("worst (first): %f (second) %f\n", fit_1, fit_2);
		);
	CUR(
			printw("\nworst (first): %f (second) %f\n", fit_1, fit_2);
	   );
	*/

	/// or

	// the average fitness
	fit_1 = 0;
	fit_2 = 0;
	for ( i = 0 ; i < samples ; i++ )
	{
		fit_1 += fitness_1[i];
		fit_2 += fitness_2[i];
	}
	fit_1 /= samples;
	fit_2 /= samples;
	fit_1 = fabs(fit_1);
	fit_2 = fabs(fit_2);
	INFO(
			printf("avg (first): %f (second) %f\n", fit_1, fit_2);
		);
	CUR(
			printw("\navg (first): %f (second) %f\n", fit_1, fit_2);
	   );

	///

	if ( fit_1 <= THRESHOLD )
		return 1;	// first wins
	if ( fit_2 <= THRESHOLD )
		return 2;	// second wins

	// no clear winner; squash weights and biases
	flatten_vector(MAX_B, MIN_B, size, B1);
	flatten_vector(MAX_B, MIN_B, size, B2);
	flatten_matrix(MAX_W, MIN_W, size, W1);
	flatten_matrix(MAX_W, MIN_W, size, W2);
	normalize_matrix(ZERO_W, size, W1);
	normalize_matrix(ZERO_W, size, W2);
	normalize_vector(ZERO_B, size, B1);
	normalize_vector(ZERO_B, size, B2);

	///

	// select best
	if ( fit_1 < fit_2 )
		// first is fitter
		cross_and_breed_all(size, W1, B1, tau1, W2, B2, tau2);
	else
		// second is fitter
		cross_and_breed_all(size, W2, B2, tau2, W1, B1, tau1);

	return 0;	// no winner anyway, so back to main loop
}


/**
 * CTRNN derivative of state over time + euler integration
 */
void derive_et_integra(int size, double timestep, double I[size], 
		double W[size][size], double B[size], double tau, double Y[size]) {
	int i, j;
	double d_Y[size];

	// for each neuron, calculate derivatives over time
	for ( i = 0 ; i < size ; i++ )
	{
		d_Y[i] = ( -Y[i] + big_sigma(size, 
					W, Y, B, i) 
				+ I[i] ) / tau;
		DBG(
				printf("d_y[%d] = %f\n", i, d_Y[i]);
		   );
	}

	// then update the neurons' state
	// note: has to be in a different for() loop
	// because the derivatives have to be updated at the same time
	// before the integration step. sorry, CPU.
	for  ( i = 0 ; i < size ; i++ )
	{
		// integration
		Y[i] += timestep * d_Y[i];
		DBG(
				printf("y[%d] = %f\n", i, Y[i]);
		   );

	}
}

///

int main ( int argc, char **argv ) {
	if ( argc <= 5 )
	{
		printf("Usage: %s <random seed> <number of nodes> <timestep> <number of samples> <length of simulation>\n", argv[0]);
		return -1;
	}

	CUR(
			initscr();
			getmaxyx(stdscr, scr_lastrow, scr_lastcol);
	   );

	// signal handling
	signal (SIGTERM, stats);

	// initialize the random number generators
	// linux generic rand()
	srand(atoi(argv[1]));
	// deprecated but still good for this purposes 48-bin rand()
	srand48(atoi(argv[1]));
	// a new version of the mersenne twister PRNG
	init_gen_rand(atoi(argv[1]));

	// get the size of the matrix from the user
	int size = atoi(argv[2]);
	const double timestep = atof(argv[3]);
	const int samples = atoi(argv[4]);
	//const int trials = atoi(argv[5]);
	const double simLength = atof(argv[5]);

	// the number of runs for each CTRNN
	const int trials = simLength / timestep;

	int i, j;

	///

	// the number of samples on which to test each CTRNN
	double velocities[samples];

	///

	// the final weight matrix
	double W[size][size];
	// zero because it will be overwritten
	zero_matrix(size, W);
	// the final state vector
	double Y[size];
	zero_vector(size, Y);
	// the final input vector
	double I[size];
	zero_vector(size, I);
	// the final biases
	double B[size];
	zero_vector(size, B);
	// the final learning rate
	double tau;

	/// first ANN

	// the weight matrix
	double W1[size][size];
	random_matrix(size, W1);
	//zero_matrix(size, W1);
	// the state vector
	double Y1[size];
	zero_vector(size, Y1);
	// the input vector
	double I1[size];
	zero_vector(size, I1);
	// the biases
	double B1[size];
	random_vector(size, B1);
	// the learning rate
	// XXX
	double tau1 = myrand_01() * 20;

	// the fitness for each sample
	double fitness_1[samples];
	zero_vector(samples, fitness_1);

	/// second ANN

	// the weight matrix
	double W2[size][size];
	random_matrix(size, W2);
	//zero_matrix(size, W2);
	// the state vector
	double Y2[size];
	zero_vector(size, Y2);
	// the input vector
	double I2[size];
	zero_vector(size, I2);
	// the biases
	double B2[size];
	random_vector(size, B2);
	// XXX
	double tau2 = myrand_01() * 20;

	// the fitness for each sample
	double fitness_2[samples];
	zero_vector(samples, fitness_1);

	///

	// build the array of random velocities
	// over which to test each CTRNN
	// velocities are in the interval (0,10)
	//random_vector(samples, velocities);
	for ( i = 0 ; i < samples ; i++ )
		//velocities[i] = myrand_01() * 100;
		velocities[i] = myrand_01();
	// sort them to have a better display (?)
	// XXX
	//sort_vector(samples, velocities);

	DBG(
		print_vector(samples, velocities);
	   );

	// a counter for the total number of tournaments
	int tournaments = 0;
	// a counter for the number of samples
	int num_sample;

	// the distance fed into the input neuron
	double distance;

	// the velocity at which each network is trained
	double vel;

	// counter for the number of trials each network is trained
	int count;

	// counter for debug printing
	int modcount;

	// temporary fitnesses for evaluation
	double fit1, fit2;

	do {
		CUR (
			// to adjust the last line (tournament counter)
			// when the user changes size of the terminal
			getmaxyx(stdscr, scr_lastrow, scr_lastcol);
			);

		// sample over many different velocities
		for ( num_sample = 0 ; num_sample < samples ; num_sample++ )
		{
			// the starting distance (0,1)
			//distance = myrand_01() * 100;
			distance = myrand_01();
			// the velocity we're evaluating
			vel = velocities[num_sample];

			// resetting counters
			count = 1;
			modcount = 1;
			// resetting the input vectors
			zero_vector(size, Y1);
			zero_vector(size, Y2);

			// do <trials> number of attempt for each network
			do {
				// input is variying over time
				if ( count <= simLength / (2 * timestep) )
					distance = distance - vel * timestep;
				else
					distance = 0;
				distance = distance < 0 ? 0 : distance;

				DBG(printf("after; count/trials: %d/%d, seen distance: %f\n", count, trials, distance););


				/*
				//distance += vel*TIMESTEP;
				// target is seen only until a certain point
				// and distance is scaled 0<x<1
				   if ( count < trials*4/5 )
					   distance += vel * timestep;
				   else
					   distance = 0;
			    */

				// no negative distances, please
				//distance = ( distance < 0 ? 0 : distance );
				//printf("seen distance: %f\n", distance);

				// feeding the input into the neurons
				Y1[0] = distance;
				Y2[0] = distance;

				derive_et_integra(size, timestep, 
						I1, W1, B1, tau1, Y1);
				derive_et_integra(size, timestep, 
						I2, W2, B2, tau2, Y2);

				count++;
			} while ( count < trials );

			fit1 = eval_fit(Y1[size-1], vel);
			fit2 = eval_fit(Y2[size-1], vel);
			// eval and save the fitness of each network
			fitness_1[num_sample] = fit1;
			fitness_2[num_sample] = fit2;
			INFO(
					printf("expected: %f first network out: %f second: %f", vel, Y1[size-1], Y2[size-1]);
					printf("\t\tfirst delta: %f second: %f\n",
						fit1, fit2);
					);
			CUR(
					mvprintw(num_sample,0,
						"expected: %f first network out: %f second: %f", 
						vel, Y1[size-1], Y2[size-1]);
					printw("\t\tfirst delta: %f second: %f\n",
						fit1, fit2);
					refresh();
			   );
		}


		// result is the best value given by the winning ANN
		short evaluation = eval_and_mutate( size, samples,
				Y1, B1, W1, fitness_1, tau1,
				Y2, B2, W2, fitness_2, tau2);
		if ( evaluation == 1 )
		{
			// first wins
			INFO (
					printf("first network wins\n");
				 );
			CUR(
					mvprintw(6,0,"first network wins\n");
			   );
			copy_matrix(size, W1, W);
			copy_vector(size, B1, B);
			copy_vector(size, Y1, Y);
			break;
		}
		if ( evaluation == 2 )
		{
			// second wins
			INFO (
					printf("second network wins\n");
				 );
			CUR(
					mvprintw(6,0,"second network wins\n");
			   );
			copy_matrix(size, W2, W);
			copy_vector(size, B2, B);
			copy_vector(size, Y2, Y);
			break;
		}

		if ( user_interaction )
		{
			printf("\nUser interrupt; tournaments done: %d", tournaments);
			printf("\n");
			for ( i = 0 ; i < samples ; i++ )
			{
				printf("expected: %f\t", velocities[i]);
				printf("\t\tfitness 1: %f fitness 2: %f\n",
						fitness_1[i], fitness_2[i]);
			}
			// reset switch
			user_interaction = 0;
		}

		// keep track of the number of tournaments done
		tournaments++;
		CUR(
				mvprintw(scr_lastrow-1,0,"tournaments: %d\n", tournaments);
		   );
	} while ( 1 );

	CUR(
			printw("\n\npress RETURN...\n");
			getch();
			endwin();
	   );

	printf("\nafter %d tournaments;\n", tournaments);
	printf("\n");
	printf("evolved CTRNN output: %f\n", Y[size-1]);
	printf("weights:\n");
	print_matrix(size, W);
	printf("biases\n");
	print_vector(size, B);

	printf("\n");

	// TODO testing with a new speed?

	return 0;
}

