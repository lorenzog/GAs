#ifndef _DAISIES_H
#define _DAISIES_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <strings.h>	// bzero
#include <getopt.h>
#include <float.h>		// DBL_EPSILON
#include <time.h>		// time()

// using .c for the inline functions
#include "SFMT.c"
#include "fcmp.c"

/////
// debugging

// all debugging code has been stripped from the source;
// to debug, use
// DBG(
// 		printf("debug variable here\n");
// 	);
#ifdef _WITH_DEBUG
#define DBG(x);		x
#else
#define DBG(x);
#endif

// all printouts have been put into a similar construct
// future development: use ncurses or gui or network
// instead of printf
#ifdef _WITH_INFO
#define INFO(x);		x
#else
#define INFO(x);
#endif

/////

// note-> the order is important {0,1}: 0 is false in C
typedef enum { false, true } bool;

/////
// the default parameters

// the number of simulations
#define NUM_RUNS_DEFAULT		1
// the mutation rate
#define MUT_RATE_DEFAULT		0.1
// the timesteps
#define TIMESTEPS_DEFAULT		100000
// the maximum and minimum perturbation values
#define PERT_MIN_DEFAULT		0
#define PERT_MAX_DEFAULT		100
// the resource values
// at which individuals best adapt
#define A_START_DEFAULT			15
#define A_END_DEFAULT			85
// the standard deviation for the allele mutation
#define MUT_STDEV_DEFAULT		0.05
// intervals for the averaging of the rate of change
// of the resource
#define INTERVALS_DEFAULT		200

// version-dependent parameters

// ECAL version

// the population size
#define POP_SIZE_DEFAULT_ECAL		2000
// the death rate
#define DEATH_RATE_DEFAULT_ECAL		0.01
// the span of the parabola
// (width of the fitness function)
#define LAMBDA_DEFAULT_ECAL			0.04
// population effect and resource back-effect
#define ALPHA_DEFAULT_ECAL			0.0025
#define BETA_DEFAULT_ECAL			0.1
// max and min resource values
#define RESOURCE_START_DEFAULT_ECAL	-50
#define RESOURCE_END_DEFAULT_ECAL	150

// Journal of theoretical biology version

// the population size
#define POP_SIZE_DEFAULT_BIOL		1000
// the death rate
#define DEATH_RATE_DEFAULT_BIOL		0.005
// the span of the parabola
// (width of the fitness function)
#define LAMBDA_DEFAULT_BIOL			5
// population effect and resource back-effect
#define ALPHA_DEFAULT_BIOL			0.0005
#define BETA_DEFAULT_BIOL			0.01
// max and min resource values
#define RESOURCE_START_DEFAULT_BIOL	-50
#define RESOURCE_END_DEFAULT_BIOL	150


#define VERSION_ECAL			1
#define VERSION_BIOL			2

#define VERSION_DEFAULT			1


/////

// used for command-line parameter parsing

static int NUM_RUNS;
static int POP_SIZE;
static float MUT_RATE;
static float DEATH_RATE;
static float LAMBDA;
static float ALPHA;
static float BETA;
static int TIMESTEPS;
static float RESOURCE_START;
static float RESOURCE_END;
static float PERT_MIN;
static float PERT_MAX;
static float A_START;
static float A_END;
static float MUT_STDEV;
static int INTERVALS;
static int RAND_SEED;

static int VERSION = VERSION_DEFAULT;


// the individuals
typedef struct individual_t {
	// the genotype
	// effect on the resource
	double theta;
	// resource at which the indiv is best adapted
	double A;
} individual_t;


// parameters used for the fitness function
// in the ECAL paper
double over_sqrt_lambda;
// in the journal of th. biol. paper
double over_lambda2;

// the rate of change of the perturbation (linear)
double PERT_ROC;

// the rate of change of the resource
double dRdT;

// the number of unique A types for a single run
int unique_A;
// the average number of A alleles
double avg_num_A;

/////
// threading

// the simulation thread
pthread_t simulation_thread;

// to use other thread, use semaphoring for reading
// from global variables


#endif
