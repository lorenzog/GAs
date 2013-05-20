/**
 */

#include "daisies.h"


///////////////////////// BIOL journal version ///////////////////////


/**
 * initializes the population with random individuals
 */
void init_pop2(individual_t *population) {
	int i;
	DBG(
			double initial_theta = 0;
			double initial_A = 0;
	   );
	
	double r;
	for ( i = 0 ; i < POP_SIZE ; i++ )
	{
		/////
		// the theta allele
		// for theta either -1 or 1
		// random
		//population[i].theta = ( genrand_real1() >= 0.5 ? 1.0 : -1.0 );
		// half and half
		//population[i].theta = ( i < POP_SIZE/2 ? -1 : 1 );

		// for theta [0,1]
		population[i].theta = genrand_real1();

		DBG(initial_theta += population[i].theta);

		//////
		// the A allele
		// use random values between A_START and A_END
		population[i].A = genrand_real2() * (A_END-A_START+1) + A_START;

		DBG(initial_A += population[i].A;);
	}
	DBG(printf("initial theta: %f, initial A: %f\n", initial_theta*ALPHA, initial_A/POP_SIZE););
}

/**
 * evaluates the resource at a certain timestep
 */
void evaluate_resource2(int timestep, double *resource, individual_t *population) {
	int i;
	DBG(int index;);

	// used to calculate how many unique A values are present
	int num_uniques = (int)A_END - A_START + 1;
	int uniques[num_uniques];
	bzero((void*)uniques, num_uniques*sizeof(int));

	// zero the number of unique A types
	unique_A = 0;

	double global_theta = 0;
	double global_A = 0;
	// for each individual
	// calculate how they affect the resource
	for ( i = 0 ; i < POP_SIZE ; i++ )
	{
		if ( population[i].theta < 0.5 )
			global_theta--;
		if ( population[i].theta > 0.5 )
			global_theta++;
		
		global_A += population[i].A;

		DBG(
		//index = (int)rint(population[i].A) - (int)A_START;
		index = (int)floor(population[i].A) - (int)A_START;
		printf("pop[%d].A = %f\n", i, population[i].A);
		printf("uniques[%d]: %d\n", index, uniques[index]);
		);

		// also round and count how many uniques A values are there
		//if ( (uniques[(int)rint(population[i].A) - (int)A_START]++) == 0 )
		if ( (uniques[(int)floor(population[i].A) - (int)A_START]++) == 0 )
			unique_A++;

		DBG( printf("tot uniques: %d\n", unique_A););
	}

	// calculate the perturbation of the resource
	// (grows linearly)
	double P = PERT_MIN + PERT_ROC * timestep;

	INFO(
			printf("%f\t", *resource);
			printf("%f\t", P);
			printf("%f\t", global_A/(POP_SIZE*1.0));
			printf("%f\t", global_theta*ALPHA);
			printf("%f\t", global_theta);
			printf("%d\t", unique_A);
			printf("\n");
		);

	// updates the resource
	// note: because timesteps are discrete,
	// there's no need to use numerical integration
	// (which would be res += timestep * new_res)
	dRdT = ( ALPHA * global_theta + BETA * ( P - *resource ) ); 
	*resource += dRdT;
	DBG(
			printf("time: %d, global_theta: %f, dRdT: %f, P: %f, new res: %f\n", timestep, global_theta, dRdT, P, *resource);
	   );
}

/**
 * a second version of the mutation routine
 */
void mutate2 (struct individual_t *individual) {
	double r;

	// mutate the A allele
	//if ( genrand_real1() < MUT_RATE )
	if ( fcmp(genrand_real1(), MUT_RATE, DBL_EPSILON) < 0 )
	{
		// a random number between [0,1)
		r = genrand_real2();
		r *= MUT_STDEV;
		individual->A += ( gen_rand32() % 2 ? r : -r );

		// check boundaries
		// use of reflection!
		//if ( individual->A > A_END )
		if ( fcmp(individual->A, A_END, DBL_EPSILON) > 0 )
			individual->A = A_END - (individual->A - A_END);
		//if ( individual->A < A_START )
		if ( fcmp(individual->A, A_START, DBL_EPSILON) < 0 )
			individual->A = A_START + (A_START - individual->A);

		r = genrand_real2();
		r *= MUT_STDEV;

		// reflection
		individual->theta += ( gen_rand32() % 2 ? r : -r );
		//if ( individual->theta > 1.0 )
		if ( fcmp(individual->theta, 1.0, DBL_EPSILON) > 0 )
			individual->theta = 1.0 - ( individual->theta - 1.0 );
		//if ( individual->theta < 0 )
		if ( fcmp(individual->theta, 0, DBL_EPSILON) < 0 )
			individual->theta = 0 + ( - individual->theta );

		DBG(printf("(new theta): %f\n", individual->theta););
	}
}

/**
 * using "environmental regulation..." paper
 */
int evaluate_fitness2(struct individual_t *first, 
		struct individual_t *second, double* resource) {

	double fitness_first, fitness_second;

	fitness_first = fmax((1 - over_lambda2 * pow((first->A - *resource), 2)), 0.0);
	fitness_second = fmax((1 - over_lambda2 * pow((second->A - *resource), 2)), 0.0);
	//fprintf(stderr, "first: %f, second: %f\n", fitness_first, fitness_second);

	// if both have zero fitness, return to calling function
	// which will create a random genotype 
	// for the replacement
	//if ( fitness_first == 0 && fitness_second == 0 )
	if ( fcmp(fitness_first, 0, DBL_EPSILON) == 0
			&& fcmp(fitness_second, 0, DBL_EPSILON) == 0)
		return 0;
	
	//if ( fitness_first > fitness_second )
	if ( fcmp(fitness_first, fitness_second, DBL_EPSILON) > 0 )
		return 1;
	//if ( fitness_first < fitness_second )
	if ( fcmp(fitness_first, fitness_second, DBL_EPSILON) < 0 )
		return -1;

}

/**
 * using the "environmental regulation..." paper
 */
void execute_ga2(individual_t *population, 
		int *possible_choices, double *resource) {
	// pick n individuals
	// note: individuals must be different (?)
	int i;
	int candidate;
	// number of chosen individuals
	// (death list!!)
	int tourn_pool = POP_SIZE*DEATH_RATE;
	int choices[tourn_pool];
	for ( i = 0 ; i < tourn_pool; i++ )
	{
		// pick a number between 0 and pop_size-1
		do
			candidate = (int) floor(genrand_real2() * POP_SIZE);
		// unless it has been chosen already
		while ( possible_choices[candidate] == -1 );
		// once selected, remove it from list
		possible_choices[candidate] = -1;
		// this is our i-th individual
		choices[i] = candidate;
	}

	int winner; //, loser;
	int first, second;
	//int counter = 0;
	int result;
	// for each couple of individuals
	for ( i = 0 ; i < tourn_pool/2 ; i++ )
	{
		// new version: tournament, the winner gets copied over
		// another random individual

		// get two random individuals
		// different from the chosen one
		do
			first = (int) floor(genrand_real2() * POP_SIZE);
		while ( first == choices[i] );

		do 
			second = (int) floor(genrand_real2() * POP_SIZE);
		while ( second == choices[i] || second == first );

		// using new version
		result = evaluate_fitness2(&population[first], 
				&population[second], resource);

		/*
		 * if no clear winner, no replacement!
		 * -> new paper: generate random genotype
		 * and replace individual
		 */
		if ( result == 0 )
		{
			population[choices[i]].theta = genrand_real2();
			population[choices[i]].A = genrand_real2() * (A_END-A_START+1) + A_START;
			continue;
		}

		// the winner gets copied
		if ( result > 0 )
		{
			// first is winner
			winner = first;
			//loser = second;
		}
		else
			{
				// second is winner
				winner = second;
				//loser = first;
			}

		// (new version)
		// winner takes place of the dead one
		population[choices[i]] = population[winner];
		
		// eventually mutate
		mutate2(&population[choices[i]]);
	}

	// now reset the choices
	for ( i = 0 ; i < tourn_pool; i++ )
		possible_choices[choices[i]] = 0;
}




////////////////////////// ECAL version ////////////////////////////



/**
 * initializes the population with random individuals
 */
void init_pop1(individual_t *population) {
	int i;
	DBG(
			double initial_theta = 0;
			double initial_A = 0;
	   );
	
	double r;
	for ( i = 0 ; i < POP_SIZE ; i++ )
	{
		/////
		// the theta allele
		// for theta either -1 or 1
		// random
		//population[i].theta = ( genrand_real1() >= 0.5 ? 1.0 : -1.0 );
		// half and half
		//population[i].theta = ( i < POP_SIZE/2 ? -1 : 1 );

		// for theta [-1, 1]
		//population[i].theta = genrand_real2() * 2.0 - 1.0;
		population[i].theta = genrand_real1() * 2.0 - 1.0;


		DBG(initial_theta += population[i].theta);

		//////
		// the A allele
		// use random values between A_START and A_END?
		//population[i].A = genrand_real1() * (A_END-A_START) + A_START;
		// XXX using [0,1), see MT faqs
		population[i].A = genrand_real2() * (A_END-A_START+1) + A_START;

		DBG(initial_A += population[i].A;);
	}
	DBG(printf("initial theta: %f, initial A: %f\n", initial_theta*ALPHA, initial_A/POP_SIZE););
}

/**
 * evaluates the resource at a certain timestep
 */
void evaluate_resource1(int timestep, double *resource, individual_t *population) {
	int i;
	DBG(int index;);

	// used to calculate how many unique A values are present
	int num_uniques = (int)A_END - A_START + 1;
	int uniques[num_uniques];
	bzero((void*)uniques, num_uniques*sizeof(int));

	// zero the number of unique A types
	unique_A = 0;

	double global_theta = 0.0;
	double global_A = 0.0;
	// for each individual
	// calculate how they affect the resource
	for ( i = 0 ; i < POP_SIZE ; i++ )
	{
		global_theta += population[i].theta;
		global_A += population[i].A;

		/*
		DBG(
		index = (int)floor(population[i].A) - (int)A_START;
		printf("pop[%d].A = %f\n", i, population[i].A);
		printf("uniques[%d]: %d\n", index, uniques[index]);
		);
		*/

		// also round and count how many uniques A values are there
		// use floor!! otherwise it will overflow the array
		//if ( (uniques[(int)rint(population[i].A) - (int)A_START]++) == 0 )
		if ( (uniques[(int)floor(population[i].A) - (int)A_START]++) == 0 )
			unique_A++;

		//DBG( printf("tot uniques: %d\n", unique_A););
	}

	// calculate the perturbation of the resource
	// (grows linearly)
	double P = PERT_MIN + PERT_ROC * timestep;

	INFO(
			printf("%f\t", *resource);
			printf("%f\t", P);
			printf("%f\t", global_A/POP_SIZE*1.0);
			printf("%f\t", global_theta*ALPHA);
			printf("%f\t", global_theta);
			printf("%d\t", unique_A);
			printf("\n");
		);

	// updates the resource
	// note: because timesteps are discrete,
	// there's no need to use numerical integration
	// (which would be res += timestep * new_res)
	//dRdT = ( ALPHA * global_theta + BETA * ( P - *resource ) ); 
	dRdT = ( ALPHA * global_theta + BETA * ( P - *resource ) ); 
	*resource += dRdT;
	DBG(
			printf("time: %d, global_theta: %f, dRdT: %f, P: %f, new res: %f\n", timestep, global_theta, dRdT, P, *resource);
	   );
}

/**
 * the mutation routine
 */
void mutate1 (struct individual_t *individual) {
	double r;

	// mutate the A allele
	//if ( genrand_real1() < MUT_RATE )
	if ( fcmp(genrand_real1(), MUT_RATE, DBL_EPSILON) < 0  )
	{
		// do mutate each gene (or just A?)
		// a random number between [0,1]
		//r = genrand_real1();
		// XXX using [0,1)
		r = genrand_real2();
		/*
		 * XXX
		// scaled down to (80-15), 0.05
		r *= 0.05 * (A_END-A_START);
		// add or subtract the value
		DBG(printf("mutation (old A): %f\t", individual->A););
		individual->A += ( gen_rand32() % 2 ? rint(r) : -rint(r) );
		*/
		r *= MUT_STDEV;
		individual->A += ( gen_rand32() % 2 ? r : -r );

		// check boundaries
		//if ( individual->A > A_END )
		if ( fcmp(individual->A, A_END, DBL_EPSILON) > 0 )
			individual->A = A_END;
			// with reflection
			//individual->A = A_END - (individual->A - A_END);
		//if ( individual->A < A_START )
		if ( fcmp(individual->A, A_START, DBL_EPSILON) < 0 )
			individual->A = A_START;
			// with reflection
			//individual->A = A_START + (A_START - individual->A);
		DBG(printf("(new A): %f\n", individual->A););


		/*
	}

	// mutate the theta allele
	if ( genrand_real1() < MUT_RATE )
	{
		// mutate theta
		DBG(printf("mutation (old theta): %f\t", individual->theta););
		*/
		//r = genrand_real1();
		// XXX using [0,1)
		r = genrand_real2();
		r *= MUT_STDEV;
		individual->theta += ( gen_rand32() % 2 ? r : -r );
		//if ( individual->theta > 1.0 )
		if ( fcmp(individual->theta, 1.0, DBL_EPSILON) > 0 )
			individual->theta = 1.0;
			// reflection
			//individual->theta = 1.0 - ( individual->theta - 1.0 );
		//if ( individual->theta < -1.0 )
		if ( fcmp(individual->theta, -1.0, DBL_EPSILON) < 0 )
			individual->theta = -1.0;
			// reflection
			//individual->theta = -1.0 + ( -1.0 - individual->theta );

		DBG(printf("(new theta): %f\n", individual->theta););
	}
}

/**
 * evaluates fitness of individuals
 */
int evaluate_fitness1(struct individual_t *first, 
		struct individual_t *second, double* resource) {

	// A is the only individual-dependent parameter
	/*
	if ( first->A == second->A )
		return 0;
		*/

	double fitness_first, fitness_second;
	//if ( fabs(first->A - *resource) < over_sqrt_lambda )
	if ( fcmp(fabs(first->A - *resource), over_sqrt_lambda, DBL_EPSILON) < 0 )
		fitness_first = 1 - LAMBDA * pow((first->A - *resource), 2);
	else
		fitness_first = 0;

	//if ( fabs(second->A - *resource) < over_sqrt_lambda )
	if ( fcmp(fabs(second->A - *resource), over_sqrt_lambda, DBL_EPSILON) < 0 )
		fitness_second = 1 - LAMBDA * pow((second->A - *resource), 2);
	else
		fitness_second = 0;
	DBG(printf("A1: %f, fit1: %f, A2: %f, fit2: %f\n", first->A, fitness_first, second->A, fitness_second););

	return fcmp(fitness_first, fitness_second, DBL_EPSILON);
	
	/*
	//if ( fitness_first == fitness_second )
	if ( fcmp(fitness_first, fitness_second, DBL_EPSILON) == 0 )
	{
		DBG(printf("no win\n"););
		return 0;
	}

	//if ( fitness_first > fitness_second )
	if ( fcmp(fitness_first, fitness_second, DBL_EPSILON) > 0 )
	{
		DBG(printf("first wins\n"););
		return 1;
	}
	else
	{
		DBG(printf("second wins\n"););
		return -1;
	}
	*/
}

/**
 * genetic algorithm
 */
void execute_ga1(individual_t *population, 
		int *possible_choices, double *resource) {
	// pick n individuals
	// note: individuals must be different (?)
	int i;
	int candidate;
	// number of chosen individuals
	// (death list!!)
	int tourn_pool = POP_SIZE*DEATH_RATE;
	int choices[tourn_pool];
	DBG(printf("Dead list: "););
	for ( i = 0 ; i < tourn_pool; i++ )
	{
		// pick a number between 0 and pop_size-1
		do
		{
			// using suggested method from SFMT FAQ's
			candidate = (int) floor(genrand_real2() * POP_SIZE);
		}
		// unless it has been chosen already
		while ( possible_choices[candidate] == -1 );
		// once selected, remove it from list
		possible_choices[candidate] = -1;
		// this is our i-th individual
		choices[i] = candidate;
		DBG(printf("%d ", candidate););
	}
	DBG(printf("\n"););

	int winner;//, loser;
	int first, second;
	//int counter = 0;
	int result;
	// for each couple of individuals
	//for ( i = 0 ; i < tourn_pool/2 ; i++ )
	for ( i = 0 ; i < tourn_pool ; i++ )
	{
		/*
		 * XXX
		 * old version: tournament, the loser disappears
		first = counter++;
		second = counter++;
		// evaluate fitness
		result = evaluate_fitness(&population[choices[first]], &population[choices[second]], resource);
		*/
		
		// new version: tournament, the winner gets copied over
		// another random individual

		// get two random individuals
		do
			//first = gen_rand32() % POP_SIZE;
			first = (int) floor(genrand_real2() * POP_SIZE);
		while ( first == choices[i] );
		do 
			//second = gen_rand32() % POP_SIZE;
			second = (int) floor(genrand_real2() * POP_SIZE);
		while ( second == choices[i] || second == first );
		DBG(printf("first: %d, second: %d\t", first, second););

		result = evaluate_fitness1(&population[first], 
				&population[second], resource);

		/*
		 * if no clear winner, no replacement!
		 * (is there a mutation?)
		 */
		if ( result == 0 )
		{
			//mutate1(&population[choices[i]]);
			DBG(printf("no winner\n"););
			continue;
		}

		// the winner gets copied
		if ( result > 0 )
		{
			// first is winner
			winner = first;
			//loser = second;
		}
		else
			{
				// second is winner
				winner = second;
				//loser = first;
			}
		// if no winner, keep going
		// winner takes the place of the loser
		// XXX (old version)
		//population[loser] = population[winner];

		// (new version)
		// winner takes place of the dead one
		population[choices[i]] = population[winner];
		//population[choices[i]].A = population[winner].A;
		//population[choices[i]].theta = population[winner].theta;
		
		// eventually mutate
		//mutate(&population[loser]);
		mutate1(&population[choices[i]]);
	}

	// now reset the choices
	for ( i = 0 ; i < tourn_pool; i++ )
		possible_choices[choices[i]] = 0;
}

////////////////////////// common routines /////////////////////////

/**
 * the simulation thread
 */
void *simulation() {
	// the population
	individual_t population[POP_SIZE];

	double mean_A = 0;
	double mean_time_regulated = 0;

	// for every requested run
	int i;
	for ( i = 0 ; i < NUM_RUNS; i++ )
	{
		INFO(
				fprintf(stderr, "run %i/%d\n", i+1, NUM_RUNS);
			);
		// initializes the population
		if ( VERSION == VERSION_ECAL )
		{
			DBG(printf("running ECAL version\n"););
			init_pop1(population);
		}
		else 
			if ( VERSION == VERSION_BIOL )
			{
				DBG(printf("running Journal of Theoretical Biology version\n"););
				init_pop2(population);
			}
			else
				pthread_exit(NULL);

		// initializes the resource
		double resource = RESOURCE_START * 1.0;

		// initializes the perturbation
		// note: keep the * 1.0 otherwise it will be cast to int
		PERT_ROC = (PERT_MAX * 1.0) / TIMESTEPS;

		// reset the average number of A alleles
		avg_num_A = 0;
		// reset the time counter
		int timestep = 0;

		// the array of possible choices 
		// used for the random selection of individuals
		// in the genetic algorithm
		// (runs in o(number of deaths)
		int possible_choices[POP_SIZE];
		// zero the array of possible choices
		bzero((void*)possible_choices, POP_SIZE*sizeof(int));

		double dRdT_averaged = 0;
		int time_regulated = 0;

		// for every instant in time
		while ( timestep < TIMESTEPS )
		{
			DBG(printf("\n"););

			// evaluate the resource
			if ( VERSION == VERSION_ECAL )
				evaluate_resource1(timestep, &resource, population);
			else
				evaluate_resource2(timestep, &resource, population);
				

			// save the rate of INCREASE of the resource
			// to be averaged
			// (using absolute values?!)
			dRdT_averaged += dRdT;

			// if the interval has elapsed
			if ( timestep % INTERVALS == 0 )
			{
				// compute the average of the growth rate
				dRdT_averaged /= INTERVALS;
				DBG( printf("avg: %f, PERT_ROC: %f\n", dRdT_averaged, PERT_ROC););
				// if the average is lesser than 
				// the rate of change of the perturbation
				//if ( fabs(dRdT_averaged) < PERT_ROC )
				if ( fcmp(fabs(dRdT_averaged), PERT_ROC, DBL_EPSILON) < 0 )
				{
					// we assume that the resource was regulated
					// for the latest interval
					time_regulated += INTERVALS;
				}
			}

			// add up the number of unique A alleles
			avg_num_A += unique_A;

			// execute the genetic algorithm
			if ( VERSION == VERSION_ECAL )
				execute_ga1(population, possible_choices, &resource);
			else
				execute_ga2(population, possible_choices, &resource);

			// time passes for everybody
			timestep++;
		}

		// average the number of A alleles over the time interval
		avg_num_A/=TIMESTEPS;
		// add up the number of average alleles
		mean_A+=avg_num_A;
		mean_time_regulated+=time_regulated;

		INFO(
				// print on standard error for easiness of collecting data
				// the time the resource was regulated
				fprintf(stderr, "%d\t", time_regulated);
				fprintf(stderr, "%f\n", avg_num_A);
			);
	}

	// average the mean
	mean_A/=NUM_RUNS;
	mean_time_regulated/=NUM_RUNS;
	INFO(
			fprintf(stderr, "Mean number of unique A traits: %f\n", mean_A);
			fprintf(stderr, "Mean time regulated: %f\n", mean_time_regulated);
		);


	pthread_exit(NULL);
}


////////////////////////// GUI and user interaction ////////////////


inline void print_usage(char *prog_name) {
	printf("Usage: %s [--OPTIONS (value)]\n", prog_name);
	printf("options\t\t(type)\t\tdefault value:\n\n");
	printf("--num_runs\t(int)\t1\n");
	printf("--pop_size\t(int)\t\t2000\n--mut_rate\t(double)\t0.1\n");
	printf("--death_rate\t(double)\t0.01\n--lambda\t(double)\t0.04\n");
	printf("--alpha\t\t(double)\t0.0024\n--beta\t\t(double)\t0.1\n");
	printf("--timesteps\t(int)\t\t100000\n");
	printf("--res_start\t(double)\t-50\n--res_end\t(double)\t150\n");
	printf("--pert_min\t(double)\t0\n--pert_max\t(double)\t100\n");
	printf("--A_start\t(double)\t15\n--A_end  \t(double)\t85\n");
	printf("--mut_stdev\t(double)\t0.05\n--random\t(int)\t\t42\n");
	printf("--intervals\t(int)\t200\n");
	printf("\n--version\t(int)\t1\n");
}

/**
 * command-line option parsing
 */
bool parse_options(int argc, char** argv) {

	///
	///taken from: getopt(3)

	int c;
	int digit_optind = 0;
	// do not show error messages
	opterr = 0;		

	while (1) {
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		/*
		 * struct option {
		 * const char *name;	// name of the long option
		 * int has_arg;			// 0 (no arg) 
		 * 						// 1 (mandatory arg) 
		 *						// 2 (optional arg)
		 * int *flag;			// NULL or a variable set to val 
		 * 						// if opt is found
		 * int val;				// value to return or to load 
		 * 						// into variable flag
		 * };
		 */
		static struct option long_options[] = {
			{"num_runs", 1, 0, 'n'},
			{"pop_size", 1, 0, 's'},
			{"mut_rate", 1, 0, 'm'},
			{"death_rate", 1, 0, 'd'},
			{"lambda", 1, 0, 'l'},
			{"alpha", 1, 0, 'a'},
			{"beta", 1, 0, 'b'},
			{"timesteps", 1, 0, 't'},
			{"res_start", 1, 0, 'r'},
			{"res_end", 1, 0, 'R'},
			{"pert_min", 1, 0, 'p'},
			{"pert_max", 1, 0, 'P'},
			{"A_start", 1, 0, '1'},
			{"A_end", 1, 0, '2'},
			{"mut_stdev", 1, 0, 'M'},
			{"random", 1, 0, 'x'},
			{"intervals", 1, 0, 'i'},
			{"version", 1, 0, 'v'},
		};

		// IMPORTANT keep the '-' at the beginning of the optstring
		// to avoid a permutation of the command-line arguments
		c = getopt_long(argc, argv, "-n:s:m:d:l:a:b:t:r:R:p:P:1:2:M:x:i:v:",
				long_options, &option_index);
		if (c == -1)
			// no more options, move on
			return true;

		//use case '1' (one) to handle non-recognized options,
		switch (c) {
			case'n':
				DBG(printf("num runs: %d\n", atoi(optarg)););
				NUM_RUNS = atoi(optarg);
				break;
			case 's':
				DBG(printf("population size: %d\n", atoi(optarg)););
				POP_SIZE = atoi(optarg);
				break;
			case 'm':
				DBG(printf("mutation rate: %f", atof(optarg)););
				MUT_RATE = atof(optarg);
				break;
			case 'd':
				DBG(printf("death rate:"););
				DEATH_RATE = atof(optarg);
				break;
			case 'l':
				DBG(printf("lambda:"););
				LAMBDA = atof(optarg);
				break;
			case 'a':
				DBG(printf("alpha:"););
				ALPHA = atof(optarg);
				break;
			case 'b':
				DBG(printf("beta:"););
				BETA = atof(optarg);
				break;
			case 't':
				DBG(printf("timesteps:"););
				TIMESTEPS = atoi(optarg);
				break;
			case 'r':
				DBG(printf("resource starting value:"););
				RESOURCE_START = atof(optarg);
				break;
			case 'R':
				DBG(printf("resource ending value:"););
				RESOURCE_END = atof(optarg);
				break;
			case 'p':
				DBG(printf("perturbation starting value:"););
				PERT_MIN = atof(optarg);
				break;
			case 'P':
				DBG(printf("perturbation ending value:"););
				PERT_MAX = atof(optarg);
				break;
			case '1':
				DBG(printf("individual minimum resource level:"););
				A_START = atof(optarg);
				break;
			case '2':
				DBG(printf("individual maximum resource level:"););
				A_END = atof(optarg);
				break;
			case 'M':
				DBG(printf("mutation rate standard deviation:"););
				MUT_STDEV = atof(optarg);
				break;
			case 'x':
				DBG(printf("random seed: %d\n", atoi(optarg)););
				RAND_SEED = atoi(optarg);
				break;
			case 'i':
				DBG(printf("intervals: %d\n", atoi(optarg)););
				INTERVALS = atoi(optarg);
				break;
			case 'v':
				DBG(printf("version: %d\n", atoi(optarg)););
				VERSION = atoi(optarg);
				break;
			default:
				printf("unrecognized option\n");
				print_usage(argv[0]);
				return false;
		}
	}
	return true;
}

/**
 * set up the default parameters
 * for the ECAL version
 * and the Journal of Theoretical Biology version
 */
void set_default_params(int version ) {
	// the default parameter
	NUM_RUNS = NUM_RUNS_DEFAULT;
	MUT_RATE = MUT_RATE_DEFAULT;
	TIMESTEPS = TIMESTEPS_DEFAULT;
	PERT_MIN = PERT_MIN_DEFAULT;
	PERT_MAX = PERT_MAX_DEFAULT;
	A_START = A_START_DEFAULT;
	A_END = A_END_DEFAULT;
	MUT_STDEV = MUT_STDEV_DEFAULT;
	INTERVALS = INTERVALS_DEFAULT;
	// random seed is the number of seconds elapsed since the Epoch
	RAND_SEED = time(NULL);

	// the different default parameters
	if ( version == VERSION_ECAL )
	{
		POP_SIZE = POP_SIZE_DEFAULT_ECAL;
		DEATH_RATE = DEATH_RATE_DEFAULT_ECAL;
		LAMBDA = LAMBDA_DEFAULT_ECAL;
		ALPHA = ALPHA_DEFAULT_ECAL;
		BETA = BETA_DEFAULT_ECAL;
		RESOURCE_START = RESOURCE_START_DEFAULT_ECAL;
		RESOURCE_END = RESOURCE_END_DEFAULT_ECAL;
	}
	if ( version == VERSION_BIOL )
	{
		POP_SIZE = POP_SIZE_DEFAULT_BIOL;
		DEATH_RATE = DEATH_RATE_DEFAULT_BIOL;
		LAMBDA = LAMBDA_DEFAULT_BIOL;
		ALPHA = ALPHA_DEFAULT_BIOL;
		BETA = BETA_DEFAULT_BIOL;
		RESOURCE_START = RESOURCE_START_DEFAULT_BIOL;
		RESOURCE_END = RESOURCE_END_DEFAULT_BIOL;
	}
}

int main ( int argc, char **argv ) {
	// set default values for parameters
	set_default_params(VERSION);

	if ( ! parse_options(argc, argv) )
		return -1;
	
	// initialize the random number generator
	// (a new version of the mersenne twister PRNG)
	init_gen_rand(RAND_SEED);
	INFO(
			fprintf(stderr, "Random seed: %d\n", RAND_SEED);
		);
	
	// initialize other useful things
	
	// no need to evaluate this every time
	// (initialized here anyway)
	// used for the 1st version
	over_sqrt_lambda = 1/sqrt(LAMBDA);
	// used for the 2nd version
	over_lambda2 = 1/pow(LAMBDA, 2);

	/////
	
	// starts the simulation
	if ( pthread_create(&simulation_thread, NULL, simulation, NULL) != 0 )
	{
		printf("Error in creating simulation thread\n");
		return -1;
	}
	
	// waits for the simulation thread to finish
	pthread_join(simulation_thread, NULL);
	
	return 0;
}
