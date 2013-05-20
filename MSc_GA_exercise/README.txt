Genetic algorithm exercise
coursework for Artificial Life, Sussex University 2007
lorenzo grespan
oct 16, 2007

Exercise
--------
	Chose a way to divide a set of 10 cards into two piles, such that the cards in the first pile sum to a number close to 36 and the cards in the second pile multiply to a number close to 360.

Initial analysis
----------------
	The genotype chosen is a binary array of length 10. A value of '1' in positin n indicates a card numbered n+1. For example, if the first pile contains cards 1, 3, 5, 7, 9 its encoding would be

1 0 1 0 1 0 1 0 1 0

The second pile would be complementary: cards that are in a pile cannot be in the other. Therefore, the second pile would be

0 1 0 1 0 1 0 1 0 1

this configuration can be obtained by bit-flipping the first pile.
Also, the number of possible combinations is 2^10 = 1024 - which is the size of the entire population.

The fitness function(s)
-----------------------
	It has been chosen to calculate the fitness as absolute distance from the target; so the fitness for the nth pile would be
	
abs ( f(pile_n) - TARGET_PILEn )

where f(pile_n) is the appropriate phenotype function. The phenotype of the first pile is the sum of the array index of the elements marked '1'. If a genotype is

1 0 1 0 1 0 1 0 1 0

its phenotype for the first pile would be

1 + 3 + 5 + 7 + 9 = 25

and its phenotype for the second pile would be 

3 * 5 * 7 * 9 = 945

The GA
------
	The genetic algorithm used is a variance of the "microbial", where two individuals are chosen (initially at random) and then compared for their fitness regarding the first pile. The loser gets one of its genes randomly mutated. The process is reiterated until there are no more mutations left to be tried.  After this first trial, the resulting winner and loser are bitflipped and compared for their fitness regarding the second pile. 
	This process leads first to a couple of individuals that are as close as possible to the first value (sum near to 36). Afterwards the resulting individuals are compared and evolved to be as close as possible to 360.

The results
-----------
	The simulation has been tried 1000 times on randomly selected initial individuals; the order of the fitness functions has been changed to see if testing for multiplication would lead to a faster convergence with respect of the testing for sum (see files "cards2-sumfirst.txt" and "cards2-prodfirst.txt"). The results have been plotted putting on the x-axis the number of the trial and on the y-axis the number of runs to achieve a value close to the requested results. 
	As can be seen from the graphs "prodfirst-graph.jpg" and "sumfirst-graph.jpg" there is no noticeable difference, given the fact that the y-axis scale is slightly different due to some high values appearing from time to time. 
	The result show also that, on average, all near-optimal solution are located below the boundary of 4000 mutations.

Analytical considerations
-------------------------
	Given a search space of 1024 possible individuals, and given the fact that some of them are valid only with respect to only one of the fitness functions (many have sum = 36 but very different products), a GA seems like a good solution to the problem. A standard search algorithm could take as most as n*log(n) to complete (quicksort with randomized partitioning), which should be used if an optimal solution were to be found; however, with a larger population size its complexity would probably grow more than the complexity required by the GA. 
	A more approfondite analysis is required to determine the complexity of the GA given the actual configurations (loser gets only mutated one gene away).
	Also, because the algorithm stops when the loser cannot mutate any more, this approach has a bias towards local minima.

Notes
-----
	The seed for the random number generator is given by the user; in this example a pseudo-random number generator has been used to seed the program via the BASH environment variable $RANDOM.
	Cross-mutating has not been chosen as a viable alternative because even small changes in the genotype can lead to very large changes in the phenotype, expecially if the changes are made on the leftmost genes (see below). There is some kind of cross-mutating, however, when the winner and the loser keep exchanging their roles; they become individuals with just a single difference, and the subsequent mutation of the loser is equivalent of a random cross-mutation as taken from the winner. For example, supposing a simplified environment with a target of 100,

	step 0:
		a: 	001
		b:	000 (a wins, b mutates)
	step 1: 
		a:	001
		b:	010 (b wins, a mutates)
	step 2:
		a:	011
		b:	010 (a wins, b mutates)

step 2 could be also achieved by cross-mutating gene no. 2 from B. Obviously, since the mutating gene is chosen at random, there is some resemblance to the random-selection of which gene to cross in the microbial algorithm. 
	Also it should be noted that this specific genetic encoding is not optimal for both the fitness functions: for example, given an individual like

	- - - - 5 6 - - - -
	sum: 11
	prod: 30

a small change in one of its genes can lead to

	1 - - - 5 6 - - - -
	sum: 12
	prod: 30

or

	- - - - 5 6 - - - 10
	sum: 22
	prod: 300

As can be seen, small changes can both lead to small or no changes in both functions. However, it can be said that mutating genes on the left hand side of the genotype leads to smaller changes than mutating on the right hand side.

