A Simple Genetic Algorithm
==========================
This code is a Python3 implementation of the Simple Genetic Algorithm presented in book 
Genetic Algorithms + Data Structures = Evolution Programs (3ed) from Michalewicz Z.

Options
-------

**Usage:** sga [options]

**Optional arguments:**

  `-h, --help`                                          Show help about this program
  
   `-p POPSIZE, --popsize POPSIZE`                      GA Population Size
  
  `-g MAXGENS, --maxgens MAXGENS`                       GA Generations
  
  `-n NVARS, --nvars NVARS`                             GA N Vars
  
  `-a ANTIGENPOPSIZE, --antigenpopsize ANTIGENPOPSIZE`  Antigens Population Size
  
  `-r RUNTIMES, --runtimes RUNTIMES`                    Number of times GA will be called. If this option is set above 1 
                                                        then GA will be called multipple times and a mean fitness 
                                                        is returned
  
  `-x CROSSOVER, --crossover CROSSOVER`                 Crossover Probability
  
  `-m MUTATION, --mutation MUTATION`                    Mutation Probability
  
  `-v, --verbose`                                       Print verbose info
 
*this is a work in progress*