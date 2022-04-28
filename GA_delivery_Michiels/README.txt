===== Genetic Algorithm Delivery - Eline Michiels =====

We give all commands assuming the working directory of the terminal is the GA_delivery subfolder.

== FILE STRUCTURE ==
We split our functions over multiple files to make the code more easy to understand.

GA.c:		Main program, calls all functions to execute the genetic algorithm. Also contains help functions necessary for the genetic algorithm.
RKF78.c:	Functions and definitions to integrate the ODE. These are used in the funcion Generate_EDO_Prediction.
randombits.c:	Auxiliary functions that garanty a true randomness when this is needed. 


- Executing the algorithm -
(I used the Cygwin terminal to compile and run my files)
The program in GA.c calls all necessary functions to execute the algorithm.
It also needs the auxiliary definitions in RKF78.c and randombits.c, so the command for building the program is:
gcc -Wall -o GA GA.c RKF78.c randombits.c

To run, we simply put:
./GA

This will run the GA-algorithm 1 time with population size 10000 and 60 generations. Alterations can be made in the main in GA.c.
The output ( the predicted values for the population of Andouin's gulls and the parameters) will be printed in the terminal. 