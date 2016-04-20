# ld_simulator


REQUIREMENTS

Python 2.7x.


INSTALLATION

Copy the file ld_simulator.py to a folder. Either make it executable or run by calling "python ld_simulator.py".


DESCRIPTION

A simple model of neutral evolution in a population of diploid individuals. You can use ld_simulator.py when teaching a class to illustrate how linkage disequilibrium is created.

The program creates a population of 500 organisms with no genetic variation at the beginning. The genome of the organisms in the population is very small: only one chromosome with 150 base pairs.

For each generation, 125 individuals are randomly selected to mate and create the next generation of 500 offspring and so forth. During each mating, genetic recombination and point mutations may occur and change the genetic makeup of the next generation. Both the mutation and recombination rates of these organisms are greater than e.g. in humans, so their evolution is pretty quick! At the end of each generation, the program randomly samples 24 individuals from the population and shows their entire genomes on the screen, so that you can see what's going on with your organisms.

Launch the program by typing

> ld_simulator.py

and pressing enter. This makes the program print the genomes of the first generation on the screen. Then, you can move forward in time by holding down the Enter key and see how the genetic makeup of the population changes. If you stop pressing Enter, time stops as well. You can quit by typing 

> exit

on your keyboard and then pressing the Enter key. 	


Launching the program by typing

> ld_simulator.py --realistic

and pressing enter makes the program does the same as above but it now uses the known human mutation and recombination rates, to give an indication of the true speed of evolution. 

Launching the program by typing

> ld_simulator.py --sample-disease

makes the program automatically fast-forward 500 generations into future before printing anything. Now, mutation occasionally introduces disease alleles (D) into your population. On the first line the program prints the disease mutation as it happened before the genomes of your organisms (the ancestral haplotype). It also marks the polymorphic site(s) most strongly correlated with the disease allele with an asterisk ("*"), if at least one of the sites is correlated.

When printing the genomes, the program now first prints the genomes of 12 randomly selected sick individuals (if there are at least 12 sick organisms in the entire population. If there are less, the program prints as many as there are).

You can save a snapshot of the program's output at this point by typing

> save

on the keyboard and pressing Enter.

