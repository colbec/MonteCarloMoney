# MonteCarloMoney

Examining the role of money in a small economy using a Monte Carlo approach

In this code I am trying to examine a small economy and how personal money in the form of promissory notes 
can perform a time shifting function. Small farmers grow beans - some farmers produce a surplus and some
have a deficit and need seeds. Assuming that some farmers have varying confidence in others they may or may
not be willing to lend surplus seed in exchange for promissory notes. Where trade is possible, the GDP of the 
economy improves.

The code uses a network graph approach. To use the code, we decide on a number of nodes (farmers) and 
connect a random number of these nodes with edges carrying random weights. Initial endowment of bean seeds
is distributed with (0,1) normal distribution.

To use the code, open an instance of Julia (initially tested using version 1.8.4), include() the code in sim1.jl
and issue the command:

monte = main(20,16,100,10);

This asks for a network of 20 nodes, 16 of which are connected with edges. Edge weights are randomly applied.
The Monte Carlo run repeats 100 times and at every 10 produces a plot and a report on matches found.
The output is a 100x4 matrix containing the results of each run in each row.

Sample output from a Monte Carlo run can be viewed at https://www.youtube.com/watch?v=s88b38EcXfc
