#Hierarchical Bin Buffering C++ Library

Daniel Lemire and Owen Kaser


This is the C++ source code for the paper: 

Daniel Lemire and Owen Kaser, Hierarchical Bin Buffering: Online Local Moments for Dynamic External Memory Arrays, ACM Transactions on Algorithms 4 (1), pages 1-31, 2008. 
http://arxiv.org/abs/cs.DS/0610128

This paper includes a survey as well as new methods to precompute 
polynomial range queries, as they are used in polynomial curve fitting
and statistics.

Examples of queries that can be written as polynomial range queries
include: the sum of the array between index i and index j, 
the center of mass of the array between the index i and the index j, and 
so on.

To build the software, do (gcc toolchain is required):

1) Go to the dubuc directory and type "make", this builds some of the essential stuff;
2) Go to lemurcore and type "make";
3) Still in the dubuc director, type "make benchmark1";

To check the code, type "make test"


