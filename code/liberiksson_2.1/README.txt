liberiksson:

This is C++ code which classifies the minimally non-Deodhar elements
of various (finite) Coxeter groups by embedded factor containment.  It
also produces candidates for embedded factor patterns in type D, and
verifies that all of the mu coefficients for Deodhar elements are 0 
or 1.

The CoxeterSystem and CoxeterElement files define a shared library
(liberiksson) which implements a "numbers game" of Kimmo Eriksson,
described in his PhD thesis, and in Bjorner-Brenti's "Combinatorics of
Coxeter groups."  The game is used to enumerate the elements of
Coxeter groups, and perform basic operations on Coxeter elements such
as multiplication, and the determination of (a cannonical) reduced
expression.  

The client program verify.cpp contains the main() function, performs 
the classification, and generates enumerative data.

To compile the code with gcc on unix, use "make unix".
To compile the code with gcc on cygwin, use "make cygwin".
Once compiled, set LD_LIBRARY_PATH=. so that the shared library can be
loaded.
Then, run "./verify".

It takes about 4 minutes to run the classification on the modern
multi-processor machine at our university.


Here are #fully-commutative elements in ranks 2 - 14:  (See
Stembridge's paper on "The enumeration of fully commutative elements
of Coxeter groups" for details.)

Type A:  {5,14,42,132,429,1430,4862,16796,58786,208012,742900,2674440,9694845}
Type B:  {7,24,83,293,1055,3860,14299,53481,201551,764217,2912167,11143499,42791039}
Type D:  {4,14,48,167,593,2144,7864,29171,109173,411501,1560089,5943199,22732739}
Type E:
 E6:  662
 E7:  2670
 E8:  10846


The total number of elements in type D is 2^n-1 n!:
  n = 6:  23040
  n = 7:  322560
  n = 8:  5160960
  n = 9:  92897280
  n = 10: 1857945600  (btw. 12! and 13!)
  n = 11: 40874803200 (btw. 13! and 14!)


----------------------------------------------------------------------
verify.cpp output:
----------------------------------------------------------------------

Type G2: 
0 6 
6 0 

  (evaluating length 0 elements, with 0 elements left to process...) 
  (evaluating length 1 elements, with 1 elements left to process...) 
  (evaluating length 2 elements, with 1 elements left to process...) 
Finished:  found 5 Deodhar elements (out of 5 short-braid-avoiding elements processed, used 0 non-Deodhar patterns).  Also verified mu in {0,1} property for these Deodhar elements (check stdout for errors).
  (Consistency:  check that #short-braid-avoiding elts = 5 = total non-Deodhar elts + Deodhar elts.)

G2 has 5 short-braid-avoiding elements.


Type F4: 
0 3 2 2 
3 0 4 2 
2 4 0 3 
2 2 3 0 

  (evaluating length 0 elements, with 0 elements left to process...) 
  (evaluating length 1 elements, with 3 elements left to process...) 
  (evaluating length 2 elements, with 11 elements left to process...) 
  (evaluating length 3 elements, with 17 elements left to process...) 
  (evaluating length 4 elements, with 13 elements left to process...) 
  (evaluating length 5 elements, with 5 elements left to process...) 
  (evaluating length 6 elements, with 1 elements left to process...) 
Finished:  found 42 Deodhar elements (out of 42 short-braid-avoiding elements processed, used 0 non-Deodhar patterns).  Also verified mu in {0,1} property for these Deodhar elements (check stdout for errors).
  (Consistency:  check that #short-braid-avoiding elts = 42 = total non-Deodhar elts + Deodhar elts.)

F4 has 42 short-braid-avoiding elements.



Minimal non-Deodhar embedded factor patterns: 
Generating minimal patterns for A7 (includes linear types BC, F4, G2, H3, H4): 
0 3 2 2 2 2 2 
3 0 3 2 2 2 2 
2 3 0 3 2 2 2 
2 2 3 0 3 2 2 
2 2 2 3 0 3 2 
2 2 2 2 3 0 3 
2 2 2 2 2 3 0 

  (evaluating length 0 elements, with 0 elements left to process...) 
  (evaluating length 1 elements, with 6 elements left to process...) 
  (evaluating length 2 elements, with 41 elements left to process...) 
  (evaluating length 3 elements, with 134 elements left to process...) 
  (evaluating length 4 elements, with 284 elements left to process...) 
  (evaluating length 5 elements, with 434 elements left to process...) 
  (evaluating length 6 elements, with 502 elements left to process...) 
  (evaluating length 7 elements, with 471 elements left to process...) 
  (evaluating length 8 elements, with 386 elements left to process...) 
  (evaluating length 9 elements, with 291 elements left to process...) 
  (evaluating length 10 elements, with 198 elements left to process...) 
  (evaluating length 11 elements, with 123 elements left to process...) 
  (evaluating length 12 elements, with 67 elements left to process...) 
  (evaluating length 13 elements, with 33 elements left to process...) 
  (evaluating length 14 elements, with 13 elements left to process...) 
Found MINIMAL PATTERN of rank 7:  [ 2 1 -6 7 -6 1 2 ] { 4 6 7 1 8 2 3 5 } ( 4 5 6 2 3 4 5 1 2 3 4 0 1 2 )
  (evaluating length 15 elements, with 3 elements left to process...) 
Finished:  found 1426 Deodhar elements (out of 1428 short-braid-avoiding elements processed, used 4 non-Deodhar patterns).  Also verified mu in {0,1} property for these Deodhar elements (check stdout for errors).
  (Consistency:  check that #short-braid-avoiding elts = 1430 = total non-Deodhar elts + Deodhar elts.)

A7 has 1430 short-braid-avoiding elements.


Generating minimal patterns for D8 (excluding 1-line pattern {-1, 6, 7, 8, -5, 2, 3, 4, 9}): 
0 2 3 2 2 2 2 2 
2 0 3 2 2 2 2 2 
3 3 0 3 2 2 2 2 
2 2 3 0 3 2 2 2 
2 2 2 3 0 3 2 2 
2 2 2 2 3 0 3 2 
2 2 2 2 2 3 0 3 
2 2 2 2 2 2 3 0 

  (evaluating length 0 elements, with 0 elements left to process...) 
  (evaluating length 1 elements, with 7 elements left to process...) 
  (evaluating length 2 elements, with 55 elements left to process...) 
  (evaluating length 3 elements, with 209 elements left to process...) 
  (evaluating length 4 elements, with 531 elements left to process...) 
  (evaluating length 5 elements, with 1000 elements left to process...) 
  (evaluating length 6 elements, with 1477 elements left to process...) 
  (evaluating length 7 elements, with 1799 elements left to process...) 
  (evaluating length 8 elements, with 1906 elements left to process...) 
  (evaluating length 9 elements, with 1835 elements left to process...) 
  (evaluating length 10 elements, with 1659 elements left to process...) 
  (evaluating length 11 elements, with 1422 elements left to process...) 
Found MINIMAL PATTERN of rank 6:  [ -7 1 8 -6 1 2 4 1 ] { -5 -4 6 -2 -1 3 7 8 9 } ( 3 4 5 0 2 3 4 1 2 3 0 )
  (evaluating length 12 elements, with 1174 elements left to process...) 
  (evaluating length 13 elements, with 923 elements left to process...) 
  (evaluating length 14 elements, with 689 elements left to process...) 
Found MINIMAL PATTERN of rank 7:  [ 1 9 -7 8 -7 1 3 4 ] { -5 6 -3 7 -2 -1 4 8 9 } ( 4 5 6 0 2 3 4 5 1 2 3 4 0 2 )
Found MINIMAL PATTERN of rank 7:  [ -9 1 2 8 -6 1 2 4 ] { -6 -5 -3 7 -1 2 4 8 9 } ( 4 5 6 2 3 4 5 0 2 3 4 1 2 0 )
Found MINIMAL PATTERN of rank 7:  [ 2 8 1 -6 7 -6 1 2 ] { -4 6 7 -1 8 2 3 5 9 } ( 5 6 7 3 4 5 6 2 3 4 5 0 2 3 )
  (evaluating length 15 elements, with 465 elements left to process...) 
Found MINIMAL PATTERN of rank 7:  [ 9 1 1 -9 3 1 1 5 ] { 5 6 7 -4 -1 2 3 8 9 } ( 3 4 5 6 2 3 4 5 0 2 3 4 1 2 3 )
  (evaluating length 16 elements, with 262 elements left to process...) 
  (evaluating length 17 elements, with 108 elements left to process...) 
  (D8 1-line pattern found in [ 5 5 1 1 -11 5 1 1 ] { -1 6 7 8 -5 2 3 4 9 }( 4 5 6 7 3 4 5 6 2 3 4 5 1 0 2 3 4 ) of rank 8)
  (evaluating length 18 elements, with 21 elements left to process...) 
  (evaluating length 19 elements, with 1 elements left to process...) 
Finished:  found 6791 Deodhar elements (out of 7064 short-braid-avoiding elements processed, used 1073 non-Deodhar patterns).  Also verified mu in {0,1} property for these Deodhar elements (check stdout for errors).
  (Consistency:  check that #short-braid-avoiding elts = 7864 = total non-Deodhar elts + Deodhar elts.)

D8 has 7864 short-braid-avoiding elements.


Generating minimal patterns for E7: 
0 3 2 2 2 2 2 
3 0 3 2 2 2 2 
2 3 0 3 2 3 2 
2 2 3 0 3 2 2 
2 2 2 3 0 2 3 
2 2 3 2 2 0 2 
2 2 2 2 3 2 0 

  (evaluating length 0 elements, with 0 elements left to process...) 
  (evaluating length 1 elements, with 6 elements left to process...) 
  (evaluating length 2 elements, with 41 elements left to process...) 
  (evaluating length 3 elements, with 134 elements left to process...) 
  (evaluating length 4 elements, with 289 elements left to process...) 
  (evaluating length 5 elements, with 460 elements left to process...) 
  (evaluating length 6 elements, with 575 elements left to process...) 
  (evaluating length 7 elements, with 607 elements left to process...) 
  (evaluating length 8 elements, with 579 elements left to process...) 
  (evaluating length 9 elements, with 510 elements left to process...) 
  (evaluating length 10 elements, with 418 elements left to process...) 
  (evaluating length 11 elements, with 332 elements left to process...) 
Found MINIMAL PATTERN of rank 6:  [ 9 -7 8 -6 1 1 2 ] { 1 6 3 7 5 2 8 4 } ( 3 4 6 1 2 5 3 4 2 3 1 )
Found MINIMAL PATTERN of rank 6:  [ 5 1 8 -6 1 -7 2 ] { 1 5 6 4 2 3 8 7 } ( 3 4 6 5 2 3 1 2 5 4 3 )
  (evaluating length 12 elements, with 253 elements left to process...) 
  (evaluating length 13 elements, with 177 elements left to process...) 
Found MINIMAL PATTERN of rank 6:  [ 1 -8 9 1 1 -8 6 ] { 5 6 2 3 4 1 7 8 } ( 0 1 2 5 3 4 2 3 1 2 5 0 1 )
Found MINIMAL PATTERN of rank 6:  [ -11 2 1 2 1 1 6 ] { 7 3 4 5 1 6 2 8 } ( 5 1 2 3 0 1 2 5 4 3 2 1 0 )
Found MINIMAL PATTERN of rank 6:  [ 1 1 -9 10 1 2 5 ] { 5 6 3 1 4 2 7 8 } ( 1 2 5 3 4 2 3 1 2 5 0 1 2 )
Found MINIMAL PATTERN of rank 6:  [ 11 -10 2 1 1 1 7 ] { 4 7 3 5 1 6 2 8 } ( 2 5 1 2 3 0 1 2 5 4 3 2 1 )
Found MINIMAL PATTERN of rank 6:  [ 1 1 2 -10 11 1 4 ] { 5 6 2 4 1 3 7 8 } ( 2 5 3 4 2 3 1 2 5 0 1 2 3 )
Found MINIMAL PATTERN of rank 6:  [ 1 10 -9 1 1 2 8 ] { 5 3 7 4 1 6 2 8 } ( 3 2 5 1 2 3 0 1 2 5 4 3 2 )
Found MINIMAL PATTERN of rank 6:  [ 1 2 1 2 -11 1 14 ] { 5 6 2 3 4 1 7 8 } ( 5 3 4 2 3 1 2 5 0 1 2 3 4 )
Found MINIMAL PATTERN of rank 6:  [ 1 1 9 -8 1 -8 9 ] { 6 3 4 7 1 5 2 8 } ( 4 3 2 5 1 2 3 0 1 2 5 4 3 )
  (evaluating length 14 elements, with 99 elements left to process...) 
  (evaluating length 15 elements, with 39 elements left to process...) 
Found MINIMAL PATTERN of rank 7:  [ 1 -8 9 -8 1 8 1 ] { 5 6 2 7 3 4 8 1 } ( 0 1 2 3 4 6 5 2 3 4 1 2 3 0 1 )
Found MINIMAL PATTERN of rank 7:  [ -13 2 1 4 1 1 2 ] { 7 3 5 6 1 2 8 4 } ( 3 4 6 1 2 3 0 1 2 5 4 3 2 1 0 )
Found MINIMAL PATTERN of rank 7:  [ 1 1 -9 2 1 16 1 ] { 5 6 7 1 3 4 8 2 } ( 1 2 3 4 6 5 2 3 4 1 2 3 0 1 2 )
Found MINIMAL PATTERN of rank 7:  [ 13 -12 1 5 1 1 1 ] { 4 7 5 6 1 2 8 3 } ( 2 3 4 6 1 2 3 0 1 2 5 4 3 2 1 )
  (evaluating length 16 elements, with 5 elements left to process...) 
Found MINIMAL PATTERN of rank 7:  [ 7 1 1 5 1 -13 1 ] { 4 7 6 1 5 8 2 3 } ( 5 2 3 4 6 1 2 5 3 4 2 3 0 1 2 5 )
Finished:  found 2341 Deodhar elements (out of 2419 short-braid-avoiding elements processed, used 329 non-Deodhar patterns).  Also verified mu in {0,1} property for these Deodhar elements (check stdout for errors).
  (Consistency:  check that #short-braid-avoiding elts = 2670 = total non-Deodhar elts + Deodhar elts.)

E7 has 2670 short-braid-avoiding elements.


