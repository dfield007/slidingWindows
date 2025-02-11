# Sliding window genome scan (SlidingWindows_v1.4.py)

### What does this script do?

This script calculates various diversity (pi_w) and divergence (Dxy) and differentiation (Fst) statistics for individual or pooled sequence data for two or more populations. Two types of output are generated: (i) individual site statistics, (ii) sliding window averages of the same statistics (with adjustable overlaps). Different data manipulations with regards to data quality can be set by the user including, minimum and maximum depth and the minimum number of populations with required depth range for site inclusion. The script can also work with individual based data, as long as this is in the required format.

## Quick guide

### How to install

Written for Python 3+. Simply place the script and all required input files in a common directory. Python modules required include: re, sys, getopt, itertools, os, math.

Ideally place this in a directory included in your PATH. See online how to edit your PATH 

Some Python modules that may not come with base installation are required. If some modules are missing, download and install the required modules with, for example:

	> sudo apt-get install python-numpy python-scipy

### Input files

•	Genomic data: in this version, the format is fixed on *.sync files as used for pooled data in the program Popoolation2 (Kofler et al., 2010).
•	Scaffold data: a file listing the set of scaffolds or chromosomes to analyse. This will assume *.sync files exist with the same names.
•	Population details: a file listing population names and sample sizes.

### How to run the program

Assuming you have the script and the required input files in your PATH directory, you can run the program as follows:

> python SlidingWindows_v1.14.py scaffoldList_example.txt popDetails_example.txt example_out 15 200 2 2 10000 0 1 0 0 1

where the system arguments for the script include (corresponding values in example):

to run, you need the following arguments: 
to run, you need the following arguments: 
(1) list of scaffolds to analyse, 
(2) file with population names and sample size, 
(3) output file name,
(4) the min depth to include a site, 
(5) the max depth to include a site, 
(6) min copies for an allele call (i.e. > 1, no singletons)
(7) min number of populations (pools) with min depth to include a site, 
(8) window size (bp), 
(9) overlap between windows (bp),
(10) keep site specific file (1 = Y,0 = N), 
(11) window analysis only (1 = Y,0 = N),
(12) site analysis only (1 = Y, 0 = N) 
(13) keep all positions (1 = Y,0 = N),

min depth: if running on individual based data, make minimum depth the size of the smallest number of haploid genomes (i.e. 4x2 = 8),
min copies: if running on individual based data, make minimum copies = 1 (i.e. allow for singletons),
input: scaffold input list must have line breaks saved as Unix (LF). Use text wrangler or similar text editor,
window analysis only: only works if a site specific file(s) already created and present in the same folder, min scaffold size: currently skips over any scaffold < 15kB. Ensures at least two windows minimum.


## News/updates

11/02/2025 - v1.4 released. Updated for Python 3+

14/05/2018 - v1.10 released. Now deals with non-consecutive position errors in the input file (if exists), whilst accounting for real sequence gaps.

31/05/2017 - v0.4 realeased. Now allows for any window size and degree of window overlap. 

For more details see *SlidingWindows_readme.pdf*
