# CUTAL: Cut Alignment for genetic sequences

## Compiling

run `make -f Makefile.gcc`


## EXAMPLES

# Run CUTAL on file test.txt and write output to testoutput.out

```sh
./CUTAL -s test.txt -o testoutput.out
```


# Run CUTAL on file test.txt solving only for problem Max Homoplasy Score, with up to 3 blocks allowed in the partition

```sh
./CUTAL -s test.txt -P 1 -b 3
```


# Run CUTAL on file test.txt solving only for problem Max Homoplasy Ratio, with up to 3 blocks allowed in the partition, and in addition find the block partition with fewest blocks (up to 3) that achieves a max homoplasy ratio of at most 0.05

```sh
./CUTAL -s test.txt -P 4 -b 3 -r 0.05
```

# Run CUTAL on file test.txt with randomization seed 345 and give a parsimony tree for each block in each block partition found

```sh
./CUTAL -s test.txt -p 345 -v 2
```


# Run CUTAL on file test.txt with up to 4 blocks allowed in the partition, and compare the returned solutions on 3 blocks with the block partition [0---3], [6---15], [17---20] 
```sh
./CUTAL -s test.txt -B [0---3],[6---15],[17---20] -b 4
```

./CUTAL -s filename [Options]



## ARGUMENTS AND OPTIONS

-s filename	Specifies the file from which to read sequence data. Files should be fully aligned
		and in PHYLIP format: 
		http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html
		This is the only argument that must be specified.


[-b int]	(Default 2) Sets the maximum number of block partitions to attempt. 
		CUTAL will attempt to find the optimal block partition (with respect to the 
		problem(s) under consideration) with at most this many blocks.


[-P int]	(Default 15) Sets which problems to solve. The last 4 bits of the binary expression
		of this number determine which of the problems {Max Homoplasy Score, Total 
		Homoplasy Score, Max Homoplasy Ratio, Total Homoplasy Ratio} to solve.
		Thus, for example:
		-P 1	 (= 0001)  means solve Max Homoplasy Score
		-P 2	 (= 0010)  means solve Total Homoplasy Score
		-P 4	 (= 0100)  means solve Max Homoplasy Ratio
		-P 8	 (= 1000)  means solve Total Homoplasy Ratio

		-P 5	 (= 0101)  means solve and Max Homoplasy Score Max Homoplasy Ratio

		By default, this value is set to 15 (= 1111), so all four problems will be solved.

[-h int]	If this value is specified and Max Homoplasy Score is one of the problems to be 
		solved, CUTAL will also return the block partition with fewest blocks that has 
		maximum homoplasy score at most int, if such a partition exists with at most
		the requested number of blocks).


[-H int]	If this value is specified and Total Homoplasy Score is one of the problems to be 
		solved, CUTAL will also return the block partition with fewest blocks that has 
		total homoplasy score at most int, if such a partition exists with at most
		the requested number of blocks).

[-r int]	If this value is specified and Max Homoplasy Ratio is one of the problems to be 
		solved, CUTAL will also return the block partition with fewest blocks that has 
		maximum homoplasy ratio at most int, if such a partition exists with at most
		the requested number of blocks).

[-R int]	If this value is specified and Total Homoplasy Ratio is one of the problems to be 
		solved, CUTAL will also return the block partition with fewest blocks that has 
		total homoplasy ratio at most int, if such a partition exists with at most
		the requested number of blocks).

[-B string]	Specifies an optional block partition to compare with the block partition(s) found
		by CUTAL.
		Must be a string of the form 
			[a---b],[c---d], ... [y---z] 
		Where a >= 0, z is at most the number of sites/characters -1, 
		any block [a--b] satisfies a <= b, 
		and any pair of consecutive blocks [a---b],[c---d] satisfies c > b.
		If this value is specified, CUTAL will return the max homoplasy score, total 
		homoplasy score, max homoplasy ratio and total homoplasy ratio of the given block
		partition; furthemore it will compare this block partition with an optimal block 
		partition (with respect to the problem(s) under consideration) that has the same 
		number of blocks.
		(The comparision includes the 'average internal boundary error' of the two 
		partitions, a measure of how much the two partitions disagree.)
		A warning will be given if the specified block partition contains a block
		with no informative sites, or if there exists an informative site not contained
		in any block.

[-v int]	(Default 1) Specifies how 'verbose' CUTAL should be.
		If set to 0, CUTAL will return the optimal solution to each problem under 
		consideration (as well as the solutions acheiving a certain desired value and/or
		comparisions to an input block partition, as specified by -h,-H,-r,-R and -B).
		If set to 1, CUTAL will also return the optimal solution(s) for each possible 
		number of of blocks, up to the specified maximum.
		If  set to 2, CUTAL will also return the homoplasy score and a parsimonious tree
		for each block in each partition.

[-o filename]	(Default output.out) Specifies a file to which all output information should be 
		written. Note that CUTAL appends to this file rather than overwriting; thus,
		calling CUTAL several times with the same output file will lead to that file 
		containing reports from multiple experiments.

[-p int]	(Default 1) Sets the randomization seed to be used by CUTAL.

-n, -N, -t	Undocumented parameters used by Parsimonator. We recommend not touching these in
		case something breaks.

