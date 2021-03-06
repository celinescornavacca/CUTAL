Here is the underlying data from our experiments. Please refer to the experimental section of the paper, and its appendix, for the full
context.

##########  2-block experiment ##########

* The EXPO files, in the folder "final2block" are for the 2-block experiment. They have the following naming convention:

EXPO_numberOfTaxa_branchLength_LengthOfSmallBlock_replicateNumber . fileType

where filetype is:

dawg = the file we input to DAWG to tell it how to create the alignment
phylip = the alignment generated by DAWG (which we then feed to CUTAL)
runtime = the running time for CUTAL on this input (this includes the time to run Parsimonator)
output = the output of CUTAL on this input.

(The _MAXHOMSCORE, _TOTHOMSCORE, _MAXHOMRATIO and _TOTHOMRATIO files summarize the boundary error, and optimal number of blocks, for the
corresponding objective function, ranging over the 20 replicates.)

In the .output files, about six lines down there is a line "Input block partition" which is telling us the experimentally specified boundaries: we tell
CUTAL this explicitly, so CUTAL can compute the breakpoint errors itself. (Above all that is a vector explaining which sites are informative). Then
the most important bit is at the end, "Comparsion to input partition" where CUTAL computes the breakpoint error for each of the 4 objective functions.
Note that this is for the -specified- number of blocks (in this case: 2), even if for a given objective function the optimum was achieved by fewer
blocks.

* The excel file "final2block" summarizes the results of the 2-block experiment. This contains slightly more information than the corresponding
table in the paper. The extra information is contained in the "# replicates too few blocks?" columns and the rightmost column,
"Number of discarded Dawg alignments (due to at least one block being uninformative)". There is also a tab summarizing running times
(time in seconds / name of replicate).

The "# replicates too few blocks?" counts the number of replicates, for that parameter combination and that objective function, in which
the algorithm for that objective function concluded that the optimal number of blocks was -lower- than the experimentally specified
number of blocks.

In the "Number of discarded Dawg alignments (due to at least one block being uninformative)" column we list the the number of extra times that we
had to call Dawg in order to obtain an aliginment that does not contain any uninformative blocks (i.e. how many Dawg calls -- above the baseline
of 20 -- were in total necessary to produce 20 replicates, none of which contained uninformative blocks). 

########## multi-block experiment ##########

* The EGGSPO files, in the folder "finalMultiblock", are for the multi-block experiment.

The content of the files is the same as above, but the naming convention is slightly different:

EGGSPO_numberOfTaxa_branchLength_numberOfBlocks_replicateNumber_{vector of lengths of the blocks} . fileType

* The excel file "finalMultiblock" summarizes the results of the multiple-block experiment. It relates to Table 2, in the same way
the file "final2block" relates to Table 1, so please see above for the explanation of the meaning of the extra columns.


The authors, September 2019.




