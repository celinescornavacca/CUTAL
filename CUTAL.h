/*  CUTAL (version 0.1) a program for cutting alignement into non-recombinant blocks under 
 *  a parsimony framework 
 *  Copyright May 2018 by Celine Scornavacca and Mark Jones
 * 
 *  The program is heavily based on Parsimonator-1.0.2 by Alexandros Stamatakis, which was partially 
 *  derived from fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an email to celine.scornavacca@umontpellier.fr and 
 *  markelliotlloyd@gmail.com
 *
 */


#include <assert.h>
#include <stdint.h>



#define nmlngth        256         /* number of characters in species name */



#define NUM_BRANCHES 10


#define TRUE             1
#define FALSE            0


#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))

#define programName        "Parsimonator"
#define programVersion     "1.0.3"
#define programDate        "June 2013"




#define M_GTRCAT         1




#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2


#define DNA_DATA         1










typedef  int boolean;




typedef unsigned int parsimonyNumber;



#define PCF 32






typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;




struct noderec;













typedef  struct noderec
{
  
  struct noderec  *next;
  struct noderec  *back;
  int              number;
  char             x;
}
  node, *nodeptr;





typedef  struct
{
  int              numsp;
  int              sites;
  unsigned char             **y;
  unsigned char             *y0;
  unsigned char             *yBUF;
  int              *wgt;
} rawdata;

typedef  struct {
  int             *alias;       /* site representing a pattern */
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory;
  int              endsite;     /* # of sequence patterns */
  double          *patrat;      /* rates per pattern */
  double          *patratStored; 
} cruncheddata;








typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;

struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};

typedef struct stringEnt stringEntry;

typedef unsigned int hashNumberType;

typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;

typedef  struct  {
 
 

  parsimonyNumber **parsimonyState_A;
  parsimonyNumber **parsimonyState_C;
  parsimonyNumber **parsimonyState_G;
  parsimonyNumber **parsimonyState_T;
  unsigned int *parsimonyScore; 
  int *ti;
  unsigned int compressedWidth;

  unsigned char             **yVector;
  

 
 
  stringHashtable  *nameHash;

 
  node           **nodep;
  node            *start;
  int              mxtips;
   
  int              ntips;
  int              nextnode;
  
 
  rawdata         *rdta;
  //cruncheddata    *cdta; //not needed anymore 
  int             *minParsimonyPerSite;  //minimum possible parsimony scores for each site (used to calculate homoplasy)
  int             *recordOfInformativeSites;  //information about which sites are informative (1 for informative, 0 for uninformative)

  char **nameList;
  char *tree_string;
  int *nodesInTree;

  int treeStringLength;
  unsigned int bestParsimony;
  
  nodeptr removeNode;
  nodeptr insertNode;

} tree;


/***************************************************************/



/**************************************************************/






typedef  struct {


  long          parsimonySeed;
  boolean       restart;
  int           numberOfTrees;
  int 			numberOfBlocks;
  float 		maxHomoplasyRatio;
  int	 		maxHomoplasyScore;
  float 		totalHomoplasyRatio;
  int	 		totalHomoplasyScore;
  int 			verbose;
  int 			problems;
  char  		trueBlockPartitionString[1024];
  int 			truePartitionBlockCount;
  int* 			truePartitionStartPoints; 
  int* 			truePartitionEndPoints;
} analdef;


/****************************** FUNCTIONS ****************************************************/


extern int makeParsimonyTreeFastDNA(tree *tr, analdef *adef, int startSite, int endSite);
extern int getMinParsimonyScoreForSite(tree *tr, int site);
extern boolean isInformative(tree *tr, int site);
extern void determineUninformativeSites(tree *tr, int *informative);
extern void printBothOpen(const char* format, ... );
extern double gettime(void);
extern unsigned int precomputed16_bitcount (unsigned int n);
extern FILE *myfopen(const char *path, const char *mode);
extern void treeReadLen (FILE *fp, tree *tr);
