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



#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>



#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>

#endif

#include "CUTAL.h"



FILE   *INFILE;





char run_id[128] = "",  
  seq_file[1024] = "", 
  tree_file[1024]="",   
  resultFileName[1024] = "",    
  infoFileName[1024] = "", 
  randomFileName[1024] = "",
  partitionInfoFileName[1024] = "",
  trueBlockPartitionString[1024] = ""; 


static boolean whitechar (int ch);
static void  treeEchoContext (FILE *fp1, FILE *fp2, int n);

static stringHashtable *initStringHashTable(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)calloc(tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}


static hashNumberType  hashString(char *p, hashNumberType tableSize)
{
  hashNumberType h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}

 

static void addword(char *s, stringHashtable *h, int nodeNumber)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)malloc((size_t)(strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}

static int lookupWord(char *s, stringHashtable *h)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return p->nodeNumber;	  	
    }

  return -1;
}

static int treeFinishCom (FILE *fp, char **strp)
{
  int  ch;
  
  while ((ch = getc(fp)) != EOF && ch != ']') {
    if (strp != NULL) *(*strp)++ = ch;    /* save character  */
    if (ch == '[') {                      /* nested comment; find its end */
      if ((ch = treeFinishCom(fp, strp)) == EOF)  break;
      if (strp != NULL) *(*strp)++ = ch;  /* save closing ]  */
    }
  }
  
  if (strp != NULL) **strp = '\0';        /* terminate string  */
  return  ch;
} /* treeFinishCom */


static int treeGetCh (FILE *fp)         /* get next nonblank, noncomment character */
{ /* treeGetCh */
  int  ch;

  while ((ch = getc(fp)) != EOF) {
    if (whitechar(ch)) ;
    else if (ch == '[') {                   /* comment; find its end */
      if ((ch = treeFinishCom(fp, (char **) NULL)) == EOF)  break;
    }
    else  break;
  }
  
  return  ch;
} /* treeGetCh */


static boolean treeLabelEnd (int ch)
{
  switch (ch) 
    {
    case EOF:  
    case '\0':  
    case '\t':  
    case '\n':  
    case '\r': 
    case ' ':
    case ':':  
    case ',':   
    case '(':   
    case ')':  
    case ';':
      return TRUE;
    default:
      break;
    }
  return FALSE;
} 


static boolean  treeGetLabel (FILE *fp, char *lblPtr, int maxlen)
{
  int      ch;
  boolean  done, quoted, lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *) NULL; 
  else 
    if (lblPtr == NULL) 
      maxlen = 0;

  ch = getc(fp);
  done = treeLabelEnd(ch);

  lblfound = ! done;
  quoted = (ch == '\'');
  if (quoted && ! done) 
    {
      ch = getc(fp); 
      done = (ch == EOF);
    }

  while (! done) 
    {
      if (quoted) 
	{
	  if (ch == '\'') 
	    {
	      ch = getc(fp); 
	      if (ch != '\'') 
		break;
	    }
        }
      else 
	if (treeLabelEnd(ch)) break;     

      if (--maxlen >= 0) *lblPtr++ = ch;
      ch = getc(fp);
      if (ch == EOF) break;
    }

  if (ch != EOF)  (void) ungetc(ch, fp);

  if (lblPtr != NULL) *lblPtr = '\0';

  return lblfound;
}


static boolean  treeFlushLabel (FILE *fp)
{ 
  return  treeGetLabel(fp, (char *) NULL, (int) 0);
} 




static int treeFindTipByLabelString(char  *str, tree *tr)                    
{
  int lookup = lookupWord(str, tr->nameHash);

  if(lookup > 0)
    {
      assert(! tr->nodep[lookup]->back);
      return lookup;
    }
  else
    { 
      printf("ERROR: Cannot find tree species: %s\n", str);
      return  0;
    }
}


static int treeFindTipName(FILE *fp, tree *tr)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabel(fp, str, nmlngth+2))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   

  return  n;
} 

static boolean treeProcessLength (FILE *fp, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(fp)) == EOF)  return FALSE;    /*  Skip comments */
  (void) ungetc(ch, fp);
  
  if (fscanf(fp, "%lf", dptr) != 1) {
    printf("ERROR: treeProcessLength: Problem reading branch length\n");
    treeEchoContext(fp, stdout, 40);
    printf("\n");
    return  FALSE;
  }
  
  return  TRUE;
}

static int treeFlushLen (FILE  *fp)
{
  double  dummy;  
  int     ch;
  
  ch = treeGetCh(fp);
  
  if (ch == ':') 
    {
      ch = treeGetCh(fp);
      
      ungetc(ch, fp);
      if(!treeProcessLength(fp, & dummy)) return 0;
      return 1;	  
    }
  
  
  
  if (ch != EOF) (void) ungetc(ch, fp);
  return 1;
} 

static void  treeEchoContext (FILE *fp1, FILE *fp2, int n)
{ /* treeEchoContext */
  int      ch;
  boolean  waswhite;
  
  waswhite = TRUE;
  
  while (n > 0 && ((ch = getc(fp1)) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = TRUE;
    }
    else {
      waswhite = FALSE;
    }
    
    if (ch > '\0') {putc(ch, fp2); n--;}
  }
} /* treeEchoContext */

static boolean treeNeedCh (FILE *fp, int c1, char *where)
{
  int  c2;
  
  if ((c2 = treeGetCh(fp)) == c1)  return TRUE;
  
  printf("ERROR: Expecting '%c' %s tree; found:", c1, where);
  if (c2 == EOF) 
    {
      printf("End-of-File");
    }
  else 
    {      	
      ungetc(c2, fp);
      treeEchoContext(fp, stdout, 40);
    }
  putchar('\n');

  if(c1 == ':')    
    printf("RAxML may be expecting to read a tree that contains branch lengths\n");

  return FALSE;
} 



static void addElementLen (FILE *fp, tree *tr, nodeptr p)
{   
  nodeptr  q;
  int      n, ch, fres;
  
  if ((ch = treeGetCh(fp)) == '(') 
    { 
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      assert(0);	      
	    }
	  else 	    
	    assert(0);	    	   
	}
      
      q = tr->nodep[n];

      addElementLen(fp, tr, q->next);
      if (! treeNeedCh(fp, ',', "in"))             
	assert(0);
      addElementLen(fp, tr, q->next->next);
      if (! treeNeedCh(fp, ')', "in"))            
	assert(0);
          	
      (void) treeFlushLabel(fp);
    }
  else 
    {   
      ungetc(ch, fp);
      if ((n = treeFindTipName(fp, tr)) <= 0)
	assert(0);
      q = tr->nodep[n];
      tr->nodesInTree[n] = 1;

      if (tr->start->number > n)  
	tr->start = q;
      (tr->ntips)++;
    }
  
 
  fres = treeFlushLen(fp);
  if(!fres)
    assert(0);
      

  p->back = q;
  q->back = p;          
} 






void treeReadLen (FILE *fp, tree *tr)
{
  nodeptr  
    p;
  
  int      
    i, 
    ch; 

  for (i = 1; i <= tr->mxtips; i++)   
    tr->nodep[i]->back = (node *) NULL; 
 
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;     
    }

  tr->start = tr->nodep[tr->mxtips + 1];

  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;         

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');
        
  addElementLen(fp, tr, p);
  
  if (!treeNeedCh(fp, ',', "in"))                
    assert(0);
  
  addElementLen(fp, tr, p->next);
  
 
  if ((ch = treeGetCh(fp)) == ',') 
    addElementLen(fp, tr, p->next->next);
  else
    assert(0);
  
  if (! treeNeedCh(fp, ')', "in"))                
    assert(0);
  

  treeFlushLabel(fp);
  
  if (! treeFlushLen(fp))                         
    assert(0);
 
  if (! treeNeedCh(fp, ';', "at end of"))       
    assert(0);
  
    
}










extern FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  
	  printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  exit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	 
	  printf("The file %s RAxML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
	  exit(-1);
	  return (FILE *)NULL;
	}
    }


}


void printBothOpen(const char* format, ... )
{
  //FILE *f = myfopen(infoFileName, "ab");
  FILE *f = myfopen(partitionInfoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}















double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}




static int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}










static void getnums (rawdata *rdta)
{
  if (fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2)
    {
      
      printf("ERROR: Problem reading number of species and sites\n");
      exit(-1);
    }

  if (rdta->numsp < 4)
    {
      
      printf("TOO FEW SPECIES\n");
      exit(-1);
    }

  if (rdta->sites < 1)
    {      
      printf("TOO FEW SITES\n");
      exit(-1);
    }

  return;
}





static boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


static void uppercase (int *chptr)
{
  int  ch;

  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
}



// Comments are by Mark Jones trying to figure out what is happening
static void getyspace (rawdata *rdta)
{
  size_t size = 4 * ((size_t)(rdta->sites / 4 + 1));	// memory size of a 'row' in the 2d array . somehow 4 things are enough to represent 4 elements I guess?
  int    i;
  unsigned char *y0;

  rdta->y = (unsigned char **) malloc((size_t)(rdta->numsp + 1) * sizeof(unsigned char *));   // create space for the top-level array
  assert(rdta->y);   //we found space right?

  y0 = (unsigned char *) malloc(((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char)); //this a pointer address, after which is the requested amount of free space...
  assert(y0);   

  // the location of y and its array may be different from the location of the y_0 and the associated actual 2d array.
  // each element of y points to the corresponding point in the 2d array space corresponding to its row.
  // this isn't going to make sense to anyone but me I guess but that's ok, I'm going to delete this  all before pushing.
  // SO... what I'm going to want to do is something like this, but for 3 dimensions???

  rdta->y0 = y0;  // pointer of first element of 2d  array

  for (i = 0; i <= rdta->numsp; i++)
    {
      rdta->y[i] = y0;		// row i of rdta is assigned to pointer (original y_0) + i*(row size)
      y0 += size;
    }

  return;
}





// Returns a 2d int array of dimensions [d1][d2], and allocates memory for it
static int **get2dIntArray(int d1, int d2)
  {
    // "Data array"
    // Reserve a block of memory large enough to store d1*d2 ints
    // Set y0 to be the start position of this block
    // Once everything is set up, the int corresponding to arr[i][j] will live at position  y0 + (i*d2 + j*) * sizeof(int)
    int *y0 = (int*) malloc(sizeof(int) * d1 * d2);
    assert(y0);

    // "Top level array"
    // Reserve a block of memory large enough to store d1 int-pointers
    // Set x0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i] will live at position x0 + i*sizeof(int*).
    // It will point to  y0 + i*d2 *  sizeof(int),  i.e. the location of the int corresponding to arr[i][0]
    // (It follows that arr[i][j] will return the int at y0 + i*d2*sizeof(int*) + j*sizeof(int*)
    //   = y0 + (i*d2 + j) * sizeof(int), as required.)
    int **x0 = (int**) malloc(sizeof(int*) * d1);
    assert(x0);


    // Now we have allocated memory space, we can set up our pointers.
    // Define arr as a pointer to x0 
    // (so that arr[0] will be the pointer living at position x0 and arr[i] will be the pointer living at position 
    //    x0 + i*sizeof(int*)  )
    int **arr = x0;

    // Define the pointer living at arr[i] to be a pointer to y0 + i*d2*sizeof(int)
    // (so that arr[i][j] will be the int living at position y0 + (i*d2 + j)*sizeof(int) )
    int *y = y0;
    int i;
    for (i = 0; i < d1; i++)
      {
        arr[i] = y;
        y += d2;
      } 

    // All the pointers are set up, so return the pointer to the start of the array.
    return arr;

  }

// Returns a 3d float array of dimensions [d1][d2][d3], and allocates memory for it
static float ***get3dFloatArray(int d1, int d2, int d3)
  {
    // "Data array"
    // Reserve a block of memory large enough to store d1*d2*d3 floats
    // Set z0 to be the start position of this block
    // Once everything is set up, the float corresponding to arr[i][j][k] will live at position  z0 + (i*(d2*d3) + j*d3 + k) * sizeof(float)
    float *z0 = (float*) malloc(sizeof(float) * d1 * d2 * d3);
    assert(z0);

    // "Intermediate array"
    // Reserve a block of memory large enough to store d1*d2 float-pointers
    // Set y0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i][j] will live at position y0 + (i*d2 + j) * sizeof(float*).
    // It will point to z0 + (i*(d2*d3) + j*d3)  * sizeof(float),  i.e. the location of the float corresponding to arr[i][j][0]
    // (It follows that arr[i][j][k]  will return the float at z0 + (i*(d2*d3) + j*d3)*sizeof(float) + k*sizeof(float)  
    //   =  z0 + (i*(d2*d3) + j*d3 + k) * sizeof(float), as required.)
    float **y0 = (float**) malloc(sizeof(float*) * d1 * d2);
    assert(y0);

    // "Top level array"
    // Reserve a block of memory large enough to store d1 float-pointer-pointers
    // Set x0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i] will live at position x0 + i*sizeof(float**).
    // It will point to  y0 + i*d2 *  sizeof(float*),  i.e. the location of the pointer corresponding to arr[i][0]
    // (It follows that arr[i][j] will return the float-pointer at y0 + i*d2*sizeof(float*) + j*sizeof(float*)
    //   = y0 + (i*d2 + j) * sizeof(float*), as required.)
    float ***x0 = (float***) malloc(sizeof(float**) * d1);
    assert(x0);


    // Now we have allocated memory space, we can set up our pointers.

    // Define arr as a pointer to x0 
    // (so that arr[0] will be the pointer living at position x0 and arr[i] will be the pointer living at position 
    //    x0 + i*sizeof(float**)  )
    float ***arr = x0;

    // Define the pointer living at arr[i] to be a pointer to y0 + i*d2*sizeof(float*)
    // (so that arr[i][j] will be the pointer living at position y0 + (i*d2 + j)*sizeof(float*) )
    float **y = y0;
    int i;
    for (i = 0; i < d1; i++)
      {
        arr[i] = y;
        y += d2;
      } 

    float *z = z0;
    int j;
    // Define the pointer living at arr[i][j] to be a pointer to z0 + (i*(d2*d3) + j*d3) * sizeof(float)
    // (so that arr[i][j][k] will be the float living at position z0 + (i*(d2*d3) + j*d3 + k) * sizeof(float) )
    for (i = 0; i < d1; i++)
    {
      for (j = 0; j < d2; j++)
        {
           arr[i][j] = z;
           z += d3;
        } 
    }

    // All the pointers are set up, so return the pointer to the start of the array.
    return arr;
  }


// Returns a 3d int array of dimensions [d1][d2][d3], and allocates memory for it
static int ***get3dIntArray(int d1, int d2, int d3)
  {
    // "Data array"
    // Reserve a block of memory large enough to store d1*d2*d3 ints
    // Set z0 to be the start position of this block
    // Once everything is set up, the int corresponding to arr[i][j][k] will live at position  z0 + (i*(d2*d3) + j*d3 + k) * sizeof(int)
    int *z0 = (int*) malloc(sizeof(int) * d1 * d2 * d3);
    assert(z0);

    // "Intermediate array"
    // Reserve a block of memory large enough to store d1*d2 int-pointers
    // Set y0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i][j] will live at position y0 + (i*d2 + j) * sizeof(int*).
    // It will point to z0 + (i*(d2*d3) + j*d3)  * sizeof(int),  i.e. the location of the int corresponding to arr[i][j][0]
    // (It follows that arr[i][j][k]  will return the int at z0 + (i*(d2*d3) + j*d3)*sizeof(int) + k*sizeof(int)  
    //   =  z0 + (i*(d2*d3) + j*d3 + k) * sizeof(int), as required.)
    int **y0 = (int**) malloc(sizeof(int*) * d1 * d2);
    assert(y0);

    // "Top level array"
    // Reserve a block of memory large enough to store d1 int-pointer-pointers
    // Set x0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i] will live at position x0 + i*sizeof(int**).
    // It will point to  y0 + i*d2 *  sizeof(int*),  i.e. the location of the pointer corresponding to arr[i][0]
    // (It follows that arr[i][j] will return the int-pointer at y0 + i*d2*sizeof(int*) + j*sizeof(int*)
    //   = y0 + (i*d2 + j) * sizeof(int*), as required.)
    int ***x0 = (int***) malloc(sizeof(int**) * d1);
    assert(x0);


    // Now we have allocated memory space, we can set up our pointers.

    // Define arr as a pointer to x0 
    // (so that arr[0] will be the pointer living at position x0 and arr[i] will be the pointer living at position 
    //    x0 + i*sizeof(int**)  )
    int ***arr = x0;

    // Define the pointer living at arr[i] to be a pointer to y0 + i*d2*sizeof(int*)
    // (so that arr[i][j] will be the pointer living at position y0 + (i*d2 + j)*sizeof(int*) )
    int **y = y0;
    int i;
    for (i = 0; i < d1; i++)
      {
        arr[i] = y;
        y += d2;
      } 

    int *z = z0;
    int j;
    // Define the pointer living at arr[i][j] to be a pointer to z0 + (i*(d2*d3) + j*d3) * sizeof(int)
    // (so that arr[i][j][k] will be the int living at position z0 + (i*(d2*d3) + j*d3 + k) * sizeof(int) )
    for (i = 0; i < d1; i++)
    {
      for (j = 0; j < d2; j++)
        {
           arr[i][j] = z;
           z += d3;
        } 
    }

    // All the pointers are set up, so return the pointer to the start of the array.
    return arr;
  }



// Returns a 3d char array of dimensions [d1][d2][d3], (i.e; a 2D array of length-d3 strings) and allocates memory for it
static char ***get3dCharArray(int d1, int d2, int d3)
  {
    // "Data array"
    // Reserve a block of memory large enough to store d1*d2*d3 chars
    // Set z0 to be the start position of this block
    // Once everything is set up, the char corresponding to arr[i][j][k] will live at position  z0 + (i*(d2*d3) + j*d3 + k) * sizeof(char)
    char *z0 = (char*) malloc(sizeof(char) * d1 * d2 * d3);
    assert(z0);

    // "Intermediate array"
    // Reserve a block of memory large enough to store d1*d2 char-pointers
    // Set y0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i][j] will live at position y0 + (i*d2 + j) * sizeof(char*).
    // It will point to z0 + (i*(d2*d3) + j*d3)  * sizeof(int),  i.e. the location of the char corresponding to arr[i][j][0]
    // (It follows that arr[i][j][k]  will return the int at z0 + (i*(d2*d3) + j*d3)*sizeof(char) + k*sizeof(char)  
    //   =  z0 + (i*(d2*d3) + j*d3 + k) * sizeof(int), as required.)
    char **y0 = (char**) malloc(sizeof(char*) * d1 * d2);
    assert(y0);

    // "Top level array"
    // Reserve a block of memory large enough to store d1 char-pointer-pointers
    // Set x0 to be the start position of this block
    // Once everything is set up, the pointer corresponding to arr[i] will live at position x0 + i*sizeof(char**).
    // It will point to  y0 + i*d2 *  sizeof(char*),  i.e. the location of the pointer corresponding to arr[i][0]
    // (It follows that arr[i][j] will return the char-pointer at y0 + i*d2*sizeof(char*) + j*sizeof(char*)
    //   = y0 + (i*d2 + j) * sizeof(char*), as required.)
    char ***x0 = (char***) malloc(sizeof(char**) * d1);
    assert(x0);


    // Now we have allocated memory space, we can set up our pointers.

    // Define arr as a pointer to x0 
    // (so that arr[0] will be the pointer living at position x0 and arr[i] will be the pointer living at position 
    //    x0 + i*sizeof(char**)  )
    char ***arr = x0;

    // Define the pointer living at arr[i] to be a pointer to y0 + i*d2*sizeof(char*)
    // (so that arr[i][j] will be the pointer living at position y0 + (i*d2 + j)*sizeof(char*) )
    char **y = y0;
    int i;
    for (i = 0; i < d1; i++)
      {
        arr[i] = y;
        y += d2;
      } 

    char *z = z0;
    int j;
    // Define the pointer living at arr[i][j] to be a pointer to z0 + (i*(d2*d3) + j*d3) * sizeof(char)
    // (so that arr[i][j][k] will be the char living at position z0 + (i*(d2*d3) + j*d3 + k) * sizeof(char) )
    for (i = 0; i < d1; i++)
    {
      for (j = 0; j < d2; j++)
        {
           arr[i][j] = z;
           z += d3;
        } 
    }

    // All the pointers are set up, so return the pointer to the start of the array.
    return arr;
  }





static boolean setupTree (tree *tr)
{
  nodeptr  p0, p, q;
  int
    i,
    j,   
    tips,
    inter; 

 

  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  
  tr->yVector      = (unsigned char **)  malloc((size_t)(tr->mxtips + 1) * sizeof(unsigned char *));  

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 

  /*TODO, must that be so long ?*/

 

  tr->nameList = (char **)malloc(sizeof(char *) * (size_t)(tips + 1));    

  if (!(p0 = (nodeptr) malloc((size_t)(tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (!(tr->nodep = (nodeptr *) malloc((size_t)(2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

    
      p->x      =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;
     

      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    p->x = 1;
	  else
	    p->x =  0;
	  p->number = i;
	  p->next   = q;
	 
	  p->back   = (node *) NULL;
	 


	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  
  tr->start       = (node *) NULL;

  

  tr->ntips       = 0;
  tr->nextnode    = 0;

  

  return TRUE;
}


static void checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
    {
      boolean valid;

      switch(buffer[i])
	{
	case '\0':
	case '\t':
	case '\n':
	case '\r':
	case ' ':
	case ':':
	case ',':
	case '(':
	case ')':
	case ';':
	case '[':
	case ']': 
	case '\'':
	  valid = FALSE;
	  break;
	default:
	  valid = TRUE;
	}

      if(!valid)
	{
	  printf("ERROR: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n", buffer, i, buffer[i]);
	  printf("Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\", \"\'\"\n");
	  printf("Exiting\n");
	  exit(-1);
	}

    }
  assert(buffer[len - 1] == '\0');
}



static boolean getdata(rawdata *rdta, tree *tr)
{
  size_t
    i, 
    j, 
    basesread, 
    basesnew;

  int
    ch, my_i, meaning,
    len,
    meaningDNA[256];
  
  boolean  
    allread, 
    firstpass;
  
  char 
    buffer[nmlngth + 2];
  
  for(i = 0; i < 256; i++)     
    meaningDNA[i] = -1;

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;  
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9; 
  meaningDNA['Y'] = 10;

  meaningDNA['N'] = 
    meaningDNA['O'] = 
    meaningDNA['X'] = 
    meaningDNA['-'] = 
    meaningDNA['?'] = 
    15;

 

  /*******************************************************************/

  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread)
    {
      for (i = 1; i <= (size_t)tr->mxtips; i++)
	{
	  if (firstpass)
	    {
	      ch = getc(INFILE);
	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);

	      my_i = 0;

	      do
		{
		  buffer[my_i] = ch;
		  ch = getc(INFILE);
		  my_i++;
		  if(my_i >= nmlngth)
		    {		     
		      printf("Taxon Name to long at taxon %zd, adapt constant nmlngth in\n", i);
		      printf("axml.h, current setting %d\n", nmlngth);
		      exit(0);
		    }
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');

	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);
	      
	      ungetc(ch, INFILE);

	      buffer[my_i] = '\0';
	      len = strlen(buffer) + 1;
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * (size_t)len);
	      strcpy(tr->nameList[i], buffer);
	    }

	  j = basesread;

	  while ((j < (size_t)rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
	    {
	      uppercase(& ch);
	      meaning = meaningDNA[ch];	     

	      if (meaning != -1)
		{
		  j++;
		  rdta->y[i][j - 1] = meaning;	
		}
	      else
		{
		  if(!whitechar(ch))
		    {
		      printf("ERROR: Bad base (%c) at site %zd of sequence %zd\n",
			     ch, j + 1, i);
		      return FALSE;
		    }
		}
	    }
	    


	  if (ch == EOF)
	    {
	      printf("ERROR: End-of-file at site %zd of sequence %zd\n", j + 1, i);
	      return  FALSE;
	    }

	  if (! firstpass && (j == basesread))
	    i--;
	  else
	    {
	      if (i == 1)
		basesnew = j;
	      else
		if (j != basesnew)
		  {
		    printf("ERROR: Sequences out of alignment\n");
		    printf("%zd (instead of %zd) residues read in sequence %zd %s\n",
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= (size_t)rdta->sites);
    }

 


  return  TRUE;
}











static void getinput(analdef *adef, rawdata *rdta, tree *tr) //, cruncheddata *cdta //not needed anymore
{
  int i;

  
  INFILE = myfopen(seq_file, "rb");
  
  getnums(rdta);

  tr->mxtips            = rdta->numsp;
  
          

  getyspace(rdta);
    

  setupTree(tr);

  
  if(!getdata(rdta, tr))
    {
      printf("Problem reading alignment file \n");
      exit(1);
    }
 // cdta->endsite = rdta->sites; // not needed anymore

  if(adef->restart)
    {
      tr->nameHash = initStringHashTable(10 * tr->mxtips);
      for(i = 1; i <= tr->mxtips; i++)
	addword(tr->nameList[i], tr->nameHash, i);
    }

  fclose(INFILE);
}


















static boolean makevalues(rawdata *rdta, tree *tr) //, cruncheddata *cdta not needed anymore
{
  int  
    i;

  tr->rdta       = rdta;
  //tr->cdta       = cdta; //not needed anymore
  
  for(i = 0; i <= rdta->numsp; i++)
    tr->yVector[i] = rdta->y[i];

  return TRUE;
}














static void initAdef(analdef *adef)
{     
  adef->restart                = FALSE;
  adef->parsimonySeed          = 0;  
  adef->numberOfTrees          = 1;
}






static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }





static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("Error character %c not allowed in run ID\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("Error: please provide a string for the run id after \"-n\" \n");
      assert(0);
    }

}

static void get_args(int argc, char *argv[], analdef *adef)
{
 
  int  
    optind = 1,  
    c;
  
  boolean
    bad_opt    =FALSE;

  char *optarg;

  run_id[0] = 0; 
  seq_file[0] = 0;
  strcpy(partitionInfoFileName,        "output.");
  
  adef->parsimonySeed=1;
  adef->numberOfBlocks=2;
  adef->maxHomoplasyRatio= FLT_MAX;
  adef->verbose= 1;
  
  while(!bad_opt && ((c = mygetopt(argc,argv,"p:n:s:o:t:N:b:B:r:v:", &optind, &optarg))!=-1))
    {      
      switch(c)
	{     
	case 'N':
	  sscanf(optarg,"%d", &(adef->numberOfTrees));	
	  if(adef->numberOfTrees <= 0)
	    {
	      printf("number of trees can't be smaller than 1\n");
	      exit(-1);
	    }
	  break;
	case 'b':       
	  sscanf(optarg,"%d", &(adef->numberOfBlocks));	
	  break;
	case 't':       
	  strcpy(tree_file, optarg);
	  adef->restart = TRUE;
	  break;  
	case 's':		  
	  strcpy(seq_file, optarg);      
	  break;
	case 'o':		  
	  strcpy(partitionInfoFileName, optarg);      
	  break;	
	case 'p':
	  sscanf(optarg,"%ld", &(adef->parsimonySeed));	
	  if(adef->parsimonySeed <= 0)
	    {
	      printf("Parsimony seed specified via -p must be greater than zero\n");
	      exit(-1);
	    }
	  break;      
	case 'n':
	  strcpy(run_id, optarg);
	  analyzeRunId(run_id);
	  break;    
	case 'B':
	  strcpy(trueBlockPartitionString, optarg);
	  break;        
	case 'r':
	  sscanf(optarg,"%f", &(adef->maxHomoplasyRatio));	// desired maximum homoplasy ratio - we're interested in shortest block partition that achieves this ratio
	  break;         
	case 'v':
	  sscanf(optarg,"%d", &(adef->verbose));	// desired verbosity - 0: just return summary information. 1: return block partition (and max homoplasy ratio) for each possible block count. 2: return homoplasy ratio and parsimonious tree for each block in each partition. (set to 1 by default)
	  break;         
	default:
	  assert(0);
	}
    }
  
  return;
}







static void makeFileNames(void)
{ 
  strcpy(resultFileName,       "RAxML_parsimonyTree.");  
  strcpy(infoFileName,         "RAxML_info.");
  //strcpy(infoFileName,         userInputString);
   
  strcat(resultFileName,       run_id); 
  strcat(infoFileName,         run_id);

  if(filexists(infoFileName))
    {
      printf("RAxML output files with the run ID <%s> already exist \n", run_id);     
      
      exit(-1);
    }
}




 










/************************************************************************************/




static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

static char bits_in_16bits [0x1u << 16];

static void compute_bits_in_16bits(void)
{
    unsigned int i;    
    
    assert(sizeof(unsigned int) == 4);

    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);
    
    return ;
}

unsigned int precomputed16_bitcount (unsigned int n)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}



 /************************/

static float getHomoplasyRatio(int i, int j, int homoplasyValue) 
{
  assert  (i <= j);  // program breaks if called with i greater than j
  return ((float) homoplasyValue) / ((float)(j-i+1));

}

// Given a number of blocks and two arrays describing the start and end points of each block, 
// and a lookup table of homoplasy scores for each possible block,
// returns the maximum homoplasy ratio of the given set of blocks
static float getMaxHomoplasyRatioFromPartition(int blockCount, int* blockStartPoints, int* blockEndPoints, int** blockScoreArray)
  {
    // okay, first step: look up how i did similarish stuff back when i was doing the dyamic programming.
    // actually no. First, set up the bits you know  about
    float maxRatioSoFar = 0;
    float tempHomoplasyRatio;
    int tempBlockStart;
    int tempBlockEnd;
    int tempHompolasyScore;
    int b;
    for (b = 0; b < blockCount; b++)
      {      
        tempBlockStart = blockStartPoints[b];
        tempBlockEnd = blockEndPoints[b];

        tempHompolasyScore = blockScoreArray[tempBlockStart][tempBlockEnd];
        tempHomoplasyRatio = getHomoplasyRatio(tempBlockStart,tempBlockEnd,tempHompolasyScore);

        if (tempHomoplasyRatio > maxRatioSoFar)
          {
            maxRatioSoFar = tempHomoplasyRatio;
          }
      }
    return maxRatioSoFar;

  }

// Converts a block count and two arrays into a string describing a block partition
static char* getBlockPartitionString(int blockCount, int* blockStartPoints, int* blockEndPoints)
{

  char* outstr = (char*) malloc(sizeof(char) * 1000);
  strcpy(outstr, "");
  char tempstr[1000];
  int tempBlockStart;
  int tempBlockEnd;
  for (int h = 0; h < blockCount; h++)
    {   

      tempBlockStart = blockStartPoints[h];
      tempBlockEnd = blockEndPoints[h];
      sprintf(tempstr, "[%d---%d]", tempBlockStart, tempBlockEnd); // make formatted string discribing block h
      strcat(outstr, tempstr);  // append block to output string
      if (h < blockCount - 1)
        {
          strcat(outstr, ", "); // append a comma if there are more blocks to come
        } 
      //else 
      //  {
      //    strcat(outstr, "\n"); // otherwise, end line
      //  }
    }
  return outstr;
}


// Returns the average internal boundary error for a block partition with blockCount blocks and endpoints defined by testBlockStarts, testBlockEnds, when compared to the block partition defined by truePartitionStartPoints and truePartitionEndPoints
// For two block partitions [x_1 -- y_1],...,[x_b---y_b] and [w_1---z_1],...,[w_b---z_b],
// And assuming that x_{i+1} = y_i +1, w_{i+1} = z_i + 1 for all  1 <= i < b,
// the average internal boundary error is defined as 
//   ( abs(y_1-z_1) + abs(y_2-z_2) + ... + abs(y_{b-1}-z_{b-1}) ) / (b-1)
// ERMMM
// So hey I guess it turned out we didn't need the startpoint arrays here after all but whatever.
// TODO go back and put in a proper check that the arrays are well-formated, later.
static float averageInternalBoundaryError(int blockCount, int* truePartitionStartPoints, int* truePartitionEndPoints, int* testBlockStarts, int* testBlockEnds)
{
  int totalBoundaryError = 0;
  int i;
  for (i = 0; i < blockCount - 1; i++)  // We don't test the endpoints for the final block, as these are assumed to both be taxaCount-1
    {
      int localBoundaryError = abs(truePartitionEndPoints[i] - testBlockEnds[i]);
      // printf("local error: %d\n", localBoundaryError);
      totalBoundaryError += localBoundaryError;
      // printf("total error so far: %d\n", totalBoundaryError);
    }
  // printf("Total error %d / internal boundary count %d\n", totalBoundaryError, blockCount -1);
  float averageError = (float) totalBoundaryError / (blockCount -1);
  // printf("Average error %f\n", averageError);
  return averageError;
}


static int getOptHomoplasyOfContiguousBlock(tree *tr, analdef *adef, int startCharacter, int endCharacter)
{



  //printf("Calculating [%d -- %d ]:", startCharacter, endCharacter);
  assert  (startCharacter <= endCharacter); 

  int parsimonyScore = makeParsimonyTreeFastDNA(tr, adef, startCharacter, endCharacter);
  
  // Calculate the lower bound on parsimony, based on how many states are needed to cover each site.
  int parsimonyScoreLowerBound = 0;  
  int i;
  for (i = startCharacter; i <= endCharacter; i++)
    {
      // Parsimonator treats characters in which only one state appears more than once as 'uninformative', and so ignores them.
      // Uninformative characters can still have a non-0 minimum parsimony score, and so can lead to an innacurate homoplasy score
      // (e.g. if site i is uninformative with 3 states (one A, one C, and the rest T's), Parsimonator returns a score of 0 for site i despite the fact that the minimum parsimony score is technically 2. Thus we would get a homoplasy score of -2.)
      // To avoid this problem, we too must ignore uninformative sites when calculating the parsimony score lower bound.
      if (isInformative(tr, i))
        {
          parsimonyScoreLowerBound = parsimonyScoreLowerBound + tr->minParsimonyPerSite[i];
        }
    }

  int returnVal = parsimonyScore - parsimonyScoreLowerBound;
  //printf("%d - %d = %d\n", parsimonyScore, parsimonyScoreLowerBound, returnVal);
  return  returnVal;
  //return 0;
}



// Given a string trueBlockPartitionString of the form [0---x],[x+1---y],.... ,[w+1---z] representing a block partition
// Sets truePartitionBlockCount to be the number of blocks in the partition,
// And sets truePartitionStartPoints (respectively, truePartitionEndPoint) to be an array of length truePartitionBlockCount containing 
// the start points (respectively, end points) of all blocks in the partition.
// TODO: spell out the format required for trueBlockPartitionString: no spaces, blocks separated by commas, blocks have the form [x---y] where x and y are integers. TODO other constraints?
static void getTrueBlockPartitionData( char* trueBlockPartitionString, int *truePartitionBlockCount, int **truePartitionStartPoints, int **truePartitionEndPoints, int sites)
{
  // First calculate number of blocks in true partition

  // Go through thr string and see how many appearances of ',' there are. For each ',' we have an additional block
  int tempBlockCount = 1;
  int i = 0;
  while (trueBlockPartitionString[i] != '\0')
    {
      if (trueBlockPartitionString[i] == ',')
        {
          tempBlockCount++;
        }
      i++;
    }
  *truePartitionBlockCount = tempBlockCount;

  // make tempStartArray and tempEndArray be arrays of length truePartitionBlockCount
  int* tempStartArray = (int*) malloc(sizeof(int) * *truePartitionBlockCount);
  int* tempEndArray = (int*) malloc(sizeof(int) * *truePartitionBlockCount);

  // now set truePartitionStartPoints and truePartitionEndPoints to be the two arrays we just constructed.
  free(*truePartitionStartPoints);
  *truePartitionStartPoints = tempStartArray;
  free(*truePartitionEndPoints);
  *truePartitionEndPoints = tempEndArray;

  // populate tempStartArray and tempEndArray with integers from the input string

  int nextStartInt = -1;
  int nextEndInt = -1;
  char* currentString = trueBlockPartitionString;
  char* nextString = trueBlockPartitionString;

  int b = 0;

  // populate the first item of each array
  if (sscanf(currentString, "[%d%s", &nextStartInt, nextString) > 0)
    {
      tempStartArray[b] = nextStartInt;
      currentString = nextString;
      if (sscanf(currentString, "---%d]%s", &nextEndInt, nextString) > 0)
        {
          tempEndArray[b] = nextEndInt;
          currentString = nextString;
        }
      else
        {
          nextEndInt = -1;
          tempEndArray[b] = nextEndInt;
          currentString = "";
          printf("Error1: badly formatted block input (blocks must be of the form [i---j] for integers i <= j, separated by ',' without spaces)\n");
          exit(-1);
        }
      b++;
    }
  // populate the subsequent items of each array
  while (sscanf(currentString, ",[%d%s)", &nextStartInt, nextString) > 0)
    {
      tempStartArray[b] = nextStartInt;
      currentString = nextString;
      if (sscanf(currentString, "---%d]%s", &nextEndInt, nextString) > 0)
        {
          tempEndArray[b] = nextEndInt;
          currentString = nextString;
        }
      else
        {
          nextEndInt = -1;
          tempEndArray[b] = nextEndInt;
          currentString = "";
          printf("Error2: badly formatted block input (blocks must be of the form [i---j] for integers i <= j, separated by ',' without spaces)\n");
          exit(-1);
        }
      b++;
    }
  //// commenting out - intended to throw an error if there is extra garbage after a well-formatted block partition string. 
  //// Was not working, threw errors for correctly formatted instances.
  //// Currently no test for this type of formatting error
  //if(currentString[0] != '\0')
  //  {
  //    printf("Error3: badly formatted block input %s (blocks must be of the form [i---j] for integers i <= j, separated by ',' without spaces)\n", currentString);
  //    exit(-1);
  //  }


  // check that the block partition is well formatted
  if (tempStartArray[0] != 0)
    {
      printf("Error4: badly formatted block input (first block must begin with 0)\n");
      exit(-1);
    }
  // check that the block partition is well formatted
  if (tempEndArray[*truePartitionBlockCount -1] != sites -1)
    {
      printf("Error5: badly formatted block input (last block [%d---%d] must end with %d (sequence length -1))\n", tempStartArray[*truePartitionBlockCount -1], tempEndArray[*truePartitionBlockCount -1], sites-1);
      exit(-1);
    }
  // TODO:check that last endpoint is equal to number of species - 1
  //if (tempEndArray[*truePartitionBlockCount -1] )

  for (b = 0; b < *truePartitionBlockCount; b++)
    {
       if (tempStartArray[b] > tempEndArray[b])
         {
           printf("Error6: badly formatted block input [%d---%d](block start cannot be greater than block end)\n", tempStartArray[b], tempEndArray[b]);
           exit(-1);
         }
    }
  for (b = 0; b < *truePartitionBlockCount -1; b++)
    {
       if (tempStartArray[b+1] != tempEndArray[b] + 1)
         {
           printf("Error7: badly formatted block input [%d---%d],[%d---%d] (block start must equal previous block end + 1)\n", tempStartArray[b], tempEndArray[b], tempStartArray[b+1], tempEndArray[b+1]);
           exit(-1);
         }
    }


  //printf ("ANOTHER TEST [%d --- %d]\n", nextStartInt, nextEndInt);
//  tempStartArray[0] = nextStartInt;
//  tempEndArray[0] = nextEndInt;
  //printf ("ANOTHER TEST [%d --- %d]\n", tempStartArray[0], tempEndArray[0]);
  //printf ("ANOTHER TEST [%d --- %d]\n", *truePartitionStartPoints[0], *truePartitionEndPoints[0]);

}

static void getOptBlockPartition(tree *tr, analdef *adef) 
{

  int debugOutput = 0;

  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;
  float desiredHomoplasyRatio = adef->maxHomoplasyRatio;
  int verbosity = adef->verbose;// controls how much information is printed to terminal. 
                                // 0: Just return summary info
                                // 1: Also return opt block partition + homoplasy ratio for each possible block count
                                // 2: Also return homoplasy ratio + parsimonious tree for each block in each partition

  int truePartitionBlockCount;
  int* truePartitionStartPoints;
  int* truePartitionEndPoints; 

  // parse and store data on the True Block Partition
  // NOTE / TODO: putting this construction here was the simplest thing to do, but not necessarily the most hygenic.
  // We might want to move this into get_args or similar.
  // Maybe we want to store the information in adef ?
  if (trueBlockPartitionString[0] != '\0')
    {
      char tempBlockPartitionString[1024] = ""; 
      strcpy(tempBlockPartitionString,trueBlockPartitionString);
      getTrueBlockPartitionData(tempBlockPartitionString, &truePartitionBlockCount, &truePartitionStartPoints, &truePartitionEndPoints, sites);
    }


  
  printBothOpen("---------------------------------\n"),
  printBothOpen("Input file: %s | Output file: %s | Maximum number of blocks: %d | Randomization seed: %d\n", seq_file, partitionInfoFileName, adef->numberOfBlocks, adef->parsimonySeed);  
  printBothOpen("---------------------------------\n");

  //int **BP_blockScoreArray  = get2dIntArray(sites, sites);

  int BP_minParsimonyScorePerSiteArray [sites];  // BP_minParsimonyScorePerSiteArray[i] stores the theoretical minimum parsimony score of  any character on site i.  (i.e. the minimum number of states used in site i, minus 1)

  // BP_blockScoreArray[i][j]  stores the optimum homoplasy score of the block on characters i to j (inclusive)
  int **BP_blockScoreArray = get2dIntArray(sites, sites);


  // BP_blockTreeStringArray[i][j] is the newick string corresponding to an optimum tree for the block on characters i to j (inclusive)
  char ***BP_blockTreeStringArray = get3dCharArray(sites, sites, tr->treeStringLength);

  int impossiblyHighHomoplasyScore = numsp*sites;

// Entry BP_DPHomoplasyRatioLookupTable[b][i][j] for i<= j is the minimum value of 
// the homoplasy ratio of any block partition P of character 0 to j, 
// for which the last block is from i to j and there are b+1 blocks in total;
// where "homoplasy ratio of P" means the maximum homoplasy ratio of 
// any block in P. A value of Float.POSITIVE_INFINITY means there is no solution 
// for this set of values (e.g. if i > j)
   float ***BP_DPHomoplasyRatioLookupTable = get3dFloatArray(desiredBlockCount, sites, sites);

   //float BP_DPHomoplasyRatioLookupTable[desiredBlockCount][sites][sites];

// Used to reconstruct an optimal solution.
// BP_backtrackTable[b][i][j]  is the value i' such that an 
// optimal block partition with b+1 blocks and last block [i,j]
// has [i',j-1] as its second-to-last block.
   int ***BP_backtrackTable = get3dIntArray(desiredBlockCount, sites, sites);

   //int BP_backtrackTable[desiredBlockCount][sites][sites];




  //printf("populating BP_minParsimonyScorePerSiteArray\n");
  // Populate BP_minParsimonyScorePerSiteArray, which stores the theoretical minimum parsimony score for each individal site
  for (int i = 0; i < sites; i++)
    {
      BP_minParsimonyScorePerSiteArray[i] = getMinParsimonyScoreForSite(tr, i);
    }
  // Add min Parsimony data to tr, so that we can look it up easily later on.
  tr->minParsimonyPerSite = BP_minParsimonyScorePerSiteArray;


  //printf("populating BP_blockScoreArray\n");
  // Populate BP_blockScoreArray, which stores the optimum homoplasy score every contiguous block of characters
  // Also populate BP_blockTreeStringArray, which stores an optimum treee for every contiguous block of characters

  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++) 
        {
          if (i > j)
            BP_blockScoreArray[i][j] = INT_MAX;
          else
            BP_blockScoreArray[i][j] = getOptHomoplasyOfContiguousBlock(tr, adef, i, j);
            // Also populate BP_blockTreeStringArray with the tree we just found

//            BP_blockTreeStringArray[i][j] = tr->tree_string;
            // populate BP_blockTreeStringArray one character at a time (not 100% sure why this works and the commented-out line above does not, probably something to do with pointers)
            for (int k = 0; k < tr->treeStringLength; k++)
            {
              BP_blockTreeStringArray[i][j][k] = tr->tree_string[k];
            }
        }  
      if (debugOutput > 0)
        { 
          printf("Calculated parsimony score for blocks [%d -- j] for all j\n",i);
        }

    } 

  //printf("populating lookup table.\n");

  for (int i = 0; i < sites; i++)   
    {

      for (int j = 0; j < sites; j++)
        {
          for (int b = 0; b < desiredBlockCount; b++)
            {
              //printf("%f\n   ",BP_DPHomoplasyRatioLookupTable[0][0][0]);
              float val = FLT_MAX;
              int bestPrevI = -1;

              if (b == 0 && i > 0)
                {
                  val = FLT_MAX;		// if first block, i should be 0
                }
              if (b > 0 && i == 0)
                {
                  val = FLT_MAX;		// if not first block, i should be bigger than 0
                }

              if (b == 0 && i == 0)			// if first block, take homoplasy ratio of that block.
                {
                   val =  getHomoplasyRatio(i,j, BP_blockScoreArray[i][j]);	
                   //printf("%d %d %d = %f\n   ",b, i,j,val);
                }
              if (b > 0 && i > 0)
                {
                  if (i > j)			// i should be at most j.
                    {
                      val = FLT_MAX;	
                    } 
                  else 
                    {

                      //find the minimum possible homoplasy score for any b-1 block partition up to this point
                      float currentMinimum = FLT_MAX;
                      int lastj = i-1;
                      for(int lasti = 0; lasti <= lastj; lasti++)
                        {
                          float tempVal = BP_DPHomoplasyRatioLookupTable[b-1][lasti][lastj];
 	
                           if (tempVal < currentMinimum )
                            {
                              currentMinimum = tempVal;
                              // update the bestPrevI value because we have a new best block partition
                              bestPrevI = lasti;
                            }
                        }
                      // currentMinimum is now the optimum homoplasy ratio for any block partition of 0 to i-1
                      // using b-1 blocks
                      float valCurrentBlock =  getHomoplasyRatio(i,j, BP_blockScoreArray[i][j]);
                      // If the current block has worse homoplasy ratio than what we could 
                      // get up to this point, that ratio is now the best we can do for a block partition of this type.
                      if (currentMinimum < valCurrentBlock )
                        {
                          val = valCurrentBlock;
                        } 
                      // Otherwise, the best we could do up to this point is still the best we can do
                      else 
                        {
                          val = currentMinimum;
                        }
                    }
                }

              // populate the lookup table with the calculated value
              BP_DPHomoplasyRatioLookupTable[b][i][j] = val;

              // also populate the backtrack table so we can reconstruct a solution later
              BP_backtrackTable[b][i][j] = bestPrevI;

              if (debugOutput > 0)
                {
                  printf("BP_DPHomoplasyRatioLookupTable[%d]",b);
                  printf("[%d]",i);
                  printf("[%d]",j);
                  printf("= %f   ",val);
                  printf("BP_backtrackTable[%d]",b);
                  printf("[%d]",i);
                  printf("[%d]",j);
                  printf("= %d\n",bestPrevI);
                }
            } 
        }
    }  

  // Okay, now we have the DP homoplasy lookup table and also the backtrack table. it remains to just find the best block partition using these tabes!

  if (debugOutput > 0)
    { 
      printf("Determining optimal block partition.\n");
    }
  

  // Given DPHomoplasyRatioLookupTable, it remains to find the value corresponding a block partition of the full set of characters (i.e. j = characterCount-1) that has minimum homoplasy ratio.
  // Find the minimum value over all entries DPHomoplasyRatioLookupTable[b][i][j] with j = characterCount-1.
  // (Note that b may be smaller than desiredBlockCount-1, if splitting into fewer than desiredBlockcount blocks 
  // actually gives better homoplasy ratio.

  // optFinalIForBlockCount[b] will return the value i for which BP_DPHomoplasyRatioLookupTable[b][i][sites-1] is minimal,
  // i.e. the start point of the end block in an optimal partition with b blocks.
  int optFinalIForBlockCount [desiredBlockCount]; 
  float optimumRatioForBlockCount [desiredBlockCount]; 


  //int optFinalI = -1;
  for(int b = 0; b < desiredBlockCount; b++)
    {
      optimumRatioForBlockCount[b] = FLT_MAX;
      optFinalIForBlockCount[b] = -1;
      for(int i = 0; i < sites; i++)
        {
          float tempVal = BP_DPHomoplasyRatioLookupTable[b][i][sites-1];
          if (tempVal < optimumRatioForBlockCount[b] )
            {
              optimumRatioForBlockCount[b] = tempVal;
              optFinalIForBlockCount[b] = i;
            }
        }
    }
  // Store the minimum value found
  // Also store the minimum block number for which ratio at most desiredHomoplasyRatio was achieved (if this was requested), and the achieved ratio
  int optFinalB = -1;
  float currentOptimum = FLT_MAX;
  int minimumBForRequestedRatio = -1;
  float ratioForMinimumBForRequestedRatio = FLT_MAX;
  float ratioForCorrectBlockCount = FLT_MAX;
  float internalBoundaryError = FLT_MAX;


  char* optimumBlockPartitionString = (char*) malloc(sizeof(char) * 1000);
  char* fewestBlockForRatioPartitionString = (char*) malloc(sizeof(char) * 1000);
  char* blockPartitionOfCorrectSizeString = (char*) malloc(sizeof(char) * 1000);


  for(int b = 0; b < desiredBlockCount; b++)
    {
      if (optimumRatioForBlockCount[b] < currentOptimum)
      {
        currentOptimum = optimumRatioForBlockCount[b];
        optFinalB = b;
      }
      // if this is the first block count for which we achieve desired ratio, record the block count and acheived ratio
      if (optimumRatioForBlockCount[b] <= desiredHomoplasyRatio && minimumBForRequestedRatio == -1)
      {
        minimumBForRequestedRatio = b;
        ratioForMinimumBForRequestedRatio = optimumRatioForBlockCount[b];
      }
    }


  float BP_optHomoplasyRatio = currentOptimum;
  int optBlockCount = optFinalB + 1;


  // Tada! We have found the optimal block partition. Hopefully. Time to report it!




  // Exciting new feature: report back an optimal block partition for EACH block count in the desired range, not just the optimum
  for (int tempB = 0; tempB < desiredBlockCount; tempB++)
    {
      // Now use backtracking to reconstruct the details of an optimal solution
  
      int BP_optBlockStarts[tempB]; // List of the first elements of each block in an optimal block partition
      int BP_optBlockEnds[tempB];	// List of the last elements of each block in an optimal block partition
      //int BP_optTrees[optFinalB + 1];	// List of (the indices of?????) optimal trees for each block in an optimal block partition
  
      // The first block we store information on is actually the last block, which we know ends with the last character.
      int currentBlockStart = optFinalIForBlockCount[tempB];
      int currentBlockEnd = sites - 1;
      //int currentBlockCount = optFinalB;


 
      // Work backwards, finding the previous block and recording its details, until we reach the first block
      for (int b = tempB; b >= 0; b--)  
        {
         // printf("heeey%d\n", b);
          if (b >0)
            {
              assert (currentBlockStart <= currentBlockEnd);
              assert (currentBlockStart >= 0);
              assert (currentBlockEnd >= 0);
            }
          // append the details of the current block to their respective lists.
          BP_optBlockStarts[b] = currentBlockStart;
          BP_optBlockEnds[b] = currentBlockEnd;
          //printf("hooo %d\n", currentBlockStart);
          //printf("hrrrrrro %d\n", currentBlockEnd);
          //BP_optTrees.add(BP_optTreePerBlock[currentBlockStart][currentBlockEnd]);
  
          int previousBlockStart = BP_backtrackTable[b][currentBlockStart][currentBlockEnd];
          int previousBlockEnd = currentBlockStart - 1;
          //int previousBlockCount= currentBlockCount - 1;
          //printf("haaaaa\n");


          // previous block becomes the new current block
          currentBlockStart = previousBlockStart;
          currentBlockEnd = previousBlockEnd;
          //currentBlockCount = previousBlockCount;
  
        }

      // Now we have the optimal block partition for this block count, we can do things with it.
      // First, we calculate and store certain values / strings associated with block partitions of interest.

      // calculate the block partition string for this partition
      char* blockPartitionString = getBlockPartitionString(tempB + 1, BP_optBlockStarts, BP_optBlockEnds);
      // If this block partition is the optimal one, store the block partition string (we have already stored the optimal homoplasy ratio)
      if (tempB == optFinalB)
        {
          strcpy(optimumBlockPartitionString, blockPartitionString);
        }
      // If we were given a desired ratio and this block count is the smallest one that acheives that ratio, 
      // store the block partition string (we have already stored the acheived homoplasy ratio)
      if (desiredHomoplasyRatio != FLT_MAX && tempB == minimumBForRequestedRatio)
        {
          strcpy(fewestBlockForRatioPartitionString, blockPartitionString);
        }
      // If we were given a "true" block partition and this block partition has the same number of blocks,
      // store the block partition string, the ratio acheived and the internal boundary error
      if (tempB == truePartitionBlockCount - 1)
        {
          strcpy(blockPartitionOfCorrectSizeString, blockPartitionString);
          ratioForCorrectBlockCount = optimumRatioForBlockCount[tempB];
          internalBoundaryError = averageInternalBoundaryError(truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, BP_optBlockStarts, BP_optBlockEnds);
        }


      // If verbosity is at least 1, report the optimum block partition for this block count.
      if (verbosity >= 1)
        {
          // Introduce the block partition (and any interesting properties)
          printBothOpen("Best %d-block partition:\n", tempB+1); 
          if (tempB == optFinalB)
            {
              printBothOpen("  **MINIMUM HOMOPLASY RATIO**\n", tempB+1);
            }
          if (desiredHomoplasyRatio != FLT_MAX && tempB == minimumBForRequestedRatio)
            {
              printBothOpen("  **FEWEST BLOCKS ACHIEVING HOMOPLASY RATIO AT MOST %f**\n", desiredHomoplasyRatio);
            }
          if (tempB == truePartitionBlockCount - 1)
            {
              printBothOpen("  **SAME BLOCK COUNT AS INPUT BLOCK PARTITION**\n");
            }

         // Print the block partition
         printBothOpen("  %s\n", blockPartitionString);

         // If verbosity is at least 2, report tree for each block in the optimum partition  TODO and ratio!
         if (verbosity >= 2)
           {
             int tempBlockStart;
             int tempBlockEnd;
             int blockHomoplasyScore;
             float blockHomoplasyRatio;
             char *blockTreeNewickString;
             for (int h = 0; h <= tempB; h++)
               {
                 tempBlockStart = BP_optBlockStarts[h];
                 tempBlockEnd = BP_optBlockEnds[h];
                 blockHomoplasyScore = BP_blockScoreArray[tempBlockStart][tempBlockEnd];
                 blockHomoplasyRatio = getHomoplasyRatio(tempBlockStart, tempBlockEnd, blockHomoplasyScore);
                 blockTreeNewickString = BP_blockTreeStringArray[tempBlockStart][tempBlockEnd];
                 printBothOpen("    [%d---%d] Homoplasy score: %d Homoplasy Ratio: %f\n",  tempBlockStart, tempBlockEnd, blockHomoplasyScore, blockHomoplasyRatio);
                 printBothOpen("    [%d---%d] Parsimony tree:\n",  tempBlockStart, tempBlockEnd);
                 printBothOpen("      %s", blockTreeNewickString);
               }
           }




         // Print the max homoplasy ratio and nulber of blocks
         printBothOpen("  Max homoplasy ratio %f, Block count %d\n", optimumRatioForBlockCount[tempB], tempB+1);

          // if the  partition we're considering has the same number of blocks as the 'true' partition, then print its average internal boundary error
          if (tempB == truePartitionBlockCount - 1)
            {
              printBothOpen("  **Average internal boundary error: %f**\n",  internalBoundaryError);
            }

        }
    }


 // Print a summary of the interesting stuff

  printBothOpen("---------------------------------\n");
  printBothOpen("Summary: \n");
  printBothOpen("  Optimal homoplasy ratio: \n");
  printBothOpen("    %s\n", optimumBlockPartitionString);
  printBothOpen("    Max homoplasy ratio %f, Block count %d\n", currentOptimum , optFinalB+1);
  // Report the shortest block partition acheiving desiredMaxHomoplasyRatio, if desiredMaxHomoplasyRatio was specified.
  if (desiredHomoplasyRatio != FLT_MAX)
    {
      printBothOpen("  Fewest blocks achieving homoplasy ratio at most %f: \n", desiredHomoplasyRatio );
      if (minimumBForRequestedRatio != -1)
        {
          printBothOpen("    %s\n", fewestBlockForRatioPartitionString);

          printBothOpen("    Max homoplasy ratio %f, Block count %d\n", ratioForMinimumBForRequestedRatio , minimumBForRequestedRatio+1);
        }
      else
        {
          printBothOpen("    No solution found.\n");
        }
    }
 

  // If we were given info on the "true" block partition, report back on how our partition (for same number of blocks) compares
  
  if (trueBlockPartitionString[0] != '\0')
    {
      printBothOpen("Comparision to input partition: \n");
      printBothOpen("  Input block partition: \n");
      printBothOpen("    %s\n", getBlockPartitionString(truePartitionBlockCount, truePartitionStartPoints,  truePartitionEndPoints));
      float trueRatio = getMaxHomoplasyRatioFromPartition(truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, BP_blockScoreArray);
      printBothOpen("    Max homoplasy ratio %f, Block count %d\n", trueRatio, truePartitionBlockCount);
      printBothOpen("  Optimal %d-block partition: \n", truePartitionBlockCount);
      if (truePartitionBlockCount <= desiredBlockCount)
        {
          printBothOpen("    %s\n", blockPartitionOfCorrectSizeString);
          printBothOpen("    Max homoplasy ratio %f, Block count %d\n", ratioForCorrectBlockCount ,  truePartitionBlockCount);
          printBothOpen("    Average internal boundary error %f\n", internalBoundaryError);
        }
      else
        {
          printBothOpen("    No solution found (input block partition has more blocks than specified maximum)\n");
        }
    }


}



int main (int argc, char *argv[])
{
  rawdata      *rdta;
  //cruncheddata *cdta;  //not needed anymore
  tree         *tr;
  analdef      *adef;

  
 
  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  // cdta = (cruncheddata *)malloc(sizeof(cruncheddata));  //not needed anymore
  tr   = (tree *)malloc(sizeof(tree));

  /* initialize lookup table for fast bit counter */

  compute_bits_in_16bits();
  
  initAdef(adef);

  
  get_args(argc,argv, adef);    
 

  getinput(adef, rdta, tr);  // , cdta //not needed anymore
  

  //makeFileNames();     //we do not want to print for now 
 

  makevalues(rdta, tr);   //, cdta //not needed anymore     
  
   //we do not want to print for now 
  //printBothOpen("\n%s version %s a fast open-source code for builiding DNA parsimony start trees\n\n", programName, programVersion); 
  //printBothOpen("Released under GNU GPL in %s by Alexandros Stamatakis\n\n", programDate);
  //printBothOpen("Alignment has %d sites and %d taxa\n\n", tr->cdta->endsite, tr->mxtips);
  //printBothOpen("%d randomized stepwise addition order parsimony trees with a couple of SPR moves will be computed\n\n", adef->numberOfTrees);



  //int minParsimonyScoreOverRuns = makeParsimonyTreeFastDNA(tr, adef, 0, rdta->sites); //rdta->sites);
  //printf ("doing a test\n");
  //int minParsimonyScoreOverRuns = makeParsimonyTreeFastDNA(tr, adef, 0, rdta->sites - 1); //rdta->sites);
  //printf ("test done\n");
  //minParsimonyScoreOverRuns = makeParsimonyTreeFastDNA(tr, adef, 0, 0); //rdta->sites);


  getOptBlockPartition(tr, adef);


  //printf("Best parsimony score %u\n", minParsimonyScoreOverRuns); 


  return 0;
}


