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
  partitionInfoFileName[1024] = "";


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


// Frees memory associated with a 2d int array arr
static void free2dIntArray(int **arr)
  {
    free(arr[0]);
    free(arr);
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

// Frees memory associated with a 3d float array arr
static void free3dFloatArray(float ***arr)
  {
    free(arr[0][0]);
    free(arr[0]);
    free(arr);
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

// Frees memory associated with a 3d int array arr
static void free3dIntArray(int ***arr)
  {
    free(arr[0][0]);
    free(arr[0]);
    free(arr);
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


// Frees memory associated with a 3d char array arr
static void free3dCharArray(char ***arr)
  {
    free(arr[0][0]);
    free(arr[0]);
    free(arr);
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




// Given the start and end points of a block and a list of which sites are uninformative, returns the informative restriction of that block,
// i.e. the start and end points of that block with uninformative sites trimmed off the beginning and end.
// e.g. if the input is [3---10] and sites 3, 4, 7 and 10 are uninformative but sites 5,6,8,9 are informative, the returned block will be [5---9].
// outputBlockStart is set to the start point of the informative restriction.
// outputBlockEnd is set to the end point of the informative restriction.
// If all sites in the block are uninformative, outputBlockStart and outputBlockEnd are both set to -1
static void getInformativeRestrictionOfBlock(int* recordOfInformativeSites, int blockStart, int blockEnd, int *outputBlockStart, int *outputBlockEnd)
{

  // if  blockStart == -1 or blockEnd == -1: we've been given a 'null block' / a representation of an uninformative block, and we should return the same thing
  if (blockStart == -1 || blockEnd == -1)
  {
    *outputBlockStart = -1;
    *outputBlockEnd = -1;
  } else {

    assert(blockStart <= blockEnd);
    int trueBlockStart = -1;
    // the true block start is the first site between blockStart and blockEnd that is informative.
    for (int i = blockStart; i <= blockEnd; i++) 
    {
      if (recordOfInformativeSites[i] != 0)
      {
         trueBlockStart = i;
         break;
      }
    }
  
    int trueBlockEnd = -1;
    // the true block end is the *last* site between blockStart and blockEnd that is informative.
    for (int i = blockEnd; i >= blockStart; i--) 
    {
      if (recordOfInformativeSites[i] != 0)
      {
         trueBlockEnd = i;
         break;
      }
    }
    //printBothOpen("Informative restriction of [%d---%d]: [%d---%d]: \n",  blockStart, blockEnd, trueBlockStart, trueBlockEnd);
    *outputBlockStart = trueBlockStart;
    *outputBlockEnd = trueBlockEnd;
  }
}

// Given the start and end points for a block partition, returns the informative restriction of  each block - i;e. trims uninformative sites from the ends of blocks.
// e.g. If numsites = 10 and the informative sites are 2,3,4, 7,9, and the input partition is [0---5],[6---9], then the return block partiton is [2---4],[7---9].
// If some blocks are completely uninformative then those blocks will be replaced with [-1----1] (i.e. -1 at both ends)
static void getInformativeRestrictionOfPartition (int numBlocks, int* recordOfInformativeSites, int *blockStarts, int *blockEnds, int *outputBlockStarts, int *outputBlockEnds)
{
  for (int h = 0; h < numBlocks; h++)
    {
      getInformativeRestrictionOfBlock(recordOfInformativeSites, blockStarts[h], blockEnds[h], &outputBlockStarts[h], &outputBlockEnds[h]);
    }
}


// Given a block partition that is already in restricted form, remove any blocks that have start and end points -1, and reduce the block count accordingly
static void removeBadBlocksFromRestrictedPartition (int *numBlocks, int* recordOfInformativeSites, int *blockStarts, int *blockEnds)
{

  int i = 0; // i keeps track of the number of actual informative blocks we've found so far
  for (int h = 0; h < *numBlocks; h++)
    {
      // if block h is informative, copy its details to the next available position   
      // (this does nothing if all blocks so far have been informative, but otherwise copies details to an earlier position
      if (blockStarts[h] != -1 && blockEnds[h] != -1)
      {
        blockStarts[i] = blockStarts[h];
        blockEnds[i] = blockEnds[h];
        i += 1;
      }
    }
  int numInformativeBlocks = i;
  // For the sake of tidying up, remove redundant copies of data from later in the arrays.
  for (int h = numInformativeBlocks; h < *numBlocks; h++)
  {
    blockStarts[h] = -1;
    blockEnds[h] = -1;
  }

  // Finally, update the block count
  *numBlocks = numInformativeBlocks;
}



// Given a string trueBlockPartitionString of the form [0---x],[x+1---y],.... ,[w+1---z] representing a block partition
// Sets truePartitionBlockCount to be the number of blocks in the partition,
// And sets truePartitionStartPoints (respectively, truePartitionEndPoint) to be an array of length truePartitionBlockCount containing 
// the start points (respectively, end points) of all blocks in the partition.
// TODO: spell out the format required for trueBlockPartitionString: no spaces, blocks separated by commas, blocks have the form [x---y] where x and y are integers. TODO other constraints?
static void getTrueBlockPartitionData( char* trueBlockPartitionString, int *truePartitionBlockCount, int **truePartitionStartPoints, int **truePartitionEndPoints, int sites,  int* recordOfInformativeSites)
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
  // UPDATE: taking out this check as we now allow 'partial' block partitions as input
  //if (tempStartArray[0] != 0)
  //  {
  //    printf("Error4: badly formatted block input (first block must begin with 0)\n");
  //    exit(-1);
  //  }

  // check that the block partition is well formatted
  // UPDATE: taking out this check as we now allow 'partial' block partitions as input
  //if (tempEndArray[*truePartitionBlockCount -1] != sites -1)
  //  {
  //    printf("Error5: badly formatted block input (last block [%d---%d] must end with %d (sequence length -1))\n", tempStartArray[*truePartitionBlockCount -1], tempEndArray[*truePartitionBlockCount -1], sites-1);
  //    exit(-1);
  //  }

  for (b = 0; b < *truePartitionBlockCount; b++)
    {
       if (tempStartArray[b] > tempEndArray[b])
         {
           printf("Error6: badly formatted block input [%d---%d](block start cannot be greater than block end)\n", tempStartArray[b], tempEndArray[b]);
           exit(-1);
         }
    }


  // UPDATE: rewriting this check as we now allow 'partial' block partitions as input
  for (b = 0; b < *truePartitionBlockCount -1; b++)
    {
  //     if (tempStartArray[b+1] != tempEndArray[b] + 1)
  //       {
  //         printf("Error7: badly formatted block input [%d---%d],[%d---%d] (block start must equal previous block end + 1)\n", tempStartArray[b], tempEndArray[b], tempStartArray[b+1], tempEndArray[b+1]);
  //         exit(-1);
  //       }
       if (tempStartArray[b+1] <= tempEndArray[b])
         {
           printf("Error7: badly formatted block input [%d---%d],[%d---%d] (block start must be greater than previous block end)\n", tempStartArray[b], tempEndArray[b], tempStartArray[b+1], tempEndArray[b+1]);
           exit(-1);
         }
    }



  char* receivedBlockPartitionString = getBlockPartitionString(*truePartitionBlockCount, tempStartArray, tempEndArray);
  printBothOpen("Input block partition: %s\n", receivedBlockPartitionString);   


  // now get the informative restriction of this block partition

  // make restrictedStartArray and restrictedEndArray be arrays of length truePartitionBlockCount
  int* restrictedStartArray = (int*) malloc(sizeof(int) * *truePartitionBlockCount);
  int* restrictedEndArray = (int*) malloc(sizeof(int) * *truePartitionBlockCount);
  getInformativeRestrictionOfPartition(*truePartitionBlockCount, recordOfInformativeSites, tempStartArray, tempEndArray, restrictedStartArray, restrictedEndArray);


  //char* restrictedBlockPartitionString = getBlockPartitionString(*truePartitionBlockCount, restrictedStartArray, restrictedEndArray);
  //printBothOpen("Restricted block partition: %s\n", restrictedBlockPartitionString);   

  int newBlockCount = *truePartitionBlockCount;
  removeBadBlocksFromRestrictedPartition (&newBlockCount, recordOfInformativeSites,  restrictedStartArray, restrictedEndArray);
  char* lessBadBlockPartitionString = getBlockPartitionString(newBlockCount, restrictedStartArray, restrictedEndArray);

  // check that number of blocks has not changed i.e. input partition did not contain uninformative blocks
  if (newBlockCount !=  *truePartitionBlockCount)
    {
      printf("WARNING 1: input block partition contains some blocks with no informative sites.\n");
    }

  *truePartitionBlockCount = newBlockCount;


  // Another check we want to do: Is the input block partition skipping over any sites that *are* informative?
  boolean skippedInformativeSites = FALSE;
  int h = 0;
  // check all the sites before the first block
  for (i = 0; i <  restrictedStartArray[0]; i++)
    {
      if (recordOfInformativeSites[i] != 0)
        {
            skippedInformativeSites = TRUE;
            break;
        }
    }
  // check all sites between end of block h and start of block h+1 (notinclusive), for all h < newblockCount - 1
  for (h = 0; h < newBlockCount - 1; h++)
    {
      for (i = restrictedEndArray[h]+1; i <  restrictedStartArray[h+1]; i++)
        {
          if (recordOfInformativeSites[i] != 0)
            {
                skippedInformativeSites = TRUE;
                break;
            }
        }

      if (skippedInformativeSites)
        {
          break;
        }
    }

  // check all sites after end of the first block
  for (i = restrictedEndArray[newBlockCount-1]+1; i < sites; i++)
    {
      if (recordOfInformativeSites[i] != 0)
        {
            skippedInformativeSites = TRUE;
            break;
        }
    }

  if (skippedInformativeSites)
    {
      printf("WARNING 2: input block partition skips some informative sites.\n");
    }



  // report informative restriction of input partition
  printBothOpen("Informative restriction of input block partition: %s\n", lessBadBlockPartitionString);   

  // now set truePartitionStartPoints and truePartitionEndPoints to be the two *restricted* arrays we just constructed.
  free(*truePartitionStartPoints);
  *truePartitionStartPoints = restrictedStartArray;
  free(*truePartitionEndPoints);
  *truePartitionEndPoints = restrictedEndArray;

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


  // parse and store data on the True Block Partition  
 // char* trueBlockPartitionString = adef->trueBlockPartitionString; 
 // int truePartitionBlockCount;
 // int* truePartitionStartPoints;
 // int* truePartitionEndPoints; 
  //int sites = rdta->sites;
  //if (adef->trueBlockPartitionString[0] != '\0')
  //  {
  //    char tempBlockPartitionString[1024] = ""; 
  //    strcpy(tempBlockPartitionString,adef->trueBlockPartitionString);
  //    getTrueBlockPartitionData(tempBlockPartitionString, &(adef->truePartitionBlockCount), &(adef->truePartitionStartPoints), &(adef->truePartitionEndPoints), sites, tr->recordOfInformativeSites);
  //  }
  

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
  strcpy(partitionInfoFileName,        "output.out");
  
  adef->parsimonySeed=1;
  adef->numberOfBlocks=2;
  adef->maxHomoplasyRatio= FLT_MAX;
  adef->maxHomoplasyScore= INT_MAX;	
  adef->totalHomoplasyRatio= FLT_MAX;   
  adef->totalHomoplasyScore= INT_MAX;   
  strcpy(adef->trueBlockPartitionString, "\0"); // 27.11.2018 // trying to fix an "unitialized variable" warning further down the road - MJ
  adef->verbose= 1;
  adef->problems= 15; // solve all problems by default
  
  while(!bad_opt && ((c = mygetopt(argc,argv,"p:n:s:o:t:N:b:B:r:R:h:H:v:P:", &optind, &optarg))!=-1))
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
	  strcpy(adef->trueBlockPartitionString, optarg);
	  break;        
	case 'r':
	  sscanf(optarg,"%f", &(adef->maxHomoplasyRatio));	// desired maximum homoplasy ratio - we're interested in shortest block partition that achieves this ratio
	  break;           
	case 'R':
	  sscanf(optarg,"%f", &(adef->totalHomoplasyRatio));	// desired total homoplasy ratio - we're interested in shortest block partition that achieves this ratio
	  break;   
	case 'h':
	  sscanf(optarg,"%d", &(adef->maxHomoplasyScore));	// desired maximum homoplasy score - we're interested in shortest block partition that achieves this score
	  break;           
	case 'H':
	  sscanf(optarg,"%d", &(adef->totalHomoplasyScore));	// desired total homoplasy score - we're interested in shortest block partition that achieves this score
	  break;    
	case 'v':
	  sscanf(optarg,"%d", &(adef->verbose));	// desired verbosity - 0: just return summary information. 1: return block partition (and max homoplasy ratio) for each possible block count. 2: return homoplasy ratio and parsimonious tree for each block in each partition. (set to 1 by default)
	  break;           
	case 'P':
	  sscanf(optarg,"%d", &(adef->problems));	// description of  which problems to solve
							// Decided using bitwise operators
							// If (P & 1) > 0 :  solve MaxScore
							// If (P & 2) > 0 :  solve TotalScore
							// If (P & 4) > 0 :  solve MaxRatio
							// If (P & 8) > 0 :  solve TotalRatio
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



// Returns the homoplasy ratio of the *informative restriction* of the block [i---j], assuming the block has the specified  homoplasy value.
// Returns FLT_MAX if the block is completely uninformative
static float getHomoplasyRatio(tree *tr, int i, int j, int homoplasyValue) 
{
  assert  (i <= j);  // program breaks if called with i greater than j


  int informativeBlockStart;
  int informativeBlockEnd;
  getInformativeRestrictionOfBlock(tr->recordOfInformativeSites, i, j, &informativeBlockStart, &informativeBlockEnd);
  if (informativeBlockStart == -1 || informativeBlockEnd == -1)
  { 
    return FLT_MAX;
  } else {  
    return ((float) homoplasyValue) / ((float)(informativeBlockEnd-informativeBlockStart+1));
  }
}

// Given a number of blocks and two arrays describing the start and end points of each block, 
// and a lookup table of homoplasy scores for each possible block,
// returns the maximum homoplasy ratio of the given set of blocks
static float getMaxHomoplasyRatioFromPartition(tree *tr, int blockCount, int* blockStartPoints, int* blockEndPoints, int** blockScoreArray)
  {
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
        tempHomoplasyRatio = getHomoplasyRatio(tr, tempBlockStart,tempBlockEnd,tempHompolasyScore);
        if (tempHomoplasyRatio > maxRatioSoFar)
          {
            maxRatioSoFar = tempHomoplasyRatio;
          }
      }
    return maxRatioSoFar;
  }

// Given a number of blocks and two arrays describing the start and end points of each block, 
// and a lookup table of homoplasy scores for each possible block,
// returns the total homoplasy ratio of the given set of blocks
static float getTotalHomoplasyRatioFromPartition(tree *tr, int blockCount, int* blockStartPoints, int* blockEndPoints, int** blockScoreArray)
  {
    float totalRatioSoFar = 0;
    int tempBlockStart;
    int tempBlockEnd;
    int tempHomoplasyScore;
    float tempHomoplasyRatio;
    int b;
    for (b = 0; b < blockCount; b++)
      {      
        tempBlockStart = blockStartPoints[b];
        tempBlockEnd = blockEndPoints[b];
        tempHomoplasyScore = blockScoreArray[tempBlockStart][tempBlockEnd];
        tempHomoplasyRatio = getHomoplasyRatio(tr, tempBlockStart,tempBlockEnd,tempHomoplasyScore);
        totalRatioSoFar += tempHomoplasyRatio;
      }
    return totalRatioSoFar;
  }

// Given a number of blocks and two arrays describing the start and end points of each block, 
// and a lookup table of homoplasy scores for each possible block,
// returns the maximum homoplasy score of the given set of blocks
static int getMaxHomoplasyScoreFromPartition(int blockCount, int* blockStartPoints, int* blockEndPoints, int** blockScoreArray)
  {
    int maxScoreSoFar = 0;
    int tempHomoplasyScore;
    int tempBlockStart;
    int tempBlockEnd;
    int b;
    for (b = 0; b < blockCount; b++)
      {      
        tempBlockStart = blockStartPoints[b];
        tempBlockEnd = blockEndPoints[b];
        tempHomoplasyScore = blockScoreArray[tempBlockStart][tempBlockEnd];
        if (tempHomoplasyScore > maxScoreSoFar)
          {
            maxScoreSoFar = tempHomoplasyScore;
          }
      }
    return maxScoreSoFar;
  }

// Given a number of blocks and two arrays describing the start and end points of each block, 
// and a lookup table of homoplasy scores for each possible block,
// returns the total homoplasy score of the given set of blocks
static int getTotalHomoplasyScoreFromPartition(int blockCount, int* blockStartPoints, int* blockEndPoints, int** blockScoreArray)
  {
    float totalScoreSoFar = 0;
    int tempBlockStart;
    int tempBlockEnd;
    int tempHomoplasyScore;
    int b;
    for (b = 0; b < blockCount; b++)
      {      
        tempBlockStart = blockStartPoints[b];
        tempBlockEnd = blockEndPoints[b];
        tempHomoplasyScore = blockScoreArray[tempBlockStart][tempBlockEnd];
        totalScoreSoFar += tempHomoplasyScore;
      }
    return totalScoreSoFar;
  }



// Returns the number of informative sites difference between site i and site j
// Assumes that both i and j are start or end points of the informative restrictions of some blocks;
// in particular requires that both sites i and j are informative sites
static int getInformativeSiteDifference(int* recordOfInformativeSites, int i, int j)
{

  // check that both i and j are informative sites; if this is not the case the informative site difference depends on information we don't have access to  (namely whether we're testing start points or endpoints)
  assert(recordOfInformativeSites[i] == 1);
  assert(recordOfInformativeSites[j] == 1);

  int minSite;
  int maxSite;
  if (i <= j)
    {
      minSite = i;
      maxSite = j;
    } else {
      minSite = j;
      maxSite = i;
    }
  // Calculate the number of informative sites between i and j inclusive.
  // If there are none, the informative site difference is 1 less than this  (assuming i and j are both informative)
  int numInformativeSites = 0;
  int h;
  for (h = minSite; h<= maxSite; h++)
    {
      if (recordOfInformativeSites[h] > 0)
        {
          numInformativeSites++;
        }
    }
  return numInformativeSites - 1; //
}

// TODO update this description
// Returns the average internal boundary error for a block partition with blockCount blocks and endpoints defined by testBlockStarts, testBlockEnds, when compared to the block partition defined by truePartitionStartPoints and truePartitionEndPoints
// For two block partitions [x_1 -- y_1],...,[x_b---y_b] and [w_1---z_1],...,[w_b---z_b],
// And assuming that x_{i+1} = y_i +1, w_{i+1} = z_i + 1 for all  1 <= i < b,
// the average internal boundary error is defined as 
//   ( abs(y_1-z_1) + abs(y_2-z_2) + ... + abs(y_{b-1}-z_{b-1}) ) / (b-1)
static float averageInternalBoundaryError(int* recordOfInformativeSites, int blockCount, int* truePartitionStartPoints, int* truePartitionEndPoints, int* testBlockStarts, int* testBlockEnds)
{
  int totalBoundaryError = 0;
  int i;
  int restrictedTrueBlockStarts[blockCount];
  int restrictedTrueBlockEnds[blockCount];
  int restrictedTestBlockStarts[blockCount];
  int restrictedTestBlockEnds[blockCount];



  char* tempStr = getBlockPartitionString(blockCount, truePartitionStartPoints, truePartitionEndPoints);
  //printBothOpen("AIBE true block partition: %s\n", tempStr);

  // get the informative restriction of both partitions
  getInformativeRestrictionOfPartition(blockCount, recordOfInformativeSites, truePartitionStartPoints, truePartitionEndPoints, restrictedTrueBlockStarts, restrictedTrueBlockEnds);

  getInformativeRestrictionOfPartition(blockCount, recordOfInformativeSites, testBlockStarts, testBlockEnds, restrictedTestBlockStarts, restrictedTestBlockEnds);



  for (i = 0; i < blockCount - 1; i++)  // We don't test the endpoints for the final block, as these are assumed to both be taxaCount-1
    {      
      // local boundary error is the average of the difference in end points plus difference in start points
      // (these two numbers should be the same, assuming both partitions cover all informative sites)
      float localBoundaryError =  ( getInformativeSiteDifference(recordOfInformativeSites, restrictedTrueBlockEnds[i], restrictedTestBlockEnds[i]) + getInformativeSiteDifference(recordOfInformativeSites, restrictedTrueBlockStarts[i+1], restrictedTestBlockStarts[i+1])) / 2;
      // printf("local error: %d\n", localBoundaryError);
      totalBoundaryError += localBoundaryError;
      // printf("total error so far: %d\n", totalBoundaryError);
    }
  // TODO: there's also a question of whether we should do any measurement on the start of the first block and end of the last block,
  // Which again should be 0 if the partitions cover all informative sites, but maybe it would good for future-proofing.
  // And if we did want to  count the start end end differences, would they also be divided by 2?

  // printf("Total error %d / internal boundary count %d\n", totalBoundaryError, blockCount -1);
  float averageError = (float) totalBoundaryError / (blockCount -1);
  // printf("Average error %f\n", averageError);
  return averageError;
}


static int getOptHomoplasyOfContiguousBlock(tree *tr, analdef *adef, int startCharacter, int endCharacter)
{
  //printf("Calculating [%d -- %d ]:", startCharacter, endCharacter);
  assert  (startCharacter <= endCharacter); 

  // first check that the block is informative: if it's not, we treat it as having score \infty and don't need to run Parsimonator
  int tempStart;
  int tempEnd;
  getInformativeRestrictionOfBlock(tr->recordOfInformativeSites, startCharacter,endCharacter, &tempStart, &tempEnd);
  if (tempStart == -1 || tempEnd == -1)
   {
     return INT_MAX;
   } else {


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
}



// Populate each element of maxScoreDPTable and maxScoreBacktrackTable.
//
// For a block partition P, let maxScore(P) denote the maximum homoplasy score of any block in P.
//
// Then maxScoreDPTable[b][i][j] for i<= j is the minimum value of maxScore(P)
// over any block partition P of characters 0 to j, 
// for which the last block is from i to j and there are b+1 blocks in total.
// A value of INT_MAX means there is no solution 
// for this set of values (e.g. if i > j)
// 
// maxScorePacktrackTable[b][i][j] is the value i' such that
// there exists an optimum block partition P satisfying  maxScoreDPTable[b][i][j],
// in which the second-last block is on characters i' to i-1.
// A value of -1 means there is no such i' for this set of values.  (e.g. if i = 0, or maxScorePacktrackTable[b][i][j] = INT_MAX.)
static void  doMaxScoreDP(tree *tr, analdef *adef, int **BP_blockScoreArray, int ***maxScoreDPTable, int ***maxScoreBacktrackTable)

{
  // initialize basic variables
  int debugOutput = 0;
  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;

  // Populate each element of maxScoreDPTable and maxScoreBacktrackTable in turn, starting with the lower-indexed values of i, j and b
  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++)
        {
          for (int b = 0; b < desiredBlockCount; b++)
            {
              // val and bestPrevI will be the values assigned to maxScoreDPTable[b][i][j] and maxScorePacktrackTable[b][i][j]
              int val =  INT_MAX;
              int bestPrevI = -1;

              // First handle the trivial cases where there is no solution.

              // If only one block, i should be 0
              if (b == 0 && i > 0)
                {
                  val = INT_MAX;		
                }

              // If more than one block, i should be bigger than 0
              if (b > 0 && i == 0)
                {
                  val = INT_MAX;		
                }

              // Now process the case when we have a single block.
              // In this case, maxScore(P) is just the homoplasy score of that block.
              if (b == 0 && i == 0)			
                {
                   val =  BP_blockScoreArray[i][j];
                }

              // Now process the case when we have more than one block
              if (b > 0 && i > 0)
                {
                  // For i > j, no block partition is defined.
                  if (i > j)			
                    {
                      val = INT_MAX;	
                    } 

                  // In all other cases, maxScore is based on the homoplasy score of block [i---j], 
                  // and the optimum maxScore of a partition with one fewer block on characters 0 to i-1. 
                  else 
                    {
                      //find the minimum possible homoplasy score for any b-1 (+1) block partition up to this point
                      int currentMinimum = INT_MAX;
                      int lastj = i-1;
                      for(int lasti = 0; lasti <= lastj; lasti++)
                        {
                          int tempVal = maxScoreDPTable[b-1][lasti][lastj];
 	
                           if (tempVal < currentMinimum )
                            {
                              currentMinimum = tempVal;
                              // update the bestPrevI value because we have a new best block partition
                              bestPrevI = lasti;
                            }
                        }
                      // currentMinimum is now the optimum homoplasy score for any block partition of 0 to i-1
                      // using b-1 (+1) blocks
                      int valCurrentBlock =  BP_blockScoreArray[i][j];
                      // If the current block has worse homoplasy score than what we could 
                      // get up to this point, that score is now the best we can do for a block partition of this type.
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
              maxScoreDPTable[b][i][j] = val;

              // also populate the backtrack table so we can reconstruct a solution later
              maxScoreBacktrackTable[b][i][j] = bestPrevI;

              if (debugOutput > 0)
                {
                  printf("maxScoreDPTable[%d][%d][%d]= %d   maxScoreBacktrackTable[%d][%d][%d]=%d\n",b,i,j, val,b,i,j,bestPrevI);
                }
            } 
        }
    }   
}








// Populate each element of totalScoreDPTable and totalScoreBacktrackTable.
//
// For a block partition P, let totalScore(P) denote the total homoplasy score of all blocks in P.
//
// Then totalScoreDPTable[b][i][j] for i<= j is the minimum value of totalScore(P)
// over any block partition P of characters 0 to j, 
// for which the last block is from i to j and there are b+1 blocks in total.
// A value of INT_MAX means there is no solution 
// for this set of values (e.g. if i > j)
// 
// totalScorePacktrackTable[b][i][j] is the value i' such that
// there exists an optimum block partition P satisfying  totalScoreDPTable[b][i][j],
// in which the second-last block is on characters i' to i-1.
// A value of -1 means there is no such i' for this set of values.  (e.g. if i = 0, or totalScorePacktrackTable[b][i][j] = INT_MAX.)
static void  doTotalScoreDP(tree *tr, analdef *adef, int **BP_blockScoreArray, int ***totalScoreDPTable, int ***totalScoreBacktrackTable)

{
  // initialize basic variables
  int debugOutput = 0;
  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;

  // Populate each element of totalScoreDPTable and totalScoreBacktrackTable in turn, starting with the lower-indexed values of i, j and b
  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++)
        {
          for (int b = 0; b < desiredBlockCount; b++)
            {
              // val and bestPrevI will be the values assigned to totalScoreDPTable[b][i][j] and totalScorePacktrackTable[b][i][j]
              int val =  INT_MAX;
              int bestPrevI = -1;

              // First handle the trivial cases where there is no solution.

              // If only one block, i should be 0
              if (b == 0 && i > 0)
                {
                  val = INT_MAX;		
                }

              // If more than one block, i should be bigger than 0
              if (b > 0 && i == 0)
                {
                  val = INT_MAX;		
                }

              // Now process the case when we have a single block.
              // In this case, totalScore(P) is just the homoplasy score of that block.
              if (b == 0 && i == 0)			
                {
                   val =  BP_blockScoreArray[i][j];
                }

              // Now process the case when we have more than one block
              if (b > 0 && i > 0)
                {
                  // For i > j, no block partition is defined.
                  if (i > j)			
                    {
                      val = INT_MAX;	
                    } 

                  // In all other cases, totalScore is based on the homoplasy score of block [i---j], 
                  // and the optimum totalScore of a partition with one fewer block on characters 0 to i-1. 
                  // UPDATE (11 March 2019) as with doTotalRatioDP, we are editing things quite a bit here, and merging the i==j and i <j cases into one. This is because the existence of infinite scores on uninformative blocks means we can no longer assume block [i---i] has score 0 in the i===j case, and the 'subtract score of [i---(j-1)]  and add score of [i---j]' trick no longer works in the i<j case.
                  else // if (i==j)
                    {
                      //find the minimum possible homoplasy score for any b-1 (+1) block partition up to this point
                      int currentMinimum = INT_MAX;
                      int lastj = i-1;
                      for(int lasti = 0; lasti <= lastj; lasti++)
                        {
                          int tempVal = totalScoreDPTable[b-1][lasti][lastj];
 	
                           if (tempVal < currentMinimum )
                            {
                              currentMinimum = tempVal;
                              // update the bestPrevI value because we have a new best block partition
                              bestPrevI = lasti;
                            }
                        }
                 //     // currentMinimum is now the optimum total homoplasy score for any block partition of 0 to i-1
                 //     // using b-1 (+1) blocks
                 //     int valCurrentBlock =  BP_blockScoreArray[i][j];

                      // The optimum total homoplasy score for a block partition of this type
                      // is equal to the best possible total up to this point, plus the score of the current block 
                      if (currentMinimum != INT_MAX)
                        {

                          float valCurrentBlock = BP_blockScoreArray[i][j];
                          if (valCurrentBlock != INT_MAX)
                            {
                              val = currentMinimum + valCurrentBlock;
                            } else {
                              val = INT_MAX;
                            }
                        }
                    }
                }

              // populate the lookup table with the calculated value
              totalScoreDPTable[b][i][j] = val;

              // also populate the backtrack table so we can reconstruct a solution later
              totalScoreBacktrackTable[b][i][j] = bestPrevI;

              if (debugOutput > 0)
                {
                  printf("totalScoreDPTable[%d][%d][%d]= %d   totalScoreBacktrackTable[%d][%d][%d]=%d   BP_blockScoreArray[%d][%d] = %d\n",b,i,j, val,b,i,j,bestPrevI, i, j, BP_blockScoreArray[i][j]);
                }
            } 
        }
    }   
}




// Populate each element of totalRatioDPTable and totalRatioBacktrackTable.
//
// For a block partition P, let totalRatio(P) denote the total of all homoplasy ratios of blocks in P.
//
// Then totalRatioDPTable[b][i][j] for i<= j is the minimum value of totalRatio(P)
// over any block partition P of characters 0 to j, 
// for which the last block is from i to j and there are b+1 blocks in total.
// A value of FLT_MAX means there is no solution 
// for this set of values (e.g. if i > j)
// 
// totalRatioPacktrackTable[b][i][j] is the value i' such that
// there exists an optimum block partition P satisfying  totalRatioDPTable[b][i][j],
// in which the second-last block is on characters i' to i-1.
// A value of -1 means there is no such i' for this set of values.  (e.g. if i = 0, or totalRatioPacktrackTable[b][i][j] = FLT_MAX.)
static void  doTotalRatioDP(tree *tr, analdef *adef, int **BP_blockScoreArray, float ***totalRatioDPTable, int ***totalRatioBacktrackTable)

{
  // initialize basic variables
  int debugOutput = 0;
  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;

  // Populate each element of totalRatioDPTable and totalRatioBacktrackTable in turn, starting with the lower-indexed values of i, j and b
  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++)
        {
          for (int b = 0; b < desiredBlockCount; b++)
            {
              // val and bestPrevI will be the values assigned to totalRatioDPTable[b][i][j] and totalRatioPacktrackTable[b][i][j]
              float val =  FLT_MAX;
              int bestPrevI = -1;

              // First handle the trivial cases where there is no solution.

              // If only one block, i should be 0
              if (b == 0 && i > 0)
                {
                  val = FLT_MAX;		
                }

              // If more than one block, i should be bigger than 0
              if (b > 0 && i == 0)
                {
                  val = FLT_MAX;		
                }

              // Now process the case when we have a single block.
              // In this case, totalRatio(P) is just the homoplasy ratio of that block.
              if (b == 0 && i == 0)			
                {
                   val = getHomoplasyRatio(tr, i,j, BP_blockScoreArray[i][j]);

                }

             // Now process the case when we have more than one block
              if (b > 0 && i > 0)
                {
                  // For i > j, no block partition is defined.
                  if (i > j)			
                    {
                      val = FLT_MAX;	
                    } 

                  // In all other cases, totalRatio is based on the homoplasy ratioof block [i---j], 
                  // and the optimum totalRatio of a partition with one fewer block on characters 0 to i-1. 
                  // BIG UPDATE (11 March 2019): no longer treating the i==j and i < j cases separately,
                  // as the intorduction of FLT_MAX ratios on uninformative blocks means that we always have to find to tot ratio
                  // by checking the total for fewer blocks and then adding the ratio for the current block.
                  else // if (i==j)
                    {
                      //find the minimum possible homoplasy ratio for any b-1 (+1) block partition up to this point
                      float currentMinimum = FLT_MAX;
                      int lastj = i-1;
                      for(int lasti = 0; lasti <= lastj; lasti++)
                        {
                          float tempVal = totalRatioDPTable[b-1][lasti][lastj];
 	
                           if (tempVal < currentMinimum )
                            {
                              currentMinimum = tempVal;
                              // update the bestPrevI value because we have a new best block partition
                              bestPrevI = lasti;
                            }
                        }
                      // currentMinimum is now the optimum total homoplasy ratio for any block partition of 0 to i-1
                      // using b-1 (+1) blocks

                      // The optimum total homoplasy ratio for a block partition of this type
                      // is equal to the best possible total up to this point, plus the ratio of the current block
                      if (currentMinimum != FLT_MAX)
                        {
                          float valCurrentBlock = getHomoplasyRatio(tr,i,j, BP_blockScoreArray[i][j]);
                          if (valCurrentBlock != FLT_MAX)
                            {
                              val = currentMinimum + valCurrentBlock;
                            } else {
                              val = FLT_MAX;
                            }
                        }
                    }
                }

              // populate the lookup table with the calculated value
              totalRatioDPTable[b][i][j] = val;

              // also populate the backtrack table so we can reconstruct a solution later
              totalRatioBacktrackTable[b][i][j] = bestPrevI;

              if (debugOutput > 0)
                {
                  printf("totalRatioDPTable[%d][%d][%d]= %f   totalRatioBacktrackTable[%d][%d][%d]=%d\n",b,i,j, val,b,i,j,bestPrevI);
                }
            } 
        }
    }   
}







// Populate each element of maxRatioDPTable and maxRatioBacktrackTable.
//
// For a block partition P, let maxRatio(P) denote the maximum homoplasy ratio of any block in P.
//
// Then maxRatioDPTable[b][i][j] for i<= j is the minimum value of maxRatio(P)
// over any block partition P of characters 0 to j, 
// for which the last block is from i to j and there are b+1 blocks in total.
// A value of FLT_MAX means there is no solution 
// for this set of values (e.g. if i > j)
// 
// maxRatioPacktrackTable[b][i][j] is the value i' such that
// there exists an optimum block partition P satisfying  maxRatioDPTable[b][i][j],
// in which the second-last block is on characters i' to i-1.
// A value of -1 means there is no such i' for this set of values.  (e.g. if i = 0, or maxRatioPacktrackTable[b][i][j] = FLT_MAX.)
static void  doMaxRatioDP(tree *tr, analdef *adef, int **BP_blockScoreArray, float ***maxRatioDPTable, int ***maxRatioBacktrackTable)

{
  // initialize basic variables
  int debugOutput = 0;
  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;

  // Populate each element of maxRatioDPTable and maxRatioBacktrackTable in turn, starting with the lower-indexed values of i, j and b
  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++)
        {
          for (int b = 0; b < desiredBlockCount; b++)
            {
              // val and bestPrevI will be the values assigned to maxRatioDPTable[b][i][j] and maxRatioPacktrackTable[b][i][j]
              float val = FLT_MAX;
              int bestPrevI = -1;

              // First handle the trivial cases where there is no solution.

              // If only one block, i should be 0
              if (b == 0 && i > 0)
                {
                  val = FLT_MAX;		
                }

              // If more than one block, i should be bigger than 0
              if (b > 0 && i == 0)
                {
                  val = FLT_MAX;		
                }

              // Now process the case when we have a single block.
              // In this case, maxRatio(P) is just the homoplasy ratio of that block.
              if (b == 0 && i == 0)			
                {
                   val =  getHomoplasyRatio(tr, i,j, BP_blockScoreArray[i][j]);
                }

              // Now process the case when we have more than one block
              if (b > 0 && i > 0)
                {
                  // For i > j, no block partition is defined.
                  if (i > j)			
                    {
                      val = FLT_MAX;	
                    } 

                  // In all other cases, maxRatio is based on the homoplasy ratio of block [i---j], 
                  // and the optimum maxRatio of a partition with one fewer block on characters 0 to i-1. 
                  else 
                    {
                      //find the minimum possible homoplasy score for any b-1 (+1) block partition up to this point
                      float currentMinimum = FLT_MAX;
                      int lastj = i-1;
                      for(int lasti = 0; lasti <= lastj; lasti++)
                        {
                          float tempVal = maxRatioDPTable[b-1][lasti][lastj];
 	
                           if (tempVal < currentMinimum )
                            {
                              currentMinimum = tempVal;
                              // update the bestPrevI value because we have a new best block partition
                              bestPrevI = lasti;
                            }
                        }
                      // currentMinimum is now the optimum homoplasy ratio for any block partition of 0 to i-1
                      // using b-1 (+1) blocks
                      float valCurrentBlock =  getHomoplasyRatio(tr, i,j, BP_blockScoreArray[i][j]);
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
              maxRatioDPTable[b][i][j] = val;

              // also populate the backtrack table so we can reconstruct a solution later
              maxRatioBacktrackTable[b][i][j] = bestPrevI;

              if (debugOutput > 0)
                {
                  printf("maxRatioDPTable[%d][%d][%d]= %f   maxRatioBacktrackTable[%d][%d][%d]=%d\n",b,i,j, val,b,i,j,bestPrevI);
                }
            } 
        }
    }   
}


// Get block partition from backtrack table
// Should this be recursive? probably not...

// TODO NOTE: previously blockCount was treated as though it was the maximum index for an array of 0-indexedblocks, i.e. one less than you'd expect from the name of the variable. This has now been changed, but I should be on the look out for wierd bugs.
static void getBlockPartitionViaBacktrack(int blockCount, int finalBlockStart, int finalBlockEnd, int*** backtrackTable, int* outputStartPoints, int* outputEndpoints)
{

  // Now use backtracking to reconstruct the details of an optimal solution
  
  //int BP_optBlockStarts[blockCount]; // List of the first elements of each block in an optimal block partition
  //int BP_optBlockEnds[blockCount];	// List of the last elements of each block in an optimal block partition
  
  // The first block we store information on is actually the last block, which we know ends with the last character.
  int currentBlockStart =  finalBlockStart;  // = optFinalIForBlockCount[blockCount];
  int currentBlockEnd =   finalBlockEnd;  // = sites - 1;

 
  // Work backwards, finding the previous block and recording its details, until we reach the first block
  for (int b = blockCount - 1; b >= 0; b--)  
    {
      if (b >0)
        {
          assert (currentBlockStart <= currentBlockEnd);
          assert (currentBlockStart >= 0);
          assert (currentBlockEnd >= 0);
        }

      // append the details of the current block to their respective lists.
      outputStartPoints[b] = currentBlockStart;
      outputEndpoints[b] = currentBlockEnd;
  
      int previousBlockStart = backtrackTable[b][currentBlockStart][currentBlockEnd];
      int previousBlockEnd = currentBlockStart - 1;
      if (b >1)
        {
          //assert (currentBlockStart <= currentBlockEnd);
          assert (previousBlockStart >= 0);
          //assert (currentBlockEnd >= 0);
        }


      // previous block becomes the new current block
      currentBlockStart = previousBlockStart;
      currentBlockEnd = previousBlockEnd;
  
    }
}


static void getRestrictedBlockPartitionViaBacktrack(int blockCount, int *recordOfInformativeSites, int finalBlockStart, int finalBlockEnd, int*** backtrackTable, int* outputStartPoints, int* outputEndpoints)
{ 
  // first get the 'full' block partition
  int unrestrictedBlockStarts[blockCount]; 
  int unrestrictedBlockEnds[blockCount];
  getBlockPartitionViaBacktrack(blockCount, finalBlockStart, finalBlockEnd, backtrackTable,  unrestrictedBlockStarts, unrestrictedBlockEnds);
  // then get the informative restriction of this partition

  getInformativeRestrictionOfPartition(blockCount, recordOfInformativeSites, unrestrictedBlockStarts, unrestrictedBlockEnds, outputStartPoints, outputEndpoints);
}





// Find the optimal block partition (with respect to each problem of interest)
// for each possible number of blocks.
// Then report back the optimal block partitions,
// and a comparison to the an input 'true' block partition, if one was specfied 
static void constructAndReportOptimalPartitions(tree *tr, analdef *adef, int **BP_blockScoreArray,  char ***BP_blockTreeStringArray, int ***maxScoreDPTable, int ***maxScoreBacktrackTable, int ***totalScoreDPTable, int ***totalScoreBacktrackTable,  float ***maxRatioDPTable, int ***maxRatioBacktrackTable, float ***totalRatioDPTable, int ***totalRatioBacktrackTable)
{

  // initialise basic variables
  int debugOutput = 0;
  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;
  float desiredMaxRatio = adef->maxHomoplasyRatio;
  int desiredMaxScore = adef->maxHomoplasyScore;
  float desiredTotalRatio = adef->totalHomoplasyRatio;
  int desiredTotalScore = adef->totalHomoplasyScore;
  int verbosity = adef->verbose;// controls how much information is printed to terminal. 
                                // 0: Just return summary info
                                // 1: Also return opt block partition + homoplasy ratio for each possible block count
                                // 2: Also return homoplasy ratio + parsimonious tree for each block in each partition

  char* trueBlockPartitionString = adef->trueBlockPartitionString; 
  int truePartitionBlockCount = adef->truePartitionBlockCount;
  int* truePartitionStartPoints = adef->truePartitionStartPoints;
  int* truePartitionEndPoints = adef->truePartitionEndPoints; 

  // variables to do with solutions to maxRatio
  float overallOptMaxRatio = FLT_MAX;
  int overallOptMaxRatioBlockCount;
  int overallOptMaxRatioStarts[desiredBlockCount];
  int overallOptMaxRatioEnds[desiredBlockCount];
  //
  int goodMaxRatioMinBlockCount = -1;
  int goodMaxRatioStarts[desiredBlockCount];
  int goodMaxRatioEnds[desiredBlockCount];
  //
  int correctCountMaxRatioStarts[desiredBlockCount];
  int correctCountMaxRatioEnds[desiredBlockCount];

  // variables to do with solutions to maxScore
  int overallOptMaxScore = INT_MAX;
  int overallOptMaxScoreBlockCount;
  int overallOptMaxScoreStarts[desiredBlockCount];
  int overallOptMaxScoreEnds[desiredBlockCount];
  //
  int goodMaxScoreMinBlockCount = -1;
  int goodMaxScoreStarts[desiredBlockCount];
  int goodMaxScoreEnds[desiredBlockCount];
  //
  int correctCountMaxScoreStarts[desiredBlockCount];
  int correctCountMaxScoreEnds[desiredBlockCount];

  // variables to do with solutions to totalScore
  int overallOptTotalScore = INT_MAX;
  int overallOptTotalScoreBlockCount;
  int overallOptTotalScoreStarts[desiredBlockCount];
  int overallOptTotalScoreEnds[desiredBlockCount];
  //
  int goodTotalScoreMinBlockCount = -1;
  int goodTotalScoreStarts[desiredBlockCount];
  int goodTotalScoreEnds[desiredBlockCount];
  //
  int correctCountTotalScoreStarts[desiredBlockCount];
  int correctCountTotalScoreEnds[desiredBlockCount];

  // variables to do with solutions to totalRatio
  // this method has too many variables
  float overallOptTotalRatio = FLT_MAX;
  int overallOptTotalRatioBlockCount;
  int overallOptTotalRatioStarts[desiredBlockCount];
  int overallOptTotalRatioEnds[desiredBlockCount];
  //
  int goodTotalRatioMinBlockCount = -1;
  int goodTotalRatioStarts[desiredBlockCount];
  int goodTotalRatioEnds[desiredBlockCount];
  //
  int correctCountTotalRatioStarts[desiredBlockCount];
  int correctCountTotalRatioEnds[desiredBlockCount];


  // determine which problems to report on
  boolean solveMaxScore = FALSE;
  boolean solveTotalScore = FALSE;
  boolean solveMaxRatio = FALSE;
  boolean solveTotalRatio = FALSE;
  if ((adef->problems & 1) > 0)
    {
       solveMaxScore = TRUE;
    }
  if ((adef->problems & 2) > 0)
    {
       solveTotalScore = TRUE;
    }
  if ((adef->problems & 4) > 0)
    {
       solveMaxRatio = TRUE;
    }
  if ((adef->problems & 8) > 0)
    {
       solveTotalRatio = TRUE;
    }


  if (debugOutput > 0)
    { 
      printf("Determining optimal block partition.\n");
    }
  
  
  // For each possible block count, find (and optionally report) the best partition on that many blocks, for each problem of interest.
  // Also for each problem of interest, make a note of the best partition overall,
  // and (optionally) the partition on fewest blocks that achieves a certain desired quality.
  for(int b = 0; b < desiredBlockCount; b++)
    {
  
      if (verbosity > 0)
        {
          printBothOpen("Best %d-block partition(s):\n", b+1);
        }

      //
      // Find and report the optimum block partition in terms of max homoplasy score
      // for this number of blocks
      //
      if (solveMaxScore)
        {
          // To find the optimum max score,
          // we need to check the optimum max score for every possible final block
          int optMaxScore = INT_MAX;
          int   optMaxScoreI = -1;
          for(int i = 0; i < sites; i++)
            {
              int tempMaxScore = maxScoreDPTable[b][i][sites - 1];
              if (tempMaxScore < optMaxScore)
                {
                  optMaxScore = tempMaxScore;
                  optMaxScoreI = i;
                }
            }
          // We now have the optimum partition on (b+1) blocks has maximum homoplasy score optMaxScore 
          // and the last block starts with character optMaxScoreI.
          // We now reconstruct this partition using the backtracking table.    
          //int unrestrictedOptBlockStarts[b+1]; 
          //int unrestrictedOptBlockEnds[b+1];
          int optBlockStarts[b+1]; 
          int optBlockEnds[b+1];
          getRestrictedBlockPartitionViaBacktrack(b+1,  tr->recordOfInformativeSites, optMaxScoreI, sites - 1, maxScoreBacktrackTable, optBlockStarts, optBlockEnds);
          char* blockPartitionString = getBlockPartitionString(b + 1, optBlockStarts, optBlockEnds);


          // Report optimum partition for this block count, if requested
          if (verbosity >= 1)
            {            
              printBothOpen("  Optimum max homoplasy score:\n");
              printBothOpen("    %s\n", blockPartitionString);  
              printBothOpen("    Max homoplasy score: %d\n", optMaxScore);  


              // Give additional information about the homoplasy score, ratio and parsimonious tree of each in the partiton, if requested
              if (verbosity >= 2)
                {
                  int tempBlockStart;
                  int tempBlockEnd;
                  int blockHomoplasyScore;
                  float blockHomoplasyRatio;
                  char *blockTreeNewickString;
                  for (int h = 0; h <= b; h++)
                    {
                      tempBlockStart = optBlockStarts[h];
                      tempBlockEnd = optBlockEnds[h];
                      blockHomoplasyScore = BP_blockScoreArray[tempBlockStart][tempBlockEnd];
                      blockHomoplasyRatio = getHomoplasyRatio(tr, tempBlockStart, tempBlockEnd, blockHomoplasyScore);
                      blockTreeNewickString = BP_blockTreeStringArray[tempBlockStart][tempBlockEnd];
                      printBothOpen("      [%d---%d] Homoplasy score: %d Homoplasy Ratio: %f\n",  tempBlockStart, tempBlockEnd, blockHomoplasyScore, blockHomoplasyRatio);
                      printBothOpen("        Parsimony tree: %s\n", blockTreeNewickString);
                    }
                }            
            }

          // Record this block partition if it gets a better max score than the previous best
          if(optMaxScore < overallOptMaxScore)
            {
              overallOptMaxScore = optMaxScore;
              overallOptMaxScoreBlockCount = b;
              memcpy(overallOptMaxScoreStarts, optBlockStarts, sizeof(overallOptMaxScoreStarts));
              memcpy(overallOptMaxScoreEnds, optBlockEnds, sizeof(overallOptMaxScoreEnds));
            }

          // Record this block partition if it is the shortest one (i.e. fewest blocks) to achieve the requested score
          if (optMaxScore <= desiredMaxScore && goodMaxScoreMinBlockCount == -1)
            {
              goodMaxScoreMinBlockCount = b;
              memcpy(goodMaxScoreStarts, optBlockStarts, sizeof(goodMaxScoreStarts));
              memcpy(goodMaxScoreEnds, optBlockEnds, sizeof(goodMaxScoreEnds));
            }

          // Record this block partition if it has the same length as the "true" block partition (if one was specified)
          if (b == truePartitionBlockCount - 1)
            {
              memcpy(correctCountMaxScoreStarts, optBlockStarts, sizeof(correctCountMaxScoreStarts));
              memcpy(correctCountMaxScoreEnds, optBlockEnds, sizeof(correctCountMaxScoreEnds));
            }

          // Free memory   // 27.11.2018 Trying to fix memory leaks - MJ
          free(blockPartitionString);
        }


      //
      // Find and report the optimum block partition in terms of total homoplasy score
      // for this number of blocks
      //
      if (solveTotalScore)
        {
          // To find the optimum total score,
          // we need to check the optimum total score for every possible final block
          int optTotalScore = INT_MAX;
          int   optTotalScoreI = -1;
          for(int i = 0; i < sites; i++)
            {
              int tempTotalScore = totalScoreDPTable[b][i][sites - 1];
              if (tempTotalScore < optTotalScore)
                {
                  optTotalScore = tempTotalScore;
                  optTotalScoreI = i;
                }
            }
          // We now have the optimum partition on (b+1) blocks has total homoplasy score optTotalScore 
          // and the last block starts with character optTotalScoreI.
          // We now reconstruct this partition using the backtracking table.    
          int optBlockStarts[b+1]; 
          int optBlockEnds[b+1];
          getRestrictedBlockPartitionViaBacktrack(b+1,  tr->recordOfInformativeSites, optTotalScoreI, sites - 1, totalScoreBacktrackTable, optBlockStarts, optBlockEnds);
          char* blockPartitionString = getBlockPartitionString(b + 1, optBlockStarts, optBlockEnds);


          // Report optimum partition for this block count, if requested
          if (verbosity >= 1)
            {            
              printBothOpen("  Optimum total homoplasy score:\n");
              printBothOpen("    %s\n", blockPartitionString);  
              printBothOpen("    Total homoplasy score: %d\n", optTotalScore);  


              // Give additional information about the homoplasy score, ratio and parsimonious tree of each in the partiton, if requested
              if (verbosity >= 2)
                {
                  int tempBlockStart;
                  int tempBlockEnd;
                  int blockHomoplasyScore;
                  float blockHomoplasyRatio;
                  char *blockTreeNewickString;
                  for (int h = 0; h <= b; h++)
                    {
                      tempBlockStart = optBlockStarts[h];
                      tempBlockEnd = optBlockEnds[h];
                      blockHomoplasyScore = BP_blockScoreArray[tempBlockStart][tempBlockEnd];
                      blockHomoplasyRatio = getHomoplasyRatio(tr, tempBlockStart, tempBlockEnd, blockHomoplasyScore);
                      blockTreeNewickString = BP_blockTreeStringArray[tempBlockStart][tempBlockEnd];
                      printBothOpen("      [%d---%d] Homoplasy score: %d Homoplasy Ratio: %f\n",  tempBlockStart, tempBlockEnd, blockHomoplasyScore, blockHomoplasyRatio);
                      printBothOpen("        Parsimony tree: %s\n", blockTreeNewickString);
                    }
                }            
            }

          // Record this block partition if it gets a better total score than the previous best
          if(optTotalScore < overallOptTotalScore)
            {
              overallOptTotalScore = optTotalScore;
              overallOptTotalScoreBlockCount = b;
              memcpy(overallOptTotalScoreStarts, optBlockStarts, sizeof(overallOptTotalScoreStarts));
              memcpy(overallOptTotalScoreEnds, optBlockEnds, sizeof(overallOptTotalScoreEnds));
            }

          // Record this block partition if it is the shortest one (i.e. fewest blocks) to achieve the requested score
          if (optTotalScore <= desiredTotalScore && goodTotalScoreMinBlockCount == -1)
            {
              goodTotalScoreMinBlockCount = b;
              memcpy(goodTotalScoreStarts, optBlockStarts, sizeof(goodTotalScoreStarts));
              memcpy(goodTotalScoreEnds, optBlockEnds, sizeof(goodTotalScoreEnds));
            }

          // Record this block partition if it has the same length as the "true" block partition (if one was specified)
          if (b == truePartitionBlockCount - 1)
            {
              memcpy(correctCountTotalScoreStarts, optBlockStarts, sizeof(correctCountTotalScoreStarts));
              memcpy(correctCountTotalScoreEnds, optBlockEnds, sizeof(correctCountTotalScoreEnds));
            }
          // Free memory   // 27.11.2018 Trying to fix memory leaks - MJ
          free(blockPartitionString);
        }

       

    
      // Find and report the optimum block partition in terms of max homoplasy ratio
      // for this number of blocks
      if (solveMaxRatio)
        {
          // To find the optimum max ratio,
          // we need to check the optimum max ratio for every possible final block
          float optMaxRatio = FLT_MAX;
          // char* optMaxRatioString = (char*) malloc(sizeof(char) * 1000);
          int   optMaxRatioI = -1;
          for(int i = 0; i < sites; i++)
            {
              float tempMaxRatio = maxRatioDPTable[b][i][sites - 1];
              if (tempMaxRatio < optMaxRatio)
                {
                  optMaxRatio = tempMaxRatio;
                  optMaxRatioI = i;
                }
            }
          // We now have the optimum partition on (b+1) blocks has maximum homoplasy ratio optMaxRatio 
          // and the last block starts with character optMaxRatioI.
          // We now reconstruct this partition using the backtracking table.    
          int optBlockStarts[b+1]; 
          int optBlockEnds[b+1];
          getRestrictedBlockPartitionViaBacktrack(b+1,  tr->recordOfInformativeSites, optMaxRatioI, sites - 1, maxRatioBacktrackTable, optBlockStarts, optBlockEnds);
          char* blockPartitionString = getBlockPartitionString(b + 1, optBlockStarts, optBlockEnds);


          // Report optimum partition for this block count, if requested
          if (verbosity >= 1)
            {            
              printBothOpen("  Optimum max homoplasy ratio:\n");
              printBothOpen("    %s\n", blockPartitionString);  
              printBothOpen("    Max homoplasy ratio: %f\n", optMaxRatio);  


              // Give additional information about the homoplasy score, ratio and parsimonious tree of each in the partiton, if requested
              if (verbosity >= 2)
                {
                  int tempBlockStart;
                  int tempBlockEnd;
                  int blockHomoplasyScore;
                  float blockHomoplasyRatio;
                  char *blockTreeNewickString;
                  for (int h = 0; h <= b; h++)
                    {
                      tempBlockStart = optBlockStarts[h];
                      tempBlockEnd = optBlockEnds[h];
                      blockHomoplasyScore = BP_blockScoreArray[tempBlockStart][tempBlockEnd];
                      blockHomoplasyRatio = getHomoplasyRatio(tr, tempBlockStart, tempBlockEnd, blockHomoplasyScore);
                      blockTreeNewickString = BP_blockTreeStringArray[tempBlockStart][tempBlockEnd];
                      printBothOpen("      [%d---%d] Homoplasy score: %d Homoplasy Ratio: %f\n",  tempBlockStart, tempBlockEnd, blockHomoplasyScore, blockHomoplasyRatio);
                      printBothOpen("        Parsimony tree: %s\n", blockTreeNewickString);
                    }
                }            
            }

          // Record this block partition if it gets a better max ratio than the previous best
          if(optMaxRatio < overallOptMaxRatio)
            {
              overallOptMaxRatio = optMaxRatio;
              overallOptMaxRatioBlockCount = b;
              memcpy(overallOptMaxRatioStarts, optBlockStarts, sizeof(overallOptMaxRatioStarts));
              memcpy(overallOptMaxRatioEnds, optBlockEnds, sizeof(overallOptMaxRatioEnds));
            }

          // Record this block partition if it is the shortest one (i.e. fewest blocks) to achieve the requested ratio
          if (optMaxRatio <= desiredMaxRatio && goodMaxRatioMinBlockCount == -1)
            {
              goodMaxRatioMinBlockCount = b;
              memcpy(goodMaxRatioStarts, optBlockStarts, sizeof(goodMaxRatioStarts));
              memcpy(goodMaxRatioEnds, optBlockEnds, sizeof(goodMaxRatioEnds));
            }

          // Record this block partition if it has the same length as the "true" block partition (if one was specified)
          if (b == truePartitionBlockCount - 1)
            {
              memcpy(correctCountMaxRatioStarts, optBlockStarts, sizeof(correctCountMaxRatioStarts));
              memcpy(correctCountMaxRatioEnds, optBlockEnds, sizeof(correctCountMaxRatioEnds));
            }
          // Free memory   // 27.11.2018 Trying to fix memory leaks - MJ
          free(blockPartitionString);
        }

      //
      // Find and report the optimum block partition in terms of total homoplasy ratio
      // for this number of blocks
      //
      if (solveTotalRatio)
        {
          // To find the optimum total ratio,
          // we need to check the optimum total ratio for every possible final block
          float optTotalRatio = FLT_MAX;
          int   optTotalRatioI = -1;
          for(int i = 0; i < sites; i++)
            {
              float tempTotalRatio = totalRatioDPTable[b][i][sites - 1];
              if (tempTotalRatio < optTotalRatio)
                {
                  optTotalRatio = tempTotalRatio;
                  optTotalRatioI = i;
                }
            }
          // We now have the optimum partition on (b+1) blocks has total homoplasy ratio optTotalRatio 
          // and the last block starts with character optTotalRatioI.
          // We now reconstruct this partition using the backtracking table.    
          int optBlockStarts[b+1]; 
          int optBlockEnds[b+1];
          getRestrictedBlockPartitionViaBacktrack(b+1,  tr->recordOfInformativeSites, optTotalRatioI, sites - 1, totalRatioBacktrackTable, optBlockStarts, optBlockEnds);
          char* blockPartitionString = getBlockPartitionString(b + 1, optBlockStarts, optBlockEnds);


          // Report optimum partition for this block count, if requested
          if (verbosity >= 1)
            {            
              printBothOpen("  Optimum total homoplasy ratio:\n");
              printBothOpen("    %s\n", blockPartitionString);  
              printBothOpen("    Total homoplasy ratio: %f\n", optTotalRatio);  

              // Give additional information about the homoplasy score, ratio and parsimonious tree of each in the partiton, if requested
              if (verbosity >= 2)
                {
                  int tempBlockStart;
                  int tempBlockEnd;
                  int blockHomoplasyScore;
                  float blockHomoplasyRatio;
                  char *blockTreeNewickString;
                  for (int h = 0; h <= b; h++)
                    {
                      tempBlockStart = optBlockStarts[h];
                      tempBlockEnd = optBlockEnds[h];
                      blockHomoplasyScore = BP_blockScoreArray[tempBlockStart][tempBlockEnd];
                      blockHomoplasyRatio = getHomoplasyRatio(tr, tempBlockStart, tempBlockEnd, blockHomoplasyScore);
                      blockTreeNewickString = BP_blockTreeStringArray[tempBlockStart][tempBlockEnd];
                      printBothOpen("      [%d---%d] Homoplasy score: %d Homoplasy Ratio: %f\n",  tempBlockStart, tempBlockEnd, blockHomoplasyScore, blockHomoplasyRatio);
                      printBothOpen("        Parsimony tree: %s\n", blockTreeNewickString);
                    }
                }            
            }

          // Record this block partition if it gets a better total ratio than the previous best
          if(optTotalRatio < overallOptTotalRatio)
            {
              overallOptTotalRatio = optTotalRatio;
              overallOptTotalRatioBlockCount = b;
              memcpy(overallOptTotalRatioStarts, optBlockStarts, sizeof(overallOptTotalRatioStarts));
              memcpy(overallOptTotalRatioEnds, optBlockEnds, sizeof(overallOptTotalRatioEnds));
            }

          // Record this block partition if it is the shortest one (i.e. fewest blocks) to achieve the requested total ratio
          if (optTotalRatio <= desiredTotalRatio && goodTotalRatioMinBlockCount == -1)
            {
              goodTotalRatioMinBlockCount = b;
              memcpy(goodTotalRatioStarts, optBlockStarts, sizeof(goodTotalRatioStarts));
              memcpy(goodTotalRatioEnds, optBlockEnds, sizeof(goodTotalRatioEnds));
            }

          // Record this block partition if it has the same length as the "true" block partition (if one was specified)
          if (b == truePartitionBlockCount - 1)
            {
              memcpy(correctCountTotalRatioStarts, optBlockStarts, sizeof(correctCountTotalRatioStarts));
              memcpy(correctCountTotalRatioEnds, optBlockEnds, sizeof(correctCountTotalRatioEnds));
            }
          // Free memory   // 27.11.2018 Trying to fix memory leaks - MJ
          free(blockPartitionString);
        }





    }


  // Print a summary of the interesting stuff

  printBothOpen("---------------------------------\n");
  printBothOpen("Summary: \n");

  if (solveMaxScore)
    {
      printBothOpen("  HOMOPLASY PER BLOCK: \n");
      printBothOpen("    Optimal max homoplasy score: \n");
      char *str = getBlockPartitionString(overallOptMaxScoreBlockCount + 1, overallOptMaxScoreStarts, overallOptMaxScoreEnds);
      printBothOpen("      %s\n", str);
      free(str);
      printBothOpen("      Max homoplasy score %d, Block count %d\n", overallOptMaxScore ,  overallOptMaxScoreBlockCount+1);
     // Report the shortest block partition acheiving desiredMaxHomoplasyScore, if desiredMaxHomoplasyScore was specified.
      if (desiredMaxScore != INT_MAX)
        {
          printBothOpen("    Fewest blocks achieving homoplasy score at most %d: \n", desiredMaxScore );
          if (goodMaxScoreMinBlockCount != -1)
            {
              char *str = getBlockPartitionString(goodMaxScoreMinBlockCount + 1, goodMaxScoreStarts, goodMaxScoreEnds);
              printBothOpen("    %s\n", str);
              free(str);
              int tempScore = getMaxHomoplasyScoreFromPartition(goodMaxScoreMinBlockCount+ 1, goodMaxScoreStarts, goodMaxScoreEnds, BP_blockScoreArray);
              printBothOpen("      Max homoplasy score %d, Block count %d\n", tempScore, goodMaxScoreMinBlockCount + 1);
            }
          else
            {
              printBothOpen("      No solution found.\n");
            }
        } 
    }
 
  if (solveTotalScore)
    {
      printBothOpen("  TOTAL HOMOPLASY: \n");
      printBothOpen("    Optimal total homoplasy score: \n");
      char *str = getBlockPartitionString(overallOptTotalScoreBlockCount + 1, overallOptTotalScoreStarts, overallOptTotalScoreEnds);
      printBothOpen("      %s\n", str);
      free(str);
      printBothOpen("      Total homoplasy score %d, Block count %d\n", overallOptTotalScore ,  overallOptTotalScoreBlockCount+1);
     // Report the shortest block partition acheiving desiredTotalHomoplasyScore, if desiredTotalHomoplasyScore was specified.
      if (desiredTotalScore != INT_MAX)
        {
          printBothOpen("    Fewest blocks achieving homoplasy score at most %d: \n", desiredTotalScore );
          if (goodTotalScoreMinBlockCount != -1)
            {
              char *str = getBlockPartitionString(goodTotalScoreMinBlockCount + 1, goodTotalScoreStarts, goodTotalScoreEnds);
              printBothOpen("    %s\n", str);
              free(str);
              int tempScore = getTotalHomoplasyScoreFromPartition(goodTotalScoreMinBlockCount+ 1, goodTotalScoreStarts, goodTotalScoreEnds, BP_blockScoreArray);
              printBothOpen("      Total homoplasy score %d, Block count %d\n", tempScore, goodTotalScoreMinBlockCount + 1);
            }
          else
            {
              printBothOpen("      No solution found.\n");
            }
        } 
    }


  if (solveMaxRatio)
    {
      printBothOpen("  HOMOPLASY RATE PER BLOCK: \n");
      printBothOpen("    Optimal max homoplasy ratio: \n");
      char *str =  getBlockPartitionString(overallOptMaxRatioBlockCount + 1, overallOptMaxRatioStarts, overallOptMaxRatioEnds);
      printBothOpen("      %s\n", str);
      free(str);
      printBothOpen("      Max homoplasy ratio %f, Block count %d\n", overallOptMaxRatio ,  overallOptMaxRatioBlockCount+1);
     // Report the shortest block partition acheiving desiredMaxHomoplasyRatio, if desiredMaxHomoplasyRatio was specified.
      if (desiredMaxRatio != FLT_MAX)
        {
          printBothOpen("    Fewest blocks achieving homoplasy ratio at most %f: \n", desiredMaxRatio );
          if (goodMaxRatioMinBlockCount != -1)
            {
              char *str = getBlockPartitionString(goodMaxRatioMinBlockCount + 1, goodMaxRatioStarts, goodMaxRatioEnds);
              printBothOpen("    %s\n", str);
              free(str);
              float tempRatio = getMaxHomoplasyRatioFromPartition(tr,goodMaxRatioMinBlockCount+ 1, goodMaxRatioStarts, goodMaxRatioEnds, BP_blockScoreArray);
              printBothOpen("      Max homoplasy ratio %f, Block count %d\n", tempRatio, goodMaxRatioMinBlockCount + 1);
            }
          else
            {
              printBothOpen("      No solution found.\n");
            }
        } 
    }


  if (solveTotalRatio)
    {
      printBothOpen("  TOTAL HOMOPLASY RATE: \n");
      printBothOpen("    Optimal total homoplasy ratio: \n");
      char *str = getBlockPartitionString(overallOptTotalRatioBlockCount + 1, overallOptTotalRatioStarts, overallOptTotalRatioEnds);
      printBothOpen("      %s\n", str);
      free(str);
      printBothOpen("      Total homoplasy ratio %f, Block count %d\n", overallOptTotalRatio ,  overallOptTotalRatioBlockCount+1);

     // Report the shortest block partition acheiving desiredTotalHomoplasyRatio, if desiredTotalHomoplasyRatio was specified.
      if (desiredTotalRatio != FLT_MAX)
        {
          printBothOpen("    Fewest blocks achieving homoplasy ratio at most %f: \n", desiredTotalRatio);
          if (goodTotalRatioMinBlockCount != -1)
            {
              char *str = getBlockPartitionString(goodTotalRatioMinBlockCount + 1, goodTotalRatioStarts, goodTotalRatioEnds);
              printBothOpen("    %s\n", str);
              free(str);
              float tempRatio = getTotalHomoplasyRatioFromPartition(tr, goodTotalRatioMinBlockCount+ 1, goodTotalRatioStarts, goodTotalRatioEnds, BP_blockScoreArray);
              printBothOpen("      Total homoplasy ratio %f, Block count %d\n", tempRatio, goodTotalRatioMinBlockCount + 1);
            }
          else
            {
              printBothOpen("      No solution found.\n");
            }
        } 
    }
   
   


  // TODO tidy up this bit, then start repeating for all the different problem types, and finally do implementation??
  // If we were given info on the "true" block partition, report back on how our partition (for same number of blocks) compares
  
  if (trueBlockPartitionString[0] != '\0')
    {
      printBothOpen("Comparision to input partition: \n");
      printBothOpen("  INPUT BLOCK PARTITION: \n");
      char *str = getBlockPartitionString(truePartitionBlockCount, truePartitionStartPoints,  truePartitionEndPoints);
      printBothOpen("    %s\n", str);
      free(str);
      int trueMaxScore = getMaxHomoplasyScoreFromPartition(truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, BP_blockScoreArray);
      int trueTotalScore = getTotalHomoplasyScoreFromPartition(truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, BP_blockScoreArray);
      float trueMaxRatio = getMaxHomoplasyRatioFromPartition(tr, truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, BP_blockScoreArray);
      float trueTotalRatio = getTotalHomoplasyRatioFromPartition(tr, truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, BP_blockScoreArray);
      printBothOpen("    Block count %d, Max homoplasy score %d, Total homoplasy score %d, Max homoplasy ratio %f, Total homoplasy ratio %f, \n", truePartitionBlockCount, trueMaxScore, trueTotalScore, trueMaxRatio, trueTotalRatio);
      printBothOpen("  Optimal %d-block partition(s): \n", truePartitionBlockCount);
      if (truePartitionBlockCount <= desiredBlockCount)
        {

          if (solveMaxScore)
            {
              char *str = getBlockPartitionString(truePartitionBlockCount, correctCountMaxScoreStarts,  correctCountMaxScoreEnds);
              printBothOpen("    MAX SCORE: %s\n", str);
              free(str);
              int correctCountMaxScore = getMaxHomoplasyScoreFromPartition(truePartitionBlockCount, correctCountMaxScoreStarts, correctCountMaxScoreEnds, BP_blockScoreArray);
              float correctCountMaxScoreAverageBoundaryError =  averageInternalBoundaryError(tr->recordOfInformativeSites, truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, correctCountMaxScoreStarts, correctCountMaxScoreEnds);
              printBothOpen("      Max homoplasy score %d, Average internal boundary error %f\n", correctCountMaxScore, correctCountMaxScoreAverageBoundaryError);
            }

          if (solveTotalScore)
            {
              char *str = getBlockPartitionString(truePartitionBlockCount, correctCountTotalScoreStarts,  correctCountTotalScoreEnds);
              printBothOpen("    TOTAL SCORE: %s\n", str);
              free(str);
              int correctCountTotalScore = getTotalHomoplasyScoreFromPartition(truePartitionBlockCount, correctCountTotalScoreStarts, correctCountTotalScoreEnds, BP_blockScoreArray);
              float correctCountTotalScoreAverageBoundaryError =  averageInternalBoundaryError(tr->recordOfInformativeSites, truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, correctCountTotalScoreStarts, correctCountTotalScoreEnds);
              printBothOpen("      Total homoplasy score %d, Average internal boundary error %f\n", correctCountTotalScore, correctCountTotalScoreAverageBoundaryError);
            }

          if (solveMaxRatio)
            {
              char *str = getBlockPartitionString(truePartitionBlockCount, correctCountMaxRatioStarts,  correctCountMaxRatioEnds);
              printBothOpen("    MAX RATIO: %s\n", str);
              free(str);
              float correctCountMaxRatio = getMaxHomoplasyRatioFromPartition(tr, truePartitionBlockCount, correctCountMaxRatioStarts, correctCountMaxRatioEnds, BP_blockScoreArray);
              float correctCountMaxRatioAverageBoundaryError =  averageInternalBoundaryError(tr->recordOfInformativeSites, truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, correctCountMaxRatioStarts, correctCountMaxRatioEnds);
              printBothOpen("      Max homoplasy ratio %f, Average internal boundary error %f\n", correctCountMaxRatio, correctCountMaxRatioAverageBoundaryError);
            }

          if (solveTotalRatio)
            {
              char *str = getBlockPartitionString(truePartitionBlockCount, correctCountTotalRatioStarts,  correctCountTotalRatioEnds);
              printBothOpen("    TOTAL RATIO: %s\n", str);
              free(str);
              float correctCountTotalRatio = getTotalHomoplasyRatioFromPartition(tr, truePartitionBlockCount, correctCountTotalRatioStarts, correctCountTotalRatioEnds, BP_blockScoreArray);
              float correctCountTotalRatioAverageBoundaryError =  averageInternalBoundaryError(tr->recordOfInformativeSites, truePartitionBlockCount, truePartitionStartPoints, truePartitionEndPoints, correctCountTotalRatioStarts, correctCountTotalRatioEnds);
              printBothOpen("      Total homoplasy ratio %f, Average internal boundary error %f\n", correctCountTotalRatio, correctCountTotalRatioAverageBoundaryError);
            }

        }
      else
        {
          printBothOpen("    No solution found (input block partition has more blocks than specified maximum)\n");
        }
    }
}







static void getOptBlockPartition(tree *tr, analdef *adef) 
{

  // Step 1. Initialise basic variables and declare arrays

  // initialise basic variables
  int debugOutput = 0;
  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;
  float desiredHomoplasyRatio = adef->maxHomoplasyRatio;

  // determine which problems are to be solved
  boolean solveMaxScore = FALSE;
  boolean solveTotalScore = FALSE;
  boolean solveMaxRatio = FALSE;
  boolean solveTotalRatio = FALSE;
  if ((adef->problems & 1) > 0)
    {
       solveMaxScore = TRUE;
    }
  if ((adef->problems & 2) > 0)
    {
       solveTotalScore = TRUE;
    }
  if ((adef->problems & 4) > 0)
    {
       solveMaxRatio = TRUE;
    }
  if ((adef->problems & 8) > 0)
    {
       solveTotalRatio = TRUE;
    }

  // declare arrays that may be used in the dynnamic programming
  int BP_minParsimonyScorePerSiteArray [sites];  
  int **BP_blockScoreArray = 0;
  char ***BP_blockTreeStringArray = 0;
  int *recordOfInformativeSites = 0;


  int ***maxScoreDPTable = 0;
  int ***maxScoreBacktrackTable = 0;

  int ***totalScoreDPTable = 0;
  int ***totalScoreBacktrackTable = 0;

  float ***maxRatioDPTable = 0;
  int ***maxRatioBacktrackTable = 0;

  float ***totalRatioDPTable = 0;
  int ***totalRatioBacktrackTable = 0;

  // float ***BP_DPHomoplasyRatioLookupTable = 0;
  // int ***BP_backtrackTable = 0;

  // announce start of program and basic variables
  printBothOpen("---------------------------------\n"),
  printBothOpen("Input file: %s | Output file: %s | Maximum number of blocks: %d | Randomization seed: %d\n", seq_file, partitionInfoFileName, adef->numberOfBlocks, adef->parsimonySeed);  
  printBothOpen("Problems to solve: | ");
  if (solveMaxScore)
    {
      printBothOpen("Max Homoplasy Score | ");
    }
  if (solveTotalScore)
    {
      printBothOpen("Total Homoplasy Score | ");
    }
  if (solveMaxRatio)
    {
      printBothOpen("Max Homoplasy Ratio | ");
    }
  if (solveTotalRatio)
    {
      printBothOpen("Total Homoplasy Ratio | ");
    }
  printBothOpen("\n");
  printBothOpen("---------------------------------\n");


  // Step 2. Set up arrays for the algorithms we want to run
  // (This reserves a lot of memory, so we only do it for the arrays we'll need to use)
  BP_blockScoreArray = get2dIntArray(sites, sites);
  BP_blockTreeStringArray = get3dCharArray(sites, sites, tr->treeStringLength);
  tr->recordOfInformativeSites = (int *)malloc(sizeof(int) * (size_t)tr->rdta->sites);

  if (solveMaxScore)
    {
      maxScoreDPTable = get3dIntArray(desiredBlockCount, sites, sites);
      maxScoreBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  if (solveTotalScore)
    {
      totalScoreDPTable = get3dIntArray(desiredBlockCount, sites, sites);
      totalScoreBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  if (solveMaxRatio)
    {
      maxRatioDPTable = get3dFloatArray(desiredBlockCount, sites, sites);
      maxRatioBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  if (solveTotalRatio)
    {
      totalRatioDPTable = get3dFloatArray(desiredBlockCount, sites, sites);
      totalRatioBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  

  // Step 3. Calculate the optimum homoplasy score and parsimony tree for each block

  // First populate BP_minParsimonyScorePerSiteArray, which stores the theoretical minimum parsimony score for each individal site
  for (int i = 0; i < sites; i++)
    {
      BP_minParsimonyScorePerSiteArray[i] = getMinParsimonyScoreForSite(tr, i);
    }
  // Add min Parsimony data to tr, so that we can look it up easily later on.
  tr->minParsimonyPerSite = BP_minParsimonyScorePerSiteArray;

  // Store a record of which sites are 'informative' and which are not. This will get used in TODO FILL IN
  determineUninformativeSites(tr, tr->recordOfInformativeSites);

  // test code to make sure recordOfInformativeSites is working correctly
  char informativeSitesString[3*sites];
  strcpy(informativeSitesString, " ");
  char tempstr[1000];
  for(int i = 0; i < sites; i++)
    {
      sprintf(tempstr, "%d",  tr->recordOfInformativeSites[i]); // make formatted string discribing informative status of site i
      strcat(informativeSitesString, tempstr);  // append informative status to output string
      if (i < sites - 1)
        {
          strcat(informativeSitesString, ", "); // append a comma if there are more sites to come
        } 
    }
  printBothOpen("Informative sites:\n");
  printBothOpen("    %s\n", informativeSitesString);  


  // parse and store data on the True Block Partition  
  //int sites = rdta->sites;
  if (adef->trueBlockPartitionString[0] != '\0')
    {
      char tempBlockPartitionString[1024] = ""; 
      strcpy(tempBlockPartitionString,adef->trueBlockPartitionString);
      getTrueBlockPartitionData(tempBlockPartitionString, &(adef->truePartitionBlockCount), &(adef->truePartitionStartPoints), &(adef->truePartitionEndPoints), sites, tr->recordOfInformativeSites);
    }



  // Then populate BP_blockScoreArray, which stores the optimum homoplasy score every contiguous block of characters
  // Also populate BP_blockTreeStringArray, which stores an optimum treee for every contiguous block of characters
  // (This is the part of the code that we expect to take the most time)
  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++) 
        {
          // We only store info on blocks [i--j] where i <= j
          if (i > j)

            BP_blockScoreArray[i][j] = INT_MAX;
          else
            // Calculate and store the optimum homoplasy score of block [i---j]
            // NOTE 11 MARCH 2019: We are now trying to focus on the 'informative restriction' of each block,
            // So really what we want is the optimum homoplasy score of the informative restriction of [i---j].
            // But this is the same as the optimum homoplasy score of [i---j] (since uninformative sites make no difference to homoplasy score)
            // We will have to be a bit more careful when it comes to calculating the homoplasy *ratio* of a block
            BP_blockScoreArray[i][j] = getOptHomoplasyOfContiguousBlock(tr, adef, i, j);
            // Also populate BP_blockTreeStringArray with the tree we just found
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



  // Step 4. Populate the Dynamic Programming and Backtracking tables for each algorithm we care about
  // (This is the bit where the "mathsy" stuff happens)
  if (solveMaxScore)
    {
      doMaxScoreDP(tr, adef, BP_blockScoreArray,  maxScoreDPTable, maxScoreBacktrackTable);
    }
  if (solveTotalScore)
    {
      doTotalScoreDP(tr, adef, BP_blockScoreArray,  totalScoreDPTable, totalScoreBacktrackTable);
    }
  if (solveMaxRatio)
    {
      doMaxRatioDP(tr, adef, BP_blockScoreArray,  maxRatioDPTable, maxRatioBacktrackTable);
    }
  if (solveTotalRatio)
    {
      doTotalRatioDP(tr, adef, BP_blockScoreArray,  totalRatioDPTable, totalRatioBacktrackTable);
    }
  

  // Step 5. Find optimal solutions and report them
  constructAndReportOptimalPartitions(tr, adef, BP_blockScoreArray, BP_blockTreeStringArray, maxScoreDPTable, maxScoreBacktrackTable, totalScoreDPTable, totalScoreBacktrackTable,  maxRatioDPTable, maxRatioBacktrackTable, totalRatioDPTable, totalRatioBacktrackTable);

 // constructAndReportOptimalPartitions(tr, adef, BP_blockScoreArray, BP_blockTreeStringArray, BP_DPHomoplasyRatioLookupTable, BP_backtrackTable);


  // Free memory before ending program   // 27.11.2018 Trying to fix memory leaks - MJ
  free3dCharArray(BP_blockTreeStringArray);  
  free2dIntArray(BP_blockScoreArray);

  if (solveMaxScore)
    {
      free3dIntArray(maxScoreDPTable);
      free3dIntArray(maxScoreBacktrackTable);
      // maxScoreDPTable = get3dIntArray(desiredBlockCount, sites, sites);
      // maxScoreBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  if (solveTotalScore)
    {
      free3dIntArray(totalScoreDPTable);
      free3dIntArray(totalScoreBacktrackTable);
      // totalScoreDPTable = get3dIntArray(desiredBlockCount, sites, sites);
      // totalScoreBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  if (solveMaxRatio)
    {
      free3dFloatArray(maxRatioDPTable);
      free3dIntArray(maxRatioBacktrackTable);
      // maxRatioDPTable = get3dFloatArray(desiredBlockCount, sites, sites);
      // maxRatioBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  if (solveTotalRatio)
    {
      free3dFloatArray(totalRatioDPTable);
      free3dIntArray(totalRatioBacktrackTable);
      // totalRatioDPTable = get3dFloatArray(desiredBlockCount, sites, sites);
      // totalRatioBacktrackTable = get3dIntArray(desiredBlockCount, sites, sites);
    }
  
  free(tr->recordOfInformativeSites); 
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


  // Free memory   // 27.11.2018 Trying to fix memory leaks - MJ
  // Update: actually this leaves a lot of indirect leaks unfixed. If I'm going to do this I shuld do it properly; for now I'm leaving it as is. - MJ
  //free(adef);
  //free(rdta);
  //free(tr);

  return 0;
}


