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
  randomFileName[1024] = "";


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
  FILE *f = myfopen(infoFileName, "ab");

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
    //printf("creating 2d int array... ");
    int **arr;
    arr  = (int **)malloc(sizeof(int *) * d1);
    arr[0] = (int *)malloc(sizeof(int) * d2 * d1);
    int i;  
    for(i = 0; i < d1; i++)
        arr[i] = (*arr + d2 * i);
    //printf("done \n");
   return arr;

  }

// Returns a 3d int array of dimensions [d1][d2][d3], and allocates memory for it
static int ***get3dIntArray(int d1, int d2, int d3)
  {
    //printf("creating 3d int array... ");
    int b,i;
    int ***arr;
    arr  = (int ***)malloc(sizeof(int **) * d1);
    arr[0] = (int **)malloc(sizeof(int *) * d2 * d1);
    arr[0][0] = (int *)malloc(sizeof(int) * d3 * d2 * d1);
 
    for(b = 0; b < d1; b++)
      {
        arr[b] = (*arr + d2 * b);
        for(i = 0; i < d2; i++)
          {
            arr[b][i] = (**arr + d3 * i);
          }
      }
    //printf("done \n");
    return arr;
  }


// Returns a 3d int array of dimensions [d1][d2][d3], and allocates memory for it
static float ***get3dFloatArray(int d1, int d2, int d3)
  {
    //printf("creating 3d float array... ");
    int b,i;
    float ***arr;
    arr  = (float ***)malloc(sizeof(float **) * d1);
    arr[0] = (float **)malloc(sizeof(float *) * d2 * d1);
    arr[0][0] = (float *)malloc(sizeof(float) * d3 * d2 * d1);
 
    for(b = 0; b < d1; b++)
      {
        arr[b] = (*arr + d2 * b);
        for(i = 0; i < d2; i++)
          {
            arr[b][i] = (**arr + d3 * i);
            for(int j = 0; j < d3; j++){
            	arr[b][i][j]=FLT_MAX;	
            }
            
          }
      }
    //printf("done \n");
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
  
  adef->parsimonySeed=1;
  adef->numberOfBlocks=2;
  
  while(!bad_opt && ((c = mygetopt(argc,argv,"p:n:s:t:N:b:", &optind, &optarg))!=-1))
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





static float getHomoplasyRatio(int i, int j, int homoplasyValue) 
{
  assert  (i <= j);  // program breaks if called with i greater than j
  return ((float) homoplasyValue) / ((float)(j-i+1));

}






static int getOptHomoplasyOfContiguousBlock(tree *tr, analdef *adef, int startCharacter, int endCharacter)
{

  assert  (startCharacter <= endCharacter); 

  int parsimonyScore = makeParsimonyTreeFastDNA(tr, adef, startCharacter, endCharacter);
  // TODO: makeParsimonyTreeFastDNA changes the tree, right?
  // Which means that after the first time we call makeParsimonyTreeFastDNA, subsequent calls won't be starting with a 'clean' tree.
  // This might be useful later on, but to begin with we'd probably like to use a fresh tree (with the same ->rdata )
  
  //int minPossibleParsimonyScore = 3*(endCharacter - startCharacter + 1);
  // TODO calculate actual minimum homoplasy score  
  //(the above calculation assumes that A,C,G,T appear in each character, and so the min parsimony score of each character is 3)
  int minPossibleParsimonyScore = 0;   // simplest way to ensure no negative homoplasy ratios while the min parsimony score is unknown

  int returnVal = parsimonyScore - minPossibleParsimonyScore;
  return  returnVal;
  //return 0;
}



static void getOptBlockPartition(tree *tr, analdef *adef) 
{

  int debugOutput = 0;

  int numsp = tr->rdta->numsp;  // number of species
  int sites = tr->rdta->sites;  // number of sites
  int desiredBlockCount = adef->numberOfBlocks;


  // BP_blockScoreArray[i][j]  stores the optimum homoplasy score of the block on characters i to j (inclusive)
  //int **BP_blockScoreArray  = get2dIntArray(sites, sites);

  int BP_blockScoreArray[sites][sites];//  = get2dIntArray(sites, sites);


  int impossiblyHighHomoplasyScore = numsp*sites;

// Entry BP_DPHomoplasyRatioLookupTable[b][i][j] for i<= j is the minimum value of 
// the homoplasy ratio of any block partition P of character 0 to j, 
// for which the last block is from i to j and there are b+1 blocks in total;
// where "homoplasy ratio of P" means the maximum homoplasy ratio of 
// any block in P. A value of Float.POSITIVE_INFINITY means there is no solution 
// for this set of values (e.g. if i > j)
  //float ***BP_DPHomoplasyRatioLookupTable = get3dFloatArray(desiredBlockCount, sites, sites);

   float BP_DPHomoplasyRatioLookupTable[desiredBlockCount][sites][sites];

// Used to reconstruct an optimal solution.
// BP_backtrackTable[b][i][j]  is the value i' such that an 
// optimal block partition with b+1 blocks and last block [i,j]
// has [i',j-1] as its second-to-last block.
  //int ***BP_backtrackTable = get3dIntArray(desiredBlockCount, sites, sites);

   int BP_backtrackTable[desiredBlockCount][sites][sites];

  //printf("populating BP_blockScoreArray\n");
  // Populate BP_blockScoreArray, which stores the optimum homoplasy score every contiguous block of characters
  // Todo there is stuff to do with declaring data structures and assigning memory space and stuff like that that I am super not doing at the moment.



  for (int i = 0; i < sites; i++)   
    {
      for (int j = 0; j < sites; j++) 
        {
          if (i > j)
            BP_blockScoreArray[i][j] = INT_MAX;
          else
           BP_blockScoreArray[i][j] = getOptHomoplasyOfContiguousBlock(tr, adef, i, j);

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


  printf("determining optimal block partition.\n");

  // Given DPHomoplasyRatioLookupTable, it remains to find the value corresponding a block partition of the full set of characters (i.e. j = characterCount-1) that has minimum homoplasy ratio.
  float currentOptimum = FLT_MAX;
  // Find the minimum value over all entries DPHomoplasyRatioLookupTable[b][i][j] with j = characterCount-1.
  // (Note that b may be smaller than desiredBlockCount-1, if splitting into fewer than desiredBlockcount blocks 
  // actually gives better homoplasy ratio.
  int optFinalI = -1;
  int optFinalB = -1;
  for(int b = 0; b < desiredBlockCount; b++)
    {
      for(int i = 0; i < sites; i++)
        {
          float tempVal = BP_DPHomoplasyRatioLookupTable[b][i][sites-1];
          if (tempVal < currentOptimum )
            {
              currentOptimum = tempVal;
              optFinalI = i;
              optFinalB = b;
            }
        }
    }
  // Store the minimum value found
  float BP_optHomoplasyRatio = currentOptimum;
  int optBlockCount = optFinalB + 1;
  // Now use backtracking to reconstruct the details of an optimal solution



  int BP_optBlockStarts[optBlockCount]; // List of the first elements of each block in an optimal block partition
  int BP_optBlockEnds[optBlockCount];	// List of the last elements of each block in an optimal block partition
  //int BP_optTrees[optFinalB + 1];	// List of (the indices of?????) optimal trees for each block in an optimal block partition

  // The first block we store information on is actually the last block, which we know ends with the last character.
  int currentBlockStart = optFinalI;
  int currentBlockEnd = sites - 1;
  //int currentBlockCount = optFinalB;


  // TEMP COMMENTING OUT TO GET AROUND ANOTHER BUG
 
  // Work backwards, finding the previous block and recording its details, until we reach the first block
  for (int b = optBlockCount-1; b >= 0; b--)  
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


  // Tada! We have found the optimal block partition. Hopefully. Time to report it!

  printf("Optimal block partition: %d blocks", optFinalB+1);


  for (int h = 0; h <= optFinalB; h++)
    {

      int tempBlockStart = BP_optBlockStarts[h];
      int tempBlockEnd = BP_optBlockEnds[h];
      printf("[");
      printf("%d", tempBlockStart); 
      printf("---");
      printf("%d", tempBlockEnd);
      printf("]");
      if (h < optFinalB)
        {
          printf(", ");
        } 
      else 
        {
          printf("\n");
        }
    }


  printf("Maximum homoplasy ratio: %f\n", BP_optHomoplasyRatio);
  

  // TODO: report back more details, like the homoplasy score / optimal tree for each block?


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


