/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id: util.c,v 1.4 2003/04/13 04:45:12 karypis Exp $
 */

#include <metislib.h>


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  abort();
}



/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *imalloc(size_t n, char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)GKmalloc(sizeof(int)*n, msg);
}


/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *idxmalloc(size_t n, char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)GKmalloc(sizeof(idxtype)*n, msg);
}


/*************************************************************************
* The following function allocates an array of float 
**************************************************************************/
float *fmalloc(size_t n, char *msg)
{
  if (n == 0)
    return NULL;

  return (float *)GKmalloc(sizeof(float)*n, msg);
}


/*************************************************************************
* The following function allocates an array of float 
**************************************************************************/
double *dmalloc(size_t n, char *msg)
{
  if (n == 0)
    return NULL;

  return (double *)GKmalloc(sizeof(double)*n, msg);
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *ismalloc(size_t n, idxtype ival, char *msg)
{
  if (n == 0)
    return NULL;

  return iset(n, ival, (idxtype *)GKmalloc(sizeof(int)*n, msg));
}

/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
double *dsmalloc(size_t n, double ival, char *msg)
{
  if (n == 0)
    return NULL;

  return dset(n, ival, (double *)GKmalloc(sizeof(double)*n, msg));
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *idxsmalloc(size_t n, idxtype ival, char *msg)
{
  if (n == 0)
    return NULL;

  return idxset(n, ival, (idxtype *)GKmalloc(sizeof(idxtype)*n, msg));
}


/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(size_t nbytes, char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL) 
    errexit("***Memory allocation failed for %s. Requested size: %zd bytes", msg, nbytes);

  return ptr;
}


/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void GKfree(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    free(*ptr1);
  *ptr1 = NULL;

  va_start(plist, ptr1);

  /* while ((int)(ptr = va_arg(plist, void **)) != -1) { */
  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
      free(*ptr);
    *ptr = NULL;
  }

  va_end(plist);
}            


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *iset(size_t n, idxtype val, idxtype *x)
{
  size_t i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *idxset(size_t n, idxtype val, idxtype *x)
{
  size_t i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
float *sset(size_t n, float val, float *x)
{
  size_t i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
double *dset(size_t n, double val, double *x)
{
  size_t i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
idxtype iamax(size_t n, idxtype *x)
{
  size_t i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
idxtype idxamax(size_t n, idxtype *x)
{
  size_t i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
idxtype idxamax_strd(size_t n, idxtype *x, idxtype incx)
{
  size_t i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
idxtype samax(size_t n, float *x)
{
  size_t i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
idxtype samax2(size_t n, float *x)
{
  size_t i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
idxtype idxamin(size_t n, idxtype *x)
{
  size_t i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
idxtype samin(size_t n, float *x)
{
  size_t i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
idxtype idxsum(size_t n, idxtype *x)
{
  size_t i;
  idxtype sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
idxtype idxsum_strd(size_t n, idxtype *x, idxtype incx)
{
  size_t i;
  idxtype sum = 0;

  for (i=0; i<n; i++, x+=incx) {
    sum += *x;
  }

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void idxadd(size_t n, idxtype *x, idxtype *y)
{
  for (n--; n>=0; n--)
    y[n] += x[n];
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
idxtype charsum(size_t n, char *x)
{
  idxtype i;
  idxtype sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
idxtype isum(size_t n, idxtype *x)
{
  size_t i;
  idxtype sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum(size_t n, float *x)
{
  size_t i;
  float sum = 0.0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum_strd(size_t n, float *x, idxtype incx)
{
  size_t i;
  float sum = 0.0;

  for (i=0; i<n; i++, x+=incx)
    sum += *x;

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void sscale(size_t n, float alpha, float *x)
{
  size_t i;

  for (i=0; i<n; i++)
    x[i] *= alpha;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float snorm2(size_t n, float *v)
{
  size_t i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += v[i] * v[i];

  return sqrt(partial);
}



/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float sdot(size_t n, float *x, float *y)
{
  size_t i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += x[i] * y[i];

  return partial;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
void saxpy(size_t n, float alpha, float *x, idxtype incx, float *y, idxtype incy)
{
  size_t i;
 
  for (i=0; i<n; i++, x+=incx, y+=incy) 
    *y += alpha*(*x);
}




/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(size_t n, idxtype *p, idxtype flag)
{
  size_t i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n <= 4)
    return;

  for (i=0; i<n; i+=16) {
    u = RandomInRange(n-4);
    v = RandomInRange(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
}



/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
idxtype ispow2(idxtype a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(idxtype seed)
{
  srand((seed == -1 ? 4321 : seed)); 
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
idxtype log2i(idxtype a)
{
  idxtype i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}


/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *GKfopen(char *fname, char *mode, char *msg)
{
  FILE *fp;
  char errmsg[256];

  fp = fopen(fname, mode);
  if (fp != NULL)
    return fp;

  sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
  perror(msg);
  errexit("Failed on GKfopen\n");
}


/*************************************************************************
* This function closes a file
**************************************************************************/
void GKfclose(FILE *fp)
{
  fclose(fp);
}


