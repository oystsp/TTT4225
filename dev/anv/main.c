#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sndfile.h>

#define ELEMENT_COUNT(X) (sizeof(X) / sizeof((X)[0]))
#define TRUE 1
#define FALSE 0

/* Problem 1.1 Efficient implementations of LP-filters*/
void pb1(){


// Decaration of specs
/*
int D=4, Fs = 16000, Fs1 = Fs/D, Fs2 = Fs1/D;
int Fg1 = Fs1/D, Fc1 = 2000;
int Fg2 = 150,   Fc2 = 270;
double rg  = 0.001, rc = 0.001;
double rg1 = rg/D;
int xin = 1;
*/










}



int rand_gauss (float *x, int N) {
/* Create Gaussian N(0,1) distributed numbers from uniformly distributed numbers using */
  float v1,v2,s;
  int i, j, M;

  M=N/2;

  /* Initialize uniform number generator */
  srand(time(NULL));

  /* Loop - each iteration generates two Gaussian numbers*/
  for (i=0; i<M; i++){
  j=2*i;


  	do {
    	v1 = 2.0 * ((float) rand()/RAND_MAX) - 1;
    	v2 = 2.0 * ((float) rand()/RAND_MAX) - 1;

    	s = v1*v1 + v2*v2;
  	  } while ( s >= 1.0 );

  	if (s == 0.0)
    	i=i-1;
  	  else {
    	x[j]=(v1*sqrt(-2.0 * log(s) / s));
		if (j+1<N)
			x[j+1]=(v2*sqrt(-2.0 * log(s) / s));
	  }
  	}
	return 0;
}



void convolve(const double Signal[/* SignalLen */], size_t SignalLen,
              const double Kernel[/* KernelLen */], size_t KernelLen,
              double Result[/* SignalLen + KernelLen - 1 */])
{
  size_t n;

  for (n = 0; n < SignalLen + KernelLen - 1; n++)
  {
    size_t kmin, kmax, k;
    Result[n] = 0;
    kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
    kmax = (n < SignalLen - 1) ? n : SignalLen - 1;
    for (k = kmin; k <= kmax; k++)
    {
      Result[n] += Signal[k] * Kernel[n - k];
    }
  }
}

void printSignal(const char* Name,
                 double Signal[/* SignalLen */], size_t SignalLen)
{
  size_t i;
  for (i = 0; i < SignalLen; i++)
  {
    printf("%s[%zu] = %f\n", Name, i, Signal[i]);
  }
  printf("\n");
}







int main()
{
  /*
  double signal[] = { 1, 1, 0, 1, 1 };
  double kernel[] = { 1, 1, 1, 1, 1 };
  double result[ELEMENT_COUNT(signal) + ELEMENT_COUNT(kernel) - 1];
  convolve(signal, ELEMENT_COUNT(signal),
           kernel, ELEMENT_COUNT(kernel),
           result);
  printSignal("signal", signal, ELEMENT_COUNT(signal));
  printSignal("kernel", kernel, ELEMENT_COUNT(kernel));
  printSignal("result", result, ELEMENT_COUNT(result));

*/

  FILE *f = fopen("filterCoeff_fullrate.txt", "r");
  if (f == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }

    float *h, *mem;
    int N,i;
    /* first line gives filter order, use to allocate memory */
	fscanf(f,"%d",&N);
    if ((h=(double *)calloc(N,sizeof(double))) == NULL) exit(1);
	if ((mem=(double *)calloc(N,sizeof(double))) == NULL) exit(1);
	printf("%d\n",N);
    for (i=0; i<N; i++){
		fscanf(f,"%f",&h[i]);
		printf("%7.6f\n",h[i]);
		mem[i]=0;        /* initialize filter memory at the same time */
	}








    fclose(f);

    free(h);
    free(mem);


  return 0;
}
