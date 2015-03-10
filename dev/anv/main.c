#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sndfile.h>
#include <complex.h>

#define ELEMENT_COUNT(X) (sizeof(X) / sizeof((X)[0]))
#define TRUE 1
#define FALSE 0
#define BUFSIZE 1024


/* Problem 1.1 Efficient implementations of LP-filters*/
void pb1(){


// Declaration of specs
/*
int D=4, Fs = 16000, Fs1 = Fs/D, Fs2 = Fs1/D;
int Fg1 = Fs1/D, Fc1 = 2000;
int Fg2 = 150,   Fc2 = 270;
double rg  = 0.001, rc = 0.001;
double rg1 = rg/D;
int xin = 1;
*/










}
/*
void pb2(){
    double f1 = 0.00875;
    double f2 = 0.0175;
    double pi = 3.1415;
    int Fs = 16000;

    double x[10];
    int n;
    for (n = 0; n < length(x); n++){
        x = sin(2*pi*f1*n) + sin(2*pi*f2*n);
    }
}*/


double downsampler(double *xin[], int D){
    /*
    double xout[length(xin)/D];
    j = 0;
    for (i = 0; i <= length(xin); i += D){
        xout[j] = &xin[i];
        j++;
    }
    return xout;
    */
}



double upsampler(double xin[], int I){
    /*
    double xout[length(xin)*I];
    int j = 0;
    for (i=0; i<= length(xout);i++){
        xout[i] = 0;
    }

    for (i= 0; i<= length(xin)*I; i+I)){
        xout[i] = &xin[j];
        j++;
    }
    return xout;
    */
}

double flip(double in[]){
    /*
    double out[48];
    int i;
    for (i = 0; i < length(out); i++){
        out[i] = in[length(in) - i];
    }*/
}

void CMF_coeffs(double QMF_h0[]){
    /*
    int size = 48;
    double g0[size] = flip(QMF_h0);
    double h1[size];
    int i;
    for (i=0; i < length(g0); i++){
        g0[i] = 2*g0[i];
    }
    for(i = 0; i < length(g0); i+=2){
        g1[i] = -2*QMF_h0[i];
        h1[i] = -g0[i+1]
    }*/
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

static void fir_filter (double xin, char * coeff ){
    double buffer[BUFSIZE], obuf[BUFSIZE];
    sf_count_t count, cnt;
    int i, n, nt=0, m, N, Nuse; /* N = number of coeff in filter*/
    FILE *f;

    /* Read filter coeff from input file */
    if ((f = fopen (coeff,"r")) == NULL ) exit(1);

    fscanf(f,"%d", &N);
    printf("%d",N);





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

/*

  FILE *f = fopen("filterCoeff_fullrate.txt", "r");
  if (f == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }

    float *h, *mem;
    int N,i;
    /* first line gives filter order, use to allocate memory
	fscanf(f,"%d",&N);
    if ((h=(float *)calloc(N,sizeof(float))) == NULL) exit(1);
	if ((mem=(float *)calloc(N,sizeof(float))) == NULL) exit(1);
	printf("%d\n",N);
    for (i=0; i<N; i++){
		fscanf(f,"%f",&h[i]);
		printf("%7.6f\n",h[i]);
		mem[i]=0;        /* initialize filter memory at the same time */








/*
    fclose(f);

    free(h);
    free(mem);
    h=NULL;
    mem=NULL;

//    COMPLEX a;



     double complex numb = 33 +222*I;

     printf("real part %f, im part %f\n", creal(numb), cimag(numb));




*/

  int xin[] = {1};
  return 0;
}
