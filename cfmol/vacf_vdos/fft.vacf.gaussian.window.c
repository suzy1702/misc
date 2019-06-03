/*
 *  fft.vacf.gaussian.window.c
 *
 *  Created by Luiz Felipe Pereira on 18/Jan/2012.
 *  Mainz
 *
 *  Cosine Fourier Transform with gaussian window function. FFTW3.
 *
 *  gcc fft.vacf.gaussian.window.c -lfftw3 -lm -fopenmp -O3 -Wall -o fft.vacf
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>

#define twopi 2.0*M_PI

/* Gaussian window */
void window_Gaussian(double *in, int N, double width)
{
  int n;
  for( n=0; n<N; n++){
    in[n] *= exp( -0.5*pow(0.5*width*n/N,2.0) );
    //in[n] *= exp( -0.5*pow(n/width/2.0/N, 2.0) ); /*width=1/width*/
  }
}

int main( int argc, char *argv[] )
{
  int n;
  int N;
  double dt, width;
  double *in, *out;
  
  fftw_plan p;
  
  if( argc < 3 ){
    printf("\n\tUsage: %s <N_frames> <dt> <window width>\n", argv[0] );
    printf("\t<dt> is time between frames in femtoseconds\n\n" );
    exit(1);
  }

  N = atoi( argv[1] );
  dt = atof( argv[2] )*0.001; /* dt = <dt>*0.001 for t in pico seconds */
  width = atof( argv[3] );    /* width in picoseconds ? */

  double Tmax=(N-1)*dt;
  double df=1/Tmax;
  
  in = calloc( N, sizeof(double) );
  out = calloc( N, sizeof(double) );

  double t,f;
  n=0;
  while( scanf("%lf %lf",&t,&f)==2 ){
    in[n]=f;
    n++;
  }
/* Print original data */
  for( n=0; n<N; n++ )
    printf( "%g %g\n", n*dt,in[n] );
  printf( "\n" );

/* Apply window function */
  window_Gaussian( in, N, width );
/* Print windowed data */
  for( n=0; n<N; n++ )
    printf( "%g %g\n", n*dt,in[n] );
  printf( "\n" );

/* Calculate the FFT (-i*2*pi) */
  p = fftw_plan_r2r_1d(N, in, out, FFTW_REDFT00, FFTW_ESTIMATE);
  fftw_execute(p);
  FILE *f1;
  f1 = fopen("vdos.dat","w");
  for( n=0; n<N; n++ ){
    fprintf( f1, "%g %g\n", n*df/2.0,out[n] );
  }
  fclose( f1 );
  printf( "\n" );
  
  return 0;
}

