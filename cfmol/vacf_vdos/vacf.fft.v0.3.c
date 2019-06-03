/*
 *  vacf.fft.v0.3.c  
 *
 *  Created by Luiz Felipe Pereira on 24/Jan/2012.
 *  MPI-P
 *
 *  Calculates velocity autocorrelation function using Fourier transforms.
 * 
 *  v0.2 prints VACF x-,y-,and z-components independently.
 *  v0.3 calculates <T>=m/3k_B <v^2> from velocity trajectory.
 *
 *  gcc vacf.fft.v0.3.c -lfftw3 -lm -fopenmp -O3 -Wall -o vacf.fft.v03
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>

int main( int argc, char *argv[] )
{
  int i, icount, atom;
  int N, Natoms, NN;
  int itemp, step;
  double xtemp, ytemp, ztemp;
  double dt;
  double vx2avg,vy2avg,vz2avg;
  double Tavg;
  double const m=12*1.660538921e-27; //C12 mass in Kg.
  double const kB=1.3806503e-23;  //SI units
  double const convert=1e-20/1e-24; //A/ps to m/s
  
  if( argc < 3 ){
    printf("\n\tUsage: %s <N_atoms> <N_frames> <dt>\n", argv[0] );
    printf("\t<dt> is time between frames in femtoseconds\n\n" );
    exit(1);
  }

  Natoms = atoi( argv[1] );
  N = atoi( argv[2] );
  dt = atof( argv[3] )*0.001; /* dt = <dt>*0.001 for t in pico seconds */
	
  fftw_complex **vx, **vy, **vz;
  fftw_plan p;

  NN = 2*N;
  vx = fftw_malloc(sizeof(fftw_complex) * Natoms);
  vy = fftw_malloc(sizeof(fftw_complex) * Natoms);
  vz = fftw_malloc(sizeof(fftw_complex) * Natoms);
  for( atom=0; atom<Natoms; atom++ ){
    vx[atom] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NN);
    vy[atom] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NN);
    vz[atom] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NN);
  }
	
  for( atom=0; atom<Natoms; atom++ ){
    #pragma omp parallel for private(i)
    for( i=0; i<NN; i++ ){
      vx[atom][i] = 0.0+I*0.0;
      vy[atom][i] = 0.0+I*0.0;
      vz[atom][i] = 0.0+I*0.0;
    }
  }

  /* Read data */
  icount=0;
  step=0;
  while( scanf("%d %lf %lf %lf \n",&itemp,&xtemp,&ytemp,&ztemp ) == 4 ){
    if( icount == Natoms ){
      step++; 
      icount=0;
      //if( 100*step%(10*N) == 0 ) printf("# %3g %% read\n", 100.0*step/N);
    }
    atom = itemp%Natoms;
    vx[atom][step]=xtemp;
    vy[atom][step]=ytemp; 
    vz[atom][step]=ztemp;  
    icount++;
  }

  double v2avg[Natoms];
  for( atom=0; atom<Natoms; atom++ ){
    v2avg[atom]=0.0;
    for( i=0; i<N; i++ ){
      v2avg[atom] += creal(vx[atom][i])*creal(vx[atom][i])
                   + creal(vy[atom][i])*creal(vy[atom][i])
                   + creal(vz[atom][i])*creal(vz[atom][i]);
    }
    v2avg[atom] /= (1.0*N);
  }
//  v2avg = v2avg[0]*convert;
  Tavg = m*v2avg[10]*convert/3.0/kB;
  //printf( " <T> = %g K \n", Tavg );

  /* Calculate the FFT (-i*2*pi) */
  for( atom=0; atom<Natoms; atom++ ){
    p = fftw_plan_dft_1d(NN, vx[atom], vx[atom], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    p = fftw_plan_dft_1d(NN, vy[atom], vy[atom], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    p = fftw_plan_dft_1d(NN, vz[atom], vz[atom], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
  }

  /* Multiply FFT by FFT* */
  for( atom=0; atom<Natoms; atom++ ){
    #pragma omp parallel for private(i)
    for( i=0; i<NN; i++ ){
      vx[atom][i] = vx[atom][i]*conj( vx[atom][i] );
      vy[atom][i] = vy[atom][i]*conj( vy[atom][i] );
      vz[atom][i] = vz[atom][i]*conj( vz[atom][i] ); 	    
    }
  }
  
  /* Calculate IFFT (+i*2*pi) */
  for( atom=0; atom<Natoms; atom++ ){
    p = fftw_plan_dft_1d(NN, vx[atom], vx[atom], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    p = fftw_plan_dft_1d(NN, vy[atom], vy[atom], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    p = fftw_plan_dft_1d(NN, vz[atom], vz[atom], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
  }

  /* Normalize FFT */
  for( atom=0; atom<Natoms; atom++ ){
    #pragma omp parallel for private(i)
    for( i=0; i<NN; i++ ){
      vx[atom][i] /= NN;
      vy[atom][i] /= NN;
      vz[atom][i] /= NN;
    }
  }

  /* Normalize ACF */
  for( atom=0; atom<Natoms; atom++ ){
    #pragma omp parallel for private(i)
    for( i=0; i<N; i++ ){
      vx[atom][i] /= (N-i);
      vy[atom][i] /= (N-i);
      vz[atom][i] /= (N-i);
    }
  }

  /* Calculate VACF averaging over all atoms */
  double vacf, vacfx, vacfy, vacfz;
  FILE *f0;
  f0 = fopen( "vacf.dat","w" );
  for( i=0; i<N; i++ ){
    vacf=vacfx=vacfy=vacfz = 0.0;
    for( atom=0; atom<Natoms; atom++ ){
      vacfx += creal(vx[atom][i])/creal(vx[atom][0]);
      vacfy += creal(vy[atom][i])/creal(vy[atom][0]);
      vacfz += creal(vz[atom][i])/creal(vz[atom][0]);
    }
    vacfx /= Natoms;
    vacfy /= Natoms;
    vacfz /= Natoms;
    vacf = (vacfx+vacfy+vacfz)/3.0;
    /* Print ACF data */
    fprintf( f0, "%g %g %g %g %g\n", i*dt,vacf,vacfx,vacfy,vacfz );
  }
  fclose( f0 );

  /* dealloc memory */
  for( atom=0; atom<Natoms; atom++ ){ 
    fftw_free(vx[atom]);
    fftw_free(vy[atom]);  
    fftw_free(vz[atom]);
  }
  fftw_free(vx);
  fftw_free(vy);  
  fftw_free(vz);
  fftw_destroy_plan(p);
  return 0;
}

