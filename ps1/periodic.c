#include<stdio.h>
#include<math.h>
#include <gsl/gsl_rng.h>


int main() {

  const gsl_rng_type * ranT;
  gsl_rng * ranr;

  gsl_rng_env_setup();
  ranT = gsl_rng_default;
  ranr = gsl_rng_alloc (ranT);

  double y,L,x[100+1],xtrial,disp_size,dx;

  int i,j,k,N,n,hist[1000+1]={0},count=0;
  int sweep,nsweeps,xindex,step;

  nsweeps=1000000;
  disp_size=2;
  dx = 0.1;

  N=10;
  L=15;

  for (i=1; i<=N; i++) {
    x[i]=i*1.0001;
  }
  x[0] = -1;
  x[N+1] = L+1;

  int imove,iup,idown,overlap;
  double deltax;

  for (sweep=1; sweep<=nsweeps; sweep++) {
    for (step=1; step<=N; step++) {
      overlap=0;

      y = gsl_rng_uniform(ranr);
      imove = ceil(y*N);

      y = gsl_rng_uniform(ranr);
      xtrial = x[imove] + disp_size*(y-0.5);
      
      iup = imove+1;
      idown = imove-1;

      if (imove==1) {
		idown = N;
      }
      if (imove==N) {
		iup = 1;
      }

      deltax = xtrial - x[idown];
      deltax -= L*round(deltax/L);
      if (deltax*deltax < 1) {
		overlap=1;
      }

      deltax = x[iup] - xtrial;
      deltax -= L*round(deltax/L);
      if (deltax*deltax < 1) {
		overlap=1;
      }

      if (overlap==0) {
		x[imove] = xtrial;
	
      }
      //            printf("%d %f\n",sweep,x[1]);

    }


    for (i=1; i<=N; i++) {
      x[i] -= L*floor(x[i]/L);
      
      xindex = ceil(x[i]/dx);
      hist[xindex]++;
    }
    count++;


    //    printf("%d\n",n);
  }


  printf("0 0\n");
  for (i=1; i<=1000; i++) {
    if (hist[i]>0) {
      printf("%f %f\n",(i-0.5)*dx,hist[i]/(count*dx));
    }
  }


}
