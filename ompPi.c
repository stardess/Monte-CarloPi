#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

double  get_clock() {
  struct timeval tv;
  int ok;
  ok = gettimeofday(&tv, (void *) 0);
  if (ok<0) {
    printf("gettimeofdayerror");
  }
  return (tv.tv_sec* 1.0 + tv.tv_usec* 1.0E-6);
}


int main(int argc, char* argv[])
{
  //  int N = 2;                                                                                                                        
  double t0= get_clock();
  int niter = 1000000;            //number of iterations per FOR loop                                                                   
  double x,y;                     //x,y value for the random coordinate                                                                 
  int i;                          //loop counter                                                                                        
  int count=0;                //Count holds all the number of how many good coordinates                                                 
  double z;                       //Used to check if x^2+y^2<=1                                                                         
  double pi;                      //holds approx value of pi                                                                            
  int numthreads = 16;

#pragma omp parallel firstprivate(x, y, z, i) shared(count) num_threads(numthreads)
  {
    srandom((int)time(NULL) ^ omp_get_thread_num());    //Give random() a seed value                                                    
    for (i=0; i<niter; ++i)              //main loop                                                                                    
      {
        x = (double)random()/RAND_MAX;      //gets a random x coordinate                                                                
        y = (double)random()/RAND_MAX;      //gets a random y coordinate                                                                
        z = sqrt((x*x)+(y*y));          //Checks to see if number is inside unit circle                                                 
        if (z<=1)
          {
            ++count;            //if it is, consider it a valid random point                                                            
          }
      }
    //print the value of each thread/rank                                                                                               
  }
  pi = ((double)count/(double)(niter*numthreads))*4.0;
  printf("Pi: %f\n", pi);
  double t1 = get_clock();
  printf("time per call: %f ns\n", (1000000000.0*(t1-t0)/niter));
  return 0;
}
