#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define MAX_THREADS 4

int main() {
    int niter = 1000000; // Set your desired number of iterations
    double x, y, z;
    double pi_sum = 0.0;

    for (int numthreads = 1; numthreads <= MAX_THREADS; numthreads++) {
        double start_time = omp_get_wtime();
        int total_count = 0;

        #pragma omp parallel num_threads(numthreads > MAX_THREADS ? MAX_THREADS : numthreads) private(x, y, z)
        {
            unsigned int seed = (unsigned int)(time(NULL) + omp_get_thread_num());
            int thread_count = 0;

            #pragma omp for
            for (int i = 0; i < niter; i++) {
                x = (double)rand_r(&seed) / RAND_MAX;
                y = (double)rand_r(&seed) / RAND_MAX;
                z = sqrt((x * x) + (y * y));
                if (z <= 1) {
                    thread_count++;
                }
            }

            #pragma omp atomic
            total_count += thread_count;
        }

        double end_time = omp_get_wtime();
        double pi = 4.0 * total_count / niter;
        pi_sum += pi;

        printf("Threads: %d, Pi: %f, Time: %f seconds\n", numthreads, pi, end_time - start_time);
    }

    double final_pi = pi_sum / MAX_THREADS;
    printf("Average Pi Estimation: %f\n", final_pi);

    return 0;
}
