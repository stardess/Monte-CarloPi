#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main() {
    int numthreads = 16;
    int niter = 1000000; // Set your desired number of iterations
    double x, y, z;
    int count = 0;
    int i;

    double start_time = omp_get_wtime();

    #pragma omp parallel firstprivate(x, y, z, i) shared(count) num_threads(numthreads)
    {
        srandom((int)time(NULL) ^ omp_get_thread_num());
        #pragma omp for reduction(+:count)
        for (i = 0; i < niter; i++) {
            x = (double)random() / RAND_MAX;
            y = (double)random() / RAND_MAX;
            z = sqrt((x * x) + (y * y));

            if (z <= 1) {
                count++;
            }
        }
    }

    double pi = 4.0 * count / niter;
    double end_time = omp_get_wtime();
    printf("Estimated pi: %f\n", pi);
    printf("Time taken: %f seconds\n", end_time - start_time);

    return 0;
}