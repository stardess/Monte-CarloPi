#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>

#define NITER_TOTAL 10000000 // Total number of iterations                                                                                                                                
#define MAX_THREADS 4

int main(int argc, char *argv[]) {
    int numprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start_time, end_time;

    if (rank == 0) {
        start_time = MPI_Wtime(); // Record start time only in rank 0                                                                                                                     
    }

    int niter_per_process = NITER_TOTAL / numprocs; // Split iterations among processes                                                                                                   

    double local_pi_sum = 0.0;

    // Perform local calculations for each process                                                                                                                                        
    srand((unsigned int)time(NULL) + rank); // Seed the random number generator                                                                                                           

    for (int numthreads = 1; numthreads <= MAX_THREADS; numthreads++) {
        int total_count = 0;

        #pragma omp parallel num_threads(numthreads)
        {
            double x, y, z;
            int thread_count = 0;
            unsigned int seed = (unsigned int)(rand());

            #pragma omp for
            for (int i = 0; i < niter_per_process; i++) {
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

        double pi = 4.0 * total_count / niter_per_process;
        local_pi_sum += pi;
    }

    double local_pi_avg = local_pi_sum / MAX_THREADS;

    // Gather all local_pi_avg values from each process                                                                                                                                   
    double *all_pi_avgs = NULL;
    if (rank == 0) {
        all_pi_avgs = (double *)malloc(numprocs * sizeof(double));
    }

    MPI_Gather(&local_pi_avg, 1, MPI_DOUBLE, all_pi_avgs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculate the final average of Pi from all gathered values                                                                                                                         
    if (rank == 0) {
        double final_pi = 0.0;
        for (int i = 0; i < numprocs; i++) {
            final_pi += all_pi_avgs[i];
        }
        final_pi /= numprocs;

        end_time = MPI_Wtime(); // Record end time only in rank 0                                                                                                                         

        printf("Average Pi Estimation: %f\n", final_pi);
        printf("Total Time: %f seconds\n", end_time - start_time);

        free(all_pi_avgs);
    }

    MPI_Finalize();
    return 0;
}


