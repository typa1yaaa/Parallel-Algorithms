#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

int main(int argc, char** argv) {
    double start_time, end_time;
    int rank, size;

    MPI_Init(&argc, &argv);

    start_time = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

    if (rank == 0) {        
        for (int i = 1; i < size; i++) {
            int process_rank;
            MPI_Status status;
            MPI_Recv(&process_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            printf("Hello from process %3d\n", process_rank);
        }
    } else {
        MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        end_time = MPI_Wtime();
        printf("%.6f\n", end_time-start_time);
    }
    
    MPI_Finalize();
    
    return 0;
}