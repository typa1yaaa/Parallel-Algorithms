#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#define TWO_PROC (0.000016+0.000012+0.000011)/3

int main(int argc, char** argv) {
    int rank, size;
    double start_time, end_time;
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank);

    int K = size;
    
    int* send_data = (int*)malloc(K * sizeof(int));
    int* recv_data = (int*)malloc(K * sizeof(int));
    

    for (int i = 0; i < K; i++) {
        send_data[i] = rand() % 100;  
    }

    for (int i = 0; i < K; i++) {
        if (rank == i) {
            printf("данные процесса %d: ", rank);
            for (int j = 0; j < K; j++) {
                printf("%d ", send_data[j]);
            }
            printf("\n");
            fflush(stdout);  
        }
        MPI_Barrier(MPI_COMM_WORLD);  
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start_time = MPI_Wtime();
    
    MPI_Alltoall(send_data, 1, MPI_INT, 
                 recv_data, 1, MPI_INT, 
                 MPI_COMM_WORLD);
    
    end_time = MPI_Wtime();
    
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        if (rank == i) {
            printf("процесс %d получил: ", rank);
            for (int j = 0; j < K; j++) {
                printf("%d ", recv_data[j]);
            }
            printf("\n");
            fflush(stdout); 
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(send_data);
    free(recv_data);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf("общее время выполнения %f секунд\n", end_time - start_time);
        printf("ускорение %f секунд\n",TWO_PROC/(end_time - start_time));

    }
    
    MPI_Finalize();
    return 0;
}