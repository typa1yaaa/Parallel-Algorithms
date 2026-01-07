#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int K;
    double *values = NULL;
    int *N = NULL;
    double recv_value = 0.0;

    if (rank == 0) {
        printf("введите K (<= %d): ", size - 1);
        fflush(stdout);
        scanf("%d", &K);

        values = (double *)malloc(K * sizeof(double));
        printf("введите %d вещественных чисел:\n", K);
        for (int i = 0; i < K; i++)
            scanf("%lf", &values[i]);

        N = (int *)malloc(size * sizeof(int));
        printf("обратите внимание, количество N = 1 должно быть %d\n", K);
        printf("введите значения N для процессов 1-%d (0 или 1):\n", size - 1);
        for (int i = 1; i < size; i++) {
            printf("процесс %d: N = ", i);
            fflush(stdout);
            scanf("%d", &N[i]);
        }
    }

    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0)
        N = (int *)malloc(size * sizeof(int));
    MPI_Bcast(N, size, MPI_INT, 0, MPI_COMM_WORLD);

    int color = (rank == 0 || N[rank] == 1) ? 0 : MPI_UNDEFINED;
    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &new_comm);

    if (color == 0) {
        int new_rank, new_size;
        MPI_Comm_rank(new_comm, &new_rank);
        MPI_Comm_size(new_comm, &new_size);

        if (new_rank == 0) {
            int index = 0;
            for (int i = 1; i < size; i++) {
                if (N[i] == 1 && index < K) {
                    MPI_Send(&values[index], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    index++;
                }
            }
        } else {
            MPI_Recv(&recv_value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("процесс %d получил число %.3f\n", rank, recv_value);
            fflush(stdout);
        }

        MPI_Comm_free(&new_comm);
    }

    if (rank == 0) {
        free(values);
        free(N);
    } else {
        free(N);
    }

    MPI_Finalize();
    return 0;
}