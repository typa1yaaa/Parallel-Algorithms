#include <mpi.h>
#include <stdio.h>

#define THREE_PROC 0.000208

int main(int argc, char *argv[]){
    int rank, size;
    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size % 3 != 0) {
        if (rank == 0)
            printf("ошибка: число процессов должно быть кратно 3\n");
        MPI_Finalize();
        return 0;
    }

    int N = size / 3;

    // создаем декартовую топологию 3xN
    int dims[2] = {3, N};
    int periods[2] = {0, 0};
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart);

    int coords[2];
    MPI_Cart_coords(cart, rank, 2, coords);

    int row = coords[0];
    int col = coords[1];

    // создаем три строки
    int remain_dims[2] = {0, 1};
    MPI_Comm row_comm;

    MPI_Cart_sub(cart, remain_dims, &row_comm);

    int row_rank;
    MPI_Comm_rank(row_comm, &row_rank);

    //задаем число в главных процессах 0, N, 2N
    int value = -1;

    if (rank == 0)  value = 100;   
    if (rank == N)  value = 200;
    if (rank == 2*N)  value = 300;

    MPI_Bcast(&value, 1, MPI_INT, 0, row_comm);

    // printf("глобальный ранг %d (строчка = %d, столбец = %d): значение = %d\n", rank, row, col, value);


    MPI_Barrier(MPI_COMM_WORLD); 
    double end_time = MPI_Wtime(); 

    if (rank == 0) {
        printf("время работы программы: %f секунд\n", end_time - start_time);
        printf("ускорение: %f секунд\n", THREE_PROC/(end_time - start_time));

    }

    MPI_Finalize();
    return 0;
}