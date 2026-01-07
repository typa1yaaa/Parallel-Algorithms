#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#define ARRAY_SIZE 10000000
#define ONE 0.09496339

int main(int argc, char** argv) {
    int rank, size;
    double start_time, end_time;
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_min = INT_MAX;
    int global_min = INT_MAX;
    int* global_array = NULL;

    int local_size = ARRAY_SIZE / size;
    
    
    if (rank == 0) {
        start_time = MPI_Wtime();

        global_array = (int*)malloc(ARRAY_SIZE * sizeof(int));
        srand(time(NULL));
        
        // printf("массив: ");
        for (int i = 0; i < ARRAY_SIZE; i++) {
            global_array[i] = rand() % 1000000;
            // printf("%d ", global_array[i]);
        }
        // printf("\n");

    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    int* local_ptr = NULL;

    if (rank == 0) {
        // процесс 0 получает данные "на месте"
        MPI_Scatter(global_array, local_size, MPI_INT,
                    MPI_IN_PLACE, local_size, MPI_INT,
                    0, MPI_COMM_WORLD);
        local_ptr = global_array;
    } else {
        local_ptr = (int*)malloc(local_size * sizeof(int));
        MPI_Scatter(NULL, local_size, MPI_INT,
                    local_ptr, local_size, MPI_INT,
                    0, MPI_COMM_WORLD);
    }
    
    if (rank == 0) {
        printf("процесс 0: раздал данные с помощью MPI_Scatter\n");
        local_min = INT_MAX;
    } else {
        for (int i = 0; i < local_size; i++) {
            if (local_ptr[i] < local_min) {
                local_min = local_ptr[i];
            }
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // if (rank != 0) {
        // printf("процесс %d: локальный минимум = %d\n", rank, local_min);
    // }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    int* all_mins = NULL;
    if (rank == 0) {
        all_mins = (int*)malloc(size * sizeof(int));
    }
    
    MPI_Gather(&local_min, 1, MPI_INT, all_mins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        global_min = INT_MAX;
        for (int i = 1; i < size; i++) {  
            if (all_mins[i] < global_min) {
                global_min = all_mins[i];
            }
        }
        
        printf("собранные минимумы: ");
        for (int i = 1; i < size; i++) {
            printf("%d ", all_mins[i]);
        }
        printf("\nглобальный минимум: %d\n", global_min);
        
        free(global_array);
        free(all_mins);
    }else{
        free(local_ptr);
    }
    
    if (rank==0){
        end_time = MPI_Wtime();
        double total_time = end_time - start_time;
        printf("время выполнения: %.8f\n", total_time);
        printf("ускорение: %.8f\n", ONE/total_time);

    }

    MPI_Finalize();
    return 0;
}