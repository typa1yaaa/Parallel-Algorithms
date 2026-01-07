#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#define FOUR_SEQ 0.000002
#define FOUR_PAR 0.000037

void seq_matmul(int *A, int *B, int *C, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int sum = 0;
            for (int k = 0; k < n; k++)
                sum += A[i*n + k] * B[k*n + j];
            C[i*n + j] = sum;
        }
}

void print_matrix(int *M, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%6d ", M[i*n + j]);
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if(rank==0) printf("Использование: mpiexec -n <число_процессов> %s <размер_матрицы>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int n = atoi(argv[1]);
    int q = (int)sqrt(size);
    if (q*q != size || n % q != 0) {
        if (rank == 0)
            printf("Число процессов должно быть квадратом целого, а размер матрицы делиться на sqrt(число процессов)\n");
        MPI_Finalize();
        return 1;
    }
    int block = n / q;

    int *A=NULL, *B=NULL, *C_seq=NULL, *C_par=NULL;
    double t_seq=0, t_par=0;

    if(rank==0){
        A = malloc(n*n*sizeof(int));
        B = malloc(n*n*sizeof(int));
        C_seq = malloc(n*n*sizeof(int));
        C_par = malloc(n*n*sizeof(int));

        srand(time(NULL));
        for(int i=0;i<n*n;i++){
            A[i] = rand()%100+1;
            B[i] = rand()%100+1;
        }

        double t0 = MPI_Wtime();
        seq_matmul(A,B,C_seq,n);
        double t1 = MPI_Wtime();
        t_seq = t1 - t0;

        printf("Время последовательного алгоритма: %f сек\n", t_seq);
    }

    int dims[2]={q,q}, periods[2]={0,0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&cart_comm);

    int coords[2]; MPI_Cart_coords(cart_comm, rank, 2, coords);
    int row_rank = coords[0], col_rank = coords[1];

    int remain_dims[2];
    remain_dims[0]=0; remain_dims[1]=1;
    MPI_Comm row_comm; MPI_Cart_sub(cart_comm, remain_dims, &row_comm);
    remain_dims[0]=1; remain_dims[1]=0;
    MPI_Comm col_comm; MPI_Cart_sub(cart_comm, remain_dims, &col_comm);

    int *A_strip = malloc(block*n*sizeof(int));
    int *B_strip = malloc(n*block*sizeof(int));
    int *C_strip = calloc(block*n,sizeof(int));

    if(rank==0){
        for(int i=0;i<q;i++){
            int *tmpA = malloc(block*n*sizeof(int));
            for(int j=0;j<block;j++)
                for(int k=0;k<n;k++)
                    tmpA[j*n+k] = A[(i*block+j)*n + k];
            for(int j=0;j<q;j++){
                int dest; MPI_Cart_rank(cart_comm,(int[]){i,j},&dest);
                if(dest==0) memcpy(A_strip,tmpA,block*n*sizeof(int));
                else MPI_Send(tmpA,block*n,MPI_INT,dest,0,cart_comm);
            }
            free(tmpA);
        }

        MPI_Datatype col_type; MPI_Type_vector(n,block,n,MPI_INT,&col_type); MPI_Type_commit(&col_type);
        for(int j=0;j<q;j++){
            int *tmpB = malloc(n*block*sizeof(int));
            for(int i0=0;i0<n;i0++)
                for(int k=0;k<block;k++)
                    tmpB[i0*block+k] = B[i0*n + (j*block+k)];
            for(int i1=0;i1<q;i1++){
                int dest; 
                MPI_Cart_rank(cart_comm,(int[]){i1,j},&dest);
                if(dest==0) memcpy(B_strip,tmpB,n*block*sizeof(int));
                else MPI_Send(tmpB,n*block,MPI_INT,dest,1,cart_comm);
            }
            free(tmpB);
        }
        MPI_Type_free(&col_type);
    }else{
        MPI_Recv(A_strip,block*n,MPI_INT,0,0,cart_comm,MPI_STATUS_IGNORE);
        MPI_Recv(B_strip,n*block,MPI_INT,0,1,cart_comm,MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    int *B_temp = malloc(n*block*sizeof(int));
    memcpy(B_temp,B_strip,n*block*sizeof(int));

    int left,right; MPI_Cart_shift(row_comm,0,-1,&left,&right);
    for(int step=0;step<q;step++){
        int current_block = (col_rank + step)%q;
        for(int i=0;i<block;i++)
            for(int j=0;j<block;j++){
                int sum=0;
                for(int k=0;k<n;k++)
                    sum += A_strip[i*n+k]*B_temp[k*block+j];
                C_strip[i*n + current_block*block+j] = sum;
            }
        if(step<q-1) MPI_Sendrecv_replace(B_temp,n*block,MPI_INT,right,0,left,0,row_comm,MPI_STATUS_IGNORE);
    }
    free(B_temp);

    if(rank==0){
        for(int i=0;i<block;i++)
            for(int j=0;j<n;j++)
                C_par[i*n+j] = C_strip[i*n+j];
        for(int i=1;i<q;i++){
            int *tmpC = malloc(block*n*sizeof(int));
            int src; MPI_Cart_rank(cart_comm,(int[]){i,0},&src);
            MPI_Recv(tmpC,block*n,MPI_INT,src,2,cart_comm,MPI_STATUS_IGNORE);
            for(int j=0;j<block;j++)
                for(int k=0;k<n;k++)
                    C_par[(i*block+j)*n + k] = tmpC[j*n+k];
            free(tmpC);
        }
    }else if(col_rank==0){
        MPI_Send(C_strip,block*n,MPI_INT,0,2,cart_comm);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    t_par = t1 - t0;

    if(rank==0){
        printf("Время параллельного алгоритма: %f сек\n", t_par);
        printf("Ускорение последовательного алгоритма: %f\n", FOUR_SEQ/t_seq);
        printf("Ускорение параллельного алгоритма: %f\n", t_seq/t_par);

        int ok=1;
        for(int i=0;i<n*n;i++) if(C_seq[i]!=C_par[i]){ok=0; break;}
        printf("%s\n", ok ? "Результаты совпадают":"Результаты не совпадают");
    }

    free(A_strip); free(B_strip); free(C_strip);
    if(rank==0) free(A), free(B), free(C_seq), free(C_par);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
