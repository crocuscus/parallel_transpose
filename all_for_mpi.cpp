#include <iostream>
#include <assert.h>
#include <vector>
#include <iomanip>
#include <mpi.h>
#include <cmath>
#include <sstream>
#include <limits.h>
#include "cstdlib"

#include "utils.h"

using namespace NMatrix;

void transpose_MPI(const int* matrix, MatrixSize matrix_size, int* transposed_matrix) {
    MatrixSize transposed_size = swapMatrixSize(matrix_size);

//    printf("*%d*", matrix[0]);
    int myrank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MatrixSizeItem len = transposed_size.first / static_cast<MatrixSizeItem>(size) + 1;
    MatrixSizeItem my_start = len * static_cast<MatrixSizeItem>(myrank);
    MatrixSizeItem my_end = std::min(
        len * static_cast<MatrixSizeItem>(myrank) + len,
        transposed_size.first
    );

//    printf("I'am %d of %d, work under %d - %d\n", myrank, size, (int)my_start, (int)my_end);

    for (MatrixSizeItem row_t = my_start; row_t < my_end; row_t++) {
        for (MatrixSizeItem column_t = 0; column_t < transposed_size.second; column_t++) {
            transposed_matrix[
                row_t * transposed_size.second + column_t
            ] = matrix[column_t * transposed_size.first + row_t];
        }
    }
}


int main(int argc, char *argv[])
{
    if (argc != 2) {
        printf("expect 1 argument : matrix size");
        exit(1);
    }
    std::size_t matrix_size;
    sscanf(argv[1], "%lu", &matrix_size);

    int err;
    if ((err = MPI_Init(&argc, &argv)))  {
        std::cerr << "MPI startup error!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, err);
        exit(2);
    }
    Matrix _mat = getRandomMatrix({matrix_size, matrix_size}, 0, 9);;
    MatrixCellType* mat;
    MatrixCellType* t_mat;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//    if (myrank == 0) {
//        mat = convertMatricToVoidPointer(_mat);
//        t_mat = new MatrixCellType[_mat.size() * _mat[0].size()];
//        MPI_Bcast(mat, matrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Bcast(t_mat, matrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
//    }

    mat = convertMatricToVoidPointer(_mat);
    t_mat = new MatrixCellType[_mat.size() * _mat[0].size()];
    int *buf = new int[matrix_size * matrix_size];

    MPI_Bcast(mat, (long long)matrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(t_mat, (long long)matrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

    double time_start = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    transpose_MPI(mat, {matrix_size, matrix_size}, t_mat);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(t_mat, buf, matrix_size * matrix_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double sm_time = MPI_Wtime() - time_start;

    double* sm_pnt_time = new double[1];
    sm_pnt_time[0] = sm_time;
    double* time = new double[1];
    MPI_Reduce(sm_pnt_time, time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    if (myrank == 0) {
//        Matrix _mat_t = convertVoidPointerToMatrix(buf, {matrix_size, matrix_size});
//        printf("\norigin:\n");
//        printMatrix(_mat);
//        printf("\n\ntransposed:\n");
//        printMatrix(_mat_t);
//        printf("time : %lf", sm_time);
        printf("%lf", (*time) * 1000);
    }


    MPI_Finalize();
    return 0;
}



