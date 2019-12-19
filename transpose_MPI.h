#ifndef TRANSPOSE_MPI_H
#define TRANSPOSE_MPI_H

#include "utils.h"
#include <mpi.h>

namespace NMatrix {

    void transpose_MPI(const Matrix& matrix, Matrix &transposed_matrix) {
        MatrixSize transposed_size = getMatrixSize(transposed_matrix);

        int myrank, size;
//        MPI_Status Status;

        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        MatrixSizeItem len = transposed_size.first / static_cast<MatrixSizeItem>(size) + 2;
        MatrixSizeItem my_start = len * static_cast<MatrixSizeItem>(myrank);
        MatrixSizeItem my_end = std::min(
            len * static_cast<MatrixSizeItem>(myrank) + len,
            transposed_size.first
        );

        printf("I'am %d of %d, work under %lx - %lx\n", myrank, size, my_start, my_end);

        for (MatrixSizeItem row_t = my_start; row_t < my_end; row_t++) {
            for (MatrixSizeItem column_t = 0; column_t < transposed_size.second; column_t++) {
                transposed_matrix[row_t][column_t] = matrix[column_t][row_t];
            }
        }

    }
}

#endif // TRANSPOSE_MPI_H
