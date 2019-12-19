/*
 * This file contain single function - openMP realization of matrix transpose.
 *
 */


#include "utils.h"
#include <omp.h>

#ifndef TRANSPOSE_OPENMP
#define TRANSPOSE_OPENMP

namespace NMatrix {

    void transpose_OpenMP(const Matrix& matrix, Matrix &transposed_matrix) {
    MatrixSize transposed_size = getMatrixSize(transposed_matrix);
    #pragma omp parallel shared(transposed_size, matrix) 
        #pragma omp for
        for (MatrixSizeItem row_t = 0; row_t < transposed_size.first; row_t++) {
            for (MatrixSizeItem column_t = 0; column_t < transposed_size.second; column_t++) {
                transposed_matrix[row_t][column_t] = matrix[column_t][row_t];
            }
        }
    }

}

#endif
