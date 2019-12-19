/*
 * This file contain utilities to generate, compare, print and transpose matrix.
 * This utils contain only dummy, not parallel transpose.
 *
 * Functions:
 * - generate:
 *     getRandomCellItem - generate random cell item for matrix (by default int)
 *     getRandomMatrix - generate random matrix, with given size and bounds for cell item
 * - compare
 *     isMatrixEqual - elementwise matrix comparison
 * - transpose
 *     dummyNotParallelTranspose - transpose given matrix, use single thread
 * - print
 *     printMatrix - just print matrix for debug
 * - time mesurment
 *     measureTimeOfTranspose - measure time of transpose, can check transpose correctness
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <omp.h>
#include <tgmath.h>

#ifndef TRANSPOSE_UTILS
#define TRANSPOSE_UTILS

namespace NMatrix {

    using std::vector;
    using std::cout;

    typedef size_t MatrixSizeItem;
    typedef std::pair<MatrixSizeItem, MatrixSizeItem> MatrixSize;
    typedef int MatrixCellType;
    typedef vector<vector<MatrixCellType> > Matrix;

    MatrixCellType getRandomCellItem(MatrixCellType low, MatrixCellType hight) {
        return (rand() % (hight - low + 1)) + low;
    }

    Matrix getRandomMatrix(const MatrixSize& size, MatrixCellType low = -99, MatrixCellType hight = 99) {
        // generate random matrix with shape = (size.first, size.second)
        Matrix matrix = Matrix();
        matrix.resize(size.first);

        for (MatrixSizeItem row = 0; row < size.first; row++) {
            matrix[row].reserve(size.second);
            for (MatrixSizeItem column = 0; column < size.second; column++) {
                matrix[row].push_back(getRandomCellItem(low, hight));
            }
        }
        return matrix;
    }

    MatrixSize getMatrixSize(const Matrix& matrix) {
        if (matrix.size() == 0) {
            return std::pair<MatrixCellType, MatrixCellType>(0, 0);
        }
        return std::pair<MatrixCellType, MatrixCellType>(matrix.size(), matrix[0].size());
    }

    MatrixSize swapMatrixSize(const MatrixSize& size) {
        // just swap matris shapes
        return std::pair<MatrixCellType, MatrixCellType>(size.second, size.first);
    }

    bool isMatrixCorrect(const Matrix& matrix) {
        if (matrix.size() == 0 || matrix.size() == 1) {
            return true;
        }
        for (MatrixSizeItem row = 1; row < matrix.size(); row++) {
            if (matrix[row - 1].size() != matrix[row].size()) {
                return false;
            }
        }
        return true;
    }

    bool isMatrixEqual(const Matrix& a, const Matrix& b) {
        assert(isMatrixCorrect(a) && isMatrixCorrect(b) && "one of matrixes is incorrect");
        if (getMatrixSize(a) != getMatrixSize(b)) {
            return false;
        }
        for (MatrixSizeItem row = 0; row < a.size(); row++) {
            for (MatrixSizeItem column = 0; column < a[row].size(); column++) {
                if (a[row][column] != b[row][column]) {
                    return false;
                }
            }
        }
        return true;
    }

    Matrix createMatrixForTranspose(const Matrix& matrix) {
        MatrixSize transposed_size = swapMatrixSize(getMatrixSize(matrix));
        Matrix transposed_matrix = Matrix();
        transposed_matrix.resize(transposed_size.first);
        for (MatrixSizeItem row = 0; row < transposed_size.first; row++) {
            transposed_matrix[row].resize(transposed_size.second);
        }
        return transposed_matrix;
    }

    Matrix dummyNotParallelTranspose(const Matrix& matrix) {
        MatrixSize transposed_matrix_size = swapMatrixSize(getMatrixSize(matrix));
        Matrix transposed_matrix = Matrix();

        transposed_matrix.resize(transposed_matrix_size.first);
        for (MatrixSizeItem row_t = 0; row_t < transposed_matrix_size.first; row_t++) {
            transposed_matrix.reserve(transposed_matrix_size.second);
            for (MatrixSizeItem column_t = 0; column_t < transposed_matrix_size.second; column_t++) {
                transposed_matrix[row_t].push_back(matrix[column_t][row_t]);
            }
        }

        return transposed_matrix;
    }

    MatrixCellType* convertMatricToVoidPointer(const Matrix& matrix) {
        MatrixSize size = getMatrixSize(matrix);
        MatrixCellType* mem = new MatrixCellType[size.first * size.second];
        for (MatrixSizeItem row = 0; row < size.first; row++) {
            for (MatrixSizeItem column = 0; column < size.second; column++) {
                mem[row * size.second + column] = matrix[row][column];
            }
        }
        return mem;
    }

    Matrix convertVoidPointerToMatrix(const MatrixCellType* pointer, const MatrixSize &size) {
        Matrix matrix;
        matrix.resize(size.first);
        for (MatrixSizeItem row = 0; row < size.first; row++) {
            matrix[row].reserve(size.second);
            for (MatrixSizeItem column = 0; column < size.second; column++) {
                matrix[row].push_back(pointer[row * size.second + column]);
            }
        }
        return matrix;
    }


    void printMatrix(const Matrix& matrix) {
        // juts print matrix to stdout, use for debug
        assert(isMatrixCorrect(matrix) && "incorrect matrix shapes");
        MatrixSize size = getMatrixSize(matrix);
        for (MatrixSizeItem row = 0; row < size.first; row++) {
            for (MatrixSizeItem column = 0; column < size.second; column++) {
                cout << matrix[row][column] << " ";
            }
            cout << "\n";
        }
    }

    long double measureTimeOfTranspose(
            void (*transpose_func)(const Matrix&, Matrix&),
            const MatrixSize& matrix_size)
    {
        // measure time of transpose function, if it specified check correctness of transposed matrix
        Matrix matrix = getRandomMatrix(matrix_size, -1e7, +1e7);
        Matrix transposed_matrix = createMatrixForTranspose(matrix);
        long double time_start = static_cast<long double>(omp_get_wtime());
        transpose_func(matrix, transposed_matrix);
        long double time_end = static_cast<long double>(omp_get_wtime());
        return time_end - time_start;
    }

    // template<class T>
    // T computeMean(vector<T> array) {
    //     // just function to compute mean of numeric array
    //     assert(array.size() != 0 && "can't compute mean of empty array");
    //     T sum = 0;
    //     for (std::size_t index = 0; index < array.size(); ++index) {
    //         sum += array[index];
    //     }
    //     return sum / array.size();
    // }

    // template<class T>
    // T computeStd(vector<T> array) {
    //     // just function to compute std of numeric array
    //     assert(array.size() > 1 && "can't compute std of empty or one-item array");
    //     T mean = computeMean(array);
    //     T sq_sum = 0;
    //     for (std::size_t index = 0; index < array.size(); ++index) {
    //         sq_sum += (array[index] - mean) * (array[index] - mean);
    //     }
    //     return sqrtl(sq_sum / (array.size() - 1));
    // }

    // std::pair<long double, long double> measureNTranspose(
    //         void (*transpose_func)(const Matrix&, Matrix&),
    //         int number_of_test,
    //         const MatrixSize& matrix_size,
    //         bool check_correctness)
    // {
    //     // measure mean and std of number_of_test tests

    //     // make tests
    //     vector<long double> times;
    //     for (int i = 0; i < number_of_test; ++i) {
    //         Matrix rand_matrix = getRandomMatrix(matrix_size, -1e7, +1e7);
    //         times.push_back(measureTimeOfTranspose(transpose_func, rand_matrix, check_correctness));
    //     }

    //     // compute statistics
    //     return std::pair<long double, long double>(
    //         computeMean(times),
    //         computeStd(times)
    //     );
    // }
}

#endif