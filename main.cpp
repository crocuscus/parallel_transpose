/*
Для выполнения практических заданий используются суперкомпьютерные вычислительные ресурсы факультета ВМК - Bluegene/P и Polus.
http://hpc.cs.msu.ru/

В каждой задаче требуется:
1) Реализовать параллельную версию предложенного алгоритма с использованием технологий OpenMP и MPI.
2) Начальные параметры для задачи подбираются таким образом, чтобы:
- Задача помещалась в оперативную память одного процессора.
- Время решения задачи было в примерном диапазоне 5 сек.-15 минут.
3) Исследовать масштабируемость полученной параллельной программы: построить графики зависимости времени исполнения от числа ядер/процессоров для различного объёма входных данных.
Для каждого набора входных данных найти количество ядер/процессоров, при котором время выполнения задачи перестаёт уменьшаться.
Оптимальным является построение трёхмерного графика: по одной из осей время работы программы, по другой - количество ядер/процессоров и по третьей - объём входных данных.
Каждый прогон программы с новыми параметрами рекомендуется выполнять несколько раз с последующим усреднением результата (для избавления от случайных выбросов).
Для замера времени рекомендуется использовать вызовы функции omp_get_wtime или MPI_Wtime, общее время работы должно определяться временем самого медленного из процессов/нитей.
Количество ядер/процессоров рекомендуется задавать в виде p=2n, n=0, 1, 2, ... , k, где k определяется доступными ресурсами.
4) Определить основные причины недостаточной масштабируемости программы при максимальном числе используемых ядер/процессоров.
5) Сравнить эффективность OpenMP и MPI-версий параллельной программы.
6) Подготовить отчет о выполнении задания, включающий: описание реализованного алгоритма, графики зависимости времени исполнения от числа ядер/процессоров для различного объёма входных данных, текст программы.
*/

#include <iostream>
#include <vector>
#include <unistd.h>
#include <mpi.h>

using std::cout;
using std::vector;


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

Matrix createMatrixForTranspose(const Matrix& matrix) {
    MatrixSize transposed_size = swapMatrixSize(getMatrixSize(matrix));
    Matrix transposed_matrix = Matrix();
    transposed_matrix.resize(transposed_size.first);
    for (MatrixSizeItem row = 0; row < transposed_size.first; row++) {
        transposed_matrix[row].resize(transposed_size.second);
    }
    return transposed_matrix;
}

void transpose_MPI(const Matrix& matrix, Matrix &transposed_matrix, MatrixSizeItem start, MatrixSizeItem step) {
    MatrixSize transposed_size = getMatrixSize(transposed_matrix);
    for (MatrixSizeItem row_t = start; row_t < transposed_size.first; row_t += step) {
        for (MatrixSizeItem column_t = 0; column_t < transposed_size.second; column_t++) {
            transposed_matrix[row_t][column_t] = matrix[column_t][row_t];
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&аrgс, &аrgv);
    int rank, size, width;

    scanf(argv[1], "%lu", &width);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Matrix matrix, transposed_matrix;
    if (!rank) {
        matrix = getRandomMatrix({width, width}, -1e7, +1e7);
        transposed_matrix = createMatrixForTranspose(matrix);
    }
    MPI_Bcast(matrix, width * width, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    long double mpi_time = MPI_Wtime();
    transpose_MPI(matrix, transposed_matrix, rank, size);
    mpi_time = MPI_Wtime() - mpi_time;
    MPI_Barrier(MPI_COMM_WORLD);
    long double res_time;
    MPI_Reduce(mpi_time, res_time, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
        FILE *fout = fopen("log", "a");
        fprintf(fout, "%Lf", res_time);
        fclose(fout);
    }
    MPI_Finalize();
    return 0;
}
