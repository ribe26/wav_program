#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

namespace {

    //複素数ベクトルの特定区間を抜く
    std::vector<std::complex<double>>limitVector(std::vector<std::complex<double>> vec, int startIdx, int endIdx) {
        std::vector<std::complex<double>> output;
        for (int i = 0; i < vec.size(); i++) {
            if (i >= startIdx && i <= endIdx) {
                output.push_back(vec[i]);
            }
        }
        return output;
    }


    //D(k)を作成する
    std::vector<std::vector<std::complex<double>>> createD(int N, int k, int startIdx, int endIdx) {
        // Initialize N x N matrix with complex zeros
        std::vector<std::vector<std::complex<double>>> D(N, std::vector<std::complex<double>>(N, std::complex<double>(0.0, 0.0)));

        // Fill the bottom-left k x k submatrix with diagonal 1s (from top-left to bottom-right)
        for (int i = 0; i < k; ++i) {
            D[N - k + i][i] = std::complex<double>(1.0, 0.0); // Set diagonal elements to 1
        }

        // Fill the top-right (N-k) x (N-k) submatrix with diagonal 1s (from top-left to bottom-right)
        for (int i = 0; i < N-k; ++i) {
            D[i][k + i] = std::complex<double>(1.0, 0.0); // Set diagonal elements to 1
        }

        return D;
    }

    // 特定要素を切り抜いたDを作成する関数
    std::vector<std::vector<std::complex<double>>> createLimitedD(int N, int k, int startIdx, int endIdx) {
        int size = endIdx - startIdx + 1;
        std::vector<std::vector<std::complex<double>>> partialMatrix(size, std::vector<std::complex<double>>(size, std::complex<double>(0, 0)));

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                int globalRow = startIdx + i;
                int globalCol = startIdx + j;

                // 左下の k×k 単位行列の範囲
                if (globalRow >= N - k && globalCol < k && globalRow - (N - k) == globalCol) {
                    partialMatrix[i][j] ={ 1, 0 }; // 1+0i
                }
                // 右上の (N-k)×(N-k) 単位行列の範囲
                else if (globalRow < N - k && globalCol >= k && globalRow == globalCol - k) {
                    partialMatrix[i][j] = {1, 0}; // 1+0i
                }
            }
        }
            return partialMatrix;
    }

    //スペクトルを対角行列に変換
    std::vector<std::vector<std::complex<double>>> toDiagonalMatrix(const std::vector<std::complex<double>>& vec) {
        int N = vec.size(); // Size of the input vector
        // Initialize N*N matrix with 0
        std::vector<std::vector<std::complex<double>>> diagMatrix(vec.size(), std::vector<std::complex<double>>(vec.size(), std::complex<double>(0.0, 0.0)));

        // Set the diagonal elements
        for (int i = 0; i < vec.size();i++) {
            diagMatrix[i][i] = vec[i];
        }

        return diagMatrix;
    }

    //スペクトルを対角行列に変換(複素共役転置版)
    std::vector<std::vector<std::complex<double>>> toDiagonalConjMatrix(const std::vector<std::complex<double>>& vec) {
        int N = vec.size(); // Size of the input vector
        // Initialize N*N matrix with 0
        std::vector<std::vector<std::complex<double>>> diagMatrix(vec.size(), std::vector<std::complex<double>>(vec.size(), std::complex<double>(0.0, 0.0)));

        // Set the diagonal elements
        for (int i = 0; i < vec.size(); i++) {
            diagMatrix[i][i] = std::conj(vec[i]);
        }
        return diagMatrix;
    }

    std::vector<std::complex<double>> vectorConj(std::vector<std::complex<double>>& vec) {
        std::vector<std::complex<double>> output(vec.size());
        for (int i = 0;i<vec.size();i++) {
            output[i]=std::conj(vec[i]);
        }
        return output;
    }

    // 行列を転置する
    std::vector<std::vector<std::complex<double>>> transposeMatrix(const std::vector<std::vector<std::complex<double>>>& mat) {
        int rows = mat.size();
        int cols = mat[0].size();
        std::vector<std::vector<std::complex<double>>> transposed(cols, std::vector<std::complex<double>>(rows));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                transposed[j][i] = mat[i][j];
            }
        }

        return transposed;
    }


    std::vector<std::complex<double>> multiply1D2D(
        const std::vector<std::complex<double>>& vec,
        const std::vector<std::vector<std::complex<double>>>& mat) {

        int N = mat.size();       // Number of rows in the matrix
        int M = mat[0].size();    // Number of columns in the matrix

        if (vec.size() != N) {
            std::cerr << "Error: Incompatible dimensions for 1D-2D multiplication." << std::endl;
            exit(1);
        }

        std::vector<std::complex<double>> result(M, std::complex<double>(0.0, 0.0));

        for (int j = 0; j < M; j++) {
            for (int i = 0; i < N; i++) {
                result[j] += vec[i] * mat[i][j];
            }
        }

        return result;
    }

    std::vector<std::complex<double>> multiply2D1D(
        const std::vector<std::vector<std::complex<double>>>& mat,
        const std::vector<std::complex<double>>& vec) {

        int N = mat.size();       // Number of rows in the matrix
        int M = mat[0].size();    // Number of columns in the matrix

        if (vec.size() != M) {
            std::cerr << "Error: Incompatible dimensions for 1D-2D multiplication." << std::endl;
            exit(1);
        }

        std::vector<std::complex<double>> result(N, std::complex<double>(0.0, 0.0));

        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < M; i++) {
                result[j] += vec[i] * mat[j][i];
            }
        }

        return result;
    }

    std::complex<double> multiply1D1D(
        const std::vector<std::complex<double>>& vec1,
        const std::vector<std::complex<double>>& vec2) {

        if (vec1.size() != vec2.size()) {
            std::cerr << "Error: Incompatible dimensions for 1D-2D multiplication." << std::endl;
            exit(1);
        }

        std::complex<double> output = { 0,0 };

        for (int j = 0; j < vec1.size(); j++) {
            output += vec1[j] * vec2[j];
        }

        return output;
    }

    template<typename T>
    std::vector<std::vector<T>> multiply2D2D(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {
        int rows1 = mat1.size();              // Number of rows in mat1
        int cols1 = mat1[0].size();           // Number of columns in mat1
        int rows2 = mat2.size();              // Number of rows in mat2
        int cols2 = mat2[0].size();           // Number of columns in mat2

        // Check if matrix multiplication is possible
        if (cols1 != rows2) {
            std::cerr << "Error: Incompatible matrix dimensions for multiplication." << std::endl;
            exit(1);
        }

        // Initialize the result matrix with zeros
        std::vector<std::vector<T>> result(rows1, std::vector<T>(cols2, T(0)));

        // Matrix multiplication
        for (int i = 0; i < rows1; ++i) {
            for (int j = 0; j < cols2; ++j) {
                for (int k = 0; k < cols1; ++k) {
                    result[i][j] += mat1[i][k] * mat2[k][j];
                }
            }
        }

        return result;
    }

    // Function to print the matrix
    void printMatrix2(const std::vector<std::vector<std::complex<double>>>& matrix) {
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }


    // Function to print the matrix
    void printMatrix(const std::vector<std::vector<int>>& matrix) {
        for (const auto& row : matrix) {
            for (int val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }

    // Function to print a 1D complex vector
    void printComplexVector(const std::vector<std::complex<double>>& vec) {
        for (const auto& val : vec) {
            std::cout << "(" << val.real() << ", " << val.imag() << "i) ";
        }
        std::cout << std::endl;
    }


    // Function to print the sine wave
    void printVector(const std::vector<double>& vec) {
        for (const auto& val : vec) {
            std::cout << val << " " << std::endl;
        }
        std::cout << std::endl;
    }

}