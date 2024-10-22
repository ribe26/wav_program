#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

namespace {
    //D(k)のを作成する
    std::vector<std::vector<std::complex<double>>> createD(int N, int k) {
        // Initialize N x N matrix with complex zeros
        std::vector<std::vector<std::complex<double>>> D(N, std::vector<std::complex<double>>(N, std::complex<double>(0.0, 0.0)));

        // Fill the bottom-left k x k submatrix with diagonal 1s (from top-left to bottom-right)
        for (int i = 0; i < k; ++i) {
            D[N - k + i][i] = std::complex<double>(1.0, 0.0); // Set diagonal elements to 1
        }

        // Fill the top-right (N-k) x (N-k) submatrix with diagonal 1s (from top-left to bottom-right)
        for (int i = 0; i < (N - k); ++i) {
            D[i][k + i] = std::complex<double>(1.0, 0.0); // Set diagonal elements to 1
        }

        return D;
    }

    //スペクトルを対角行列に変換
    std::vector<std::vector<std::complex<double>>> toDiagonalMatrix(const std::vector<std::complex<double>>& vec) {
        int N = vec.size(); // Size of the input vector
        // Initialize N*N matrix with 0
        std::vector<std::vector<std::complex<double>>> diagMatrix(N, std::vector<std::complex<double>>(N, std::complex<double>(0.0, 0.0)));

        // Set the diagonal elements
        for (int i = 0; i < N; ++i) {
            diagMatrix[i][i] = vec[i];
        }

        return diagMatrix;
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