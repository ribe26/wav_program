#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include "MatrixCalculation.h"
namespace {
	std::vector<vector<complex<double>>> calcDdashk(std::vector<std::complex<double>>c, std::vector<vector<complex<double>>> Dk) {
		std::vector<vector<complex<double>>> output;
		std::vector<vector<complex<double>>> diagC = toDiagonalMatrix(c);
		std::vector<vector<complex<double>>> diagConjC = toDiagonalConjMatrix(c);
		output = multiply2D2D(diagConjC, Dk);
		output = multiply2D2D(output, diagC);
		return output;
	}

	std::complex<double> calcLamda(std::vector<std::complex<double>> H, std::vector<std::vector<std::complex<double>>> Ddash) {
		std::complex<double> output;
		std::vector<std::complex<double>> conjH(H.size()); 
		conjH=vectorConj(H);
		
		std::vector<std::complex<double>> calc = multiply1D2D(conjH, Ddash);
		output = multiply1D1D(calc, H);
		
		//printComplexVector(conjH);
		//printMatrix2(Ddash);
		//printComplexVector(calc);
		//cout << output <<endl;
		return output;
	}

	void renewFilter(double mu, std::complex<double> lamda, std::vector<std::vector<std::complex<double>>> Ddash, std::vector<std::complex<double>>* H) {
		int N = Ddash.size();       // Number of rows in the matrix
		int M = Ddash[0].size();    // Number of columns in the matrix

		if (H->size() != M) {
			std::cerr << "Error: Incompatible dimensions for 1D-2D multiplication." << std::endl;
			exit(1);
		}


		std::vector<std::complex<double>> tempH = *H;
		

		for (int j = 0; j < N; ++j) {
			std::complex<double> calc = { 0,0 };
			for (int i = 0; i < M; i++) {
				calc+= tempH[i] * Ddash[j][i];
			}
			calc *= lamda*mu;
			(*H)[j] += calc;
		}
	}

}