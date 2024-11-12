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

	std::complex<double> calclamda(std::vector<std::complex<double>> H, std::vector<std::vector<std::complex<double>>> Ddash) {
		std::complex<double> output;
		std::vector<std::complex<double>> conjH = vectorConj(H);
		
		std::vector<std::complex<double>> calc = multiply1D2D(conjH, Ddash);
		output = multiply1D1D(calc, H);
		return output;
	}
}