#pragma once
#include <vector>
#include <complex>
#include <cmath>
namespace {
	std::complex<double> calcDdashk(std::vector<std::complex<double>>c,int k) {
		std::complex<double> output;
		int j = c.size() - k;
		if (j >= c.size()) {
			j -= c.size();
		}
		int i = 0;
		for (i;i < k;i++) {
			output+=(std::conj(c[i]) * c[j++]);
		}
		j = 0;
		while (i < c.size()) {
			//cout << i << "," << j << endl;
			output+=(std::conj(c[i++]) * c[j++]);
		}
		return output;
	}

	std::vector<std::complex<double>> calclamda(std::vector<std::complex<double>> H, std::vector<std::complex<double>> Ddash) {
		std::vector<std::complex<double>> output;
		double power = 0.0;
		for (int i = 0;i < H.size();i++) {
			output.push_back(conj(conj(H[i]) * H[i] * Ddash[i]));
			power += abs(output[i])*abs(output[i]);
		}
		
		for (int i = 0;i < H.size();i++) {
			output[i] /= power;
		}
		
		//output /= power;
		return output;
	}
}