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
		//output /= abs(output);
		output = conj(output);
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
		

		double power = 0;
		for (int i = 0; i < H->size(); i++) {
			power += abs((*H)[i]);
		}

		double ratio = power;
		for (int i = 0; i < H->size(); i++) {
			(*H)[i]/=ratio;
		}
		

	}



    void deReverbe(std::vector<double> rir, double startFs, double endFs, int FilterLength,double minFm,double RT) {
        int length = 48000;         // Length of the sine wave
        double samplingRate = 48000.0; // Sampling rate in Hz


        Signal c(rir, samplingRate);
        Signal c2 = c;

        char filenameC[100];
        sprintf_s(filenameC,sizeof(filenameC), "original%f_%f_%d_%f_%f.txt", startFs, endFs, FilterLength, minFm, RT);
        string filename(filenameC);
        c.calc_MTF(20,filename);


        double rt60original = calculateRT60(c.dataL, samplingRate);
        c2.squared();
        Spectrum C(c);
        Spectrum CT(c);
        CT.Conj();
        Spectrum C2(c2);

        //単位インデックスあたりの周波数
        double unitFs = samplingRate / C.dataL.size();
        //どの帯域を切り抜くか
        int startIdx = 0;
        int endIdx = 0;
        while (startIdx * unitFs < startFs) {
            startIdx++;
        }

        while (endIdx * unitFs < endFs) {
            endIdx++;
        }
        vector<int> startIndexes;
        vector<int> endIndexes;
        while (startIdx < endIdx) {
            startIndexes.push_back(startIdx);
            startIdx += FilterLength;
            endIndexes.push_back(startIdx - 1);
        }

        //MTFを考慮する範囲
        int maxFm = 20;
        int minMTFidx = 0;
        while (minMTFidx * unitFs < minFm) {
            minMTFidx++;
        }

        int maxMTFidx = 0;
        while (maxMTFidx * unitFs < maxFm) {
            maxMTFidx++;
        }
        


        //cout << c.dataL.size() <<"-,-" << startIndexes[0] << "," << endIndexes[0] << "," << minMTFidx << "," << maxMTFidx << endl;
        //フィルタ形成パラメータ
        double mu = 0.1;
        int l = 20;
            std::vector<std::complex<double>>Filters;

            for (int bandIdx = 0; bandIdx < startIndexes.size(); bandIdx++) {
                //cout << "band" << bandIdx << endl;
                std::vector<std::complex<double>>limitedC = limitVector(C.dataL, startIndexes[bandIdx], endIndexes[bandIdx]);

                std::vector<std::vector<std::vector<std::complex<double>>>> Ddash;
                for (int i = minMTFidx; i < maxMTFidx; i++) {
                    vector<vector<std::complex<double>>> Dk = createLimitedD(C.dataL.size(), i, startIndexes[bandIdx], endIndexes[bandIdx]);
                    vector<vector<std::complex<double>>> calc1 = calcDdashk(limitedC, Dk);
                    Ddash.push_back(calc1);
                }
                //cout << "done calc Ddash" << endl;
                std::vector<std::complex<double>> H(FilterLength, { 1.0,0.0 });
                for (int k = 0; k < l; k++) {
                    for (int i = 0; i < maxMTFidx - minMTFidx; i++) {
                        std::complex<double> lamda = calcLamda(H, Ddash[i]);
                        renewFilter(mu, lamda, Ddash[i], &H);
                    }
                }

                if (bandIdx == 0) {
                    Filters = H;
                }
                else {
                    std::copy(H.begin(), H.end(), std::back_inserter(Filters));
                }
            }


            std::vector<std::complex<double>> Purpose_Filter(C.dataL.size(), { 1.0,0.0 });
            int cnt = 0;
            for (int i = startIndexes.front(); i <= endIndexes.back(); i++) {
                Purpose_Filter[i] = Filters[cnt];
                Purpose_Filter[Purpose_Filter.size() - i - 1] = Filters[cnt];
                cnt++;
            }

            double power = 0;
            for (int i = 0; i < Purpose_Filter.size(); i++) {
                power += abs(Purpose_Filter[i]);
            }

            double ratio = power;
            for (int i = 0; i < Purpose_Filter.size(); i++) {
                Purpose_Filter[i] /= ratio;
            }



            for (int i = 0; i < C.dataL.size(); i++) {
                C.dataL[i] *= Purpose_Filter[i];
            }
            Signal out(C);
            out.normalize();

            double rt60filtered = calculateRT60(out.dataL, samplingRate);
            
            sprintf_s(filenameC, sizeof(filenameC), "filtered%f_%f_%d_%f_%f.txt", startFs, endFs, FilterLength, minFm, RT);
            string filename2(filenameC);
            out.calc_MTF(20, filename2);


            // double startFs, double endFs, double FilterLength,double minFm
            if (rt60original > rt60filtered) {
                cout << "-----------------------------" << endl;
                cout << "startFs" << startFs<<endl;
                cout << "endFs" << endFs << endl;
                cout << "FilterLength" << FilterLength << endl;
                cout << "minFm" << minFm << endl;
                cout << "RT" << RT << endl;
                cout << "-----------------------------" << endl;
            }
        }
}
