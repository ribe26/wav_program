// wav_program.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "Signal.h"
#include "math.h"
#include "FFT.h"
#include "matplotlibcpp.h"
#include <vector>
#include "Spectrum.h"
#include "MatrixCalculation.h"
#include "utility_functions.h"

int main()
{  
    // Parameters for the sine wave
    int length = 128;         // Length of the sine wave
    double amplitude = 1.0;   // Amplitude of the sine wave
    double frequency = 5.0;   // Frequency of the sine wave in Hz
    double samplingRate = 100.0; // Sampling rate in Hz


    // Generate the sine wave
    std::vector<double> sineWave = generateSineWave(length, amplitude, frequency, samplingRate);
    std::vector<double> sineWave2(sineWave.size(), 0.0);
    for (int i = 0; i < sineWave.size(); i++) {
        sineWave2[i] = sineWave[i] * sineWave[i];
    }

    Signal c(sineWave,samplingRate);
    Signal c2(sineWave2,samplingRate);
    Spectrum C(c);
    Spectrum CT(c);
    CT.Conj();
    Spectrum C2(c2);

    vector<complex<double>>G2;

    for (int i = 0; i < C.dataL.size(); i++) {
        std::vector<std::vector<std::complex<double>>> Dk = createD(C.dataL.size(), i);
        //printMatrix2(Dk);
        std::vector<std::complex<double>> calc = multiply1D2D(CT.dataL, Dk);
        std::complex<double> calc2 = multiply1D1D(calc, C.dataL);
        G2.push_back(calc2);
    }

    double abs_max_G2 = 0;
    double abs_max_C2 = 0;
    for (int i = 0; i < G2.size(); i++) {
        if (abs(G2[i]) > abs_max_G2) {
            abs_max_G2 = abs(G2[i]);
        }
        if (abs(C2.dataL[i]) > abs_max_C2) {
            abs_max_C2 = abs(C2.dataL[i]);
        }
    }

    std::vector <double> G2_abs;
    std::vector <double> C2_abs;
    for (int i = 0; i < G2.size(); i++) {
        G2_abs.push_back(abs(G2[i]));
        C2_abs.push_back(abs(C2.dataL[i]));

    }

    std::cout << G2.size() << std::endl;
    std::cout << C2.dataL.size() << std::endl;

    for (int i = 0; i < G2.size(); i++) {
        G2_abs[i] /= abs_max_G2;
        C2_abs[i] /= abs_max_C2;
        std::cout << G2_abs[i] << "," << C2_abs[i] << std::endl;
    }
}