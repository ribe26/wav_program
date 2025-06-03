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
#include "thesis_functions.h"
#include "algorithm"

int main()
{

    // Parameters for the sine wave
    int length = 65536;         // Length of the sine wave
    double amplitude = 1.0;   // Amplitude of the sine wave
    double frequency = 5.0;   // Frequency of the sine wave in Hz
    double samplingRate = 48000.0; // Sampling rate in Hz


    //std::vector<double> rir = generateReverbImpulsuse(length, 0.5, amplitude, samplingRate);
    //Signal c(rir, samplingRate);
    //Signal c2 = c;

    std::vector<double> rir = generateReverbImpulsuse(length, 1.3, amplitude, samplingRate);
    /*
    std::vector<double> sine1 = generateSineWave(length, 10.0, 5, samplingRate);
    std::vector<double> sine2 = generateSineWave(length, 1.0, 50, samplingRate);
    std::vector<double> sine3 = generateSineWave(length, 0.25, 100, samplingRate);
    std::vector<double> sine4 = generateSineWave(length, 0.125, 500, samplingRate);
    std::vector<double> rir(length);
    for (int i = 0; i < length; i++) {
        rir[i] = sine1[i] + sine2[i] + sine3[i] + sine4[i];
    }
    */

    Signal c_original("rir/usina_main_s1_p5.wav");
    Signal c("rir/usina_main_s1_p5.wav");
    
    //Signal c_original(rir,samplingRate);
    //Signal c(rir, samplingRate);
    c.show();
    c.show_MTF(20);
    c.calc_MTF(20, "original.txt");
    Spectrum C(c);
    Spectrum C_ORIGINAL(c_original);
    C.show();
    double original_energy = C.calc_energy();
    C.show_MTF(20);

    length = c.dataL.size();
    samplingRate = c.Fs;


    complex<double> init_value(1.0, 0.0);
    std::vector<complex<double>>h(C.dataL.size(), init_value);
    Spectrum H(h, samplingRate, length);
    H.show();
    std::cout << "H power:" << H.calc_power() << endl;
    std::cout << "H energy:" << H.calc_energy() << endl;


    double endMTFfreq = 20;
    double unitFs = c.Fs / (double)c.dataL.size();
    int endIdx = endMTFfreq / unitFs;


    double bandStartFreq = 2000;
    double bandEndFreq = 16000;
    int bandStartIdx = bandStartFreq / unitFs;
    int bandEndIdx = bandEndFreq / unitFs;
    /*
    for (int i = 0; i < C.dataL.size()/2; i++) {
        if (i<bandStartIdx || i>bandEndIdx) {
            C.dataL[i] = 0;
            C.dataL[C.dataL.size() - i] = 0;
        }
    }
    */
    //original_energy = C.calc_energy();




    //vector<int> target_index = {15,16,17,18,19,20,21,22,23,24,25,26,27};
    std::vector<int> target_index(endIdx + 1);  // サイズn+1のvectorを作成
    std::iota(target_index.begin(), target_index.end(), 0);


    complex<double> complex_zero(0.0, 0.0);
    std::vector<complex<double>>H_step(C.dataL.size(), complex_zero);
    std::vector<complex<double>>G_temp(C.dataL.size(), complex_zero);
    Spectrum G(G_temp, samplingRate, length);

    int iteration = 1000;
    double step = 0.000000000001;

    for (int j = 0; j < iteration; j++) {
        std::cout << "iteration:" << j << endl;
        //cout << "H power:" << H.calc_energy() << endl;
        std::fill(G_temp.begin(), G_temp.end(), complex_zero);
        for (int k = 0; k < target_index.size(); k++) {

            for (int i = 0; i < G.dataL.size(); i++) {
                G.dataL[i] = C.dataL[i] * H.dataL[i];
                //G.dataL[G.dataL.size() - i] = conj(G.dataL[i]);
            }
            G.set_energy(original_energy);

            
            for (int m = 0; m < C.dataL.size(); m++) {
                int idx = (target_index[k] - m + C.dataL.size()) % C.dataL.size();
                //cout << "(m,idx):(" << m << "," << idx << ")" << endl;
                G_temp[target_index[k]] += (G.dataL[m] * G.dataL[idx]);
            }
            

            G_temp[target_index[k]]/= (double)G.dataL.size();
            //G_temp[target_index[k]] = G.dataL[target_index[k]];
            double div = (C.dataL.size() * abs(G_temp[target_index[k]]));

            for (int p = 0; p < C.dataL.size()/2; p++) {
                int idx = (target_index[k] - p + C.dataL.size()) % C.dataL.size();

                H_step[p] += step * (G_temp[target_index[k]] * conj(C.dataL[p]) * conj(C.dataL[idx]) * conj(H.dataL[idx])) / div;
                //H_step[p] += step * conj(G_temp[k]) * C.dataL[p] * C.dataL[idx] * H.dataL[idx];
            }
        }
        for (int p = 0; p < C.dataL.size()/2; p++) {
            H.dataL[p] += H_step[p];
            H.dataL[H.dataL.size() - p] = conj(H.dataL[p]);
        }
        H.normalize_power();
        std::fill(H_step.begin(), H_step.end(), complex_zero);
    }

    std::cout << "H power:" << H.calc_power() << endl;
    std::cout << "H energy:" << H.calc_energy() << endl;
    H.show();


    cout << "original length:" << C.dataL.size() << endl;
    cout << "filter length:" << H.dataL.size() << endl;
    for (int p = 0; p < C.dataL.size(); p++) {
        C.dataL[p] = C_ORIGINAL.dataL[p] * H.dataL[p];
        if (p > C.dataL.size()/2-200 && p<C.dataL.size()/2) { cout << H.dataL[p] << endl; }
    }
    C.set_energy(C_ORIGINAL.calc_energy());
    C.show();
    C.show_MTF(20);

    Signal c_inv(C);
    c_inv.show();
    c_inv.show_MTF(20);
    c_inv.calc_MTF(20, "filtered.txt");


    double rt60original = calculateRT60(c_original.dataL, samplingRate);
    cout << "rt60 c original:" << rt60original << endl;

    double rt60filtered = calculateRT60(c_inv.dataL, samplingRate);
    cout << "rt60 c filtered" << rt60filtered << endl;

}


