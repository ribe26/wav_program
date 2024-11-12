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

int main()
{  
    
    // Parameters for the sine wave
    int length = 48000;         // Length of the sine wave
    double amplitude = 1.0;   // Amplitude of the sine wave
    double frequency = 5.0;   // Frequency of the sine wave in Hz
    double samplingRate = 48000.0; // Sampling rate in Hz


    /*
    // Generate the sine wave
    std::vector<double> sineWave = generateSineWave(length, amplitude, frequency, samplingRate);
    std::vector<double> rir = generateReverbImpulsuse(length, 0.5, amplitude, samplingRate);
    //Signal sig(rir, samplingRate);
    //sig.show();
    //cout << rir[0]<<endl;
    Signal c(sineWave,samplingRate);
    Spectrum C(c);
    c.show();
    C.show();
    //Signal c2 = c;
    */    
    std::vector<double> rir = generateReverbImpulsuse(length, 0.01, amplitude, samplingRate);
    Signal c(rir, samplingRate);
    Signal c2 = c;
    //Signal c("rir/usina_main_s1_p5.wav");
    //Signal c2("rir/usina_main_s1_p5.wav");
    //c.show();
    //cout << "start calculation" << endl;
   
    c2.squared();
    Spectrum C(c);
    Spectrum CT(c);
    CT.Conj();
    Spectrum C2(c2);
    //std::vector<std::vector<std::complex<double>>> diagCT = toDiagonalMatrix(CT.dataL);
    //std::vector<std::vector<std::complex<double>>> diagC = toDiagonalMatrix(C.dataL);

    double unitFs = samplingRate / C.dataL.size();

  
    //どの帯域を切り抜くか
    double startFs =2000;
    double endFs  = 3000;

    int startIdx = 0;
    int endIdx = 0;

    while (startIdx * unitFs < startFs) {
        startIdx++;
    }

    while (endIdx * unitFs < endFs) {
        endIdx++;
    }
    std::cout <<unitFs << "," << startIdx << "," << endIdx;

    int FilterLength = 100;
    endIdx = startIdx + FilterLength-1;


    
    double mu = 0.0001;
    int l = 20;
    int maxFm = 20;

    int maxMTFidx = 0;
    while (maxMTFidx * unitFs < maxFm) {
        maxMTFidx++;
    }

    cout << maxMTFidx;
   

    std::limitedC = limitVector(startIdx,endIdx);


    
    std::vector < std::vector < ::vector<std::complex<double>>>> Ddash;
    for (int i = 0; i < maxMTFidx; i++) {
        cout << i << endl;
        vector<vector<std::complex<double>>> Dk = createD(C.dataL.size(),i,startIdx,endIdx);
        vector<vector<std::complex<double>>> calc1 = calcDdashk(C.dataL,Dk);
        Ddash.push_back(calc1);
    }
    cout << "done calc Ddash";
    


    std::vector<std::complex<double>> H(C.dataL.size(), { 1.0,0.0 });
    
    /*
    for (int k = 0;k < l;k++) {
        for (int i = 0; i < maxF; i++) {
            std::vector<std::complex<double>> lamda = calclamda(H, Ddash);
            //if (k < 10) { cout << lamda << endl; }
            double power = 0;
            if (k == 0 && i==0) {
                for (int h = 0;h < 10;h++) {
                    cout <<lamda[h] << endl;
                }
            }
            
            for (int j = 0;j < H.size();j++) {
                H[j] += mu * lamda[j] * Ddash[j]* H[j];
                power += abs(H[j]);
            }
            for (int j = 0;j < H.size();j++) {
                H[j] /=power;
            }

        }
    }
    //H[0] = { 1.0,0 };
    Spectrum H_spec2(H, c.Fs, H.size());
    H_spec2.show();
    
    */

    /*
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
    */
    
}