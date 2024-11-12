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
    std::vector<double> rir = generateReverbImpulsuse(length, 0.5, amplitude, samplingRate);
    Signal c(rir, samplingRate);
    Signal c2 = c;
    //Signal c("rir/usina_main_s1_p5.wav");
    //Signal c2("rir/usina_main_s1_p5.wav");
    //c.show();
    //cout << "start calculation" << endl;
   
    c2.squared();
    Spectrum C(c);
    Spectrum CT(c);
    c.show();
    c.show_MTF(20);
    CT.Conj();
    Spectrum C2(c2);
    //std::vector<std::vector<std::complex<double>>> diagCT = toDiagonalMatrix(CT.dataL);
    //std::vector<std::vector<std::complex<double>>> diagC = toDiagonalMatrix(C.dataL);

    double unitFs = samplingRate / C.dataL.size();

  
    //どの帯域を切り抜くか
    double startFs =1000;
    double endFs  = 8000;

    int startIdx = 0;
    int endIdx = 0;

    int FilterLength = 128;

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
        endIndexes.push_back(startIdx-1);
    }

    cout <<"reshaping Fs area" << (endIndexes.back() - startIndexes.front()) * unitFs << endl;
    
    double mu = 0.0001;
    int l = 1000;
    int maxFm = 20;

    int maxMTFidx = 0;
    while (maxMTFidx * unitFs < maxFm) {
        maxMTFidx++;
    }

    cout << maxMTFidx << endl;
    cout << "done" << endl;
   
    cout << "num of bands" << startIndexes.size()<<endl;
    std::vector<std::complex<double>>Filters;

    for (int bandIdx = 0; bandIdx < startIndexes.size(); bandIdx++) {
        cout << "band" << bandIdx << endl;
        std::vector<std::complex<double>>limitedC = limitVector(C.dataL, startIndexes[bandIdx], endIndexes[bandIdx]);

        std::vector<std::vector<std::vector<std::complex<double>>>> Ddash;
        for (int i = 0; i < maxMTFidx; i++) {
            vector<vector<std::complex<double>>> Dk = createLimitedD(C.dataL.size(), i, startIndexes[bandIdx], endIndexes[bandIdx]);
            vector<vector<std::complex<double>>> calc1 = calcDdashk(limitedC, Dk);
            Ddash.push_back(calc1);
        }
        cout << "done calc Ddash" << endl;

        std::vector<std::complex<double>> H(FilterLength, { 1.0,0.0 });
        for (int k = 0; k < 1; k++) {
            for (int i = 0; i < 1; i++) {
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

    std::vector<double> plot_x;
    std::vector<double> plot_y;
    int count = 0;
    for (int i = startIndexes.front(); i <= endIndexes.back(); i++) {
        plot_x.push_back(i * unitFs);
        C.dataL[i] *= Filters[count++];
    }

    for (int i = 0; i < plot_x.size(); i++) {
        plot_y.push_back( 20 * log10(abs(Filters[i])));
        //cout << "x:" << plot_x[i] <<"y:" << plot_y[i]<<endl;
    }

    cout << "plot_x" << plot_x.size();
    cout << "plot_y" << plot_y.size();
    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();


    Signal out(C);
    out.show();
    out.show_MTF(maxFm);

    //Spectrum H_spec2(Filters, c.Fs, Filters.size());
    //H_spec2.show();

    /*
    std::vector<std::complex<double>>limitedC = limitVector(C.dataL, startIdx, endIdx);
    std::vector<std::vector<std::vector<std::complex<double>>>> Ddash;
    for (int i = 0; i < maxMTFidx; i++) {
        vector<vector<std::complex<double>>> Dk = createLimitedD(C.dataL.size(),i,startIdx,endIdx);
        vector<vector<std::complex<double>>> calc1 = calcDdashk(limitedC,Dk);
        Ddash.push_back(calc1);
    }
    cout << "done calc Ddash"<<endl;

    std::vector<std::complex<double>> H(FilterLength, { 1.0,0.0 });
    
    
    for (int k = 0;k < 1;k++) {
        for (int i = 0; i < 1; i++) {
            std::complex<double> lamda = calcLamda(H, Ddash[i]);
            renewFilter(mu,lamda, Ddash[i], &H);
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