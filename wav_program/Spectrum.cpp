#include "Spectrum.h"
#include "Signal.h"
#include "FFT.h"
#include "matplotlibcpp.h"
#include "math.h"
// スペクトラムクラスのコンストラクタ
Spectrum::Spectrum(Signal signal) {
    Fs = signal.Fs;
    original_length = signal.dataL.size();
    FFT(signal);
}

Spectrum::Spectrum(vector<complex<double>> data, double samplingRate, long original_len) {
    original_length = original_len;
    Fs = samplingRate;
    for (int i = 0;i < data.size();i++) {
        dataL.push_back(data[i]);
        dataR.push_back(data[i]);
        //if (i < 10) { cout << "data:" << data[i]; }
    }
}
// FFT変換を行う関数
void Spectrum::FFT(Signal signal) {
    waveFormatpcm = signal.waveFormatpcm;
    waveFileheader = signal.waveFileheader;

    data_length = signal.dataL.size();
    int* ip = new int[2 + static_cast<int>(sqrt(0.5 * data_length * 2)) + 1];
    std::vector<double> w(data_length * 2 * 5 / 4);
    

    vector<double> tempL;
    vector<double> tempR;

    

    for (int i = 0; i < signal.dataL.size();i++) {
        tempL.push_back(signal.dataL[i]);
        tempL.push_back(0.0);
        tempR.push_back(signal.dataR[i]);
        tempR.push_back(0.0);
    }

    // データ長が2のべき乗でなければ0埋めする
    int count = 1;
    if ((tempL.size() & (tempL.size() - 1)) != 0) {
        while (count < tempL.size()) {
            count *= 2;
        }
    }
    
    while (tempL.size() < count) {
        tempL.push_back(0.0);
        tempR.push_back(0.0);
    }
    
    
    ip[0] = 0;
    cdft(tempL.size(), -1, &tempL[0], ip, &w[0]);
    ip[0] = 0;
    cdft(tempR.size(), -1, &tempR[0], ip, &w[0]);

    delete[] ip;

    for (int i = 0; i < tempL.size() / 2; i++) {
        dataL.push_back({tempL[i*2],tempL[i*2+1]});
        dataR.push_back({ tempR[i * 2],tempR[i * 2 + 1] });
    }
}


// スペクトラムを表示する関数
void Spectrum::show() {
    //std::cout << "Fs: " << Fs << std::endl;
    std::vector<double> plot_x(dataL.size()/2);
    std::vector<double> plot_y(dataL.size()/2);

    for (int i = 0; i < dataL.size() / 2; i++) {
        plot_x[i] = Fs / (double)dataL.size() *i;
        plot_y[i] = 20 * log10(abs(dataL[i]));
        //cout << "x:" << plot_x[i] <<"y:" << plot_y[i]<<endl;
    }

    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();
}

// 複素共役を取る関数
void Spectrum::Conj() {

    for (long i = 0; i < data_length; i++) {
        dataL[i] = conj(dataL[i]);
        dataR[i] = conj(dataR[i]);
    }
}

// パワーを計算する関数
double Spectrum::calc_power(long freq) {
    double out = 0;
    long count = 0;

    for (long i = 0; i < dataL.size(); i++) {
        if (Fs / (dataL.size()) * i >= freq) {
            break;
        }
        count++;
        out += abs(dataL[i])* abs(dataL[i]);
    }

    return out / (count * count);
}
