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


    vector<double> tempL;
    vector<double> tempR;


    for (int i = 0; i < signal.dataL.size(); i++) {
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
        tempL.push_back(0.0);
        tempR.push_back(0.0);
    }


    data_length = tempL.size();
    int* ip = new int[2 + static_cast<int>(sqrt(0.5 * data_length * 2)) + 1];
    std::vector<double> w(data_length * 2 * 5 / 4);
   

    
    
    
    ip[0] = 0;
    cdft(tempL.size(), -1, &tempL[0], ip, &w[0]);
    ip[0] = 0;
    cdft(tempR.size(), -1, &tempR[0], ip, &w[0]);

    delete[] ip;

    for (int i = 0; i < tempL.size() / 2; i++) {
        dataL.push_back({tempL[i*2],tempL[i*2+1]});
        dataR.push_back({ tempR[i * 2],tempR[i * 2 + 1] });
    }
    data_length = dataL.size();
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
double Spectrum::calc_power() {
    double out = 0;
    for (long i = 0; i < dataL.size(); i++) {
        out += abs(dataL[i])* abs(dataL[i]);
    }

    return out / dataL.size();
}


double Spectrum::calc_energy() {
    double out = 0;
    for (long i = 0; i < dataL.size(); i++) {
        out += abs(dataL[i]) * abs(dataL[i]);
    }

    return out;
}

// パワーを正規化する関数
void Spectrum::normalize_power() {
    double power = this->calc_energy();
    double ratio = sqrt((double)dataL.size()/power);
    for (long i = 0; i < dataL.size(); i++) {
        this->dataL[i] *= ratio;
    }
}

void Spectrum::set_energy(double energy) {
    double now_energy = this->calc_energy();
    double ratio = sqrt(energy / now_energy);
    for (long i = 0; i < dataL.size(); i++) {
        this->dataL[i] *= ratio;
    }
}



void Spectrum::squared() {
    std::vector<complex<double>> tempL(data_length,complex<double>{0,0});
    std::vector<complex<double>> tempR(data_length, complex<double>{0, 0});

    for (int k = 0; k < data_length/2; k++) {
        for (int m = 0; m < data_length; m++) {
            int idx = (k - m + data_length) % data_length;
            //cout << "(m,idx):(" << m << "," << idx << ")" << endl;
            tempL[k] += (dataL[m] * dataL[idx]);
            tempR[k] += (dataR[m] * dataR[idx]);
        }
        tempL[k] /= data_length;
        tempR[k] /= data_length;
    }

    for (int k = 0; k < data_length/2; k++) {
        //cout << k << "," << (data_length - 1) - k<<endl;
        dataL[k] = tempL[k];
        dataR[k] = tempR[k];
        dataL[data_length- k] = conj(tempL[k]);
        dataR[data_length- k] = conj(tempR[k]);
    }

}


std::vector<complex<double>> Spectrum::squared_limit_freq(double freq) {
    double unitFs = (Fs / (double)data_length);
    int maxIdx = freq / unitFs;
    std::vector<complex<double>> tempL(data_length, complex<double>{0, 0});
    for (int k = 0; k <= maxIdx; k++) {
        for (int m = 0; m < data_length; m++) {
            int idx = (k - m + data_length) % data_length;
            //cout << "(m,idx):(" << m << "," << idx << ")" << endl;
            tempL[k] += (dataL[m] * dataL[idx]);
        }
        tempL[k] /= data_length;
    }
    return tempL;
}

void Spectrum::show_MTF(double freq) {
    double power = this->calc_power();
    std::vector<double> plot_x;
    std::vector<double> plot_y;
    std::vector<complex<double>> mtf=this->squared_limit_freq(freq);


    double unitFs = (Fs / (double)data_length);
    int maxIdx = freq / unitFs;
 

    cout << "freq:" << freq << endl;
    cout << "Fs:" << Fs << endl;
    cout << "unitFs:" << unitFs << endl;
    cout << "data_length:" << data_length << endl;
    cout << "maxIdx:" << maxIdx << endl;
    cout << "power:" << power << endl;


    for (int i = 0; i <= maxIdx; i++) {
        plot_x.push_back(unitFs*i);
        plot_y.push_back(abs(mtf[i]) / power);
    }
    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();
}