#include "Spectrum.h"
#include "Signal.h"
#include "FFT.h"
#include "matplotlibcpp.h"

// �X�y�N�g�����N���X�̃R���X�g���N�^
Spectrum::Spectrum(Signal signal) {
    Fs = signal.Fs;
    FFT(signal);
}

// FFT�ϊ����s���֐�
void Spectrum::FFT(Signal signal) {
    waveFormatpcm = signal.waveFormatpcm;
    waveFileheader = signal.waveFileheader;

    data_length = signal.dataL.size() / 2;
    int* ip = new int[2 + static_cast<int>(sqrt(0.5 * data_length * 2)) + 1];
    std::vector<double> w(data_length * 2 * 5 / 4);

    dataL = signal.dataL;
    dataR = signal.dataR;

    // �f�[�^����2�ׂ̂���łȂ����0���߂���
    int count = 1;
    if ((dataL.size() & (dataL.size() - 1)) != 0) {
        while (count < dataL.size()) {
            count *= 2;
        }
    }

    while (dataL.size() < count) {
        dataL.push_back(0.0);
        dataR.push_back(0.0);
    }

    ip[0] = 0;
    cdft(dataL.size(), -1, &dataL[0], ip, &w[0]);
    ip[0] = 0;
    cdft(dataR.size(), -1, &dataR[0], ip, &w[0]);

    delete[] ip;
}

// �X�y�N�g������\������֐�
void Spectrum::show() {
    std::cout << "Fs: " << Fs << std::endl;
    std::vector<double> plot_x(dataL.size() / 4);
    std::vector<double> plot_y(dataL.size() / 4);

    for (long i = 0; i < dataL.size() / 4; i++) {
        plot_x[i] = Fs / (dataL.size() / 2.0) * i;
        plot_y[i] = 20 * log10(sqrt(dataL[2 * i] * dataL[2 * i] + dataL[2 * i + 1] * dataL[2 * i + 1]));
    }

    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();
}

// ���f���������֐�
void Spectrum::Conj() {
    for (long i = 0; i < data_length; i++) {
        dataL[i * 2 + 1] *= -1;
        dataR[i * 2 + 1] *= -1;
    }
}

// �p���[���v�Z����֐�
double Spectrum::calc_power(long freq) {
    double out = 0;
    long count = 0;

    for (long i = 0; i < dataL.size() / 2.0; i++) {
        if (Fs / (dataL.size() / 2.0) * i >= freq) {
            break;
        }
        count++;
        out += dataL[2 * i] * dataL[2 * i] + dataL[2 * i + 1] * dataL[2 * i + 1];
    }

    return out / (count * count);
}
