#include "Spectrum.h"
#include "Signal.h"
#include "FFT.h"
#include "matplotlibcpp.h"


Spectrum::Spectrum(Signal signal) {
	Fs = 44100;
	FFT(signal);
}

void Spectrum::FFT(Signal signal) {

	waveFormatpcm = signal.waveFormatpcm;
	waveFileheader = signal.waveFileheader;
	
	data_length = signal.dataL.size()/2;
	int* ip = new int[2 + (int)sqrt(0.5 * data_length*2) + 1];
	vector<double> w(data_length*2 * 5 / 4);

	dataL = signal.dataL;
	dataR = signal.dataR;

	cdft(data_length*2, -1, &dataL[0], ip, &w[0]);
	cdft(data_length*2, -1, &dataR[0], ip, &w[0]);
	
}



void Spectrum::show() {
	vector<double> plot_x(data_length / 2);
	vector<double> plot_y(data_length / 2);


	for (long i = 0; i < data_length / 2; i++) {
		plot_x[i] = Fs / data_length * i;
		plot_y[i] = sqrt(dataL[2 * i] * dataL[2 * i] + dataL[2 * i + 1] * dataL[2 * i + 1]);
	}

	matplotlibcpp::plot(plot_x, plot_y);
	matplotlibcpp::show();
}

