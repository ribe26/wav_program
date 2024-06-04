#include "Spectrum.h"
#include "Signal.h"
#include "FFT.h"
#include "matplotlibcpp.h"


Spectrum::Spectrum(Signal signal) {
	Fs = signal.Fs;
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


	//ƒf[ƒ^’·‚³‚ª2‚Ì‚×‚«æ‚Å‚È‚¯‚ê‚Î‚O–„‚ß‚·‚é
	int count = 1;
	if ((dataL.size() & dataL.size() - 1) != 0) {
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
	
}



void Spectrum::show() {
	vector<double> plot_x(data_length / 2);
	vector<double> plot_y(data_length / 2);


	for (long i = 0; i < data_length / 2; i++) {
		plot_x[i] = Fs / data_length * i;
		plot_y[i] = 20*log10(sqrt(dataL[2 * i] * dataL[2 * i] + dataL[2 * i + 1] * dataL[2 * i + 1]));
	}

	matplotlibcpp::plot(plot_x, plot_y);
	matplotlibcpp::show();
}

void Spectrum::Conj() {
	for (long i = 0;i < data_length;i++) {
		dataL[i * 2 + 1] = dataL[i * 2 + 1] * -1;
		dataR[i * 2 + 1] = dataR[i * 2 + 1] * -1;
	}
}

