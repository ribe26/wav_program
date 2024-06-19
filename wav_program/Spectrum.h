#pragma once
#include <vector>
#include <complex>
#include "wave_defs.h"

class Spectrum
{
public:
	Spectrum(class Signal signal);

	unsigned long data_length;//何点のFFTか
	unsigned long Fs;
	//2n要素目：実部 2n+1要素目:虚部
	std::vector<double> dataL;
	std::vector<double> dataR;
	void FFT(class Signal signal);

	double calc_power(long freq);
	void Conj();

	void show();

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;
};