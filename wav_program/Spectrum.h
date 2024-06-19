#pragma once
#include <vector>
#include <complex>
#include "wave_defs.h"

class Spectrum
{
public:
	Spectrum(class Signal signal);

	unsigned long data_length;//���_��FFT��
	unsigned long Fs;
	//2n�v�f�ځF���� 2n+1�v�f��:����
	std::vector<double> dataL;
	std::vector<double> dataR;
	void FFT(class Signal signal);

	double calc_power(long freq);
	void Conj();

	void show();

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;
};