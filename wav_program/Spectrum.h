#pragma once
#include <vector>
#include <complex>
#include "wave_defs.h"

class Spectrum
{
public:
	unsigned int Fs;
	std::vector<complex<double>> data;
	void FFT(class Signal signal);
private:
	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;

};

