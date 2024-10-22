#pragma once
#include <vector>
#include <complex>
#include "wave_defs.h"

class Spectrum
{
public:
	Spectrum(class Signal signal);

	unsigned long original_length;//Signalから変換した場合，元のSignalの長さを記録しておく
	unsigned long data_length;//何点のFFTか
	unsigned long Fs;//サンプリング周波数
	
	//2n要素目：実部 2n+1要素目:虚部
	std::vector<complex<double>> dataL;
	std::vector<complex<double>> dataR;
	
	//信号をFFTしてdataLとdataR保存
	void FFT(class Signal signal);
	
	//パワーの計算
	double calc_power(long freq);
	
	//複素共役をとる
	void Conj();
	//振幅スペクトルをプロット
	void show();

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;
};