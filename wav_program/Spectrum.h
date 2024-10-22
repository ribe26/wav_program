#pragma once
#include <vector>
#include <complex>
#include "wave_defs.h"

class Spectrum
{
public:
	Spectrum(class Signal signal);

	unsigned long original_length;//Signal����ϊ������ꍇ�C����Signal�̒������L�^���Ă���
	unsigned long data_length;//���_��FFT��
	unsigned long Fs;//�T���v�����O���g��
	
	//2n�v�f�ځF���� 2n+1�v�f��:����
	std::vector<complex<double>> dataL;
	std::vector<complex<double>> dataR;
	
	//�M����FFT����dataL��dataR�ۑ�
	void FFT(class Signal signal);
	
	//�p���[�̌v�Z
	double calc_power(long freq);
	
	//���f�������Ƃ�
	void Conj();
	//�U���X�y�N�g�����v���b�g
	void show();

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;
};