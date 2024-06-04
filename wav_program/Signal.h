#pragma once
#include <vector>
#include <string>
#include "wave_defs.h"


class Signal
{
public:
	unsigned int Fs;
	std::vector<double> dataL;
	std::vector<double> dataR;
	

	//wav�t�@�C�����琶��
	Signal(const char* filename);

	//�X�y�N�g������IFFT�ɂ���Đ���
	Signal(class Spectrum);
	
	//�f�[�^���琶��
	Signal(vector<double>,double F);

	void show();

	void normalize();

	//wav�t�@�C���ǂݍ��݊֌W---------------------------------------------------------
	int readfmtChunk(FILE* fp, tWaveFormatPcm* waveFmtPcm);
	int wavHdrRead(const char* in_wavefile);

	int read8BitWavMonaural(FILE *fpIn);
	int read8BitWavStereo(FILE* fpIn);
	int read16BitWavMonaural(FILE* fpIn);
	int read16BitWavStereo(FILE* fpIn);

	int readDataSub(FILE* fpIn);
	int read(const char* in_wavfile);


	int showWavdata();

	//wave�t�@�C����������--------------------------------------------------------------
	long wavHeaderWrite(FILE *fpIn);
	int write16BitWav(FILE* fpOut);
	int write(const char* outFile);

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;

	unsigned long sizeOfData;
	unsigned long posOfData;
	unsigned long samples;//sizeOfData 

	unsigned short bytesPerSingleCh;





	
};

