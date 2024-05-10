#pragma once
#include <vector>
#include <string>
#include "wave_defs.h"


class Signal
{
public:
	unsigned long Fs;
	std::vector<double> dataL;
	std::vector<double> dataR;
	

	//wav�t�@�C�����琶��
	Signal(const char* filename);
	//�X�y�N�g������IFFT�ɂ���Đ���
	//Signal(class Spectrum);
	
private:
	
	//wav�t�@�C���ǂݍ��݊֌W
	int readfmtChunk(FILE* fp, tWaveFormatPcm* waveFmtPcm);
	int wavHdrRead(const char* in_wavefile);

	int read8BitWavMonaural(FILE *fpIn);
	int read8BitWavStereo(FILE* fpIn);
	int read16BitWavMonaural(FILE* fpIn);
	int read16BitWavStereo(FILE* fpIn);

	int readDataSub(FILE* fpIn);
	int read(const char* in_wavfile);

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;

	unsigned long sizeOfData;
	unsigned long posOfData;
	unsigned long samples;//sizeOfData 

	unsigned short bytesPerSingleCh;





	
};

