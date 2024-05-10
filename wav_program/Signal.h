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
	

	//wavファイルから生成
	Signal(const char* filename);
	//スペクトルからIFFTによって生成
	//Signal(class Spectrum);
	
private:
	
	//wavファイル読み込み関係
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

