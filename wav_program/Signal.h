#pragma once
#include <vector>
#include <string>
#include "wave_defs.h"

class Signal
{
public:
	double Fs;
	std::vector<double> dataL;
	std::vector<double> dataR;


	//wavファイルから生成
	Signal(const char* filename);

	//スペクトルからIFFTによって生成
	Signal(class Spectrum);

	//データから生成
	Signal(vector<double>, double F);

	//データの表示
	void show();
	//正規化
	void normalize();

	//信号を二乗する
	void squared();

	void get_after_peak();

	//信号のパワーの計算
	double calc_power();

	//MTFを表示
	void show_MTF(double freq);
	void calc_MTF(double freq, string filename);
	vector<double>get_MTF_vector(int endIdx);
	//ダウンサンプリング 引数分の１にダウンサンプリングする	
	void down_sampling(unsigned int ratio);

	//任意長さになるまでゼロパディング
	void add_zero(long length);

	//wavファイル読み込み関係---------------------------------------------------------
	int readfmtChunk(FILE* fp, tWaveFormatPcm* waveFmtPcm);
	int wavHdrRead(const char* in_wavefile);

	int read8BitWavMonaural(FILE* fpIn);
	int read8BitWavStereo(FILE* fpIn);
	int read16BitWavMonaural(FILE* fpIn);
	int read16BitWavStereo(FILE* fpIn);

	int readDataSub(FILE* fpIn);
	int read(const char* in_wavfile);


	int showWavdata();

	//waveファイル書き込み--------------------------------------------------------------
	long wavHeaderWrite(FILE* fpIn);
	int write16BitWav(FILE* fpOut);
	int write(const char* outFile);

	tWaveFormatPcm waveFormatpcm;
	SWaveFileHeader waveFileheader;

	unsigned long sizeOfData;
	unsigned long posOfData;
	unsigned long samples;//sizeOfData 

	unsigned short bytesPerSingleCh;






};