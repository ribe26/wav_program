#pragma once
#include <stdio.h>
#include <string.h>
#include <iostream>
using namespace std;

#ifndef _MAX_PATH
#define _MAX_PATH (255)
#endif

#pragma pack(push,1)
typedef struct tagSWaveFileHeader
{
    char            hdrRiff[4];         // 'RIFF'
    unsigned int    sizeOfFile;         // ファイルサイズ - 8
    char            hdrWave[4];         // 'WAVE'
} SWaveFileHeader;

typedef struct tagChank
{
    unsigned char   hdrFmtData[4];      // 'fmt ' or 'data'
    unsigned int    sizeOfFmtData;      // sizeof(PCMWAVEFORMAT) or Waveデーターサイズ
} tChank;

typedef struct tagWaveFormatPcm
{
    unsigned short  formatTag;          // WAVE_FORMAT_PCM
    unsigned short  channels;           // number of channels
    unsigned int    samplesPerSec;      // sampling rate
    unsigned int    bytesPerSec;        // samplesPerSec * channels * (bitsPerSample/8)
    unsigned short  blockAlign;         // block align
    unsigned short  bitsPerSample;      // bits per sampling
} tWaveFormatPcm;



#pragma pack(pop)


// defines
#define STR_RIFF    "RIFF"              //
#define STR_WAVE    "WAVE"              //
#define STR_fmt     "fmt "              //
#define STR_data    "data"              //

#define WAV_MONAURAL    1
#define WAV_STEREO      2
