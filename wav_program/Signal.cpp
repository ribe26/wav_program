#include "Signal.h"
#include <string>
#include <vector>
#include "Spectrum.h"




Signal::Signal(const char* filename) {
    read(filename);
    Fs = this->waveFormatpcm.samplesPerSec;
}


/*
Signal::Signal(Spectrum sperctum){
}
*/



//wavファイル関係

int Signal::readfmtChunk(FILE* fp, tWaveFormatPcm* waveFmtPcm)
{
    if (fread(waveFmtPcm, sizeof(tWaveFormatPcm), 1, fp) != 1)
        return -1;

    cout << "             formatTag:" << waveFmtPcm->formatTag << "(1 = PCM)" << endl;
    //printf( "             データ形式: %u (1 = PCM)\n", waveFmtPcm->formatTag);
    cout << "              channels:" << waveFmtPcm->channels << endl;
    //printf( "           チャンネル数: %u\n", waveFmtPcm->channels);
    cout << "         samplesPerSec:" << waveFmtPcm->samplesPerSec << "[Hz]" << endl;
    //printf( "     サンプリング周波数: %lu [Hz]\n", waveFmtPcm->samplesPerSec);
    cout << "           bytesPerSec:" << waveFmtPcm->bytesPerSec << "[bytes/sec]" << endl;
    //printf( "          バイト数 / 秒: %lu [bytes/sec]\n", waveFmtPcm->bytesPerSec);
    cout << "            blockAlign:" << waveFmtPcm->blockAlign << "[bytes]" << endl;
    //printf( " バイト数×チャンネル数: %u [bytes]\n", waveFmtPcm->blockAlign);
    cout << "         bitsPerSample:" << waveFmtPcm->bitsPerSample << "[bits/sample]" << endl;
    //printf( "    ビット数 / サンプル: %u [bits/sample]\n", waveFmtPcm->bitsPerSample);

    return 0;
}


int Signal::wavHdrRead(const char* in_wavefile)
{
    SWaveFileHeader waveFileHeader;
    tWaveFormatPcm  waveFmtPcm;
    tChank          chank;
    long fPos, len;
    FILE* fp;

    if (fopen_s(&fp,in_wavefile, "rb") != 0)
    {
        fprintf(stderr, " %sをオープンできません.\n", in_wavefile);
        return -1;                  // error
    }
    fprintf(stdout, "\n%s :\n", in_wavefile);

    // ヘッダ情報
    if (fread(&waveFileHeader, sizeof waveFileHeader, 1, fp) != 1)
    {
        fprintf(stderr, " %ld で読み込み失敗.\n", ftell(fp));
        fclose(fp);
        return -1;                  // error
    }

    if (memcmp(&waveFileHeader.hdrRiff, STR_RIFF, 4) != 0)
    {
        fprintf(stderr, "'RIFF' フォーマットでない.\n");
        fclose(fp);
        return -1;                  // error
    }

    // WAVE ヘッダ情報
    if (memcmp(waveFileHeader.hdrWave, STR_WAVE, 4) != 0)
    {
        fprintf(stderr, "'WAVE' が無い.\n");
        fclose(fp);
        return -1;                  // error
    }

    // 4Byte これ以降のバイト数 = (ファイルサイズ - 8)(Byte)
    len = waveFileHeader.sizeOfFile;

    // チャンク情報
    while (fread(&chank, sizeof chank, 1, fp) == 1)
    {
        if (memcmp(chank.hdrFmtData, STR_fmt, sizeof chank.hdrFmtData) == 0)
        {
            len = chank.sizeOfFmtData;
            printf("              fmt size: %ld [bytes]\n", len);
            fPos = ftell(fp);
            if (readfmtChunk(fp, &waveFmtPcm) != 0)
                return -1;
            this->waveFormatpcm.formatTag = waveFmtPcm.formatTag;
            this->waveFormatpcm.channels = waveFmtPcm.channels;
            this->waveFormatpcm.samplesPerSec = waveFmtPcm.samplesPerSec;
            this->waveFormatpcm.bytesPerSec = waveFmtPcm.bytesPerSec;
            this->waveFormatpcm.blockAlign = waveFmtPcm.blockAlign;
            this->waveFormatpcm.bitsPerSample = waveFmtPcm.bitsPerSample;
            fseek(fp, fPos + len, SEEK_SET);
        }
        else if (memcmp(chank.hdrFmtData, STR_data, 4) == 0)
        {
            this->sizeOfData = chank.sizeOfFmtData;
            cout << "         data size:" << this->sizeOfData << "[bytes]" << endl;
            this->posOfData = ftell(fp);
            fseek(fp, this->sizeOfData + this->posOfData - 4, SEEK_SET);
            break;
        }
        else
        {
            len = chank.sizeOfFmtData;
            printf("%c%c%c%c¥の長さ: %ld [bytes]\n\n",
                chank.hdrFmtData[0], chank.hdrFmtData[1],
                chank.hdrFmtData[2], chank.hdrFmtData[3], len);
            fPos = ftell(fp);
            fseek(fp, fPos + len, SEEK_SET);
        }

    }
    fclose(fp);

    return 0;
}


int Signal::read8BitWavMonaural(FILE* fpIn)
{
    unsigned int  i;
    unsigned char In;
    this->samples = this->sizeOfData / sizeof In;
    for (i = 0; i < this->samples; i++)
    {
        if (fread(&In, sizeof In, 1, fpIn) != 1)
            return -1;

        this->dataL.push_back(In);
    }
    return 0;
}

int Signal::read8BitWavStereo(FILE* fpIn)
{
    unsigned int  i;
    unsigned char In[2];

    this->samples = this->sizeOfData / sizeof In;
    for (i = 0; i < this->samples; i++)
    {
        if (fread(In, sizeof In, 1, fpIn) != 1)
            return -1;

        this->dataL.push_back(In[0]);
        this->dataR.push_back(In[1]);
    }
    return 0;
}


int Signal::read16BitWavMonaural(FILE* fpIn)
{
    unsigned int  i;
    short In;
    this->samples = this->sizeOfData / sizeof In;
    cout << "samples_size:" << samples << endl;
    for (i = 0; i < this->samples; i++)
    {
        if (fread(&In, sizeof(In), 1, fpIn) != 1)
            return -1;

        this->dataL.push_back(In);
    }

    return 0;
}

int Signal::read16BitWavStereo(FILE* fpIn)
{
    unsigned int  i;
    short In[2];
    this->samples = this->sizeOfData / sizeof In;
    for (i = 0; i < this->samples; i++)
    {
        if (fread(In, sizeof In, 1, fpIn) != 1)
            return -1;

        this->dataL.push_back(In[0]);
        this->dataR.push_back(In[1]);
    }
    return 0;
}


/*--------------------------------------------------------------------------
* wav データをダンプ
*/
int Signal::readDataSub(FILE * fpIn)
{
    fseek(fpIn, this->posOfData, SEEK_SET);    //元ファイルのデータ開始部分へ

    if (this->waveFormatpcm.channels == 1) {
        if (this->bytesPerSingleCh == 1)
            return read8BitWavMonaural(fpIn);
        else
            return read16BitWavMonaural(fpIn);
    }
    else {
        if (this->bytesPerSingleCh == 1)
            return read8BitWavStereo(fpIn);
        else
            return read16BitWavStereo(fpIn);
    }
}

/*--------------------------------------------------------------------------
* ファイル内容書き出し
*
* inFile:     input wav file, must be stereo wav file.
* sampRate:   sampling rate[Hz]
* sampBits:   sampling bits / sec
* posOfData:  strat position of data.
* sizeOfData: size of data.
*
*/
int Signal::read(const char* in_wavefile)
{
    //wavのヘッダを読み込む
    if (wavHdrRead(in_wavefile) != 0)
        return -1;

    FILE* fpIn;
    this->bytesPerSingleCh = this->waveFormatpcm.bitsPerSample / 8;
    if ((fopen_s(&fpIn,in_wavefile, "rb")) != 0)
    {
        printf("%s をオープンできません.\n", in_wavefile);
        return -1;
    }

    // ダンプ wav データ
    if (readDataSub(fpIn) != 0)
    {
        printf("readDataSubでエラー発生.\n");
        fclose(fpIn);
        return -1;
    }

    fclose(fpIn);

    return 0;
}


