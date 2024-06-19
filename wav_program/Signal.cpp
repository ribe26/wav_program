#include "Signal.h"
#include <string>
#include <vector>
#include "Spectrum.h"
#include "matplotlibcpp.h"
#include "wave_defs.h"
#include "FFT.h"
#include "math.h"



//wavファイルから生成
Signal::Signal(const char* filename) {
    read(filename);
    Fs = this->waveFormatpcm.samplesPerSec;
}

//vectorから生成
Signal::Signal(vector<double> data, double F) {
    Fs = F;
    for (int i = 0; i < data.size(); i++) {
        dataL.push_back(data[i]);
        dataR.push_back(data[i]);
        dataL.push_back(0.0);
        dataR.push_back(0.0);
    }
}


//IFFTで生成
Signal::Signal(Spectrum spectrum) {
    Fs = spectrum.Fs;
    dataL = spectrum.dataL;
    dataR = spectrum.dataR;

    int n = spectrum.dataL.size();
    int* ip = new int[2 + (int)sqrt(0.5 * n) + 1];
    vector<double> w(n * 5 / 4);

    ip[0] = 0;
    cdft(n, 1, &dataL[0], ip, &w[0]);
    ip[0] = 0;
    cdft(n, 1, &dataR[0], ip, &w[0]);


    for (int i = 0; i < n; i++) {
        dataL[i] /= (n / 2);
        dataR[i] /= (n / 2);
    }
}



void Signal::show() {
    unsigned long n = dataL.size();

    vector<double> plot_x(n / 2);
    vector<double> plot_y(n / 2);

    for (int i = 0; i < n / 2; i++) {
        plot_x[i] = i;
        plot_y[i] = dataR[2 * i];
    }

    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();

}

void Signal::normalize() {
    double maxL = *max_element(this->dataL.begin(), this->dataL.end());
    double minL = *min_element(this->dataL.begin(), this->dataL.end());
    double normL;
    if (abs(maxL) > abs(minL)) normL = abs(maxL);
    else normL = abs(minL);


    double maxR = *max_element(this->dataR.begin(), this->dataR.end());
    double minR = *min_element(this->dataR.begin(), this->dataR.end());
    double normR;
    if (abs(maxR) > abs(minR)) normL = abs(maxR);
    else normR = abs(minR);

    for (auto& e : this->dataL) {
        e /= normL;
    }

    for (auto& e : this->dataR) {
        e /= normR;
    }
}

void Signal::squared() {
    for (long i = 0; i < dataL.size() / 2; i++) {
        dataL[i * 2] = dataL[i * 2] * dataL[i * 2];
        dataR[i * 2] = dataR[i * 2] * dataR[i * 2];
    }
}

void Signal::show_MTF(double freq) {
    vector<double> squared_signal;
    for (long i = 0; i < dataL.size() / 2; i++) {
        squared_signal.push_back(dataL[i * 2] * dataL[i * 2]);
    }
    Signal sq_sig(squared_signal, Fs);
    //Spectrum SIG(*this);
    //double power = SIG.calc_power(freq);
    
    Spectrum G(sq_sig);
    //double power = G.calc_power(freq);
    double power = this->calc_power();

    vector<double> plot_x;
    vector<double> plot_y;
    double out = 0;
    for (long i = 0;  i < dataL.size() / 2; i++) {
        if (Fs / (G.dataL.size() / 2.0) * i >= freq) {
            break;
        }
        plot_x.push_back( Fs / (G.dataL.size() / 2.0) * i);
        plot_y.push_back(sqrt(G.dataL[2 * i] * G.dataL[2 * i] + G.dataL[2 * i + 1] * G.dataL[2 * i + 1]) / power);
        out +=sqrt(G.dataL[2 * i] * G.dataL[2 * i] + G.dataL[2 * i + 1] * G.dataL[2 * i + 1]) / power;
    }
    cout << out << endl;


    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();

}


double Signal::calc_power() {
    double out = 0;
    for (long i = 0; i < dataL.size() / 2; i++) {
        out += dataL[i * 2] * dataL[i * 2];
    }
    return out;
}

void Signal::down_sampling(unsigned int ratio) {
    Fs = Fs / ratio;
    vector<double> new_dataL;
    vector<double> new_dataR;
    for (long i = 0; i < dataL.size() / 2;i++) {
        if (i % ratio == 0) {
            new_dataL.push_back(dataL[i * 2]);
            new_dataL.push_back(dataL[i * 2 + 1]);

            new_dataR.push_back(dataR[i * 2]);
            new_dataR.push_back(dataR[i * 2 + 1]);
        }
    }
    dataL = new_dataL;
    dataR = new_dataR;
}


//wavファイル関係----------------------------------------------------------------------

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

    if (fopen_s(&fp, in_wavefile, "rb") != 0)
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
            this->Fs = waveFmtPcm.samplesPerSec;
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

        this->dataL.push_back((double)(In - 128.0));
        this->dataL.push_back(0);
        this->dataR.push_back((double)(In - 128.0));
        this->dataR.push_back(0);
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

        this->dataL.push_back((double)(In[0] - 128.0));
        this->dataL.push_back(0);
        this->dataR.push_back((double)(In[1] - 128.0));
        this->dataR.push_back(0);
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
        this->dataL.push_back(0);
        this->dataR.push_back(In);
        this->dataR.push_back(0);
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
        this->dataL.push_back(0);
        this->dataR.push_back(In[1]);
        this->dataR.push_back(0);
    }
    return 0;
}


int Signal::showWavdata() {
    cout << waveFileheader.hdrRiff << endl;
    cout << waveFileheader.sizeOfFile << endl;
    cout << waveFileheader.hdrWave << endl;

    cout << waveFormatpcm.formatTag << endl;


}

/*--------------------------------------------------------------------------
* wav データをダンプ
*/
int Signal::readDataSub(FILE* fpIn)
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
    if ((fopen_s(&fpIn, in_wavefile, "rb")) != 0)
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

/*waveファイル保存用関数*/



long Signal::wavHeaderWrite(FILE* fp)
{
    unsigned short bytes;
    WrSWaveFileHeader wrWavHdr;
    //strncpy_s(wrWavHdr.hdrRiff, STR_RIFF, sizeof wrWavHdr.hdrRiff); // RIFF ヘッダ
    wrWavHdr.hdrRiff[0] = 'R';
    wrWavHdr.hdrRiff[1] = 'I';
    wrWavHdr.hdrRiff[2] = 'F';
    wrWavHdr.hdrRiff[3] = 'F';

    wrWavHdr.sizeOfFile = sizeOfData + sizeof(wrWavHdr) - 8;      // ファイルサイズ

    //strncpy_s(wrWavHdr.hdrWave, STR_WAVE, sizeof wrWavHdr.hdrWave); // WAVE ヘッダ
    wrWavHdr.hdrWave[0] = 'W';
    wrWavHdr.hdrWave[1] = 'A';
    wrWavHdr.hdrWave[2] = 'V';
    wrWavHdr.hdrWave[3] = 'E';

    //strncpy_s(wrWavHdr.hdrFmt, STR_fmt, sizeof wrWavHdr.hdrFmt);    // fmt チャンク
    wrWavHdr.hdrFmt[0] = 'f';
    wrWavHdr.hdrFmt[1] = 'm';
    wrWavHdr.hdrFmt[2] = 't';
    wrWavHdr.hdrFmt[3] = ' ';

    wrWavHdr.sizeOfFmt = 16;           // fmt チャンク,無圧縮wav は 16

    wrWavHdr.stWaveFormat.formatTag = 1;                          // 無圧縮PCM = 1

    wrWavHdr.stWaveFormat.channels = 2;                          // ch (mono=1, stereo=2)

    wrWavHdr.stWaveFormat.samplesPerSec = Fs;               // sampleng rate(Hz)

    wrWavHdr.stWaveFormat.bytesPerSec = 4 * Fs;

    wrWavHdr.stWaveFormat.blockAlign = 4;                // byte/サンプル*チャンネル

    wrWavHdr.stWaveFormat.bitsPerSample = 16;               // bit/サンプル

    //strncpy_s(wrWavHdr.hdrData, STR_data, sizeof wrWavHdr.hdrData); // dataチャンク
    wrWavHdr.hdrData[0] = 'd';
    wrWavHdr.hdrData[1] = 'a';
    wrWavHdr.hdrData[2] = 't';
    wrWavHdr.hdrData[3] = 'a';

    wrWavHdr.sizeOfData = this->dataL.size() * 2;              // データ長 (byte)

    fwrite(&wrWavHdr, sizeof wrWavHdr, 1, fp);                  // write header

    return ftell(fp);
}


/*--------------------------------------------------------------------------
* 16 bits/sampling
*/
int Signal::write16BitWav(FILE* fpOut)
{
    cout << "called 16Bit stereo" << endl;
    unsigned long  i;
    short In[2];

    for (i = 0; i < this->dataL.size() / 2; i++)
    {
        In[0] = (short)this->dataL[i * 2];    //Left
        In[1] = (short)this->dataR[i * 2];    //Right
        //cout << In[0] << endl;
        if (fwrite(In, sizeof In, 1, fpOut) != 1)
            return -1;
    }
    return 0;
}

/*--------------------------------------------------------------------------
* ファイル内容書き出し
*
* inFile:     input wav file, must be stereo wav file.
* outFile:    output wav file, it only has left  channel.
* sampRate:   sampling rate[Hz]
* sampBits:   sampling bits / sec
* posOfData:  strat position of data.
* sizeOfData: size of data.
*/
int Signal::write(const char* outFile)
{
    FILE* fpOut;

    if ((fopen_s(&fpOut, outFile, "wb")) != 0)
    {
        fprintf(stderr, "%s をオープンできません.¥n", outFile);
        return -1;
    }

    // wav ヘッダ書き込み
    if (wavHeaderWrite(fpOut) != 44)
    {
        fprintf(stderr, "ヘッダを書き込めません: %s¥n", outFile);
        fclose(fpOut);
        return -1;
    }

    // wav データ書き込み
    if (write16BitWav(fpOut) != 0)
    {
        fprintf(stderr, "wavDataWriteでエラー発生.¥n");
        fclose(fpOut);
        return -1;
    }

    fclose(fpOut);

    return 0;
}