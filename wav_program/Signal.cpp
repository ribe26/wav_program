#pragma once
#include "Signal.h"
#include <string>
#include <vector>
#include "Spectrum.h"
#include "matplotlibcpp.h"
#include "wave_defs.h"
#include "FFT.h"
#include <cmath>
#include "utility_functions.h"

// WAVファイルから生成
Signal::Signal(const char* filename) {
    read(filename);
    Fs = this->waveFormatpcm.samplesPerSec;
}

// vectorから生成
Signal::Signal(std::vector<double> data, double F) {
    Fs = F;
    for (int i = 0; i < data.size(); i++) {
        dataL.push_back(data[i]);
        dataR.push_back(data[i]);
    }
}

// IFFTで生成
Signal::Signal(Spectrum spectrum) {
    Fs = spectrum.Fs;
    dataL = std::vector<double>(spectrum.original_length);
    dataR = std::vector<double>(spectrum.original_length);

    std:: vector<double> tempLsig;
    std:: vector<double> tempRsig;


    

    for (long i = 0; i < spectrum.dataL.size(); i++) {
        tempLsig.push_back(spectrum.dataL[i].real());
        tempLsig.push_back(spectrum.dataL[i].imag());
        tempRsig.push_back(spectrum.dataR[i].real());
        tempRsig.push_back(spectrum.dataR[i].imag());
    }

    int n = tempLsig.size()*2;
    int* ip = new int[2 + (int)std::sqrt(0.5 * n) + 1];
    std::vector<double> w(n * 5 / 4);


    ip[0] = 0;
    cdft(tempLsig.size(), 1, &tempLsig[0], ip, &w[0]);
    ip[0] = 0;
    cdft(tempRsig.size(), 1, &tempRsig[0], ip, &w[0]);
    
    cout << tempLsig.size() << endl;
    double ratio = 1.0 / (tempLsig.size() / 2.0);
    cout << ratio <<endl;

    double Max = 0;
    for (int i = 0; i < spectrum.original_length; i++) {
        dataL[i] = tempLsig[i * 2] *ratio;
        dataR[i] = tempRsig[i * 2] *ratio;
        if (abs(dataL[i]) > Max) {
            Max = abs(dataL[i]);
        }
    }
    cout <<"Max:" << Max << endl;
    delete[] ip;
   }

// 信号を表示する
void Signal::show(bool showflag,bool saveflag,string dir,string fname) {
    unsigned long n = dataL.size();
    std::vector<double> plot_x(n);
    std::vector<double> plot_y(n);

    for (size_t i = 0; i < n; i++) {
        plot_x[i] = i/Fs;
        plot_y[i] = dataL[i];
    }


    std::map<std::string, std::string> style1;
    style1["label"] = "Signal";
    style1["color"] = "blue";
    style1["linestyle"] = "-";
    matplotlibcpp::plot(plot_x, plot_y, style1);

    matplotlibcpp::title("RIR", { {"fontsize", "14"} });
    matplotlibcpp::xlabel("time[second]", { {"fontsize", "14"} });
    matplotlibcpp::ylabel("amplitude", { {"fontsize", "14"} });
    matplotlibcpp::legend({ {"fontsize", "14"} });
    if (saveflag) {
        string filename = dir + string("/") + fname + string(".png");
        matplotlibcpp::save(filename);
    }
    if (showflag) {
        matplotlibcpp::show();
    }
    matplotlibcpp::close();
}

// 正規化
void Signal::normalize() {
    double maxL = *std::max_element(this->dataL.begin(), this->dataL.end());
    double minL = *std::min_element(this->dataL.begin(), this->dataL.end());
    double normL = std::max(std::abs(maxL), std::abs(minL));

    double maxR = *std::max_element(this->dataR.begin(), this->dataR.end());
    double minR = *std::min_element(this->dataR.begin(), this->dataR.end());
    double normR = std::max(std::abs(maxR), std::abs(minR));

    for (auto& e : this->dataL) {
        e /= normL;
    }

    for (auto& e : this->dataR) {
        e /= normR;
    }

}

//直接音以降を取り出す
void Signal::get_after_peak() {
    int peakIndex = 0;
    double peakAbs = std::fabs(dataL[0]);

    for (int i = 1; i < static_cast<int>(dataL.size()); ++i) {
        double a = std::fabs(dataL[i]);
        if (a > peakAbs) {
            peakAbs = a;
            peakIndex = i;
        }
    }

    // 2. ピーク以降の部分だけを切り出して新しい vector を作る
    std::vector<double> tail(dataL.begin() + peakIndex, dataL.end());

    dataL = tail;
}

void Signal::set_sec(double second) {
    int size = Fs * second;
    if (dataL.size() > size) {
        dataL.resize(size);
        dataR.resize(size);
    }
    else if (dataL.size()<size) {
        for (int i = 0;dataL.size() < size;i++) {
            dataL.push_back(0.0);
            dataR.push_back(0.0);
        }
    }
}

// 二乗する
void Signal::squared() {
    for (size_t i = 0; i < dataL.size(); i++) {
        dataL[i] = dataL[i] * dataL[i];
        dataR[i] = dataR[i] * dataR[i];
    }
}

// MTFを表示する
void Signal::show_MTF(double freq,bool showflag,bool saveflag,string dir,string fname) {
    std::vector<double> squared_signal;
    for (size_t i = 0; i < dataL.size(); i++) {
        squared_signal.push_back(dataL[i] * dataL[i]);
    }
    Signal sq_sig(squared_signal, Fs);

    Spectrum G(sq_sig);
    double power = this->calc_power();

    std::vector<double> plot_x;
    std::vector<double> plot_y;
    double out = 0;
    double unitFs = Fs / (double)G.dataL.size();
    int endIdx = freq/unitFs;
    

    for (size_t i = 0; i <= endIdx; i++) {
        plot_x.push_back(unitFs * i);
        plot_y.push_back(abs(G.dataL[i]) / power);
    }



    std::map<std::string, std::string> style1;
    style1["label"] = "MTF";
    style1["color"] = "blue";
    style1["linestyle"] = "-";
    matplotlibcpp::plot(plot_x, plot_y, style1);

    matplotlibcpp::title("RIR-MTF", { {"fontsize", "14"} });
    matplotlibcpp::xlabel("Modulation Frequency[Hz]", { {"fontsize", "14"} });
    matplotlibcpp::ylabel("MTF", { {"fontsize", "14"} });
    matplotlibcpp::legend({ {"fontsize", "14"} });
    if (saveflag) {
        string filename = dir + string("/") + fname + string(".png");
        matplotlibcpp::save(filename);
    }
    if (showflag) {
        matplotlibcpp::show();
    }
    matplotlibcpp::close();
}

// MTFを計算する
void Signal::calc_MTF(double freq, string dir, string fname) {
    std::vector<double> squared_signal;
    for (size_t i = 0; i < dataL.size(); i++) {
        squared_signal.push_back(dataL[i] * dataL[i]);
    }
    Signal sq_sig(squared_signal, Fs);
    Spectrum G(sq_sig);
    double power = this->calc_power();

    std::vector<double> plot_y;
    double out = 0;
    double unitFs = Fs / dataL.size();
    int endIdx = freq / unitFs;

    cout << "unitFs:" << unitFs <<endl;
    cout << "fileout end idx:" << endIdx << endl;

    for (size_t i = 0; i < endIdx; i++) {
        plot_y.push_back(abs(G.dataL[i]) / power);
    }
    string filename = dir + string("/") + fname + string(".txt");
    writeVectorToFile(filename, plot_y);
}

vector<double> Signal::get_MTF_vector(int endIdx) {
    std::vector<double> squared_signal;
    for (size_t i = 0; i < dataL.size(); i++) {
        squared_signal.push_back(dataL[i] * dataL[i]);
    }
    Signal sq_sig(squared_signal, Fs);
    Spectrum G(sq_sig);
    double power = this->calc_power();

    std::vector<double> plot_y;
    double out = 0;
    double unitFs = Fs / dataL.size();

    for (size_t i = 0; i < endIdx; i++) {
        plot_y.push_back(abs(G.dataL[i]) / power);
    }
    return plot_y;
}


// パワーを計算する
double Signal::calc_power() {
    double out = 0;
    for (size_t i = 0; i < dataL.size(); i++) {
        out += dataL[i] * dataL[i];
    }
    return out;
}

// ダウンサンプリング
void Signal::down_sampling(double Fs_in) {
    /*
    Fs = Fs / ratio;
    std::vector<double> new_dataL;
    std::vector<double> new_dataR;

    for (size_t i = 0; i < dataL.size(); i++) {
        if (i % ratio == 0) {
            new_dataL.push_back(dataL[i]);
            new_dataR.push_back(dataR[i]);
        }
    }

    dataL = new_dataL;
    dataR = new_dataR;
    */

    std::vector<double> new_dataL = resample(dataL, Fs, Fs_in);
    std::vector<double> new_dataR = resample(dataR, Fs, Fs_in);
    dataL = new_dataL;
    dataR = new_dataR;
    Fs = Fs_in;
}

void Signal::up_sampling(double Fs_in) {
    /*
    Fs = Fs / ratio;
    std::vector<double> new_dataL;
    std::vector<double> new_dataR;

    for (size_t i = 0; i < dataL.size(); i++) {
        if (i % ratio == 0) {
            new_dataL.push_back(dataL[i]);
            new_dataR.push_back(dataR[i]);
        }
    }

    dataL = new_dataL;
    dataR = new_dataR;
    */

    std::vector<double> new_dataL = upsample(dataL, Fs, Fs_in);
    std::vector<double> new_dataR = upsample(dataR, Fs, Fs_in);
    dataL = new_dataL;
    dataR = new_dataR;
    Fs = Fs_in;
}

void Signal::add_zero(long length) {
    long needed = 0;
    needed = length - dataL.size();
    if (needed <= 0) {
        return;
    }
    std::vector<double> zeros(needed, 0.0);
    dataL.insert(dataL.begin(), zeros.begin(), zeros.end());
}

// WAVファイル関係

// fmtチャンクを読み込む
int Signal::readfmtChunk(FILE* fp, tWaveFormatPcm* waveFmtPcm) {
    if (fread(waveFmtPcm, sizeof(tWaveFormatPcm), 1, fp) != 1)
        return -1;

    std::cout << "             formatTag:" << waveFmtPcm->formatTag << "(1 = PCM)" << std::endl;
    std::cout << "              channels:" << waveFmtPcm->channels << std::endl;
    std::cout << "         samplesPerSec:" << waveFmtPcm->samplesPerSec << "[Hz]" << std::endl;
    std::cout << "           bytesPerSec:" << waveFmtPcm->bytesPerSec << "[bytes/sec]" << std::endl;
    std::cout << "            blockAlign:" << waveFmtPcm->blockAlign << "[bytes]" << std::endl;
    std::cout << "         bitsPerSample:" << waveFmtPcm->bitsPerSample << "[bits/sample]" << std::endl;

    return 0;
}

// WAVファイルヘッダを読み込む
int Signal::wavHdrRead(const char* in_wavefile) {
    SWaveFileHeader waveFileHeader;
    tWaveFormatPcm  waveFmtPcm;
    tChank          chank;
    long fPos, len;
    FILE* fp;

    if (fopen_s(&fp, in_wavefile, "rb") != 0) {
        fprintf(stderr, " %sをオープンできません.\n", in_wavefile);
        return -1; // error
    }
    fprintf(stdout, "\n%s :\n", in_wavefile);

    // ヘッダ情報を読み込む
    if (fread(&waveFileHeader, sizeof waveFileHeader, 1, fp) != 1) {
        fprintf(stderr, " %ld で読み込み失敗.\n", ftell(fp));
        fclose(fp);
        return -1; // error
    }

    if (memcmp(&waveFileHeader.hdrRiff, STR_RIFF, 4) != 0) {
        fprintf(stderr, "'RIFF' フォーマットでない.\n");
        fclose(fp);
        return -1; // error
    }

    // WAVE ヘッダ情報を読み込む
    if (memcmp(waveFileHeader.hdrWave, STR_WAVE, 4) != 0) {
        fprintf(stderr, "'WAVE' が無い.\n");
        fclose(fp);
        return -1; // error
    }

    // 4Byte これ以降のバイト数 = (ファイルサイズ - 8)(Byte)
    len = waveFileHeader.sizeOfFile;

    // チャンク情報を読み込む
    while (fread(&chank, sizeof chank, 1, fp) == 1) {
        if (memcmp(chank.hdrFmtData, STR_fmt, sizeof chank.hdrFmtData) == 0) {
            len = chank.sizeOfFmtData;
            std::printf("              fmt size: %ld [bytes]\n", len);
            fPos = ftell(fp);
            if (readfmtChunk(fp, &waveFmtPcm) != 0)
                return -1;
            this->waveFormatpcm = waveFmtPcm;
            this->Fs = waveFmtPcm.samplesPerSec;
            fseek(fp, fPos + len, SEEK_SET);
        }
        else if (memcmp(chank.hdrFmtData, STR_data, 4) == 0) {
            this->sizeOfData = chank.sizeOfFmtData;
            std::cout << "         data size:" << this->sizeOfData << "[bytes]" << std::endl;
            this->posOfData = ftell(fp);
            fseek(fp, this->sizeOfData + this->posOfData - 4, SEEK_SET);
            break;
        }
        else {
            len = chank.sizeOfFmtData;
            std::printf("%c%c%c%cの長さ: %ld [bytes]\n\n",
                chank.hdrFmtData[0], chank.hdrFmtData[1],
                chank.hdrFmtData[2], chank.hdrFmtData[3], len);
            fPos = ftell(fp);
            fseek(fp, fPos + len, SEEK_SET);
        }
    }

    fclose(fp);
    return 0;
}

// 8ビットモノラルWAVを読み込む
int Signal::read8BitWavMonaural(FILE* fpIn) {
    unsigned int i;
    unsigned char In;
    this->samples = this->sizeOfData / sizeof In;

    for (i = 0; i < this->samples; i++) {
        if (fread(&In, sizeof In, 1, fpIn) != 1)
            return -1;

        this->dataL.push_back((double)(In - 128.0));
        this->dataL.push_back(0);
        this->dataR.push_back((double)(In - 128.0));
        this->dataR.push_back(0);
    }

    return 0;
}

// 8ビットステレオWAVを読み込む
int Signal::read8BitWavStereo(FILE* fpIn) {
    unsigned int i;
    unsigned char In[2];

    this->samples = this->sizeOfData / sizeof In;
    for (i = 0; i < this->samples; i++) {
        if (fread(In, sizeof In, 1, fpIn) != 1)
            return -1;

        this->dataL.push_back((double)(In[0] - 128.0));
        this->dataL.push_back(0);
        this->dataR.push_back((double)(In[1] - 128.0));
        this->dataR.push_back(0);
    }

    return 0;
}

// 16ビットモノラルWAVを読み込む
int Signal::read16BitWavMonaural(FILE* fpIn) {
    unsigned int i;
    short In;
    this->samples = this->sizeOfData / sizeof In;
    std::cout << "samples_size:" << samples << std::endl;

    for (i = 0; i < this->samples; i++) {
        if (fread(&In, sizeof(In), 1, fpIn) != 1)
            return -1;

        this->dataL.push_back(In);
        this->dataR.push_back(In);
    }

    return 0;
}

// 16ビットステレオWAVを読み込む
int Signal::read16BitWavStereo(FILE* fpIn) {
    unsigned int i;
    short In[2];
    this->samples = this->sizeOfData / sizeof In;

    for (i = 0; i < this->samples; i++) {
        if (fread(In, sizeof In, 1, fpIn) != 1)
            return -1;

        this->dataL.push_back(In[0]);
        this->dataR.push_back(In[1]);
    }

    return 0;
}

// WAVデータを表示する
int Signal::showWavdata() {
    std::cout << waveFileheader.hdrRiff << std::endl;
    std::cout << waveFileheader.sizeOfFile << std::endl;
    std::cout << waveFileheader.hdrWave << std::endl;
    std::cout << waveFormatpcm.formatTag << std::endl;
    return 0;
}

// WAVデータをダンプする
int Signal::readDataSub(FILE* fpIn) {
    fseek(fpIn, this->posOfData, SEEK_SET); // 元ファイルのデータ開始部分へ

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

// ファイル内容を書き出し
int Signal::read(const char* in_wavefile) {
    // WAVのヘッダを読み込む
    if (wavHdrRead(in_wavefile) != 0)
        return -1;

    FILE* fpIn;
    this->bytesPerSingleCh = this->waveFormatpcm.bitsPerSample / 8;

    if ((fopen_s(&fpIn, in_wavefile, "rb")) != 0) {
        std::printf("%s をオープンできません.\n", in_wavefile);
        return -1;
    }

    // ダンプWAVデータ
    if (readDataSub(fpIn) != 0) {
        std::printf("readDataSubでエラー発生.\n");
        fclose(fpIn);
        return -1;
    }

    fclose(fpIn);
    return 0;
}

// WAVヘッダを書き込む
long Signal::wavHeaderWrite(FILE* fp) {
    WrSWaveFileHeader wrWavHdr;

    wrWavHdr.hdrRiff[0] = 'R';
    wrWavHdr.hdrRiff[1] = 'I';
    wrWavHdr.hdrRiff[2] = 'F';
    wrWavHdr.hdrRiff[3] = 'F';

    wrWavHdr.sizeOfFile = sizeOfData + sizeof(wrWavHdr) - 8;

    wrWavHdr.hdrWave[0] = 'W';
    wrWavHdr.hdrWave[1] = 'A';
    wrWavHdr.hdrWave[2] = 'V';
    wrWavHdr.hdrWave[3] = 'E';

    wrWavHdr.hdrFmt[0] = 'f';
    wrWavHdr.hdrFmt[1] = 'm';
    wrWavHdr.hdrFmt[2] = 't';
    wrWavHdr.hdrFmt[3] = ' ';

    wrWavHdr.sizeOfFmt = 16;

    wrWavHdr.stWaveFormat.formatTag = 1;

    wrWavHdr.stWaveFormat.channels = 2;

    wrWavHdr.stWaveFormat.samplesPerSec = Fs;

    wrWavHdr.stWaveFormat.bytesPerSec = 4 * Fs;

    wrWavHdr.stWaveFormat.blockAlign = 4;

    wrWavHdr.stWaveFormat.bitsPerSample = 16;

    wrWavHdr.hdrData[0] = 'd';
    wrWavHdr.hdrData[1] = 'a';
    wrWavHdr.hdrData[2] = 't';
    wrWavHdr.hdrData[3] = 'a';

    wrWavHdr.sizeOfData = this->dataL.size() * 2;

    fwrite(&wrWavHdr, sizeof wrWavHdr, 1, fp);

    return ftell(fp);
}

// 16ビットWAVを書き込む
int Signal::write16BitWav(FILE* fpOut) {
    std::cout << "called 16Bit stereo" << std::endl;
    unsigned long i;
    short In[2];

    for (i = 0; i < this->dataL.size(); i++) {
        In[0] = (short)this->dataL[i]; // 左チャンネル
        In[1] = (short)this->dataR[i]; // 右チャンネル

        if (fwrite(In, sizeof In, 1, fpOut) != 1)
            return -1;
    }

    return 0;
}

template <typename T>
T clampValue(T v, T lo, T hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

void writeWavMono16(
    const std::string& filename,
    const std::vector<double>& data,
    int sampleRate
) {
    const uint16_t numChannels = 1;   // モノラル
    const uint16_t bitsPerSample = 16;  // 16bit PCM
    const uint32_t numSamples = static_cast<uint32_t>(data.size());

    const uint32_t byteRate = sampleRate * numChannels * bitsPerSample / 8;
    const uint16_t blockAlign = numChannels * bitsPerSample / 8;
    const uint32_t subchunk2Size = numSamples * numChannels * bitsPerSample / 8;
    const uint32_t chunkSize = 36 + subchunk2Size;

    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // ---- RIFFヘッダ ----
    ofs.write("RIFF", 4);
    ofs.write(reinterpret_cast<const char*>(&chunkSize), 4);
    ofs.write("WAVE", 4);

    // ---- fmt チャンク ----
    ofs.write("fmt ", 4);

    uint32_t subchunk1Size = 16;               // PCMなら16
    uint16_t audioFormat = 1;                // PCM=1 (リニアPCM)
    ofs.write(reinterpret_cast<const char*>(&subchunk1Size), 4);
    ofs.write(reinterpret_cast<const char*>(&audioFormat), 2);
    ofs.write(reinterpret_cast<const char*>(&numChannels), 2);
    ofs.write(reinterpret_cast<const char*>(&sampleRate), 4);
    ofs.write(reinterpret_cast<const char*>(&byteRate), 4);
    ofs.write(reinterpret_cast<const char*>(&blockAlign), 2);
    ofs.write(reinterpret_cast<const char*>(&bitsPerSample), 2);

    // ---- data チャンク ----
    ofs.write("data", 4);
    ofs.write(reinterpret_cast<const char*>(&subchunk2Size), 4);

    // ---- 波形データ本体 (16bit PCM) ----
    for (size_t i = 0; i < data.size(); ++i) {
        double d = data[i];
        // -1.0〜1.0 にクリップしてから int16_t に変換
        double clamped = clampValue(d, -1.0, 1.0);
        int16_t sample = static_cast<int16_t>(clamped * 32767.0);
        ofs.write(reinterpret_cast<const char*>(&sample), sizeof(sample));
    }

    ofs.close();
}


// ファイル内容を書き出し
int Signal::write(string filename) {
    writeWavMono16(filename, dataL, Fs);
     
    /*
    FILE* fpOut;
    if ((fopen_s(&fpOut, outFile, "wb")) != 0) {
        fprintf(stderr, "%s をオープンできません.\n", outFile);
        return -1;
    }

    // WAVヘッダ書き込み
    if (wavHeaderWrite(fpOut) != 44) {
        fprintf(stderr, "ヘッダを書き込めません: %s\n", outFile);
        fclose(fpOut);
        return -1;
    }

    // WAVデータ書き込み
    if (write16BitWav(fpOut) != 0) {
        fprintf(stderr, "wavDataWriteでエラー発生.\n");
        fclose(fpOut);
        return -1;
    }

    fclose(fpOut);
    return 0;
    */
    return 0;
    }



