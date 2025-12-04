// wav_program.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "Signal.h"
#include "math.h"
#include "FFT.h"
#include "matplotlibcpp.h"
#include <vector>
#include "Spectrum.h"
#include "MatrixCalculation.h"
#include "utility_functions.h"
#include "thesis_functions.h"
#include "algorithm"
#include <windows.h>

int main()
{

    // init Parameter
    int length = 0;         // Length of the sine wave
    double samplingRate = 0; // Sampling rate in Hz
    double MTF_MAX = 15;
    

    std::string dir = "air_type1_air_phone_stairway_hfrp";   // std::string

    if (CreateDirectoryA(dir.c_str(), NULL)) {   // ANSI版 API を使用
        std::cout << "ディレクトリ作成成功\n";
    }
    else {
        DWORD err = GetLastError();
        if (err == ERROR_ALREADY_EXISTS) {
            std::cout << "既に存在します\n";
        }
        else {
            std::cout << "作成失敗: エラーコード " << err << "\n";
        }
    }

    //std::vector<double> IR = generateReverbImpulsuse(96000, 1.5, 1, 48000);

    //残響を除去したいインパルス応答の読み込み
    //インパルス応答はSignalクラスが持つ変数のdataLというvector<double>型のベクトルで保存されている。
    std::string filename = dir + string(".wav");
    std::string rir_dir = string("rir/") + dir + string(".wav");


    //Signal c_original("usina_main/original.wav");
    //Signal c("usina_main/filtered.wav");
    Signal c_original(rir_dir.c_str());
    Signal c(rir_dir.c_str());
    //Signal c_original(IR, 48000);
    //Signal c(IR, 48000);

    long original_length = c.dataL.size();

    c.down_sampling(16000.0);
    c_original.down_sampling(16000.0);

    c_original.normalize();
    c.normalize();

    c.set_sec(1.5);
    c_original.set_sec(1.5);

    //c_original.get_after_peak();
    //c.get_after_peak();
    //long AP_length = c.dataL.size();


    //c_original.write("IR_alcuin_AP_48000Hz.wav");

    //show_two_signal(c_original, c);
    //show_two_MTF(c_original, c,MTF_MAX);


    //c.add_zero(original_length);
    show_two_signal(false,c_original, c,dir,"");

    //読み込んだインパルス応答の情報を表示
    c.show(true,dir,"original_IR");//インパルス応答をプロット
    //c.show_MTF(600);//インパルス応答のMTFを引数で指定した数値の変調周波数までプロット
    c.calc_MTF(600, "original.txt");//引数で指定した変調周波数までのMTFをテキストファイルに保存

    //元のインパルスと形成したフィルタを畳み込んだインパルスそれぞれの残響時間を表示する。


    //読み込んだインパルス応答をスペクトルに変換
    //スペクトルははSpectrumクラスが持つ変数のdataLというvector<complex<double>>型のベクトルで保存されている。
    Spectrum C(c);
    Spectrum C_ORIGINAL(c_original);

    //スペクトル情報の表示
    C.show(true, dir, "original_IR_Spectrum");//スペクトルをプロット
    C.show_MTF(MTF_MAX+5.0,true,dir,"original_IR_MTF");
    //初期状態のインパルス応答のエネルギーを計算しておく。
    double original_energy = C.calc_energy();

    //読み込んだインパルス応答の長さ、サンプリング周波数を変数に保存しておく。
    length = c.dataL.size();
    samplingRate = c.Fs;
    cout << "Fs" << c.Fs << endl;

    //元のインパルスと形成したフィルタを畳み込んだインパルスそれぞれの残響時間を表示する。
    double rt60original = computeRT20(c_original.dataL, samplingRate);
    double STIoriginal = computeSTI_fromIR_multiband(c_original.dataL, samplingRate);
    cout << "rt60 c original:" << rt60original << endl;
    cout << "STI c original:" << STIoriginal << endl;

    




    //フィルタHをここで定義、これを更新して最後に畳み込むことでインパルス応答の残響時間抑制を狙う。
    complex<double> init_value(1.0, 0.0);//フィルタHを初期化するための複素数(1.0+0.0i)　全て１で初期化(時間領域ではt=0でパルスが立つ信号にあたる。)
    std::vector<complex<double>>h(C.dataL.size(), init_value);
    Spectrum H(h, samplingRate, length);//フィルタHをspectrumクラスとして定義
    H.show(false,"","");//フィルタの初期状態を表示
    std::cout << "H power:" << H.calc_power() << endl;
    std::cout << "H energy:" << H.calc_energy() << endl;


    //最適化のパラメータを定義



    double startMTFfreq = 1.0;
    double endMTFfreq = MTF_MAX;//最適化においてMTFを考慮する上限
    double unitFs = c.Fs / (double)c.dataL.size();
    int startIdx = startMTFfreq / unitFs;
    int endIdx = endMTFfreq / unitFs;//MTFの上限周波数がvectorのインデックスの値でどこにあたるかを計算



    //0~endIdxの整数を連番で持つvectorを定義(forループでこのvectorに格納されている添え字の変調周波数において、MTFのフィルタHについての勾配を計算し、最適化を行っている。)
    //std::vector<int> target_index(endIdx + 1);  // サイズn+1のvectorを作成
    //std::iota(target_index.begin(), target_index.end(), 0);

    std::vector<int> target_index = make_range(startIdx, endIdx);




    //フィルタ計算における暫定の値を保存するための変数やクラス
    complex<double> complex_zero(0.0, 0.0);
    std::vector<complex<double>>H_step(C.dataL.size(), complex_zero);
    std::vector<complex<double>>G_temp(C.dataL.size(), complex_zero);
    Spectrum G(G_temp, samplingRate, length);


    //最適化のループ回数とステップサイズを定義
    int iteration = 20000;
    double step = 1.0;
    //double step= 0.00000000000000000001;

    //フィルタHの最適化を開始
    for (int j = 0; j < iteration; j++) {
        if (j % 100 == 0) { std::cout << "iteration:" << j << endl; }
        //cout << "H power:" << H.calc_energy() << endl;
        std::fill(G_temp.begin(), G_temp.end(), complex_zero);
        for (int k = 0; k < target_index.size(); k++) {
            
            /*
            if (k % 2 == 1) {
                continue;
            }
            */

            /*
            for (int i = 0; i < G.dataL.size(); i++) {
                G.dataL[i] = C.dataL[i] * H.dataL[i];
                //G.dataL[G.dataL.size() - i] = conj(G.dataL[i]);
            }

            G.set_energy(original_energy);
            */

            for (int m = 0; m < C.dataL.size(); m++) {
                int idx = (target_index[k] - m + C.dataL.size()) % C.dataL.size();
                //cout << "(m,idx):(" << m << "," << idx << ")" << endl;
                //G_temp[target_index[k]] += (G.dataL[m] * G.dataL[idx]);
                G_temp[target_index[k]] = C.dataL[m] *C.dataL[idx]*H.dataL[m] * H.dataL[idx];
            }


            G_temp[target_index[k]] /= (double)G.dataL.size();
            //G_temp[target_index[k]] = G.dataL[target_index[k]];
            //double div = (C.dataL.size() * abs(G_temp[target_index[k]]));
            //double div = C.dataL.size()/2.0;
            double div = 1.0;
            for (int p = 0; p < C.dataL.size() / 2; p++) {
                int idx = (target_index[k] - p + C.dataL.size()) % C.dataL.size();
                //フィルタHの各周波数についての目的関数の勾配を計算
                //H_step[p] += step * (G_temp[target_index[k]] * conj(C.dataL[p]) * conj(C.dataL[idx]) * conj(H.dataL[idx])) / div;
                H_step[p] += step * (G_temp[target_index[k]] * conj(C.dataL[p]) * conj(C.dataL[idx]) * conj(H.dataL[idx])) / div;
                //H_step[p] += step * conj(G_temp[k]) * C.dataL[p] * C.dataL[idx] * H.dataL[idx];
                //H_step[p] += step * (conj(H.dataL[p])*abs(C.dataL[p]*C.dataL[p]));

            }
        }
        //フィルタ更新
        for (int p = 0; p < C.dataL.size() / 2; p++) {
            H.dataL[p] += H_step[p];
            H.dataL[H.dataL.size() - p] = conj(H.dataL[p]);
        }
        //フィルタのエネルギーが初期状態と同様になるように調節
        H.normalize_power();
        std::fill(H_step.begin(), H_step.end(), complex_zero);
    }

    //完成したフィルタを表示
    H.normalize_power();
    std::cout << "H power:" << H.calc_power() << endl;
    std::cout << "H energy:" << H.calc_energy() << endl;
    H.show(true,dir,"filter_spectrum");
    Signal H_inv(H);
    H_inv.show(true,dir,"filter_signal");

    //フィルタをインパルス応答に畳み込む。
    for (int p = 0; p < C.dataL.size(); p++) {
        C.dataL[p] = C_ORIGINAL.dataL[p] * H.dataL[p];
        //if (p > C.dataL.size() / 2 - 200 && p < C.dataL.size() / 2) { cout << H.dataL[p] << endl; }
    }
    //畳み込んだ結果が元のインパルス応答のエネルギーと同様になるように調節。
    C.set_energy(C_ORIGINAL.calc_energy());
    //畳み込んだ結果のインパルス応答のスペクトルを表示
    //C.show();
    //C.show_MTF(20);

    //時間領域に逆変換して表示
    Signal c_inv(C);
    //c_inv.show();
    //c_inv.show_MTF(600);
    //MTFを計算してテキストに出力する
    c_inv.calc_MTF(600, "filtered.txt");


    show_two_MTF(true,c_original, c_inv, MTF_MAX+5.0,dir,"MTF_comparison");
    c_original.normalize();
    c_inv.normalize();


    show_two_signal(true,c_original, c_inv,dir,"Signal_comparison");
    
    std::string original_save_dir = dir + string("/original.wav");
    c_original.write(original_save_dir.c_str());
    std::string save_dir = dir + string("/filtered.wav");
    c_inv.write(save_dir.c_str());

    //元のインパルスと形成したフィルタを畳み込んだインパルスそれぞれの残響時間を表示する。
    rt60original = computeRT20(c_original.dataL, samplingRate);
    STIoriginal = computeSTI_fromIR_multiband(c_original.dataL, samplingRate);
    cout << "rt60 c original:" << rt60original << endl;
    cout << "STI c original:" << STIoriginal << endl;

    //元のインパルスと形成したフィルタを畳み込んだインパルスそれぞれの残響時間を表示する。
    double rt60inv = computeRT20(c_inv.dataL, samplingRate);
    double STIinv = computeSTI_fromIR_multiband(c_inv.dataL, samplingRate);
    cout << "rt60 c inv:" << rt60inv << endl;
    cout << "STI c inv:" << STIinv << endl;

    vector<double> STI_result = { STIoriginal,STIinv };
    std::string save_dir_STI = dir + string("/STI_result.txt");
    writeVectorToFile(save_dir_STI, STI_result);

}


