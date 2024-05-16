// wav_program.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "Signal.h"
#include "math.h"
#include "fft.h"
#include "matplotlibcpp.h"

int main()
{
    int n = 8192;
    vector<double> a(n);
    
    std::vector<double>y(n/2);
    std::vector<double>x(n / 2);

    double omega = 3.1415;
    double T = 0.001;
    for (int i = 0; i < n / 2; i++) {
        a[2 * i] = sin(i * omega * T) + 2*sin(i * omega * 30 * T);
        a[2 * i + 1] = 0;
        y[i] = a[2 * i];
        x[i] = n*T * i;
    }
    
    
    int* ip = new int[2+(int)sqrt(0.5*n)+1];
    vector<double> w(n * 5 / 4);
    ip[0] = 0;
    cdft(n ,-1,a,ip, w);

    for (int i = 0; i < n / 2; i++) {
        std::cout << "Re:" << a[i * 2] << "Im:" << a[2 * i + 1] << endl;
        y[i] = sqrt(a[i * 2] * a[i * 2] + a[i * 2 + 1] * a[i * 2 + 1]);
    }
    
    matplotlibcpp::plot(x, y);
    matplotlibcpp::show();
    
    
    //Signal signal("Tvtokyo.wav");
    //std::cout << "Hello World!\n";
}

// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー

// 作業を開始するためのヒント: 
//    1. ソリューション エクスプローラー ウィンドウを使用してファイルを追加/管理します 
//   2. チーム エクスプローラー ウィンドウを使用してソース管理に接続します
//   3. 出力ウィンドウを使用して、ビルド出力とその他のメッセージを表示します
//   4. エラー一覧ウィンドウを使用してエラーを表示します
//   5. [プロジェクト] > [新しい項目の追加] と移動して新しいコード ファイルを作成するか、[プロジェクト] > [既存の項目の追加] と移動して既存のコード ファイルをプロジェクトに追加します
//   6. 後ほどこのプロジェクトを再び開く場合、[ファイル] > [開く] > [プロジェクト] と移動して .sln ファイルを選択します
