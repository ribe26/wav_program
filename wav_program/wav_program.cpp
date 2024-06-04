// wav_program.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "Signal.h"
#include "math.h"
#include "FFT.h"
#include "matplotlibcpp.h"
#include <vector>
#include "Spectrum.h"
int main()
{   /*
    int n = 512;
    vector<double> a(n);
    double sinF = 4000.0;

    double omega = 2 * 3.1415 * sinF;
    double F = 44100.0;
    double T = 1.0 / F;
    for (double i = 0.0; i < n; i++) {
        a[i] = sin(i * omega / F);
    }

    Signal signal(a, F);
    signal.show();
    
    Spectrum spec(signal);
    spec.show();

    Signal signal2(spec);
    signal2.show();
    */


    /*
    Signal signal("L-30e000a.wav");
    signal.show();
    Spectrum spec(signal);
    spec.show();

    Signal signal2(spec);
    signal2.show();
    cout << signal2.dataL.size() << endl;
    */

    Signal signal("Tvtokyo.wav");
    signal.show();
    signal.normalize();
    signal.write("test.wav");
    signal.show();
    //git test
    //Spectrum spec(signal);
}