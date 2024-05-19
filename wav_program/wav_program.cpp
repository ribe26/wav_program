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
{   
    /*
    int n = 512;
    vector<double> a(n);
    vector<double> a2(n);
    double sinF = 5000.0;

    double omega = 2*3.1415*sinF;
    double F = 44100;
    double T = 1.0/F;
    for (int i = 0; i < n / 2; i++) {
        a[2 * i] = sin(i * omega * T) + sin(i * omega * 2 * T);
        a[2 * i + 1] = 0;
    }


    sinF = 2500.0;
    omega = 2 * 3.1415 * sinF;
    for (int i = 0; i < n / 2; i++) {
        a2[2 * i] = sin(i * omega * T) + sin(i * omega * 2 * T);
        a2[2 * i + 1] = 0;
    }


    std::vector<double>plot_y(n / 2);
    std::vector<double>plot_x(n / 2);

    for (int i = 0; i < n / 2; i++) {
        plot_y[i] = a[2 * i];
        plot_x[i] = i;
    }
    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();

    for (int i = 0; i < n / 2; i++) {
        plot_y[i] = a2[2 * i];
    }
    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();
    
    int* ip = new int[2+(int)sqrt(0.5*n)+1];
    vector<double> w(n * 5 / 4);
    ip[0] = 0;
    cdft(n ,-1,&a[0], ip, &w[0]);

    //ip[0] = 0;
    cdft(n, -1, &a2[0], ip, &w[0]);

    
    std::vector<double>plot_y2(n/4);
    std::vector<double>plot_x2(n/4);
    for (int i = 0; i < n / 4; i++) {
        plot_y2[i]=sqrt(a[2*i]*a[2*i]+a[2*i+1]*a[2*i+1]);
    }


    for (int i = 0; i < n / 4;i++) {
        plot_x2[i] = F/(n/2)*i;
    }


    matplotlibcpp::plot(plot_x2, plot_y2);
    matplotlibcpp::show();

    for (int i = 0; i < n / 4; i++) {
        plot_y2[i] = sqrt(a2[2 * i] * a2[2 * i] + a2[2 * i + 1] * a2[2 * i + 1]);
    }

    matplotlibcpp::plot(plot_x2, plot_y2);
    matplotlibcpp::show();


    cdft(n, 1, &a[0], ip, &w[0]);
    for (int i = 0; i < n; i++) {
        a[i] /= (n/2);
    }

  
    
    for (int i = 0; i < n / 2; i++) {
        plot_y[i] = a[2 * i];
        plot_x[i] = i;
    }
    

    matplotlibcpp::plot(plot_x, plot_y);
    matplotlibcpp::show();

    */

    /*

    int* ip2 = new int[2 + (int)sqrt(0.5 * n) + 1];
    vector<double> w2(n * 5 / 4);
    rdft(n,1,&a[0], ip2, &w2[0]);

    for (int i = 0; i < n; i++) {
        cout << a[i] << endl;
    }
    */
    
    
    Signal signal("Tvtokyo.wav");
    std::cout << "Hello World!\n";
    signal.show();

    Spectrum spec(signal);
}