// wav_program.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "Signal.h"
#include "math.h"
#include "FFT.h"
#include "matplotlibcpp.h"
#include <vector>
#include "Spectrum.h"
#include "utility_funcitons.h"	


//様々な残響時間でMTFの表示をする

int main()
{
	int samples = 96000;
	vector<double> sig(samples, 0);
	double TR = 0.1;
	double amp = 1.0;
	double freq = 48000;
	double F = 48000;
	for (int i = 0;i < samples;i++) {
		sig[i] = amp * exp(-6.9 * (i / freq) / TR) * generateGaussianNoise(0.0, 0.01);
	}

	//残響時間のリスト
	vector<double> TRs = { 0.01,0.05,0.1,0.3,0.5,1,1.5 };


	for (int i = 0;i < TRs.size();i++) {
		cout << "TR=" << TRs[i] << endl;
		TR = TRs[i];
		for (int i = 0;i < samples;i++) {
			sig[i] = amp * exp(-6.9 * (i / freq) / TR) * generateGaussianNoise(0.0, 0.01);
		}
		Signal c_signal(sig, freq);
		c_signal.show();
		c_signal.show_MTF(2000);
	}
}