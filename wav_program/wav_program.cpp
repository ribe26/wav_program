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
	const char* filename = "rir/usina_main_s1_p5.wav";
	Signal c_signal(filename);
	c_signal.show();
	Spectrum spec(c_signal);
	spec.show();
	c_signal.show_MTF();
	/*
	vector<double> imp_vec((c_signal.dataL.size()+1)/2.0, 0.0);
	Signal impulse(imp_vec, c_signal.Fs);

	Spectrum C(c_signal);
	Spectrum C_conj(c_signal);
	C_conj.Conj();

	C.show();
	C_conj.show();
	Spectrum H(impulse);
	*/

	//git _test
}