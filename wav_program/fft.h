#pragma once
//Copyright Takuya OOURA, 1996-2001
/*
元ソース：
 * 京都大学助教授の大浦拓哉氏がフリーソフトとして提供する
 * 「汎用 FFT (高速 フーリエ/コサイン/サイン 変換) パッケージ」
 * (http://www.kurims.kyoto-u.ac.jp/~ooura/fft-j.html)
 * のfft4g.cより抜粋,一部改変を行った
*/


/*
n :データ長 二の累乗長にする
a :データ群

ip:int型ポインタ
w :cos/sinテーブル
↑ip[0]=0のに初期化される
*/

void cdft(int n, int isgn, double* a, int* ip, double* w);
void makewt(int nw, int* ip, double* w);

/* -------- child routines -------- */

void bitrv2(int n, int* ip, double* a);
void bitrv2conj(int n, int* ip, double* a);
void cftfsub(int n, double* a, double* w);
void cftbsub(int n, double* a, double* w);
void cft1st(int n, double* a, double* w);
void cftmdl(int n, int l, double* a, double* w);



