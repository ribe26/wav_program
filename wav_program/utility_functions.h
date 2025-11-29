#pragma once
#include <cmath>
#include <cstdlib>  // rand(), RAND_MAX
#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>

namespace {
    const double PI = 3.14159265358979323846;
    
    void show_two_signal(Signal signal1,Signal signal2) {
        unsigned long n = signal1.dataL.size();
        double Fs = signal1.Fs;
        std::vector<double> plot_x(n);
        std::vector<double> plot_y(n);
        std::vector<double> plot_y2(n);

        for (size_t i = 0; i < n; i++) {
            plot_x[i] = i/Fs;
            plot_y[i] = signal1.dataL[i];
            plot_y2[i] = signal2.dataL[i];
        }

        std::map<std::string, std::string> style1;
        style1["label"] = "original";
        style1["color"] = "blue";
        //style1["marker"] = "+";
        //style1["markeredgecolor"] = "green";
        style1["linestyle"] = "-";

        std::map<std::string, std::string> style2;
        style2["label"] = "filtered";
        style2["color"] = "red";
        //style1["marker"] = "+";
        //style1["markeredgecolor"] = "green";
        style2["linestyle"] = "-";

        matplotlibcpp::plot(plot_x, plot_y,style1);
        matplotlibcpp::plot(plot_x, plot_y2,style2);
        matplotlibcpp::title("RIR", { {"fontsize", "14"} });
        matplotlibcpp::xlabel("time[second]", { {"fontsize", "14"} });
        matplotlibcpp::ylabel("amplitude", { {"fontsize", "14"} });
        matplotlibcpp::legend({ {"fontsize", "14"} });
        matplotlibcpp::show();
    }
    

    void show_two_MTF(Signal signal1, Signal signal2,double endFreq) {
        unsigned long n = signal1.dataL.size();
        double Fs = signal1.Fs;
        double unitFs = Fs / n;
        int endIdx = endFreq / unitFs;
        std::vector<double> plot_x(endIdx);
        std::vector<double> plot_y(n);
        std::vector<double> plot_y2(n);


        plot_y = signal1.get_MTF_vector(endIdx);
        plot_y2 = signal2.get_MTF_vector(endIdx);


        for (size_t i = 0; i < endIdx; i++) {
            plot_x[i] = i *unitFs;
        }

        /*
        cout << "plot_MTF_x_length:" << plot_x.size();
        cout << "plot_MTF_y1_length:" << plot_y.size();
        cout << "plot_MTF_y2_length:" << plot_y2.size();
        */

        std::map<std::string, std::string> style1;
        style1["label"] = "original";
        style1["color"] = "blue";
        //style1["marker"] = "+";
        //style1["markeredgecolor"] = "green";
        style1["linestyle"] = "-";

        std::map<std::string, std::string> style2;
        style2["label"] = "filtered";
        style2["color"] = "red";
        //style1["marker"] = "+";
        //style1["markeredgecolor"] = "green";
        style2["linestyle"] = "-";

        matplotlibcpp::plot(plot_x, plot_y,style1);
        matplotlibcpp::plot(plot_x, plot_y2,style2);
        matplotlibcpp::title("MTF", { {"fontsize", "14"} });
        matplotlibcpp::xlabel("modulation Frequency[Hz]", { {"fontsize", "14"} });
        matplotlibcpp::ylabel("MTF", { {"fontsize", "14"} });
        matplotlibcpp::legend({ {"fontsize", "14"} });
        matplotlibcpp::show();
    }
    
    
    
    std::vector<double> generateSineWave(int length, double amplitude, double frequency, double samplingRate) {
        std::vector<double> sineWave(length);

        // Generate sine wave
        for (int n = 0; n < length; ++n) {
            double t = n / samplingRate; // Time index
            sineWave[n] = amplitude * std::sin(2.0 * PI * frequency * t);
        }

        return sineWave;
    }

    // ボックスミュラー法を使用して正規分布に従う乱数を生成する関数
    double generateGaussianNoise(double mean, double stddev) {
        static bool hasSpare = false;
        static double spare;

        if (hasSpare) {
            hasSpare = false;
            return mean + stddev * spare;
        }

        hasSpare = true;

        double u, v, s;
        do {
            u = 2.0 * rand() / RAND_MAX - 1.0;
            v = 2.0 * rand() / RAND_MAX - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);

        s = sqrt(-2.0 * log(s) / s);
        spare = v * s;
        return mean + stddev * u * s;
    }

    std::vector<double> generateReverbImpulsuse(int length, double Tr, double amplitude, double samplingRate) {
        vector<double> output(length);
        for (int i = 0; i < length; i++) {
            output[i] = amplitude * exp(-6.9 * i / samplingRate / Tr);// * generateGaussianNoise(0, 0.1);
        }
        return output;
    }

    void writeVectorToFile(const std::string& filename, const std::vector<double>& vec) {
        // 出力ファイルストリームを開く
        std::ofstream outFile(filename);

        // ファイルが開けない場合のエラーチェック
        if (!outFile.is_open()) {
            //std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }

        // ベクトルの内容をファイルに書き込む
        for (const double& val : vec) {
            outFile << val << "\n"; // 改行で区切って書き込む
        }

        // ファイルを閉じる
        outFile.close();
        //std::cout << "Successfully wrote vector to " << filename << std::endl;
    }
    // 累積エネルギー計算 (逆累積和)
    std::vector<double> calculateCumulativeEnergy(const std::vector<double>& impulse_response) {
        std::vector<double> energy(impulse_response.size(), 0.0);
        double cumulative_sum = 0.0;

        for (int i = impulse_response.size() - 1; i >= 0; --i) {
            cumulative_sum += impulse_response[i] * impulse_response[i];
            energy[i] = cumulative_sum;
        }

        return energy;
    }

    // 線形フィッティングで傾きを計算
    std::pair<double, bool> linearFit(const std::vector<double>& x, const std::vector<double>& y, double& intercept) {
        if (x.size() < 2 || y.size() < 2) {
            std::cerr << "Error: Not enough points for linear fit." << std::endl;
            return { 0.0, false };  // 傾きゼロ、失敗フラグ
        }

        double n = x.size();
        double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
        double sum_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
        double sum_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

        double denominator = n * sum_xx - sum_x * sum_x;
        if (std::abs(denominator) < 1e-6) {  // 分母が小さい場合のエラー処理
            std::cerr << "Error: Denominator is too small in linear fit." << std::endl;
            return { 0.0, false };
        }

        double slope = (n * sum_xy - sum_x * sum_y) / denominator;
        intercept = (sum_y - slope * sum_x) / n;

        return { slope, true };  // 傾き、成功フラグ
    }

    // 残響時間を計算
#include <vector>
#include <cmath>
#include <limits>

// RIR から RT20 (秒) を求める
// 返り値：RT20 [s]（計算できない場合は NaN）
    double computeRT20(const std::vector<double>& rir, double sampleRate)
    {
        const double eps = 1e-20;

        if (rir.empty() || sampleRate <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        const int N = static_cast<int>(rir.size());

        // 1. Schroeder 積分でエネルギー減衰曲線 (EDC) を計算
        std::vector<double> edc(N);
        double sum = 0.0;
        for (int i = N - 1; i >= 0; --i) {
            sum += rir[i] * rir[i];
            edc[i] = sum;
        }

        // 全体エネルギーが 0 以下なら計算できない
        if (edc[0] <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // 2. 0 dB を基準に正規化して dB に変換
        std::vector<double> edcDb(N);
        const double edc0 = edc[0];
        for (int i = 0; i < N; ++i) {
            double norm = edc[i] / edc0;
            if (norm < eps) norm = eps;
            edcDb[i] = 10.0 * std::log10(norm);  // 0 dB から負の方向へ落ちていく
        }

        // 3. -5 dB と -25 dB を超える点を探す
        int idxStart = -1; // -5 dB になった最初のインデックス
        int idxEnd = -1; // -25 dB になった最初のインデックス

        for (int i = 0; i < N; ++i) {
            if (idxStart < 0 && edcDb[i] <= -5.0) {
                idxStart = i;
            }
            if (idxStart >= 0 && edcDb[i] <= -25.0) {
                idxEnd = i;
                break;
            }
        }

        // 適切な範囲が見つからなかったら NaN
        if (idxStart < 0 || idxEnd <= idxStart) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // 4. [-5 dB, -25 dB] の範囲で t vs dB の線形回帰
        const int M = idxEnd - idxStart + 1;
        double Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;

        for (int k = 0; k < M; ++k) {
            int i = idxStart + k;
            double t = static_cast<double>(i) / sampleRate;  // 時刻 [s]
            double y = edcDb[i];                             // レベル [dB]

            Sx += t;
            Sy += y;
            Sxx += t * t;
            Sxy += t * y;
        }

        double denom = static_cast<double>(M) * Sxx - Sx * Sx;
        if (std::fabs(denom) < 1e-20) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double a = (static_cast<double>(M) * Sxy - Sx * Sy) / denom; // 傾き [dB/s]
        double b = (Sy - a * Sx) / static_cast<double>(M);           // 切片 [dB]

        // a は負の値（減衰）になるはず
        if (a >= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // 5. RT20 (T20)：-5〜-25 dB の 20 dB 減衰から 60 dB 分を外挿
        //   20 dB 減衰時間 = t(-25) - t(-5)
        //   RT20 (T20) = 3 * (t(-25) - t(-5))
        // 線形モデル y(t) = a t + b から t(-5), t(-25) を求めて使う。
        auto t_at_level = [&](double levelDb) -> double {
            return (levelDb - b) / a;     // a < 0 なので正になる
            };

        double t_m5 = t_at_level(-5.0);
        double t_m25 = t_at_level(-25.0);
        double deltaT20 = t_m25 - t_m5;

        if (deltaT20 <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // RT20 = 3 * 20dB減衰時間（=換算 60 dB）
        double rt20 = 3.0 * deltaT20;
        return rt20;
    }


    // ハニング窓を複素ベクトルに掛ける関数
    // spectrum : 複素スペクトル
    // f        : 中心周波数（Hz）
    // fs       : サンプリング周波数（Hz）
    void applyHanningWindow(std::vector<std::complex<double>>& spectrum, double f, double fs, size_t window_size) {
        size_t N = spectrum.size();
        if (window_size > N) {
            window_size = N; // 窓サイズはスペクトル長までに制限
        }
        if (window_size == 0) return; // 窓サイズ0なら何もしない

        // 周波数軸作成
        std::vector<double> freq_axis(N);
        for (size_t i = 0; i < N; ++i) {
            freq_axis[i] = fs * i / N;
        }

        // 中心周波数fに最も近いインデックスを求める
        size_t center_idx = 0;
        double min_diff = std::abs(freq_axis[0] - f);
        for (size_t i = 1; i < N; ++i) {
            double diff = std::abs(freq_axis[i] - f);
            if (diff < min_diff) {
                min_diff = diff;
                center_idx = i;
            }
        }

        // ハニング窓を作成（窓サイズ分）
        std::vector<double> window(window_size);
        for (size_t i = 0; i < window_size; ++i) {
            window[i] = 0.5 * (1 - std::cos(2 * 3.14159265 * i / (window_size - 1)));
        }

        // 窓をスペクトルに掛ける
        // 窓を中心インデックスに合わせて左右に展開
        int half = static_cast<int>(window_size / 2);
        for (size_t i = 0; i < N; ++i) {
            int dist = static_cast<int>(i) - static_cast<int>(center_idx);
            if (dist >= -half && dist <= half) {
                // 窓の該当インデックス
                size_t w_idx = dist + half;
                spectrum[i] *= window[w_idx];
            }
            else {
                // 窓外は0にする（スペクトル成分を消す）
                spectrum[i] = 0;
            }
        }
    }

    void applyBandPassFilter(std::vector<std::complex<double>>& spectrum, double fs, double f_low, double f_high) {
        size_t N = spectrum.size();

        for (size_t i = 0; i < N; ++i) {
            double freq = fs * i / N;
            if (freq < f_low || freq > f_high) {
                spectrum[i] = 0;  // 通過帯域外はゼロにする
            }
        }
    }


    void applySmoothBandPassFilter(std::vector<std::complex<double>>& spectrum, double fs,
        double f_low, double f_high, double transition_bw) {
        size_t N = spectrum.size();

        for (size_t i = 0; i < N; ++i) {
            double freq = fs * i / N;

            if (freq < f_low - transition_bw || freq > f_high + transition_bw) {
                // 通過帯域外＋フェード領域外は0に
                spectrum[i] = 0;
            }
            else if (freq >= f_low - transition_bw && freq < f_low) {
                // 下限フェードイン部分：0 → 1にハニングで滑らかに増加
                double x = (freq - (f_low - transition_bw)) / transition_bw; // 0〜1の正規化
                double w = 0.5 * (1 - std::cos(3.14159265 * x));  // ハニング窓（半分だけ）
                spectrum[i] *= w;
            }
            else if (freq > f_high && freq <= f_high + transition_bw) {
                // 上限フェードアウト部分：1 → 0にハニングで滑らかに減少
                double x = (freq - f_high) / transition_bw; // 0〜1の正規化
                double w = 0.5 * (1 + std::cos(3.14159265 * x));  // ハニング窓（半分だけ）
                spectrum[i] *= w;
            }
            else {
                // 通過帯域はそのまま通す
                // spectrum[i] *= 1.0;
            }
        }
    }


}