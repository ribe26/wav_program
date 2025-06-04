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

    // ピーク以降を抜き出す関数
    std::vector<double> extractPeakAndBeyond(const std::vector<double>& signal) {
        // ピークのインデックスを探す
        auto max_it = std::max_element(signal.begin(), signal.end());
        size_t peak_index = std::distance(signal.begin(), max_it);

        // ピーク以降のデータを抜き出す
        std::vector<double> peak_and_beyond(signal.begin() + peak_index, signal.end());

        return peak_and_beyond;
    }
    // RT60を計算する関数
    double calculateRT60(const std::vector<double>& impulse_response, double sampling_rate) {
        size_t N = impulse_response.size();

        // エネルギー（振幅の2乗）を格納するベクター
        std::vector<double> energy(N);

        // エネルギーの積分
        energy[0] = impulse_response[0] * impulse_response[0];
        for (size_t i = 1; i < N; ++i) {
            energy[i] = energy[i - 1] + impulse_response[i] * impulse_response[i];
        }

        // エネルギーの最大値を取得
        double max_energy = energy[N - 1];

        // 減衰部分を探す（エネルギーが最大値の1%まで減少した地点）
        size_t start_idx = N - 1;
        for (size_t i = N - 1; i > 0; --i) {
            if (energy[i] <= max_energy * 0.01) {
                start_idx = i;
                break;
            }
        }

        // 減衰部分（後半部分）のエネルギーを取り出す
        std::vector<double> decay_energy(energy.begin() + start_idx, energy.end());
        size_t decay_length = decay_energy.size();

        // 最小二乗法による指数関数のフィッティング
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;

        // 時間とエネルギーの対数を計算
        for (size_t i = 0; i < decay_length; ++i) {
            double t = i / sampling_rate;  // サンプル番号から秒への変換
            double y = std::log(decay_energy[i]);  // エネルギーの対数
            double x = t;

            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }

        // 最小二乗法による傾き（α）の計算
        double denominator = (decay_length * sum_x2) - (sum_x * sum_x);
        double alpha = ((decay_length * sum_xy) - (sum_x * sum_y)) / denominator;

        // RT60の計算
        double RT60 = 60.0 / alpha;

        return RT60;
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