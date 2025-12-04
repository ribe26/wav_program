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

    //-------------------------------------
// Biquad band-pass フィルタ
//-------------------------------------
    struct Biquad {
        double b0, b1, b2;
        double a1, a2;
        double z1, z2;

        Biquad() : b0(1), b1(0), b2(0), a1(0), a2(0), z1(0), z2(0) {}

        double process(double x) {
            // Direct Form II
            double y = b0 * x + z1;
            z1 = b1 * x - a1 * y + z2;
            z2 = b2 * x - a2 * y;
            return y;
        }

        void reset() {
            z1 = z2 = 0.0;
        }
    };

    // RBJ cookbook 版 Band-pass (constant skirt gain)
    // centerFreq: 中心周波数 [Hz], Q: 品質係数, fs: サンプリング周波数
    Biquad designBandpass(double centerFreq, double Q, double fs)
    {
        Biquad biq;
        double w0 = 2.0 * PI * centerFreq / fs;
        double alpha = std::sin(w0) / (2.0 * Q);

        double b0 = Q * alpha;
        double b1 = 0.0;
        double b2 = -Q * alpha;
        double a0 = 1.0 + alpha;
        double a1 = -2.0 * std::cos(w0);
        double a2 = 1.0 - alpha;

        // 正規化
        biq.b0 = b0 / a0;
        biq.b1 = b1 / a0;
        biq.b2 = b2 / a0;
        biq.a1 = a1 / a0;
        biq.a2 = a2 / a0;

        biq.reset();
        return biq;
    }

    //-------------------------------------
    // 1 バンド分のフィルタ処理
    //-------------------------------------
    std::vector<double> applyBandpass(
        const std::vector<double>& x, double fc, double fs)
    {
        // Q は適当に 4.0 くらい（オクターブ帯域くらいの広さ）
        double Q = 4.0;
        Biquad bp = designBandpass(fc, Q, fs);

        std::vector<double> y(x.size());
        for (size_t n = 0; n < x.size(); ++n) {
            y[n] = bp.process(x[n]);
        }
        return y;
    }

    std::vector<int> make_range(int start, int end) {
        std::vector<int> result;

        if (start <= end) {
            result.reserve(end - start + 1);
            for (int i = start; i <= end; ++i) {
                result.push_back(i);
            }
        }
        else {
            result.reserve(start - end + 1);
            for (int i = start; i >= end; --i) {
                result.push_back(i);
            }
        }

        return result;
    }
    
    void show_two_signal(bool saveflag,Signal signal1, Signal signal2, string save_dir,string fname) {
        unsigned long n = signal1.dataL.size();
        double Fs = signal1.Fs;
        std::vector<double> plot_x(n);
        std::vector<double> plot_y(n);
        std::vector<double> plot_y2(n);

        for (size_t i = 0; i < n; i++) {
            plot_x[i] = i / Fs;
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

        matplotlibcpp::plot(plot_x, plot_y, style1);
        matplotlibcpp::plot(plot_x, plot_y2, style2);
        matplotlibcpp::title("RIR", { {"fontsize", "14"} });
        matplotlibcpp::xlabel("time[second]", { {"fontsize", "14"} });
        matplotlibcpp::ylabel("amplitude", { {"fontsize", "14"} });
        matplotlibcpp::legend({ {"fontsize", "14"} });
        string filename = save_dir +string("/") + fname + string(".png");
        if (saveflag)
        {
            matplotlibcpp::save(filename);
        }
        matplotlibcpp::show();
    }
    

    void show_two_MTF(bool saveflag,Signal signal1, Signal signal2,double endFreq, string save_dir, string fname) {
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
        string filename = save_dir + string("/") + fname + string(".png");
        if (saveflag)
        {
            matplotlibcpp::save(filename);
        }
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
        //vector<double> SineWave = generateSineWave(length, 1.0, 48000, samplingRate);



        for (int i = 0; i < length; i++) {
            output[i] = amplitude * exp(-6.9 * i / samplingRate / Tr)*generateGaussianNoise(0.0,1.0);
        }

        double maxL = *std::max_element(output.begin(), output.end());
        double minL = *std::min_element(output.begin(), output.end());
        double normL = std::max(std::abs(maxL), std::abs(minL));

        for (auto& e : output) {
            e /= normL;
            //e *= 100.0;
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
            if (idxStart >= 0 && edcDb[i] <= -15.0) {
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

    //--- Hann 窓 FIR ローパスフィルタ設計 ---
    std::vector<double> designLowpassFIR(int taps, double cutoff, double fs)
    {
        std::vector<double> h(taps);
        int M = taps - 1;
        double fc = cutoff / fs; // 正規化カットオフ (0〜0.5)

        for (int n = 0; n < taps; ++n) {
            if (n == M / 2) {
                h[n] = 2.0 * fc;
            }
            else {
                double x = PI * (n - M / 2);
                h[n] = sin(2.0 * PI * fc * (n - M / 2)) / x;
            }
            // Hann window
            h[n] *= 0.5 * (1 - cos(2 * PI * n / M));
        }
        return h;
    }

    //--- FIR filtering ---
    std::vector<double> firFilter(const std::vector<double>& x, const std::vector<double>& h)
    {
        int N = x.size();
        int M = h.size();
        std::vector<double> y(N, 0.0);

        for (int n = 0; n < N; ++n) {
            double acc = 0.0;
            for (int k = 0; k < M; ++k) {
                int idx = n - k;
                if (idx >= 0)
                    acc += h[k] * x[idx];
            }
            y[n] = acc;
        }
        return y;
    }

    //--- Downsampling ---
    std::vector<double> downsample(const std::vector<double>& x, int factor)
    {
        std::vector<double> y;
        for (int i = 0; i < x.size(); i += factor) {
            y.push_back(x[i]);
        }
        return y;
    }

    //--- サンプルレート変換 ---
    std::vector<double> resample(const std::vector<double>& input,
        double F1, double F2)
    {
        int factor = (int)(F1 / F2);     // ダウンサンプリング比 (整数前提)
        double cutoff = F2 / 2.0;        // ローパスカットオフ

        // FIR フィルタ設計
        int taps = 101; // フィルタ長（必要に応じて調整）
        auto h = designLowpassFIR(taps, cutoff, F1);

        // フィルタ適用
        auto filtered = firFilter(input, h);

        // 間引き
        return downsample(filtered, factor);
    }


    //STIを計算
    double computeSTI_fromIR_multiband(
        const std::vector<double>& ir,
        double fs)
    {
        if (ir.empty() || fs <= 0.0) {
            return 0.0;
        }

        const size_t N = ir.size();

        // 7 オクターブバンド中心周波数 [Hz]
        const double bandCenters[7] = {
            125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0
        };

        // 14 変調周波数 [Hz] (IEC に準拠した 1/3oct ステップ)
        const double modulationFreqs[14] = {
            0.63, 0.80, 1.00, 1.25, 1.60, 2.00, 2.50,
            3.15, 4.00, 5.00, 6.30, 8.00, 10.0, 12.5
        };

        // male 用のバンド重要度 α_k （正規化前, IEC の値）
        double alpha_raw[7] = {
            0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173
        };

        // α を合計 1 に正規化して重み w_k として使用
        double alpha_sum = 0.0;
        for (int k = 0; k < 7; ++k) alpha_sum += alpha_raw[k];
        double w[7];
        for (int k = 0; k < 7; ++k) w[k] = alpha_raw[k] / alpha_sum;

        // バンドごとの MTI
        double MTI[7] = { 0 };

        // ---- 各バンドごとに処理 ----
        for (int band = 0; band < 7; ++band) {

            // 1) オクターブバンドでフィルタリング
            std::vector<double> h_band = applyBandpass(ir, bandCenters[band], fs);

            // 2) エネルギー系列 e[n] = h^2
            std::vector<double> e(N);
            double energySum = 0.0;
            for (size_t n = 0; n < N; ++n) {
                e[n] = h_band[n] * h_band[n];
                energySum += e[n];
            }
            if (energySum <= 0.0) {
                // このバンドはエネルギーがほぼゼロ → TI=0 とみなす
                MTI[band] = 0.0;
                continue;
            }

            // 3) 14変調周波数での MTF → TI
            const int M = 14;
            double TI_band[14];

            for (int i = 0; i < M; ++i) {
                double Fm = modulationFreqs[i];

                std::complex<double> acc(0.0, 0.0);

                // 離散積分 : sum e[n] * exp(-j 2π Fm t)
                for (size_t n = 0; n < N; ++n) {
                    double t = static_cast<double>(n) / fs;
                    double phase = -2.0 * PI * Fm * t;
                    std::complex<double> ex(std::cos(phase), std::sin(phase));
                    acc += e[n] * ex;
                }

                double mag = std::abs(acc);
                double m = mag / energySum; // 0〜1

                // 数値安定性のためにクリップ
                m = std::max(1e-6, std::min(0.999999, m));

                // MTF → SNR [dB]
                double snr = 10.0 * std::log10(m / (1.0 - m));

                // SNR を [-15, +15] dB にクリップ (IEC)
                snr = std::max(-15.0, std::min(15.0, snr));

                // SNR → TI
                double ti = (snr + 15.0) / 30.0;
                ti = std::max(0.0, std::min(1.0, ti));

                TI_band[i] = ti;
            }

            // 4) このバンドの MTI = TI の平均
            double sumTI = 0.0;
            for (int i = 0; i < M; ++i) sumTI += TI_band[i];
            MTI[band] = sumTI / static_cast<double>(M);
        }

        // 5) 7バンドの MTI を重み w_k で合成 → STI
        double sti = 0.0;
        for (int k = 0; k < 7; ++k) {
            sti += w[k] * MTI[k];
        }

        // 最終的に 0〜1 にクリップ
        sti = std::max(0.0, std::min(1.0, sti));
        return sti;
    }


}