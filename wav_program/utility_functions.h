#pragma once
#include <cmath>
#include <cstdlib>  // rand(), RAND_MAX
#include <iostream>
#include <vector>

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
    for (int i = 0;i < length;i++) {
        output[i] = amplitude * exp(-6.9 * i / samplingRate / Tr) * generateGaussianNoise(0, 0.1);
    }
    return output;
}