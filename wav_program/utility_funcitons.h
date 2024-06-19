#pragma once


#include <cmath>
#include <cstdlib>  // rand(), RAND_MAX

// �{�b�N�X�~�����[�@���g�p���Đ��K���z�ɏ]�������𐶐�����֐�
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