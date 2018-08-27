/*

Copyright 2018 Lingqi Yan

This file is part of WaveOpticsBrdf.

WaveOpticsBrdf is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WaveOpticsBrdf is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WaveOpticsBrdf.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef GABOR_KERNEL_H
#define GABOR_KERNEL_H

#define SCALE_FACTOR 2.0
// #define SCALE_FACTOR 2.35482004503

#include "helpers.h"

class GaborKernel {
public:
    GaborKernel() {}
    GaborKernel(Vector2 mu_, Float sigma_, Vector2 a_, comp C_ = comp(Float(1.0), Float(0.0))) :
                mu(mu_), sigma(sigma_), a(a_), C(C_) {}
    comp eval(Vector2 s);
    comp xform(Vector2 u);
public:
    comp C;
    Vector2 mu;
    Float sigma;
    Vector2 a;
};

class GaborKernelPrime {
public:
    GaborKernelPrime() {}
    GaborKernelPrime(Vector2 mu_, Float sigma_, Float cInfo_, Vector2 aInfo_) :
                        mu(mu_), sigma(sigma_), cInfo(cInfo_), aInfo(aInfo_) {}
    GaborKernel toGaborKernel(Float lambda);
    Vector2 getFFTCenter(Float lambda);     // -a = -aInfo / lambda.
    Float getFFTWidth();                    // 1 / (2 pi sigma).
public:
    Float cInfo;
    Vector2 mu;
    Float sigma;
    Vector2 aInfo;
};

#endif
