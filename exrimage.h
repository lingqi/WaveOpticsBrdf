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

#ifndef EXR_IMAGE_H
#define EXR_IMAGE_H

#include <iostream>
#include <OpenEXR/OpenEXRConfig.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>
#include <cmath>
#include <cstdio>
#include "helpers.h"

using namespace std;
using namespace Imf;

class EXRImage {
public:
    EXRImage(const char *filename);

    void readImage(const char *filename);

    // Bicubic interpolation.
    Float getValue(Float x, Float y);

    // Bicubic interpolation. u and v are from 0 to 1.
    Float getValueUV(Float u, Float v);

    static void writeImage(const Float *image, const char *filename, int outputHeight, int outputWidth);

private:
    void computeCoeff(Float *alpha, const Float *x);

    inline int mod(int x, int y) {
        return ((x % y) + y) % y;
    }

    inline Float hp(int x, int y) {
        return image[mod(x, height)][mod(y, width)].r;
    }
     
    inline Float hpx(int x, int y) {
        return (image[mod(x + 1, height)][mod(y, width)].r - image[mod(x - 1, height)][mod(y, width)].r) / 2.0;
    }
     
    inline Float hpy(int x, int y) {
        return (image[mod(x, height)][mod(y + 1, width)].r - image[mod(x, height)][mod(y - 1, width)].r) / 2.0;
    }
     
    inline Float hpxy(int x, int y) {
        return (hp(x + 1, y + 1) - hp(x + 1, y) - hp(x, y + 1) + 2.0 * hp(x, y) - hp(x - 1, y) - hp(x, y - 1) + hp(x - 1, y - 1)) / 2.0;
    }
     
public:
    Array2D<Rgba> image;
    int width, height;
};

#endif
