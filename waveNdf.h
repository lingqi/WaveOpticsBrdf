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

#ifndef WNDF_H
#define WNDF_H

#include "waveBrdf.h"

typedef Eigen::Matrix<std::complex<float>, Eigen::Dynamic, 1> VectorXc;
typedef Eigen::Array<float, 3, 1> Color;
typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FloatImage;
typedef Eigen::Array<Color, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ColorImage;
typedef Eigen::Array<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ComplexImage;


class WaveNDF
{
    Heightfield& hf;
    int resolution, crop;
    ComplexImage img, tmp;
    FloatImage out;
    ColorImage rgb, cropped;
    VectorXc a, b;
    std::vector<float> lambdas, colors;
    std::vector<FloatImage> channels;

    void fft2();
    void fftshift();
    void mkCrop();

public:
    WaveNDF(Heightfield& h, int res, int c): hf(h), resolution(res), crop(c) {}
    void generate(const Query& query, const char* outputFilename);
    void generateSpectral(const Query& query, const char* outputFilename);
};


#endif
