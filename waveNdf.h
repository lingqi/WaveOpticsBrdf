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
    int resolution;
    ComplexImage img, tmp;
    FloatImage out;
    VectorXc a, b;

    void fft2();
    void fftshift();

public:
    WaveNDF(Heightfield& h, int res): hf(h), resolution(res) {}
    void generate(Query& query, char* outputFilename);
    void generateSpectral(Query& query, char* outputFilename);
};


#endif
