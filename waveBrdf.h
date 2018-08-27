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

#ifndef WDF_H
#define WDF_H

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "exrimage.h"
#include "gaborkernel.h"

using namespace std;
using namespace Eigen;

class EXRImage;

struct Query {
    Vector2 mu_p;
    Float sigma_p;
    Vector2 omega_i;
    Vector2 omega_o;

    Float lambda;           // in microns.
};

// Heightfield.
// Suppose it starts from (0, 0) and
// extends to positive (w * mTexelWidth, h * mTexelWidth) and
// is tiling along both axes.
class Heightfield {
public:
    Heightfield() {}
    Heightfield(EXRImage *heightfieldImage, Float texelWidth = 1.0, Float vertScale = 1.0)
        : mHeightfieldImage(heightfieldImage), mTexelWidth(texelWidth), mVertScale(vertScale) {}
    GaborKernel g(int i, int j, Float F, Float lambda);

    Vector2 n(Float i, Float j);
public:
    EXRImage *mHeightfieldImage;
    Float mTexelWidth;      // in microns.
    Float mVertScale;
};

class BrdfBase {
public:
    BrdfBase() {}
    BrdfBase(Heightfield *heightfield) : mHeightfield(heightfield) {}
public:
    virtual Float queryBrdf(const Query &query) = 0;
    virtual Float* genBrdfImage(const Query &query, int resolution) = 0;

protected:
    Heightfield *mHeightfield;
};

class GeometricBrdf: public BrdfBase {
public:
    GeometricBrdf() {}
    GeometricBrdf(Heightfield *heightfield, int sampleNum = 10000000) : mHeightfield(heightfield), mSampleNum(sampleNum) {}
    Float* genNdfImage(const Query &query, int resolution);
    virtual Float queryBrdf(const Query &query);
    virtual Float* genBrdfImage(const Query &query, int resolution);
protected:
    Heightfield *mHeightfield;
private:
    int mSampleNum;
};


class WaveBrdf: public BrdfBase {
public:
    WaveBrdf() {}
    WaveBrdf(Heightfield *heightfield) : mHeightfield(heightfield) {}
    virtual Float queryBrdf(const Query &query);
    virtual Float* genBrdfImage(const Query &query, int resolution);
protected:
    Heightfield *mHeightfield;
};

class WaveBrdfAccel: public WaveBrdf {
public:
    WaveBrdfAccel() {}

    WaveBrdfAccel(Heightfield *heightfield, string method);

    comp queryIntegral(const Query &query, int layer, int xIndex, int yIndex);
    virtual Float queryBrdf(const Query &query);
    virtual Float* genBrdfImage(const Query &query, int resolution);

protected:
    Heightfield *mHeightfield;

private:
    vector<vector<GaborKernelPrime>> gaborKernelPrime;
    vector<vector<vector<AABB>>> angularBB;
    int mTopLayer;

public:
    string mMethod;
};

#endif
