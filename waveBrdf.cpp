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

#include "waveBrdf.h"
#include "spectrum.h"

GaborKernel Heightfield::g(int i, int j, Float F, Float lambda) {
    Vector2 m_k(Float((i + 0.5) * mTexelWidth), Float((j + 0.5) * mTexelWidth));
    Float l_k = mTexelWidth;

    Vector2 mu_k = m_k;
    Float sigma_k = l_k / SCALE_FACTOR;

    Float H_mk = mHeightfieldImage->getValue(i + 0.5, j + 0.5) * mTexelWidth * mVertScale;   // Assuming mTexelWidth doesn't affect the heightfield's shape.
    Vector2 HPrime_mk(Float((mHeightfieldImage->getValue(i + 1, j) - mHeightfieldImage->getValue(i, j)) * mVertScale),
                      Float((mHeightfieldImage->getValue(i, j + 1) - mHeightfieldImage->getValue(i, j)) * mVertScale));

    comp C_k = sqrt(F) * l_k * l_k * cnis(4.0 * M_PI / lambda * (H_mk - HPrime_mk.dot(m_k)));
    Vector2 a_k = 2.0 * HPrime_mk / lambda;

    return GaborKernel(mu_k, sigma_k, a_k, C_k);
}

Vector2 Heightfield::n(Float i, Float j) {
    Vector2 HPrime(Float((mHeightfieldImage->getValue(i + 0.5f, j) - mHeightfieldImage->getValue(i - 0.5f, j)) * mVertScale),
                   Float((mHeightfieldImage->getValue(i, j + 0.5f) - mHeightfieldImage->getValue(i, j - 0.5f)) * mVertScale));

    Vector3 n(-HPrime(0), -HPrime(1), Float(1.0));
    n.normalize();

    return n.head(2);
}

WaveBrdfAccel::WaveBrdfAccel(Heightfield *heightfield, string method) {
    mHeightfield = heightfield;
    mMethod = method;

    const int &height = mHeightfield->mHeightfieldImage->height;
    const int &width = mHeightfield->mHeightfieldImage->width;
    assert(height == width);


    mTopLayer = (int) floor(log(height * 1.0) / log(2.0) + 1e-4); // 0, 1, 2, ..., mTopLayer.

    cout << "Preprocessing heightfield: layer 0..." << endl;
    angularBB.push_back(vector<vector<AABB>>());
    vector<vector<AABB>> &angularBBLayer = angularBB.back();
    angularBBLayer.reserve(height);
    for (int i = 0; i < height; i++) {
        angularBBLayer.push_back(vector<AABB>());
        angularBBLayer.back().reserve(width);
        gaborKernelPrime.push_back(vector<GaborKernelPrime>());
        gaborKernelPrime.back().reserve(width);
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < height; i++) {
        vector<GaborKernelPrime> &gaborKernelPrimeRow = gaborKernelPrime[i];
        vector<AABB> &angularBBRow = angularBBLayer[i];
        for (int j = 0; j < width; j++) {
            Vector2 m_k(Float((i + 0.5) * mHeightfield->mTexelWidth), Float((j + 0.5) * mHeightfield->mTexelWidth));
            Float l_k = mHeightfield->mTexelWidth;

            Vector2 mu_k = m_k;
            Float sigma_k = l_k / SCALE_FACTOR;

            Float H_mk = mHeightfield->mHeightfieldImage->getValue(i, j) * mHeightfield->mTexelWidth * mHeightfield->mVertScale;   // Assuming mTexelWidth doesn't affect the heightfield's shape.

            Vector2 HPrime_mk(Float((mHeightfield->mHeightfieldImage->getValue(i + 1, j) - mHeightfield->mHeightfieldImage->getValue(i - 1, j)) / 2.0 * mHeightfield->mVertScale),
                              Float((mHeightfield->mHeightfieldImage->getValue(i, j + 1) - mHeightfield->mHeightfieldImage->getValue(i, j - 1)) / 2.0 * mHeightfield->mVertScale));

            Float cInfo_k = H_mk - HPrime_mk.dot(m_k);
            Vector2 aInfo_k = 2.0 * HPrime_mk;

            gaborKernelPrimeRow.push_back(GaborKernelPrime(mu_k, sigma_k, cInfo_k, aInfo_k));
            angularBBRow.push_back(AABB(-aInfo_k(0), -aInfo_k(0), -aInfo_k(1), -aInfo_k(1)));
        }
    }

    int currentHeight = height;
    int currentWidth = width;
    for (int currentLayer = 1; currentLayer <= mTopLayer; currentLayer++) {
        cout << "Preprocessing heightfield: layer " << currentLayer << "..." << endl;
        vector<vector<AABB>> angularBBLayer;
        currentHeight >>= 1;
        currentWidth >>= 1;
        angularBBLayer.reserve(currentHeight);
        for (int i = 0; i < currentHeight; i++) {
            vector<AABB> angularBBRow;
            angularBBRow.reserve(currentWidth);
            for (int j = 0; j < currentWidth; j++) {
                AABB aabb1 = angularBB[currentLayer - 1][i * 2 + 0][j * 2 + 0];
                AABB aabb2 = angularBB[currentLayer - 1][i * 2 + 1][j * 2 + 0];
                AABB aabb3 = angularBB[currentLayer - 1][i * 2 + 1][j * 2 + 1];
                AABB aabb4 = angularBB[currentLayer - 1][i * 2 + 0][j * 2 + 1];
                angularBBRow.push_back(combineAABB(aabb1, aabb2, aabb3, aabb4));
            }
            angularBBLayer.push_back(angularBBRow);
        }
        angularBB.push_back(angularBBLayer);
    }
}

comp WaveBrdfAccel::queryIntegral(const Query &query, int layer, int xIndex, int yIndex) {
    Float C3 = 2.0f;
    if (mMethod == "GHS" || mMethod == "RGHS" || mMethod == "Kirchhoff") {
        Float oDotN = sqrt(abs(1.0f - query.omega_o.dot(query.omega_o)));
        Float iDotN = sqrt(abs(1.0f - query.omega_i.dot(query.omega_i)));
        C3 = oDotN + iDotN;
    }

    if (layer == 0) {
        Vector2 omega_a = (query.omega_i + query.omega_o) / 2.0;
        Vector2 m_k(Float((xIndex + 0.5) * mHeightfield->mTexelWidth), Float((yIndex + 0.5) * mHeightfield->mTexelWidth));

        Float period = mHeightfield->mTexelWidth * mHeightfield->mHeightfieldImage->height;
        Float pDistSqr = distSqrPeriod(m_k, query.mu_p, period);
        if (pDistSqr > (3.0 * query.sigma_p) * (3.0 * query.sigma_p))
            return comp(Float(0.0), Float(0.0));

        Vector2 uQuery = 2.0 * omega_a / query.lambda;
        AABB &aabb = angularBB[layer][xIndex][yIndex];
        Vector2 abLambda(aabb.xMin / query.lambda, aabb.yMin / query.lambda);
        
        abLambda *= C3 / 2.0f;

        Float sigma_k = mHeightfield->mTexelWidth / SCALE_FACTOR;
        Float aSigma = 1.0 / (2.0 * M_PI * sigma_k);
        Float aDistSqr = (abLambda - uQuery).squaredNorm();
        if (aDistSqr > (3.0 * aSigma) * (3.0 * aSigma))
            return comp(Float(0.0), Float(0.0));

        Float C2 = 1.0f;
        if (mMethod == "Kirchhoff") {
            Vector2 HPrime = gaborKernelPrime[xIndex][yIndex].aInfo / 2.0f;
            Vector3 m(-HPrime(0), -HPrime(1), 1.0);
            m.normalize();

            Float oDotN = sqrt(abs(1.0f - query.omega_o.dot(query.omega_o)));
            Float iDotN = sqrt(abs(1.0f - query.omega_i.dot(query.omega_i)));
            Vector3 omega_o(query.omega_o(0), query.omega_o(1), oDotN);
            Vector3 omega_i(query.omega_i(0), query.omega_i(1), iDotN);

            C2 = (omega_o + omega_i).dot(m) / m(2);
        }

        GaborKernelPrime gp = gaborKernelPrime[xIndex][yIndex];
        gp.aInfo *= C3 / 2.0f;
        gp.cInfo *= C3 / 2.0f;

        Float w_mk = Float(1.0 / (sqrt(M_PI) * query.sigma_p)) * exp(-pDistSqr / (2.0 * query.sigma_p * query.sigma_p));
        GaborKernel g = gp.toGaborKernel(query.lambda);
        return w_mk * C2 * g.xform(2.0 * omega_a / query.lambda);
    }

    // Reject the node positionally.
    int layerScale = floor(pow(2.0, layer) + 1e-4);
    AABB positionalBB((xIndex + 0.0) * layerScale * mHeightfield->mTexelWidth,
                      (xIndex + 1.0) * layerScale * mHeightfield->mTexelWidth,
                      (yIndex + 0.0) * layerScale * mHeightfield->mTexelWidth,
                      (yIndex + 1.0) * layerScale * mHeightfield->mTexelWidth);
    AABB queryBB(query.mu_p(0) - 3.0 * query.sigma_p, query.mu_p(0) + 3.0 * query.sigma_p,
                 query.mu_p(1) - 3.0 * query.sigma_p, query.mu_p(1) + 3.0 * query.sigma_p);
    Float period = mHeightfield->mTexelWidth * mHeightfield->mHeightfieldImage->height;
    if (!intersectAABBRot(positionalBB, queryBB, period))
        return comp(Float(0.0), Float(0.0));

    // Reject the node angularly.
    Vector2 omega_a = (query.omega_i + query.omega_o) / 2.0;
    Vector2 uQuery = 2.0 * omega_a / query.lambda;
    AABB &aabb = angularBB[layer][xIndex][yIndex];
    AABB aabbLambda(aabb.xMin * C3 / 2.0f / query.lambda, aabb.xMax * C3 / 2.0f / query.lambda,
                    aabb.yMin * C3 / 2.0f / query.lambda, aabb.yMax * C3 / 2.0f / query.lambda);

    Float sigma_k = mHeightfield->mTexelWidth / SCALE_FACTOR;
    Float sigma = 1.0 / (2.0 * M_PI * sigma_k);

    AABB aabbLambdaExpanded(aabbLambda.xMin - sigma * 3.0,
                            aabbLambda.xMax + sigma * 3.0,
                            aabbLambda.yMin - sigma * 3.0,
                            aabbLambda.yMax + sigma * 3.0);

    if (!insideAABB(aabbLambdaExpanded, uQuery))
        return comp(Float(0.0), Float(0.0));

    return queryIntegral(query, layer - 1, xIndex * 2 + 0, yIndex * 2 + 0) +
           queryIntegral(query, layer - 1, xIndex * 2 + 1, yIndex * 2 + 0) +
           queryIntegral(query, layer - 1, xIndex * 2 + 1, yIndex * 2 + 1) +
           queryIntegral(query, layer - 1, xIndex * 2 + 0, yIndex * 2 + 1);
}

Float WaveBrdfAccel::queryBrdf(const Query &query) {
    // Move the query center to the first HF period.
    Query q = query;
    Float hfHeightWorld = mHeightfield->mTexelWidth * mHeightfield->mHeightfieldImage->height;
    Float hfWidthWorld = mHeightfield->mTexelWidth * mHeightfield->mHeightfieldImage->width;
    q.mu_p(0) -= ((int) floor(q.mu_p(0) / hfHeightWorld)) * hfHeightWorld;
    q.mu_p(1) -= ((int) floor(q.mu_p(1) / hfWidthWorld)) * hfWidthWorld;

    // Traverse the hierarchy to calculate the inner integral.
    comp I = queryIntegral(q, mTopLayer, 0, 0);

    // Update corresponding C1, C2 and C3 based on different diffraction models.
    Float oDotN = sqrt(abs(1.0f - query.omega_o.dot(query.omega_o)));
    Float iDotN = sqrt(abs(1.0f - query.omega_i.dot(query.omega_i)));
    Float C1 = 0.0f;

    if (mMethod == "OHS" || mMethod == "GHS")
        C1 = oDotN / (query.lambda * query.lambda * iDotN);
    else if (mMethod == "ROHS" || mMethod == "RGHS")
        C1 = (iDotN + oDotN) * (iDotN + oDotN) / (query.lambda * query.lambda * 4.0f * iDotN * oDotN);
    else if (mMethod == "Kirchhoff")
        C1 = 1.0f / (query.lambda * query.lambda * 4.0f * iDotN * oDotN);

    return C1 * pow(abs(I), 2.0f);
}

Float* WaveBrdfAccel::genBrdfImage(const Query &query, int resolution) {
    Float *brdfImage = new Float[resolution * resolution * 3];

    if (query.lambda != 0.0) {
         // Single wavelength.
         #pragma omp parallel for schedule(dynamic)
         for (int i = 0; i < resolution; i++) {
             printf("Generating BRDF image: row %d...\n", i);
             for (int j = 0; j < resolution; j++) {
                 Vector2 omega_o((i + 0.5) / resolution * 2.0 - 1.0,
                                 (j + 0.5) / resolution * 2.0 - 1.0);
                 Float brdfValue;
                 if (omega_o.norm() > 1.0) {
                     brdfValue = 0.0;
                 } else {
                     Query q = query;
                     q.omega_o = omega_o;
                     brdfValue = queryBrdf(q);
                 }

                 if (std::isnan(brdfValue))
                     brdfValue = 0.0;

                 brdfImage[(i * resolution + j) * 3 + 0] = brdfValue;
                 brdfImage[(i * resolution + j) * 3 + 1] = brdfValue;
                 brdfImage[(i * resolution + j) * 3 + 2] = brdfValue;
             }
         }
    } else {
        // Multiple wavelengths.
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < resolution; i++) {
            printf("Generating BRDF image: row %d...\n", i);
            for (int j = 0; j < resolution; j++) {
                Vector2 omega_o((i + 0.5) / resolution * 2.0 - 1.0,
                                (j + 0.5) / resolution * 2.0 - 1.0);

                vector<float> spectrumSamples;

                for (int k = 0; k < SPECTRUM_SAMPLES; k++) {

                    Float brdfValue;
                    if (omega_o.norm() > 1.0) {
                        brdfValue = 0.0;
                    } else {
                        Query q = query;
                        q.omega_o = omega_o;
                        q.lambda = (k + 0.5) / SPECTRUM_SAMPLES * (0.83 - 0.36) + 0.36;
                        brdfValue = queryBrdf(q);
                    }

                    if (std::isnan(brdfValue))
                        brdfValue = 0.0;

                    spectrumSamples.push_back(brdfValue);
                }

                float r, g, b;
                SpectrumToRGB(spectrumSamples, r, g, b);

                brdfImage[(i * resolution + j) * 3 + 0] = r;
                brdfImage[(i * resolution + j) * 3 + 1] = g;
                brdfImage[(i * resolution + j) * 3 + 2] = b;
            }
        }
    }

    return brdfImage;
}

Float WaveBrdf::queryBrdf(const Query &query) {
    cout << "Not implemented" << endl;
    return 0.0;
}

Float* WaveBrdf::genBrdfImage(const Query &query, int resolution) {
    cout << "Not implemented" << endl;
    return NULL;
}

Float GeometricBrdf::queryBrdf(const Query &query) {
    cout << "Not implemented" << endl;
    return 0.0;
}

inline Vector2 sampleGauss2d(Float r1, Float r2) {
    // https://en.wikipedia.org/wiki/Box-Muller_transform
    Float tmp = std::sqrt(-2 * std::log(r1));
    Float x = tmp * std::cos(2 * Float(M_PI) * r2);
    Float y = tmp * std::sin(2 * Float(M_PI) * r2);
    return Vector2(x, y);
}

Float* GeometricBrdf::genNdfImage(const Query &query, int resolution) {
    int N = (int) std::sqrt(mSampleNum);
    const Float intrinsicRoughness = Float(1) / N;
    int *inds = new int[N * N];

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // sample query Gaussian, stratified
            Float rx = (i + randUniform<Float>()) / N;
            Float ry = (j + randUniform<Float>()) / N;
            Vector2 g = sampleGauss2d(rx, ry);

            // look up normal
            Vector2 x = g * query.sigma_p / mHeightfield->mTexelWidth;
            x += query.mu_p;
            Vector2 normal = mHeightfield->n(x[0], x[1]);

            // intrinsic roughness
            normal += intrinsicRoughness * sampleGauss2d(randUniform<Float>(), randUniform<Float>());

            int xi = (int)((1 + normal[0]) / 2 * resolution);
            int yi = (int)((1 - normal[1]) / 2 * resolution);
            if (xi < 0 || xi >= resolution || yi < 0 || yi >= resolution) continue;
            inds[i*N + j] = yi * resolution + xi;
        }
    }

    // bin the samples using indices computed above
    int npix = resolution * resolution;
    int *bins = new int[npix];
    memset(bins, 0, npix * sizeof(int));
    for (int i = 0; i < N*N; i++) bins[inds[i]]++;

    Vector3 *ndfImage = new Vector3[npix];
    Float scale = Float(npix) / (4 * N * N);
    for (int i = 0; i < npix; i++) ndfImage[i] = Vector3::Constant(scale * bins[i]);

    delete[] inds;
    delete[] bins;
    return (Float*) ndfImage;
}

Float* GeometricBrdf::genBrdfImage(const Query &query, int resolution) {
    const int ndfResolution = resolution * 2;
    Float *ndfImage = genNdfImage(query, ndfResolution);

    Float *brdfImage = new Float[resolution * resolution * 3];
    for (int i = 0; i < resolution; i++) {
        for (int j = 0; j < resolution; j++) {
            brdfImage[(i * resolution + j) * 3 + 0] = 0.0;
            brdfImage[(i * resolution + j) * 3 + 1] = 0.0;
            brdfImage[(i * resolution + j) * 3 + 2] = 0.0;

            const int numSamples = 16;
            for (int k = 0; k < numSamples; k++) {
                Vector2 omega_o((i + randUniform<Float>()) / resolution * 2.0 - 1.0,
                                (j + randUniform<Float>()) / resolution * 2.0 - 1.0);

                Vector3 omega_i_3(query.omega_i(0), query.omega_i(1), sqrt(1.0 - query.omega_i.norm()));
                Vector3 omega_o_3(omega_o(0), omega_o(1), sqrt(1.0 - omega_o.norm()));
                Vector2 omega_h = (omega_i_3 + omega_o_3).normalized().head(2);

                Vector2 xy = (omega_h + Vector2(1.0, 1.0)) / 2.0 * ndfResolution;

                int xInt = (int)(xy(0));
                int yInt = (int)(xy(1));
                if (xInt < 0 || xInt >= ndfResolution ||
                    yInt < 0 || yInt >= ndfResolution)
                    continue;

                Float D = ndfImage[(xInt * ndfResolution + yInt) * 3 + 0];
                Float F = 1.0;
                Float G = 1.0;
                Float cosThetaI = sqrt(abs(1.0 - query.omega_i.dot(query.omega_i)));
                Float cosThetaO = sqrt(abs(1.0 - omega_o.dot(omega_o)));
                Float brdfValue = D * F * G / (4.0 * cosThetaI * cosThetaO);

                if (std::isnan(brdfValue / numSamples))
                    continue;

                brdfImage[(i * resolution + j) * 3 + 0] += brdfValue / numSamples;
                brdfImage[(i * resolution + j) * 3 + 1] += brdfValue / numSamples;
                brdfImage[(i * resolution + j) * 3 + 2] += brdfValue / numSamples;
            }
        }
    }

    delete[] ndfImage;
    return brdfImage;
}
