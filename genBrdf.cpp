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

#include <iostream>
#include <unistd.h>
#include <chrono>
#include <string>
#include "waveBrdf.h"
#include "spectrum.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
    srand(time(NULL));

    char *heightfieldFilename = NULL;
    double texelWidth = 1.0;
    double vertScale = 1.0;

    double mu_x = 0.0;
    double mu_y = 0.0;
    double sigma_p = 10.0;

    string method = "Wave";
    int sampleNum = 10000000;
    string diffModel = "OHS";
    double lambda = 0.0;

    double omega_i_x = 0.0;
    double omega_i_y = 0.0;

    char *outputFilename = NULL;
    int outputResolution = 256;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:l:r:e:w:v:x:y:p:m:g:d:s:t:n:")) != -1) {
        switch (opt) {
            case 'i': heightfieldFilename = optarg; break;      // heightfield filename.
            case 'w': texelWidth = atof(optarg); break;         // The width of a texel in microns on the heightfield.
            case 'v': vertScale = atof(optarg); break;          // The vertical scaling factor of the heightfield.

            case 'x': mu_x = atof(optarg); break;               // Center x of the Gaussian footprint.
            case 'y': mu_y = atof(optarg); break;               // Center y of the Gaussian footprint.
            case 'p': sigma_p = atof(optarg); break;            // Size (1 sigma) of the Gaussian footprint.

            case 'm': method = optarg; break;                   // Method. Choose between "Geom" and "Wave".
            case 'n': sampleNum = atoi(optarg); break;          // Number of binning samples. Only valid for geometric optics.
            case 'd': diffModel = optarg; break;                // Diffraction model. Choose between "OHS", "GHS", "ROHS", "RGHS" and "Kirchhoff". And expect only subtle difference.
            case 'l': lambda = atof(optarg); break;             // Wavelength in microns. Once set, single wavelength mode is ON.

            case 's': omega_i_x = atof(optarg); break;          // Incoming light's x coordinate (assuming z = 1).
            case 't': omega_i_y = atof(optarg); break;          // Incoming light's y coordinate (assuming z = 1).

            case 'o': outputFilename = optarg; break;           // output filename.
            case 'r': outputResolution = atoi(optarg); break;   // output resolution.
        }
    }

    if (heightfieldFilename == NULL) {
        cout << "A heightfield file must be specified." << endl;
        return -1;
    }

    if (outputFilename == NULL) {
        cout << "An output file must be specified." << endl;
        return -1;
    }

    SpectrumInit();

    EXRImage heightfieldImage(heightfieldFilename);
    Heightfield heightfield(&heightfieldImage, texelWidth, vertScale);

    Query query;
    query.mu_p = Vector2(mu_x, mu_y);
    query.sigma_p = sigma_p;
    query.omega_i = Vector3(omega_i_x, omega_i_y, 1.0).normalized().head(2);
    query.lambda = lambda;

    if (method == "Geom") {
        GeometricBrdf geometricBrdf(&heightfield, sampleNum);
        Float *brdfImage = geometricBrdf.genBrdfImage(query, outputResolution);
        EXRImage::writeImage(brdfImage, outputFilename, outputResolution, outputResolution);
        delete[] brdfImage;
    } else if (method == "Wave") {
        WaveBrdfAccel waveBrdfAccel(&heightfield, diffModel);
        Float *brdfImage = waveBrdfAccel.genBrdfImage(query, outputResolution);
        EXRImage::writeImage(brdfImage, outputFilename, outputResolution, outputResolution);
        delete[] brdfImage;
    } else if (method == "GeomNdf") {
        GeometricBrdf geometricBrdf(&heightfield, sampleNum);
        Float *ndfImage = geometricBrdf.genNdfImage(query, outputResolution);
        EXRImage::writeImage(ndfImage, outputFilename, outputResolution, outputResolution);
        delete[] ndfImage;
    }

    return 0;
}
