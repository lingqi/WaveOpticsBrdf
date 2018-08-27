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

/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
using namespace std;

const int SPECTRUM_SAMPLES = 8;
const int CIE_samples = 471;
extern const float CIE_wavelengths[CIE_samples];
extern const float CIE_X_entries[CIE_samples];
extern const float CIE_Y_entries[CIE_samples];
extern const float CIE_Z_entries[CIE_samples];
extern const float CIE_D65_entries[CIE_samples];

void SpectrumInit();

void SpectrumToXYZ(const vector<float> &s, float &x, float &y, float &z);

void XYZToRGB(float x, float y, float z, float &r, float &g, float &b);

void SpectrumToRGB(const vector<float> &s, float &r, float &g, float &b);

#endif
