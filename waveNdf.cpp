#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <vector>
#include <iostream>
#include "waveNdf.h"

using namespace std;
using namespace Eigen;


void WaveNDF::fft2()
{
    FFT<float> fft;
	int n = img.rows();
	tmp.resize(n, n);
	a.resize(n); b.resize(n);

	for (int i = 0; i < n; i++)
	{
		a = img.row(i);
		fft.fwd(b, a);
		tmp.row(i) = b;
	}

	for (int i = 0; i < n; i++)
	{
		a = tmp.col(i);
		fft.fwd(b, a);
		img.col(i) = b;
	}
}


void WaveNDF::fftshift()
{
	int n = img.rows();
	tmp.resize(n, n);
	int nh = n / 2;

	tmp.rightCols(nh) = img.leftCols(nh);
	tmp.leftCols(nh) = img.rightCols(nh);
	img.topRows(nh) = tmp.bottomRows(nh);
	img.bottomRows(nh) = tmp.topRows(nh);
}


void WaveNDF::generate(const Query& query, char* outputFilename)
{
	int n = resolution;
	int nh = n / 2;
	img.resize(n, n);
	out.resize(n, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			// map window to [-3, 3]^2
			float u = float(i - nh) / nh;
			float v = float(j - nh) / nh;
			u *= 3; v *= 3;
			
			// unit Gaussian window
			Float g = std::exp(-0.5f * (u*u + v*v));

			// map to query footprint [-3 sigma, 3 sigma]^2
			u *= query.sigma_p; v *= query.sigma_p;
			u += query.mu_p[0]; v += query.mu_p[1];
			
			// lookup heightfield and compute R^*(u)
			float h = hf.eval(u, v);
			img(i, j) = g * cis(-4 * float(M_PI) * h / query.lambda);
		}

	fftshift();
	fft2();
	fftshift();

	Float mx = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			float v = norm(img(i, j));
			mx = std::max(mx, v);
			out(i, j).setConstant(v); // already squared=
		}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			out(i, j) /= mx;

	EXRImage::writeImage((float*) &out(0,0), outputFilename, n, n);
}
