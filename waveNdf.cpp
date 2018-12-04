#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <vector>
#include <iostream>
#include "spectrum.h"
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


void WaveNDF::generate(Query& query, char* outputFilename)
{
	if (query.lambda <= 0) { generateSpectral(query, outputFilename); return; }

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

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			out(i, j) = norm(img(i, j)); // already squared

	out *= 1.0f / (resolution * resolution);

	if (outputFilename != nullptr)
	{
		ColorImage rgb(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				rgb(i, j).setConstant(out(i, j));
		EXRImage::writeImage((float*) &rgb(0,0), outputFilename, n, n);
	}
}


void WaveNDF::generateSpectral(Query& query, char* outputFilename)
{
	const int s = SPECTRUM_SAMPLES;
	vector<float> lambdas(s);
	for (int i = 0; i < s; i++)
		lambdas[i] = (i + 0.5) / s * (0.7 - 0.3) + 0.3;

	vector<FloatImage> channels(s);
	for (int i = 0; i < s; i++)
	{
		query.lambda = lambdas[i];
		generate(query, nullptr);
		channels[i] = out;
	}

	vector<float> tmp(s);
	float r, g, b;
	int n = resolution;
	ColorImage rgb(n, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < s; k++) tmp[k] = channels[k](i, j);
			SpectrumToRGB(tmp, r, g, b);
			rgb(i, j) = Color(r, g, b);
		}

	EXRImage::writeImage((float*) &rgb(0,0), outputFilename, n, n);
}
