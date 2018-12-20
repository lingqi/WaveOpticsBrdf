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


void WaveNDF::mkCrop()
{
	cropped.resize(crop, crop);
	int i = (resolution - crop) / 2;
	cropped = rgb.block(i, i, crop, crop);
}


void WaveNDF::generate(const Query& query, const char* outputFilename)
{
	if (query.lambda <= 0) { generateSpectral(query, outputFilename); return; }

	int n = resolution;
	int nh = n / 2;
	img.resize(n, n);
	out.resize(n, n);

	// normalizing constants
	float Ac = M_PI * query.sigma_p * query.sigma_p; // integral of *squared* footprint
	float C = 1 / (Ac * query.lambda * query.lambda); // see Wave NDF doc
	C *= 4; // this approximates the (psi . n)^2 term for small angles
	float T = 2.0f * k * query.sigma_p / n; // sample step in primal domain

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			// map window to [-k, k]^2
			float u = (i + 0.5f - nh) / nh;
			float v = (j + 0.5f - nh) / nh;
			u *= k; v *= k;
			
			// unit Gaussian window
			Float g = std::exp(-0.5f * (u*u + v*v));

			// map to query footprint [-k sigma, k sigma]^2
			u *= query.sigma_p; v *= query.sigma_p;
			u += query.mu_p[0]; v += query.mu_p[1];
			
			// lookup heightfield and compute R^*(u)
			float h = hf.eval(u, v);
			img(i, j) = g * cis(-4 * float(M_PI) * h / query.lambda);
		}

	if (visRs)
	{
		rgb.resize(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				rgb(i, j).setConstant(img(i, j).real());
		EXRImage::writeImage((float*) &rgb(0,0), "Rs.exr", n, n);
	}

	fftshift();
	fft2();
	fftshift();
	img *= T * T; // area of one cell

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			out(i, j) = norm(img(i, j)); // "norm" is squared, per C++ complex docs

	// normalize
	out *= C;

	if (outputFilename != nullptr)
	{
		rgb.resize(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				rgb(i, j).setConstant(out(i, j));

		float* p = (float*) &rgb(0,0);
		if (crop > 0) { mkCrop(); p = (float*) &cropped(0,0); n = crop; }
		EXRImage::writeImage(p, outputFilename, n, n);
	}
}


void WaveNDF::generateSpectral(const Query& query, const char* outputFilename)
{
	const int s = SPECTRUM_SAMPLES;
	lambdas.resize(s);
	colors.resize(s);
	channels.resize(s);

	for (int i = 0; i < s; i++)
		lambdas[i] = (i + 0.5) / s * (0.7 - 0.3) + 0.3;

	for (int i = 0; i < s; i++)
	{
		Query q = query;
		q.lambda = lambdas[i];
		q.sigma_p *= lambdas[i] / 0.5;
		generate(q, nullptr);
		channels[i] = out;
	}

	float r, g, b;
	int n = resolution;
	rgb.resize(n, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < s; k++) colors[k] = channels[k](i, j);
			SpectrumToRGB(colors, r, g, b);
			rgb(i, j) = Color(r, g, b);
		}

	float* p = (float*) &rgb(0,0);
	if (crop > 0) { mkCrop(); p = (float*) &cropped(0,0); n = crop; }
	EXRImage::writeImage(p, outputFilename, n, n);
}
