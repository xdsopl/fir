/*
fir - just playing around with fir filters
Written in 2013 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "window.h"

void dft(complex float *X, complex float *x, int N)
{
	for (int k = 0; k < N; k++) {
		X[k] = 0;
		for (int n = 0; n < N; n++)
			X[k] += x[n] * cexpf(-I * 2.0f * (float)M_PI * (float)k * (float)n / (float)N);
	}
}

void lpf(complex float *b, float rate, float freq, int N, float (*window)(float, float, float), float a)
{
	float sum = 0.0f;
	for (int n = 0; n < N; n++) {
		float x = (float)n - (float)(N - 1) / 2.0f;
		float w = window(n, N, a);
		float v = 2.0f * freq / rate * sinc(2.0f * freq / rate * x);
		sum += v * w;
		b[n] = v * w;
	}
	for (int n = 0; n < N; n++)
		b[n] /= sum;
}

void hpf(complex float *b, float rate, float freq, int N, float (*window)(float, float, float), float a)
{
	float sum = 0.0f;
	for (int n = 0; n < N; n++) {
		float x = (float)n - (float)(N - 1) / 2.0f;
		float w = window(n, N, a);
		float v = 2.0f * freq / rate * sinc(2.0f * freq / rate * x);
		sum += v * w;
		b[n] = delta(x) - v * w;
	}
	for (int n = 0; n < N; n++)
		b[n] /= sum;
}
void bpf(complex float *b, float rate, float freq, float bw, int N, float (*window)(float, float, float), float a)
{
	float sum = 0.0f;
	float freq0 = freq - bw / 2.0f;
	float freq1 = freq + bw / 2.0f;
	for (int n = 0; n < N; n++) {
		float x = (float)n - (float)(N - 1) / 2.0f;
		float w = window(n, N, a);
		float v0 = 2.0f * freq0 / rate * sinc(2.0f * freq0 / rate * x);
		float v1 = 2.0f * freq1 / rate * sinc(2.0f * freq1 / rate * x);
		sum += 0.5f * (v0 + v1) * w;
		b[n] = (v1 - v0) * w;
	}
	for (int n = 0; n < N; n++)
		b[n] /= sum;
}

void single_dft(complex float *b, float rate, float freq, int N, float (*window)(float, float, float), float a)
{
	float sum = 0.0f;
	for (int n = 0; n < N; n++) {
		float x = (float)(N - 1) / 2.0f - (float)n;
		float w = window(n, N, a);
		complex float v = cexpf(-I * 2.0f * (float)M_PI * freq / rate * x);
		sum += cabsf(v) * w;
		b[n] = v * w;
	}
	for (int n = 0; n < N; n++)
		b[n] /= sum;
}

float decibel(float amp)
{
	return 10.0f * log10f(amp * amp);
}

int main(int argc, char **argv)
{
	(void)argc; (void)argv;

	float rate = 11025.0f;
	float freq = 1000.0f;
	int taps = 511;

	complex float x[taps];

	float (*window)(float, float, float) = kaiser;
	float a = 2.0f;

	if (2 == argc && !strcmp(argv[1], "rect")) {
		window = rect;
		a = 0.0f;
	} else if (3 == argc && !strcmp(argv[1], "gauss")) {
		window = gauss;
		a = atof(argv[2]);
	} else if (2 == argc && !strcmp(argv[1], "hann")) {
		window = hann;
		a = 0.0f;
	} else if (2 == argc && !strcmp(argv[1], "hamming")) {
		window = hamming;
		a = 0.0f;
	} else if (2 == argc && !strcmp(argv[1], "blackman")) {
		window = exact_blackman;
		a = 0.0f;
	} else if (2 == argc && !strcmp(argv[1], "lanczos")) {
		window = lanczos;
		a = 0.0f;
	} else if (3 == argc && !strcmp(argv[1], "kaiser")) {
		window = kaiser;
		a = atof(argv[2]);
	}

	lpf(x, rate, freq, taps, window, a);
//	hpf(x, rate, freq, taps, window, a);
//	bpf(x, rate, freq, rate / (float)taps, taps, window, a);
//	single_dft(x, rate, freq, taps, window, a);
#if 0
	for (int i = 0; i < taps; i++)
		printf("%f %f %f\n", (float)i * rate / (float)taps, crealf(x[i]), cimagf(x[i]));
#endif

#if 0
	for (int i = 0; i < taps; i++) {
		complex float sum = 0.0f;
		for (int j = 0; j <= i; j++)
			sum += x[j];
		printf("%f %f %f\n", (float)i * rate / (float)taps, crealf(sum), cimagf(sum));
	}
#endif
	complex float X[taps];
	dft(X, x, taps);
#if 1
	for (int i = 0; i < taps; i++)
		printf("%f %f\n", (float)i * rate / (float)taps, decibel(cabsf(X[i])));
#endif

#if 0
	for (int i = 0; i < (taps - 1); i++)
		printf("%f %f\n", (float)i * rate / (float)taps, -cargf(X[i+1] / X[i]) / (2.0f * (float)M_PI) * (float)taps * cabsf(X[i]));
#endif
	return 0;
}

