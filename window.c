/*
fir - just playing around with fir filters
Written in 2013 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <math.h>
#include "window.h"

float delta(float x)
{
	return 0.0f == x ? 1.0f : 0.0f;
}
float sinc(float x)
{
	return 0.0f == x ? 1.0f : sinf((float)M_PI * x) / ((float)M_PI * x);
}
float rect(float n, float N, float a)
{
	(void)a;
	return n >= 0.0f && n < N ? 1.0f : 0.0f;
}
float hann(float n, float N, float a)
{
	(void)a;
	return 0.5f * (1.0f - cosf(2.0f * (float)M_PI * n / (N - 1.0f)));
}
float hamming(float n, float N, float a)
{
	(void)a;
	return 0.54f - 0.46f * cosf(2.0f * (float)M_PI * n / (N - 1.0f));
}
float lanczos(float n, float N, float a)
{
	(void)a;
	return sinc(2.0f * n / (N - 1.0f) - 1.0f);
}
float gauss(float n, float N, float o)
{
	return expf(- 1.0f/2.0f * powf((n - (N - 1.0f) / 2.0f) / (o * (N - 1.0f) / 2.0f), 2.0f));
}
float blackman(float n, float N, float a)
{
	return 0.5f * (1.0f - a) - 0.5f * cosf(2.0f * (float)M_PI * n / (N - 1.0f)) + 0.5f * a * cosf(4.0f * (float)M_PI * n / (N - 1.0f));
}
float exact_blackman(float n, float N, float a)
{
	(void)a;
	return 7938.0f / 18608.0f - 9240.0f / 18608.0f * cosf(2.0f * (float)M_PI * n / (N - 1.0f)) + 1430.0f / 18608.0f * cosf(4.0f * (float)M_PI * n / (N - 1.0f));
}
float i0f(float x)
{
	// converges for -3*M_PI:3*M_PI in less than 20 iterations
	float sum = 1.0f, val = 1.0f, c = 0.0f;
	for (int n = 1; n < 20; n++) {
		float tmp = x / (2.0f * (float)n);
		val *= tmp * tmp;
		float y = val - c;
		float t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	return sum;
}
float kaiser(float n, float N, float a)
{
	return i0f((float)M_PI * a * sqrtf(1.0f - powf((2.0f * n) / (N - 1.0f) - 1.0f, 2.0f))) / i0f((float)M_PI * a);
}

