/*
fir - just playing around with fir filters
Written in 2013 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef WINDOW_H
#define WINDOW_H
float delta(float);
float sinc(float);
float rect(float, float, float);
float hann(float, float, float);
float hamming(float, float, float);
float lanczos(float, float, float);
float gauss(float, float, float);
float blackman(float, float, float);
float exact_blackman(float, float, float);
float kaiser(float, float, float);

#endif

