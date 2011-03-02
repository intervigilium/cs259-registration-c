#include "convolution2d.h"
#include "fft.h"

void convolution2D(float Are[IMAX][IMAX])
{
	int i, j;
	float Aim[IMAX][IMAX];
	for (i = 0; i < IMAX; i++)
		for (j = 0; j < IMAX; j++)
			Aim[i][j] = 0.0f;

	fft2D(Are, Aim, 0);

	float a, b, c, d;
	for (i = 0; i < IMAX; i++) {
		for (j = 0; j < IMAX; j++) {
			a = Are[i][j];
			b = Aim[i][j];
			c = Bre[i][j];
			d = Bim[i][j];
			Are[i][j] = a * c - b * d;
			Aim[i][j] = a * d + b * c;
		}
	}

	fft2D(Are, Aim, 1);
}
