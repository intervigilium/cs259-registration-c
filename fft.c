#include "fft.h"

static void permute_bitrev(float Are[], float Aim[], int size)
{
	int i, bri;
	int j, inp, rev = 0;
	float t_re, t_im;

	for (i = 0; i < size; i++) {
		inp = i;
		rev = 0;
		for (j = 0; j < 7; j++) {
			rev = (rev << 1) | (inp & 1);
			inp >>= 1;
		}
		bri = rev;

		/* skip already swapped elements */
		if (bri <= i)
			continue;
		t_re = Are[i];
		t_im = Aim[i];
		Are[i] = Are[bri];
		Aim[i] = Aim[bri];
		Are[bri] = t_re;
		Aim[bri] = t_im;
	}
}

static void fft(float A_re[], float A_im[], const float W_re[],
		const float W_im[], int pp)
{
	float w_re, w_im, u_re, u_im, t_re, t_im;
	int mm, g, b;
	int i, mt, k;
	int nn = IMAX;

	/* for each stage */
	for (mm = nn; mm >= 2; mm = mm >> 1) {
		/* m = n/2^s; mt = m/2; */
		mt = mm >> 1;

		/* for each group of butterfly */
		for (g = 0, k = 0; g < nn; g += mm, k++) {
			w_re = W_re[k];
			w_im = W_im[k];

			/* for each butterfly */
			for (b = g; b < (g + mt); b++) {
				t_re =
				    w_re * A_re[b + mt] - w_im * A_im[b + mt];
				t_im =
				    w_re * A_im[b + mt] + w_im * A_re[b + mt];

				u_re = A_re[b];
				u_im = A_im[b];
				A_re[b] = u_re + t_re;
				A_im[b] = u_im + t_im;
				A_re[b + mt] = u_re - t_re;
				A_im[b + mt] = u_im - t_im;
			}
		}
	}

	if (pp == 1) {
		for (i = 0; i < nn; i++) {
			A_re[i] = A_re[i] / nn;
			A_im[i] = A_im[i] / nn;
		}
	}
}

static void fft1D(float Are[], float Aim[], int pp)
{
	const float Wre[] = {
		1, 6.12323e-17, 0.707107, -0.707107,
		0.92388, -0.382683, 0.382683, -0.92388,
		0.980785, -0.19509, 0.55557, -0.83147,
		0.83147, -0.55557, 0.19509, -0.980785,
		0.995185, -0.0980171, 0.634393, -0.77301,
		0.881921, -0.471397, 0.290285, -0.95694,
		0.95694, -0.290285, 0.471397, -0.881921,
		0.77301, -0.634393, 0.0980171, -0.995185,
		0.998795, -0.0490677, 0.671559, -0.740951,
		0.903989, -0.427555, 0.33689, -0.941544,
		0.970031, -0.24298, 0.514103, -0.857729,
		0.803208, -0.595699, 0.14673, -0.989177,
		0.989177, -0.14673, 0.595699, -0.803208,
		0.857729, -0.514103, 0.24298, -0.970031,
		0.941544, -0.33689, 0.427555, -0.903989,
		0.740951, -0.671559, 0.0490677, -0.998795
	};

	const float Wim_c[] = {
		-0, -1, -0.707107, -0.707107,
		-0.382683, -0.92388, -0.92388, -0.382683,
		-0.19509, -0.980785, -0.83147, -0.55557,
		-0.55557, -0.83147, -0.980785, -0.19509,
		-0.0980171, -0.995185, -0.77301, -0.634393,
		-0.471397, -0.881921, -0.95694, -0.290285,
		-0.290285, -0.95694, -0.881921, -0.471397,
		-0.634393, -0.77301, -0.995185, -0.0980171,
		-0.0490677, -0.998795, -0.740951, -0.671559,
		-0.427555, -0.903989, -0.941544, -0.33689,
		-0.24298, -0.970031, -0.857729, -0.514103,
		-0.595699, -0.803208, -0.989177, -0.14673,
		-0.14673, -0.989177, -0.803208, -0.595699,
		-0.514103, -0.857729, -0.970031, -0.24298,
		-0.33689, -0.941544, -0.903989, -0.427555,
		-0.671559, -0.740951, -0.998795, -0.0490677
	};
	float Wim[IMAX / 2];
	int i;
	if (pp == 1) {
		for (i = 0; i < IMAX / 2; i++)
			Wim[i] = -Wim_c[i];
	} else {
		for (i = 0; i < IMAX / 2; i++)
			Wim[i] = Wim_c[i];
	}
	fft(Are, Aim, Wre, Wim, pp);
	permute_bitrev(Are, Aim, IMAX);
}

void fft2D(float Are[][], float Aim[][], int pp)
{
	int i, j;
	float A1_re[IMAX], A1_im[IMAX];

	for (i = 0; i < IMAX; i++) {
		for (j = 0; j < IMAX; j++) {
			A1_re[j] = Are[i][j];
			A1_im[j] = Aim[i][j];
		}
		fft1D(A1_re, A1_im, pp);
		for (j = 0; j < IMAX; j++) {
			Are[i][j] = A1_re[j];
			Aim[i][j] = A1_im[j];
		}
	}

	for (j = 0; j < IMAX; j++) {
		for (i = 0; i < IMAX; i++) {
			A1_re[i] = Are[i][j];
			A1_im[i] = Aim[i][j];
		}
		fft1D(A1_re, A1_im, pp);
		for (i = 0; i < IMAX; i++) {
			Are[i][j] = A1_re[i];
			Aim[i][j] = A1_im[i];
		}
	}
}
