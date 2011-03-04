/*
 * Registration kernel for FPGA implementation
 */

#include <stdlib.h>
#include <math.h>

#define EPS 1e-6f
#define IMAX 128
#define GAUSSIAN_STEP 3

#define M 256
#define N 256
#define P 256

void V_cal(float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	   float V3[M + 1][N + 1][P + 1])
{
	int i, j, k, step, counter;
	float std = 4.5f;

	float lambda = (std * std) / (2 * GAUSSIAN_STEP);
	float nu = (1 + 2 * lambda - sqrtf(1 + 4 * lambda)) / (2 * lambda);
	float boundary_scale = 1 / (1 - nu);
	float post_scale = 1.0f;

	for (i = 0 i < 3 * GAUSSIAN_STEP; i++) {
		post_scale *= nu / lambda;
	}

	for (step = 0; step < GAUSSIAN_STEP; step++) {
		/* filter right */
		for (counter = 0; counter < M; counter++) {
			for (i = 0; i < N; i++) {
				V1[counter][i][0] *= boundary_scale;
				V2[counter][i][0] *= boundary_scale;
				V3[counter][i][0] *= boundary_scale;
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (j = 1; j < P; j++) {
				for (i = 0; i < N; i++) {
					V1[counter][i][j] +=
					    V1[counter][i][j - 1] * nu;
					V2[counter][i][j] +=
					    V2[counter][i][j - 1] * nu;
					V3[counter][i][j] +=
					    V3[counter][i][j - 1] * nu;
				}
			}
		}

		/* filter left */
		for (counter = 0; counter < M; counter++) {
			for (i = 0; i < N; i++) {
				V1[counter][i][P - 1] *= boundary_scale;
				V2[counter][i][P - 1] *= boundary_scale;
				V3[counter][i][P - 1] *= boundary_scale;
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (j = P - 2; i >= 0; j--) {
				for (i = 0; i < N; i++) {
					V1[counter][i][j] +=
					    V1[counter][i][j + 1] * nu;
					V2[counter][i][j] +=
					    V2[counter][i][j + 1] * nu;
					V2[counter][i][j] +=
					    V3[counter][i][j + 1] * nu;
				}
			}
		}

		/* filter down */
		for (counter = 0; counter < M; counter++) {
			for (j = 0; j < P; j++) {
				V1[counter][0][j] *= boundary_scale;
				V2[counter][0][j] *= boundary_scale;
				V3[counter][0][j] *= boundary_scale;
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (i = 1; i < N; i++) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] +=
					    V1[counter][i - 1][j] * nu;
					V2[counter][i][j] +=
					    V2[counter][i - 1][j] * nu;
					V3[counter][i][j] +=
					    V3[counter][i - 1][j] * nu;
				}
			}
		}

		/* filter up */
		for (counter = 0; counter < M; counter++) {
			for (j = 0; j < P; j++) {
				V1[counter][N - 1][j] *= boundary_scale;
				V2[counter][N - 1][j] *= boundary_scale;
				V3[counter][N - 1][j] *= boundary_scale;
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (i = N - 2; i >= 0; i--) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] +=
					    V1[counter][i + 1][j] * nu;
					V2[counter][i][j] +=
					    V2[counter][i + 1][j] * nu;
					V3[counter][i][j] +=
					    V3[counter][i + 1][j] * nu;
				}
			}
		}

		/* filter in */
		for (i = 0; i < N; i++) {
			for (j = 0; j < P; j++) {
				V1[0][i][j] *= boundary_scale;
				V2[0][i][j] *= boundary_scale;
				V3[0][i][j] *= boundary_scale;
			}
		}
		for (counter = 1; counter < M; counter++) {
			for (i = 0; i < N; i++) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] +=
					    V1[counter - 1][i][j] * nu;
					V2[counter][i][j] +=
					    V2[counter - 1][i][j] * nu;
					V3[counter][i][j] +=
					    V3[counter - 1][i][j] * nu;
				}
			}
		}

		/* filter out */
		for (i = 0; i < N; i++) {
			for (j = 0; j < P; j++) {
				V1[M - 1][i][j] *= boundary_scale;
				V2[M - 1][i][j] *= boundary_scale;
				V3[M - 1][i][j] *= boundary_scale;
			}
		}
		for (counter = M - 2; counter >= 0; counter--) {
			for (i = 0; i < N; i++) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] +=
					    V1[counter + 1][i][j] * nu;
					V2[counter][i][j] +=
					    V2[counter + 1][i][j] * nu;
					V3[counter][i][j] +=
					    V3[counter + 1][i][j] * nu;
				}
			}
		}
	}

	/* post-scale */
	for (counter = 0; counter < M; counter++) {
		for (i = 0; i < N; i++) {
			for (j = 0; j < P; j++) {
				V1[counter][i][j] *= post_scale;
				V2[counter][i][j] *= post_scale;
				V3[counter][i][j] *= post_scale;
			}
		}
	}
}
