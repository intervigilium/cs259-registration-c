/*
 * Registration kernel for FPGA implementation
 *
 * left/right/up/down/in/out sweeps can be parallelized
 * on outer loop (counter) and middle loop (sweep row)
 * left/right/up/down/in/out sweeps can be pipelined
 * on outer loop (counter) and middle loop (sweep row)
 * left/right/up/down/in/out sweeps can be stenciled
 * on middle loop (sweep row)
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

	/* sweep row stencil */
	float v1_cache[M];
	float v2_cache[M];
	float v3_cache[M];

	for (i = 0 i < 3 * GAUSSIAN_STEP; i++) {
		post_scale *= nu / lambda;
	}

	for (step = 0; step < GAUSSIAN_STEP; step++) {
		/* filter right
		 * loop is split into boundary and main loops
		 * loop index by j first due to data dependency */
		for (counter = 0; counter < M; counter++) {
			for (i = 0; i < N; i++) {
				V1[counter][i][0] *= boundary_scale;
				v1_cache[i] = V1[counter][i][0];
				V2[counter][i][0] *= boundary_scale;
				v2_cache[i] = V2[counter][i][0];
				V3[counter][i][0] *= boundary_scale;
				v3_cache[i] = V3[counter][i][0];
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (j = 1; j < P; j++) {
				for (i = 0; i < N; i++) {
					V1[counter][i][j] += v1_cache[i] * nu;
					v1_cache[i] = V1[counter][i][j];
					V2[counter][i][j] += v2_cache[i] * nu;
					v2_cache[i] = V2[counter][i][j];
					V3[counter][i][j] += v3_cache[i] * nu;
					v3_cache[i] = V3[counter][i][j];
				}
			}
		}

		/* filter left
		 * previous filter right modifies boundary */
		for (counter = 0; counter < M; counter++) {
			for (i = 0; i < N; i++) {
				V1[counter][i][P - 1] *= boundary_scale;
				v1_cache[i] = V1[counter][i][P - 1];
				V2[counter][i][P - 1] *= boundary_scale;
				v2_cache[i] = V2[counter][i][P - 1];
				V3[counter][i][P - 1] *= boundary_scale;
				v3_cache[i] = V3[counter][i][P - 1];
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (j = P - 2; i >= 0; j--) {
				for (i = 0; i < N; i++) {
					V1[counter][i][j] += v1_cache[i] * nu;
					v1_cache[i] = V1[counter][i][j];
					V2[counter][i][j] += v2_cache[i] * nu;
					v2_cache[i] = V2[counter][i][j];
					V3[counter][i][j] += v3_cache[i] * nu;
					v3_cache[i] = V3[counter][i][j];
				}
			}
		}

		/* filter down */
		for (counter = 0; counter < M; counter++) {
			for (j = 0; j < P; j++) {
				V1[counter][0][j] *= boundary_scale;
				v1_cache[j] = V1[counter][0][j];
				V2[counter][0][j] *= boundary_scale;
				v2_cache[j] = V2[counter][0][j];
				V3[counter][0][j] *= boundary_scale;
				v3_cache[j] = V3[counter][0][j];
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (i = 1; i < N; i++) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] += v1_cache[j] * nu;
					v1_cache[j] = V1[counter][i][j];
					V2[counter][i][j] += v2_cache[j] * nu;
					v2_cache[j] = V2[counter][i][j];
					V3[counter][i][j] += v3_cache[j] * nu;
					v3_cache[j] = V3[counter][i][j];
				}
			}
		}

		/* filter up */
		for (counter = 0; counter < M; counter++) {
			for (j = 0; j < P; j++) {
				V1[counter][N - 1][j] *= boundary_scale;
				v1_cache[j] = V1[counter][N - 1][j];
				V2[counter][N - 1][j] *= boundary_scale;
				v2_cache[j] = V2[counter][N - 1][j];
				V3[counter][N - 1][j] *= boundary_scale;
				v3_cache[j] = V3[counter][N - 1][j];
			}
		}
		for (counter = 0; counter < M; counter++) {
			for (i = N - 2; i >= 0; i--) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] += v1_cache[j] * nu;
					v1_cache[j] = V1[counter][i][j];
					V2[counter][i][j] += v2_cache[j] * nu;
					v2_cache[j] = V2[counter][i][j];
					V3[counter][i][j] += v3_cache[j] * nu;
					v3_cache[j] = V3[counter][i][j];
				}
			}
		}

		/* filter in */
		for (i = 0; i < N; i++) {
			for (j = 0; j < P; j++) {
				V1[0][i][j] *= boundary_scale;
				v1_cache[j] = V1[0][i][j];
				V2[0][i][j] *= boundary_scale;
				v2_cache[j] = V2[0][i][j];
				V3[0][i][j] *= boundary_scale;
				v3_cache[j] = V3[0][i][j];
			}
		}
		for (counter = 1; counter < M; counter++) {
			for (i = 0; i < N; i++) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] += v1_cache[j] * nu;
					v1_cache[j] = V1[counter][i][j];
					V2[counter][i][j] += v2_cache[j] * nu;
					v2_cache[j] = V2[counter][i][j];
					V3[counter][i][j] += v3_cache[j] * nu;
					v3_cache[j] = V3[counter][i][j];
				}
			}
		}

		/* filter out */
		for (i = 0; i < N; i++) {
			for (j = 0; j < P; j++) {
				V1[M - 1][i][j] *= boundary_scale;
				v1_cache[j] = V1[M - 1][i][j];
				V2[M - 1][i][j] *= boundary_scale;
				v2_cache[j] = V2[M - 1][i][j];
				V3[M - 1][i][j] *= boundary_scale;
				v3_cache[j] = V3[M - 1][i][j];
			}
		}
		for (counter = M - 2; counter >= 0; counter--) {
			for (i = 0; i < N; i++) {
				for (j = 0; j < P; j++) {
					V1[counter][i][j] += v1_cache[j] * nu;
					v1_cache[j] = V1[counter][i][j];
					V2[counter][i][j] += v2_cache[j] * nu;
					v2_cache[j] = V2[counter][i][j];
					V3[counter][i][j] += v3_cache[j] * nu;
					v3_cache[j] = V3[counter][i][j];
				}
			}
		}
	}

	/* post-scale can be completely pipelined and parallelized */
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
