#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "registration.h"
#include "convolution2d.h"

void V_cal(float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	   float V3[M + 1][N + 1][P + 1])
{
	int i, j, k;
	float std = 4.5f;

	float lambda = (std * std) / (2 * GAUSSIAN_STEP);
	float nu = (1 + 2 * lambda - sqrtf(1 + 4 * lambda)) / (2 * lambda);
	float BoundaryScale = 1 / (1 - nu);
	int counter;
	float PostScale = 1;
	for (i = 0; i < 3 * GAUSSIAN_STEP; i++)
		PostScale *= nu / lambda;

	int step;

	for (step = 0; step < GAUSSIAN_STEP; step++) {
		//filter right 
		for (counter = 0; counter < M; counter++)
			for (i = 0; i < N; i++) {
				V1[counter][i][0] =
				    V1[counter][i][0] * BoundaryScale;
				V2[counter][i][0] =
				    V2[counter][i][0] * BoundaryScale;
				V3[counter][i][0] =
				    V3[counter][i][0] * BoundaryScale;
			}
		for (counter = 0; counter < M; counter++)
			for (j = 1; j < P; j++)
				for (i = 0; i < N; i++) {
					V1[counter][i][j] =
					    V1[counter][i][j] +
					    V1[counter][i][j - 1] * nu;
					V2[counter][i][j] =
					    V2[counter][i][j] +
					    V2[counter][i][j - 1] * nu;
					V3[counter][i][j] =
					    V3[counter][i][j] +
					    V3[counter][i][j - 1] * nu;
				}

		//filter left
		for (counter = 0; counter < M; counter++)
			for (i = 0; i < N; i++) {
				V1[counter][i][P - 1] =
				    V1[counter][i][P - 1] * BoundaryScale;
				V2[counter][i][P - 1] =
				    V2[counter][i][P - 1] * BoundaryScale;
				V3[counter][i][P - 1] =
				    V3[counter][i][P - 1] * BoundaryScale;
			}
		for (counter = 0; counter < M; counter++)
			for (j = P - 2; j >= 0; j--)
				for (i = 0; i < N; i++) {
					V1[counter][i][j] =
					    V1[counter][i][j] +
					    V1[counter][i][j + 1] * nu;
					V2[counter][i][j] =
					    V2[counter][i][j] +
					    V2[counter][i][j + 1] * nu;
					V3[counter][i][j] =
					    V3[counter][i][j] +
					    V3[counter][i][j + 1] * nu;
				}

		//filter down
		for (counter = 0; counter < M; counter++)
			for (j = 0; j < P; j++) {
				V1[counter][0][j] =
				    V1[counter][0][j] * BoundaryScale;
				V2[counter][0][j] =
				    V2[counter][0][j] * BoundaryScale;
				V3[counter][0][j] =
				    V3[counter][0][j] * BoundaryScale;
			}
		for (counter = 0; counter < M; counter++)
			for (i = 1; i < N; i++)
				for (j = 0; j < P; j++) {
					V1[counter][i][j] =
					    V1[counter][i][j] + V1[counter][i -
									    1]
					    [j] * nu;
					V2[counter][i][j] =
					    V2[counter][i][j] + V2[counter][i -
									    1]
					    [j] * nu;
					V3[counter][i][j] =
					    V3[counter][i][j] + V3[counter][i -
									    1]
					    [j] * nu;
				}

		//filter up 
		for (counter = 0; counter < M; counter++)
			for (j = 0; j < P; j++) {
				V1[counter][N - 1][j] =
				    V1[counter][N - 1][j] * BoundaryScale;
				V2[counter][N - 1][j] =
				    V2[counter][N - 1][j] * BoundaryScale;
				V3[counter][N - 1][j] =
				    V3[counter][N - 1][j] * BoundaryScale;
			}
		for (counter = 0; counter < M; counter++)
			for (i = N - 2; i >= 0; i--)
				for (j = 0; j < P; j++) {
					V1[counter][i][j] =
					    V1[counter][i][j] + V1[counter][i +
									    1]
					    [j] * nu;
					V2[counter][i][j] =
					    V2[counter][i][j] + V2[counter][i +
									    1]
					    [j] * nu;
					V3[counter][i][j] =
					    V3[counter][i][j] + V3[counter][i +
									    1]
					    [j] * nu;
				}

		//filter in 
		for (i = 0; i < N; i++)
			for (j = 0; j < P; j++) {
				V1[0][i][j] = V1[0][i][j] * BoundaryScale;
				V2[0][i][j] = V2[0][i][j] * BoundaryScale;
				V3[0][i][j] = V3[0][i][j] * BoundaryScale;
			}
		for (counter = 1; counter < M; counter++)
			for (i = 0; i < N; i++)
				for (j = 0; j < P; j++) {
					V1[counter][i][j] =
					    V1[counter][i][j] + V1[counter -
								   1][i][j] *
					    nu;
					V2[counter][i][j] =
					    V2[counter][i][j] + V2[counter -
								   1][i][j] *
					    nu;
					V3[counter][i][j] =
					    V3[counter][i][j] + V3[counter -
								   1][i][j] *
					    nu;
				}

		//filter out 
		for (i = 0; i < N; i++)
			for (j = 0; j < P; j++) {
				V1[M - 1][i][j] =
				    V1[M - 1][i][j] * BoundaryScale;
				V2[M - 1][i][j] =
				    V2[M - 1][i][j] * BoundaryScale;
				V3[M - 1][i][j] =
				    V3[M - 1][i][j] * BoundaryScale;
			}
		for (counter = M - 2; counter >= 0; counter--)
			for (i = 0; i < N; i++)
				for (j = 0; j < P; j++) {
					V1[counter][i][j] =
					    V1[counter][i][j] + V1[counter +
								   1][i][j] *
					    nu;
					V2[counter][i][j] =
					    V2[counter][i][j] + V2[counter +
								   1][i][j] *
					    nu;
					V3[counter][i][j] =
					    V3[counter][i][j] + V3[counter +
								   1][i][j] *
					    nu;
				}
	}
	for (counter = 0; counter < M; counter++)
		for (i = 0; i < N; i++)
			for (j = 0; j < P; j++) {
				V1[counter][i][j] *= PostScale;
				V2[counter][i][j] *= PostScale;
				V3[counter][i][j] *= PostScale;
			}

}

void interp(float U1[M][N][P], float U2[M][N][P], float U3[M][N][P],
	    float interpT[M][N][P], float T[M][N][P])
{
	int i, j, k;
	for (i = 1; i < M - 1; i++)
		for (j = 1; j < N - 1; j++)
			for (k = 1; k < P - 1; k++) {
				float itmp_f = i - U1[i][j][k];
				int itmp = (int)itmp_f;
				float jtmp_f = j - U2[i][j][k];
				int jtmp = (int)jtmp_f;
				float ktmp_f = k - U3[i][j][k];
				int ktmp = (int)ktmp_f;
				interpT[i][j][k] =
				    (1.0f + itmp - itmp_f) * (1.0f + jtmp -
							      jtmp_f) * (1.0f +
									 ktmp -
									 ktmp_f)
				    * T[itmp][jtmp][ktmp] + (1.0f + itmp -
							     itmp_f) * (1.0f +
									jtmp -
									jtmp_f)
				    * (ktmp_f - ktmp) * T[itmp][jtmp][ktmp +
								      1] +
				    (1.0f + itmp - itmp_f)
				    * (jtmp_f - jtmp) * (1.0f + ktmp -
							 ktmp_f) *
				    T[itmp][jtmp + 1][ktmp] + (1.0f + itmp -
							       itmp_f) *
				    (jtmp_f - jtmp) * (ktmp_f -
						       ktmp) * T[itmp][jtmp +
								       1][ktmp +
									  1] +
				    (itmp_f - itmp) * (1.0f + jtmp - jtmp_f)
				    * (1.0f + ktmp - ktmp_f) * T[itmp +
								 1][jtmp][ktmp]
				    + (itmp_f - itmp) * (1.0f + jtmp -
							 jtmp_f) * (ktmp_f -
								    ktmp) *
				    T[itmp + 1][jtmp][ktmp + 1] + (itmp_f -
								   itmp) *
				    (jtmp_f - jtmp) * (1.0f + ktmp - ktmp_f)
				    * T[itmp + 1][jtmp + 1][ktmp] + (itmp_f -
								     itmp) *
				    (jtmp_f - jtmp) * (ktmp_f - ktmp) * T[itmp +
									  1]
				    [jtmp + 1][ktmp + 1];
			}
}

void computeMI(float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	       float V3[M + 1][N + 1][P + 1], float interpT[M][N][P],
	       float S[M][N][P])
{
	int p1[IMAX];
	int p2[IMAX];
	int p12[IMAX][IMAX];
	float L[IMAX][IMAX];
	int i, j, k;

	for (i = 0; i < IMAX; i++)
		for (j = 0; j < IMAX; j++)
			p12[i][j] = 0;

	for (i = 0; i < IMAX; i++) {
		p1[i] = 0;
		p2[i] = 0;
	}
	for (i = 1; i < M - 1; i++)
		for (j = 1; j < N - 1; j++)
			for (k = 1; k < P - 1; k++) {
				float s_f = (0.5f * S[i][j][k]);
				float t_f = (0.5f * interpT[i][j][k]);

				int s = (int)(s_f);
				int t = (int)(t_f);
				p12[s][t] += 1;
			}
	for (j = 0; j < IMAX; j++)
		for (i = 0; i < IMAX; i++)
			p1[i] += p12[i][j];

	for (i = 0; i < IMAX; i++)
		for (j = 0; j < IMAX; j++)
			p2[j] += p12[i][j];

	float MI = 0.0f;

	for (i = 0; i < IMAX; i++)
		for (j = 0; j < IMAX; j++) {
			L[i][j] =
			    (M * N * P * (float)(p12[i][j] + EPS) /
			     ((float)(p1[i] + EPS) * (float)(p2[j] + EPS)));
		}

	for (i = 0; i < IMAX; i++)
		for (j = 0; j < IMAX; j++) {
			L[i][j] = 1.0f + logf(L[i][j]);
			MI += p12[i][j] * (L[i][j] - 1.0f);

		}

	convolution2D(L);

	for (i = 1; i < M - 1; i++)
		for (j = 1; j < N - 1; j++)
			for (k = 1; k < P - 1; k++) {
				float dTdx =
				    (interpT[i + 1][j][k] -
				     interpT[i - 1][j][k]) * 0.5f;
				float dTdy =
				    (interpT[i][j + 1][k] -
				     interpT[i][j - 1][k]) * 0.5f;
				float dTdz =
				    (interpT[i][j][k + 1] -
				     interpT[i][j][k - 1]) * 0.5f;

				float s = S[i][j][k] * 0.5f;
				float t = interpT[i][j][k] * 0.5f;
				int s_int = (int)s;
				int t_int = (int)t;
				float temp = L[s_int][t_int];
				V1[i][j][k] = -1.0f * temp * dTdx;
				V2[i][j][k] = -1.0f * temp * dTdy;
				V3[i][j][k] = -1.0f * temp * dTdz;
			}

	printf("MI is %f", MI / (M * N * P));
}

void U_cal(float U1[M][N][P], float U2[M][N][P], float U3[M][N][P],
	   float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	   float V3[M + 1][N + 1][P + 1])
{
	int i, j, k;
	float Rmax_L2 = 0.0;
#define u1(i,j,k) U1[i][j][k]
#define u2(i,j,k) U2[i][j][k]
#define u3(i,j,k) U3[i][j][k]
#define v1(i,j,k) V1[i][j][k]
#define v2(i,j,k) V2[i][j][k]
#define v3(i,j,k) V3[i][j][k]

	for (i = 1; i < M - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			for (k = 1; k < P - 1; k++) {
				float du1_dx =
				    (u1(i + 1, j, k) - u1(i - 1, j, k)) * 0.5f;
				float du2_dx =
				    (u2(i + 1, j, k) - u2(i - 1, j, k)) * 0.5f;
				float du3_dx =
				    (u3(i + 1, j, k) - u3(i - 1, j, k)) * 0.5f;

				float du1_dy =
				    (u1(i, j + 1, k) - u1(i, j - 1, k)) * 0.5f;
				float du2_dy =
				    (u2(i, j + 1, k) - u2(i, j - 1, k)) * 0.5f;
				float du3_dy =
				    (u3(i, j + 1, k) - u3(i, j - 1, k)) * 0.5f;

				float du1_dz =
				    (u1(i, j, k + 1) - u1(i, j, k - 1)) * 0.5f;
				float du2_dz =
				    (u2(i, j, k + 1) - u2(i, j, k - 1)) * 0.5f;
				float du3_dz =
				    (u3(i, j, k + 1) - u3(i, j, k - 1)) * 0.5f;

				float r1 =
				    v1(i, j, k) - v1(i, j, k) * du1_dx - v2(i,
									    j,
									    k) *
				    du1_dy - v3(i, j, k) * du1_dz;
				float r2 =
				    v2(i, j, k) - v1(i, j, k) * du2_dx - v2(i,
									    j,
									    k) *
				    du2_dy - v3(i, j, k) * du2_dz;
				float r3 =
				    v3(i, j, k) - v1(i, j, k) * du3_dx - v2(i,
									    j,
									    k) *
				    du3_dy - v3(i, j, k) * du3_dz;

				v1(i, j, k) = r1;
				v2(i, j, k) = r2;
				v3(i, j, k) = r3;

				float R_L2 = sqrtf(r1 * r1 + r2 * r2 + r3 * r3);
				Rmax_L2 = (Rmax_L2 > R_L2) ? Rmax_L2 : R_L2;
			}
		}
	}
	float dt = 0.1 / Rmax_L2;
	printf("dt = %f\n", dt);

	for (i = 1; i < M - 1; i++)
		for (j = 1; j < N - 1; j++)
			for (k = 1; k < P - 1; k++) {
				u1(i, j, k) += dt * v1(i, j, k);
				u2(i, j, k) += dt * v2(i, j, k);
				u3(i, j, k) += dt * v3(i, j, k);
			}
#undef u1
#undef u2
#undef u3
#undef v1
#undef v2
#undef v3
}
