#ifndef _REGISTRATION_H
#define _REGISTRATION_H

#define EPS 1e-6f

#define M 256
#define N 256
#define P 256

#define GAUSSIAN_STEP 3
#define IMAX 128

void V_cal(float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	   float V3[M + 1][N + 1][P + 1]);

void U_cal(float U1[M][N][P], float U2[M][N][P], float U3[M][N][P],
	   float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	   float V3[M + 1][N + 1][P + 1]);

void interp(float U1[M][N][P], float U2[M][N][P], float U3[M][N][P],
	    float interpT[M][N][P]);

void computeMI(float V1[M + 1][N + 1][P + 1], float V2[M + 1][N + 1][P + 1],
	       float V3[M + 1][N + 1][P + 1], float interpT[M][N][P],
	       float S[M][N][P]);

#endif
