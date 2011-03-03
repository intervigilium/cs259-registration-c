#ifndef _REGISTRATION_H
#define _REGISTRATION_H

#define EPS 1e-6f

#define M 256
#define N 256
#define P 256

#define GAUSSIAN_STEP 3
#define IMAX 128

void V_cal(float V1[][][], float V2[][][],
	   float V3[][][]);

void U_cal(float U1[][][], float U2[][][], float U3[][][],
	   float V1[][][], float V2[][][],
	   float V3[][][]);

void interp(float U1[][][], float U2[][][], float U3[][][],
	    float interpT[][][], float T[][][]);

void computeMI(float V1[][][], float V2[][][],
	       float V3[][][], float interpT[][][],
	       float S[][][]);

#endif
