#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include "papi.h"
#include "registration.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		printf("execute using ./a.out NumIter\r\n");
		exit(0);
	}

	float S[M][N][P];
	float T[M][N][P];

	float U1[M][N][P];
	float U2[M][N][P];
	float U3[M][N][P];

	float V1[M + 1][N + 1][P + 1];
	float V2[M + 1][N + 1][P + 1];
	float V3[M + 1][N + 1][P + 1];

	float interpT[M][N][P];

	int NumIter = atoi(argv[1]);

	int i, j, k, q;

	FILE *S_ptr;
	FILE *T_ptr;

	float Smax = 8091.0, Tmax = 8180.0;
	unsigned short num;

	S_ptr = fopen("data/S.img", "rb");
	if (S_ptr == NULL) {
		printf("reading error\n");
		exit(1);
	}

	T_ptr = fopen("data/T.img", "rb");
	if (T_ptr == NULL) {
		printf("reading error\n");
		exit(1);
	}

	float t;
	for (k = 0; k < P; k++) {
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {
				fread(&num, 2, 1, S_ptr);
				t = (float)num / Smax * 255.0f;
				S[i][j][k] = t;
			}
		}
	}

	fclose(S_ptr);

	for (k = 0; k < P; k++) {
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {
				fread(&num, 2, 1, T_ptr);
				t = (float)num / Tmax * 255.0f;
				T[i][j][k] = t;
			}
		}
	}
	fclose(T_ptr);

	float u0 = 0.0f;

	for (k = 0; k < P; k++) {
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {
				U1[i][j][k] = u0;
				U2[i][j][k] = u0;
				U3[i][j][k] = u0;
				V1[i][j][k] = u0;
				V2[i][j][k] = u0;
				V3[i][j][k] = u0;
				interpT[i][j][k] = u0;
			}
		}
	}

	printf("Begin \r\n");
	struct timeval a;
	struct timeval b;
	gettimeofday(&a, 0);
	unsigned g_batchid = 0;
	if (argc == 3) {
		g_batchid = atoi(argv[2]);
	}

	int Events[5];
	u_long_long papi_values[5];
	util_start_papi(g_batchid, Events);

	for (i = 0; i < NumIter; i++) {
		printf("iteration  %d\r\n", i);
		interp(U1, U2, U3, interpT, T);
		computeMI(V1, V2, V3, interpT, S);
		V_cal(V1, V2, V3);
		U_cal(U1, U2, U3, V1, V2, V3);
	}
	util_stop_papi(g_batchid, papi_values);
	util_print_papi(g_batchid, papi_values, (g_batchid == 0));

	gettimeofday(&b, 0);
	printf("End\r\n");
	printf("Seconds is %d\r\n", (b.tv_sec - a.tv_sec));
}
