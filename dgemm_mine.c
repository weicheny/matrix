#include <math.h>

const char* dgemm_desc = "My awesome dgemm.";

void square_dgemm(const int M, const double *A, const double *B, double *C)
{
	int j,k,i;
	for (j = 0; j < M; ++j) {
        for (k = 0; k < M; ++k) {
            for (i = 0; i < M; ++i) {
                C[j*M+i] += A[k*M+i] * B[j*M+k];
            }
        }
    }
/*						
   int i, j, k;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            double cij = C[j*M+i];
            for (k = 0; k < M; ++k)
                cij += A[k*M+i] * B[j*M+k];
            C[j*M+i] = cij;
        }
    }
  
int i, k, j;
    
    for (i = 0; i < M; ++i) {
        for (k = 0; k < M; ++k) {
            for (j = 0; j < M; ++j) {
                C[j*M+i] += A[k*M+i] * B[j*M+k];
            }
        }
    }
  
      int j, i, k;
    for (j = 0; j < M; ++j) {
        for (i = 0; i < M; ++i) {
            double cij = C[j*M+i];
            for (k = 0; k < M; ++k)
                cij += A[k*M+i] * B[j*M+k];
            C[j*M+i] = cij;
        }
    }	*/
}


