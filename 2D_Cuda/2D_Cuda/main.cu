#include<stdio.h>
#include<math.h>
#include<time.h>
#include<iostream>
using namespace std;
#include"cuda.h"
#include"cuda_runtime.h"
#include"device_launch_parameters.h"
#define pi 3.14159265359
double* Make2DDoubleArray(int arraySizeX, int arraySizeY) {
	double* theArray;
	theArray = (double*)calloc(arraySizeX*arraySizeY, sizeof(double));
	return theArray;
}
double* Makecudaarray(int arraySizeX, int arraySizeY){
	double* Array;
	cudaMalloc((void**)&Array, arraySizeX*arraySizeY*sizeof(double));
	return Array;
}

__global__ void velmodify(double *d_X, double *d_Y, double *d_Z, double *d_U, double *d_V, int xno, int yno, double dtbyrhodh, double freq, double dt, int xsource, int ysource, int m, double *d_A, double *d_x, double *d_y){
	int i = blockIdx.x;
	int j = threadIdx.x;
	if (i <= xno - 1 && j <= yno)
		d_U[(i + 1)*(yno + 2) + j] += dtbyrhodh*(d_X[((i + 1)*(yno + 2)) + j + 1] - d_X[((i + 1)*(yno + 2)) + j] + d_Z[((i + 1)*(yno + 1)) + j] - d_Z[(i*(yno + 1)) + j]);
	if (i <= xno  && j <= yno - 1)
		d_V[(i*(yno + 2)) + j + 1] += dtbyrhodh*(d_Z[(i*(yno + 1)) + j + 1] - d_Z[(i*(yno + 1)) + j] + d_Y[((i + 1)*(yno + 2)) + j + 1] - d_Y[(i*(yno + 2)) + j + 1]);
	__syncthreads();
	if ((m*dt) <= (double)3 / freq){
		d_V[(ysource*(yno + 2)) + xsource] = (1 - cosf(2 * pi*freq*dt*m / 3))*cosf(2 * pi*freq*m*dt);

	}
	__syncthreads();
	if (i == 0 && j == 0){
		*d_x += (d_U[(85 * (yno + 2)) + 96] * dt);
		*d_y += (d_V[(85 * (yno + 2)) + 96] * dt);
		d_A[m - 1] = sqrt(d_x[0] * d_x[0] + d_y[0] * d_y[0]);
	}
}

__global__ void strmodify(double *d_X, double *d_Y, double *d_Z, double *d_U, double *d_V, int xno, int yno, double lambdaplus2mudtbydh, double lambdadtbydh, double dtmubydh){
	int i = blockIdx.x;
	int j = threadIdx.x;
	if (i <= xno - 1 && j <= yno - 1){
		d_X[((i + 1)*(yno + 2)) + j + 1] += lambdaplus2mudtbydh*(d_U[((i + 1)*(yno + 2)) + j + 1] - d_U[((i + 1)*(yno + 2)) + j]) + lambdadtbydh*(d_V[((i + 1)*(yno + 2)) + j + 1] - d_V[(i*(yno + 2)) + j + 1]);
		d_Y[((i + 1)*(yno + 2)) + j + 1] += lambdaplus2mudtbydh*(d_V[((i + 1)*(yno + 2)) + j + 1] - d_V[(i*(yno + 2)) + j + 1]) + lambdadtbydh*(d_U[((i + 1)*(yno + 2)) + j + 1] - d_U[((i + 1)*(yno + 2)) + j]);

		__syncthreads();

	}
	if (i <= xno && j <= yno)
		d_Z[(i*(yno + 1)) + j] += dtmubydh*(d_V[(i*(yno + 2)) + j + 1] - d_V[(i*(yno + 2)) + j] + d_U[((i + 1)*(yno + 2)) + j] - d_U[(i*(yno + 2)) + j]);
	__syncthreads();

}

int main(){
	clock_t start1, end1, start2, end2;
	double time_taken;
	start1 = clock();
	int p = 2667, lo = 6396, sh = 3103;
	double mu = (double)sh*sh*p, lambda = (double)lo*lo*p - 2 * mu, freq = 2.25*1e6, wavelength = (double)lo / freq, dh = (double)wavelength / 30.0, dt = (double)dh / (lo*1.5);
	int length = 20, breadth = 20, m;
	int xno = floor(length / (1000.0*dh)), yno = floor(breadth / (1000.0*dh));
	xno = xno + 16 - (xno % 16);
	yno = yno + 16 - (yno % 16);
	xno = xno - 2;
	yno = yno - 2;
	int xsource = ((xno + 2) / 2) - 1, ysource = 0;
	double time_total = pow(10.0, -6.0);
	double timesteps = ceil(10e-6 / dt);
	double dtbyrhodh = dt / (p*dh), lambdaplus2mudtbydh = (lambda + 2 * mu)*dt / dh, lambdadtbydh = lambda*dt / dh, dtmubydh = dt*mu / dh;
	double *X, *Y, *U, *V, *Z, *A, *x, *y;
	X = Make2DDoubleArray(xno + 2, yno + 2);
	Y = Make2DDoubleArray(xno + 2, yno + 2);
	Z = Make2DDoubleArray(xno + 1, yno + 1);
	U = Make2DDoubleArray(xno + 2, yno + 2);
	V = Make2DDoubleArray(xno + 2, yno + 2);
	A = Make2DDoubleArray(1, timesteps);
	x = Make2DDoubleArray(1, 1);
	y = Make2DDoubleArray(1, 1);
	double* d_X = Makecudaarray(xno + 2, yno + 2);
	double* d_Y = Makecudaarray(xno + 2, yno + 2);
	double* d_Z = Makecudaarray(xno + 1, yno + 1);
	double* d_U = Makecudaarray(xno + 2, yno + 2);
	double* d_V = Makecudaarray(xno + 2, yno + 2);
	double* d_A = Makecudaarray(1, timesteps);
	double* d_x = Makecudaarray(1, 1);
	double* d_y = Makecudaarray(1, 1);
	cudaMemcpy(d_X, X, (xno + 2)*(yno + 2)* sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Y, Y, (xno + 2)*(yno + 2)* sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Z, Z, (xno + 1)*(yno + 1)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_U, U, (xno + 2)*(yno + 2)* sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_V, V, (xno + 2)*(yno + 2)* sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_A, A, (timesteps)* sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, y, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, sizeof(double), cudaMemcpyHostToDevice);
	dim3 threadsperblock(yno + 2);
	dim3 numblocks(xno + 2);
	end1 = clock();
	start2 = clock();
	for (m = 1; m <= timesteps; m++){

		velmodify << <numblocks, threadsperblock >> >(d_X, d_Y, d_Z, d_U, d_V, xno, yno, dtbyrhodh, freq, dt, xsource, ysource, m, d_A, d_x, d_y);

		strmodify << <numblocks, threadsperblock >> >(d_X, d_Y, d_Z, d_U, d_V, xno, yno, lambdaplus2mudtbydh, lambdadtbydh, dtmubydh);


	}
	end2 = clock();
	cudaDeviceSynchronize();
	cudaMemcpy(A, d_A, (timesteps)* sizeof(double), cudaMemcpyDeviceToHost);


	for (int i = 0; i<timesteps; i++)
		printf("%d\t%e\n", i + 1, A[i]);


	time_taken = (double)(end1 - start1) / CLOCKS_PER_SEC;
	printf("Time elapsed is %lfseconds\nGRID Size:%d*%d\nTime Steps Taken:%lf\nNo of blocks:%d\n", time_taken, xno + 2, yno + 2, timesteps, yno + 2);
	printf("Time elapsed per function in gpu is %lf secomds\n", (double)(end2 - start2) / CLOCKS_PER_SEC);
	free(X);
	free(Y);
	free(Z);
	free(U);
	free(V);
	cudaFree(d_X);
	cudaFree(d_Y);
	cudaFree(d_Z);
	cudaFree(d_U);
	cudaFree(d_V);
	cudaFree(d_A);
	cudaFree(d_x);
	cudaFree(d_y);
	getchar();
	return 0;
}
