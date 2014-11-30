#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
# define pi 3.14159265359
int p = 2667, lo = 6396, sh = 3103;
double mu = (double)sh*sh*p, lambda = (double)lo*lo*p - 2 * mu, freq = 2.25*1e6, wavelength = (double)lo / freq, dh = (double)wavelength / 30.0, dt = (double)dh / (lo*1.5);
int length = 20, breadth = 20;
int xno = (length / (1000 * dh)), yno = (breadth / (1000 * dh));
int xsource = 111, ysource = 0;
double time_total = pow(10.0, -6.0);
double timesteps = ceil(10e-6 / dt);
double dtbyrhodh = dt / (p*dh), lambdaplus2mudtbydh = (lambda + 2 * mu)*dt / dh, lambdadtbydh = lambda*dt / dh, dtmubydh = dt*mu / dh;
double** Make2DDoubleArray(int arraySizeX, int arraySizeY) {
	double** theArray;
	theArray = (double**)calloc(arraySizeX, sizeof(double*));
	for (int i = 0; i < arraySizeX; i++)
		theArray[i] = (double*)calloc(arraySizeY, sizeof(double));
	return theArray;
}
void velmodify(double**X, double**Y, double**Z, double**U, double**V){
	for (int i = 0; i <= xno - 1; i++)
	for (int l = 0; l <= yno; l++)
		U[i + 1][l] += dtbyrhodh*(X[i + 1][l + 1] - X[i + 1][l] + Z[i + 1][l] - Z[i][l]);
	for (int j = 0; j <= xno; j++)
	for (int k = 0; k <= yno - 1; k++)
		V[j][k + 1] += dtbyrhodh*(Z[j][k + 1] - Z[j][k] + Y[j + 1][k + 1] - Y[j][k + 1]);


}
void strmodify(double**X, double**Y, double**Z, double**U, double**V){
	for (int i = 0; i <= xno - 1; i++){
		for (int k = 0; k <= yno - 1; k++){
			X[i + 1][k + 1] += lambdaplus2mudtbydh*(U[i + 1][k + 1] - U[i + 1][k]) + lambdadtbydh*(V[i + 1][k + 1] - V[i][k + 1]);
			Y[i + 1][k + 1] += lambdaplus2mudtbydh*(V[i + 1][k + 1] - V[i][k + 1]) + lambdadtbydh*(U[i + 1][k + 1] - U[i + 1][k]);
		}
	}
	for (int j = 0; j <= xno; j++)
	for (int l = 0; l <= yno; l++)
		Z[j][l] += dtmubydh*(V[j][l + 1] - V[j][l] + U[j + 1][l] - U[j][l]);



}
int main(){
	clock_t start, end;
	double time_taken;
	start = clock();
	double **X, **Y, **U, **V, **Z, *A, x = 0.0, y = 0.0, result = 0.0; int n = 1;
	xno += 16 - (xno % 16);
	yno += 16 - (yno % 16);
	xno = xno - 2;
	yno = yno - 2;
	/* FILE *fp;
	fp=fopen("raj.dat","w"); */
	X = Make2DDoubleArray(xno + 2, yno + 2);
	Y = Make2DDoubleArray(xno + 2, yno + 2);
	Z = Make2DDoubleArray(xno + 1, yno + 1);
	U = Make2DDoubleArray(xno + 2, yno + 2);
	V = Make2DDoubleArray(xno + 2, yno + 2);
	A = (double*)calloc(timesteps, sizeof(double));
	for (n = 1; n <= timesteps; n++){
		velmodify(X, Y, Z, U, V);
		if ((double)(n*dt) <= (double)(3 / freq)){
			V[ysource][xsource] = (1 - cos(2 * pi*freq*dt*n / 3))*cos(2 * pi*freq*n*dt);
		}
		x += U[85][96] * dt;
		y += V[85][96] * dt;
		A[n - 1] = sqrt(x*x + y*y);
		strmodify(X, Y, Z, U, V);

	}
	for (int j = 0; j<timesteps; j++)
		printf("%d\t%e\n", j + 1, A[j]);
	end = clock();
	time_taken = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Time elapsed is %lf\nGRID Size:%d*%d\nTime Steps Taken:%lf\n", time_taken, xno + 2, yno + 2, timesteps);
	//fclose(fp);
	getchar();
	return 0;
}