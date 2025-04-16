#include <stdlib.h>
#include "rk4.h"

void runge_kutta(int n, int nt, int m, double h, void (*f)(double, double*, double*), double* x0, double* res)
{
	const int NEWTON_NMAX = 100000;
    int i, j;
	double t;
	double* k1;
	double* k2;
	double* k3;
	double* k4;
	double* xcur;
	double* fres;
	double* farg;
	double* xlast;
	k1 = (double*)malloc(n * sizeof(double));
	k2 = (double*)malloc(n * sizeof(double));
	k3 = (double*)malloc(n * sizeof(double));
	k4 = (double*)malloc(n * sizeof(double));
	farg = (double*)malloc(n * sizeof(double));
	fres = (double*)malloc(n * sizeof(double));
	xcur = (double*)malloc(n * sizeof(double));
	xlast = (double*)malloc(n * sizeof(double));
	for(j = 0; j < n; j++)
	{
		xcur[j] = x0[j];
	}
 	for(i = 0; i < nt; i++)
	{
		t = i * h;
		f(t, xcur, fres);
		for(j = 0; j < n; j++)
		{
			k1[j] = h * fres[j];
		}
		for(j = 0; j < n; j++)
		{
			farg[j] = xcur[j] + 0.5 * k1[j];
		}
		f(t + 0.5 * h, farg, fres);
		for(j = 0; j < n; j++)
		{
			k2[j] = h * fres[j];
		}
		for(j = 0; j < n; j++)
		{
			farg[j] = xcur[j] + 0.5 * k2[j];
		}
		f(t + 0.5 * h, farg, fres);
		for(j = 0; j < n; j++)
		{
			k3[j] = h * fres[j];
		}
		for(j = 0; j < n; j++)
		{
			farg[j] = k3[j] + xcur[j];
		}
		f(t + h, farg, fres);
		for(j = 0; j < n; j++)
		{
			k4[j] = h * fres[j];
		}
		for(j = 0; j < n; j++)
		{
			xlast[j] = xcur[j];
		}
		for(j = 0; j < n; j++)
		{
			xcur[j] = xlast[j] + (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
		}
		if(i % m == 0)
		{
			for(j = 0; j < n; j++)
			{
				res[(i / m) * n + j] = xlast[j];
			}
		}
	}
	free(xcur);
	free(xlast);
	free(farg);
	free(fres);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
}
