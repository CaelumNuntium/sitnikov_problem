#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "legendre.h"
#include "collo.h"
#include "rk4.h"
#include "sitnikov.h"

#define N_PARAMS 8
#define N 2

int write_result(const char*, int, int, double, double*);

double e;

int main()
{
    const char* params[N_PARAMS] = { "eccentricity", "initial_z", "initial_z_dot", "integrator", "periods", "delta_t", "m", "s" };
    int is_required[N_PARAMS] = { 1, 1, 1, 1, 1, 1, 0, 0 };
    char* values[N_PARAMS];
    double z0, dz0, n_periods, h;
    int i, integrator, m, nt, nm, s;
    double x0[2];
    double* res;
    double* c;
    for(i = 0; i < N_PARAMS; i++)
    {
        values[i] = (char*)malloc(100 * sizeof(char));
    }
    read_config("config.ini", N_PARAMS, params, values, is_required);
    sscanf(values[0], "%lf", &e);
    sscanf(values[1], "%lf", &z0);
    sscanf(values[2], "%lf", &dz0);
    if(!strcmp(values[3], "collo"))
    {
        integrator = 1;
    }
    else if(!strcmp(values[3], "rk4"))
    {
        integrator = 2;
    }
    else
    {
        fprintf(stderr, "Error: Invalid value: integrator = %s", values[3]);
        return 1;
    }
    sscanf(values[4], "%lf", &n_periods);
    sscanf(values[5], "%lf", &h);
    if(!strcmp(values[6], "#ND"))
    {
        m = 1;
    }
    else
    {
        sscanf(values[6], "%d", &m);
    }
    if(integrator == 1)
    {
        if(!strcmp(values[7], "#ND"))
        {
            fprintf(stderr, "Error: Parameter s is required for collocation method");
            return 2;
        }
        else
        {
            sscanf(values[7], "%d", &s);
        }
    }
    for(i = 0; i < N_PARAMS; i++)
    {
        free(values[i]);
    }
    nt = (int)(n_periods * 2.0 * M_PI / h);
    nm = nt / m;
    x0[0] = z0;
    x0[1] = dz0;
    res = (double*)malloc((nm + 1) * N * sizeof(double));
    switch(integrator)
    {
        case 1:
            c = (double*)malloc(s * sizeof(double));
            lobatto(s, c);
            collo(N, nt, m, h, sitnikov_eq, x0, s, c, res);
            free(c);
            break;
        case 2:
            runge_kutta(N, nt, m, h, sitnikov_eq, x0, res);
            break;
    }
    write_result("result.dat", N, nm, h * m, res);
    free(res);
    return 0;
}

int write_result(const char* filename, int n, int nt, double dt, double* res)
{
    FILE* out;
    int i, j;
    if(!(out = fopen(filename, "w")))
    {
        return 1;
    }
    for(i = 0; i < nt; i++)
    {
        fprintf(out, "%.10lf", dt * i);
        for(j = 0; j < n; j++)
        {
            fprintf(out, " %.16lf", res[i * n + j]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
}
