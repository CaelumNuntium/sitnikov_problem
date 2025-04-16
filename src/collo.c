#include <stdlib.h>
#include <math.h>
#include "legendre.h"
#include "collo.h"
#include "matmul.h"
#include "gauss.h"

void make_y_matrix(int, int, double*, double*, double*, double*);
void make_f_matrix(int, int, double, double, void (*)(double, double*, double*), double*, double*, double*);
double difference(int, int, double*, double*);

int collo(int n, int nt, int m, double h, void (*F)(double, double*, double*), double* x0, int s, double* c, double* res)
{
    const double DBL_EPSILON = 2.5e-16;
    const int NMAX = 1000;
    int i, ii, j, jj, k;
    double* x;
    double* y;
    double* f;
    double* a;
    double* a_last;
    double* ab;
    double* ct;
    double* ct_inv;
    double* b;
    double t_cur;
    double dif;
    ct = (double*)malloc(s * s * sizeof(double));
    b = (double*)malloc(s * s * sizeof(double));
    for(j = 0; j < s; j++)
    {
        for(jj = 1; jj <= s; jj++)
        {
            ct[(jj - 1) * s + j] = shifted_legendre_polynomial_derivative_value(jj, c[j]);
            b[(jj - 1) * s + j] = shifted_legendre_polynomial_value(jj, c[j]) - (jj % 2 ? -1 : 1);
        }
    }
    ct_inv = (double*)malloc(s * s * sizeof(double));
    invert(s, ct, ct_inv);
    free(ct);
    x = (double*)malloc(n * sizeof(double));
    for(j = 0; j < n; j++)
    {
        x[j] = x0[j];
    }
    y = (double*)malloc(n * s * sizeof(double));
    a = (double*)malloc(n * s * sizeof(double));
    a_last = (double*)malloc(n * s * sizeof(double));
    f = (double*)malloc(n * s * sizeof(double));
    for(i = 0; i < nt; i++)
    {
        t_cur = i * h;
        if(i % m == 0)
        {
            for(j = 0; j < n; j++)
            {
                res[i / m * n + j] = x[j];
            }
        }
        dif = 1.0;
        for(j = 0; j < n; j++)
        {
            for(jj = 0; jj < s; jj++)
            {
                a[j * s + jj] = 0.0;
            }
        }
        k = 0;
        while(dif > DBL_EPSILON)
        {
            make_y_matrix(n, s, x, a, b, y);
            make_f_matrix(n, s, h, t_cur, F, y, c, f);
            for(j = 0; j < n; j++)
            {
                for(jj = 0; jj < s; jj++)
                {
                    a_last[j * s + jj] = a[j * s + jj];
                }
            }
            matrix_multiplication(n, s, s, f, ct_inv, a);
            dif = difference(n, s, a, a_last);
            k++;
        }
        make_y_matrix(n, s, x, a, b, y);
        for(j = 0; j < n; j++)
        {
            x[j] = y[n * (s - 1) + j];
        }
    }
    free(ct_inv);
    free(x);
    free(a);
    free(a_last);
    free(b);
    free(y);
    free(f);
}

void make_y_matrix(int n, int s, double* y0, double* a, double* b, double* y)
{
    double* ab;
    int i, j;
    ab = (double*)malloc(n * s * sizeof(double));
    matrix_multiplication(n, s, s, a, b, ab);
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < s; j++)
        {
            y[j * n + i] = y0[i] + ab[i * s + j];
        }
    }
    free(ab);
}

void make_f_matrix(int n, int s, double h, double t0, void (*F)(double, double*, double*), double* y, double* c, double* f)
{
    int i, j;
    double* ftmp;
    ftmp = (double*)malloc(n * s * sizeof(double));
    for(j = 0; j < s; j++)
    {
        F(t0 + h * c[j], y + n * j, ftmp + n * j);
        for(i = 0; i < n; i++)
        {
            f[i * s + j] = ftmp[j * n + i] * h;
        }
    }
    free(ftmp);
}

double difference(int n, int s, double* a1, double* a2)
{
    int i, j;
    double res, cur_dif;
    res = 0.0;
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < s; j++)
        {
            cur_dif = fabs(a1[i * s + j] - a2[i * s + j]);
            if(cur_dif > res)
            {
                res = cur_dif;
            }
        }
    }
    return res;
}
