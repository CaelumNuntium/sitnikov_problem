#include <stdlib.h>
#include <math.h>

extern double e;

void sitnikov_eq(double t, double* x, double* res)
{
    const double DBL_EPS = 2.5e-16;
    double E, M, r;
    E = 0.0;
    M = t;
    while(fabs(M - E) > DBL_EPS)
    {
        E = M + e * sin(E);
    }
    r = 1.0 - e * cos(E);
    res[0] = x[1];
    res[1] = -x[0] / ((r * r + x[0] * x[0]) * sqrt(r * r + x[0] * x[0]));
}
