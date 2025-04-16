#include <stdlib.h>
#include "matmul.h"

void matrix_multiplication(int n, int m, int l, double* a, double* b, double* res)
{
	int i, j, k;
	double* acur;
	double* b_col;
	b_col = (double*)malloc(m * sizeof(double));
	//printf("n = %d; m = %d; l = %d\n");
	for(j = 0; j < l; j++)
	{
		for(k = 0; k < m; k++)
		{
			b_col[k] = b[k * l + j];
			//printf("b_col[%d] = b[%d * %d + %d] = b[%d]\n", k, k, l, j, k * l + j);
		}
		for(i = 0; i < n; i++)
		{
			acur = a + i * m;
			res[i * l + j] = cdot(m, acur, b_col);
		}
	}
	free(b_col);
}

double cdot(int n, double* a, double* b)
{
	int i;
	double res;
	res = 0.0;
	for(i = 0; i < n; i++)
	{
		res += a[i] * b[i];
	}
	return res;
}
