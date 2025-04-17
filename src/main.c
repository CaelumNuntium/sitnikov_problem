#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "legendre.h"
#include "collo.h"
#include "rk4.h"
#include "sitnikov.h"

#define CONFIG_NAME "config.ini"
#define N_PARAMS_GENERAL 6
#define N_PARAMS_PORTRAIT 3
#define N_PARAMS_POINCARE 6
#define N 2

int write_result(const char*, int, int, double, double*);
int write_poincare_map(const char*, int, int, int, int, double*);

double e;

int main()
{
    const char* general_parameters[N_PARAMS_GENERAL] = { "eccentricity", "periods", "integrator", "delta_t", "s", "poincare_map" };
    char* general_values[N_PARAMS_GENERAL];
    int general_is_req[N_PARAMS_GENERAL] = { 1, 1, 1, 1, 0, 0 };
    const char* portrait_parameters[N_PARAMS_PORTRAIT] = { "initial_z", "initial_z_dot", "m" };
    char* portrait_values[N_PARAMS_PORTRAIT];
    int portrait_is_req[N_PARAMS_PORTRAIT] = { 1, 1, 0 };
    const char* poincare_parameters[N_PARAMS_POINCARE] = { "z_min", "z_max", "n_z", "z_dot_min", "z_dot_max", "n_z_dot" };
    char* poincare_values[N_PARAMS_POINCARE];
    int poincare_is_req[N_PARAMS_POINCARE] = { 1, 1, 1, 1, 1, 1 };
    int i, integrator, nt, m, s, poincare, nm;
    int j, k, l, ll;
    double periods, h;
    double z_min, z_max, z_dot_min, z_dot_max, z_step, z_dot_step;
    double z0, dz0;
    int n_z, n_z_dot, steps_in_period, n_periods, n_all;
    double x0[N];
    double* c;
    double* res;
    double* poincare_map;
    double* poincare_map_loc;
    double* poincare_map_loc_loc;
    for(i = 0; i < N_PARAMS_GENERAL; i++)
    {
        general_values[i] = (char*)malloc(100 * sizeof(char));
    }
    read_config(CONFIG_NAME, N_PARAMS_GENERAL, general_parameters, general_values, general_is_req);
    sscanf(general_values[0], "%lf", &e);
    sscanf(general_values[1], "%lf", &periods);
    if(!strcmp(general_values[2], "collo"))
    {
        integrator = 2;
    }
    else if(!strcmp(general_values[2], "rk4"))
    {
        integrator = 1;
    }
    else
    {
        fprintf(stderr, "Error: Invalid value: integrator = %s\n", general_values[2]);
        return 1;
    }
    sscanf(general_values[3], "%lf", &h);
    if(integrator == 2)
    {
        if(!strcmp(general_values[4], "#ND"))
        {
            fprintf(stderr, "Error: s is required for collocation method!\n");
            return 2;
        }
        sscanf(general_values[4], "%d", &s);
    }
    if(!strcmp(general_values[5], "true"))
    {
        poincare = 1;
    }
    else if(!strcmp(general_values[5], "false") || !strcmp(general_values[5], "#ND"))
    {
        poincare = 0;
    }
    else
    {
        fprintf(stderr, "Error: Invalid value: poincare_map = %s\n", general_values[5]);
    }
    for(i = 0; i < N_PARAMS_GENERAL; i++)
    {
        free(general_values[i]);
    }
    if(poincare)
    {
        for(i = 0; i < N_PARAMS_POINCARE; i++)
        {
            poincare_values[i] = (char*)malloc(100 * sizeof(char));
        }
        read_config(CONFIG_NAME, N_PARAMS_POINCARE, poincare_parameters, poincare_values, poincare_is_req);
        sscanf(poincare_values[0], "%lf", &z_min);
        sscanf(poincare_values[1], "%lf", &z_max);
        sscanf(poincare_values[2], "%d", &n_z);
        sscanf(poincare_values[3], "%lf", &z_dot_min);
        sscanf(poincare_values[4], "%lf", &z_dot_max);
        sscanf(poincare_values[5], "%d", &n_z_dot);
        for(i = 0; i < N_PARAMS_POINCARE; i++)
        {
            free(poincare_values[i]);
        }
        steps_in_period = (int)(2.0 * M_PI / h) + 1;
        h = 2.0 * M_PI / steps_in_period;
        m = steps_in_period;
        n_periods = (int)periods;
        nt = n_periods * steps_in_period;
        nm = n_periods;
        if(n_z == 0)
        {
            z_step = 0.0;
        }
        else
        {
            z_step = (z_max - z_min) / n_z;
        }
        if(n_z_dot == 0)
        {
            z_dot_step = 0.0;
        }
        else
        {
            z_dot_step = (z_dot_max - z_dot_min) / n_z_dot;
        }
        if(integrator == 2)
        {
            c = (double*)malloc(s * sizeof(double));
            lobatto(s, c);
        }
        res = (double*)malloc(N * (nm + 1) * sizeof(double));
        poincare_map = (double*)malloc(N * (n_z + 1) * (n_z_dot + 1) * (nm - 1) * sizeof(double));
        printf("Progress: %2d %%", 0);
        n_all = (n_z + 1) * (n_z_dot + 1);
        for(j = 0; j <= n_z; j++)
        {
            x0[0] = z_min + z_step * j;
            poincare_map_loc = poincare_map + N * nm * (n_z_dot + 1) * j;
            for(k = 0; k <= n_z_dot; k++)
            {
                x0[1] = z_dot_min + z_dot_step * k;
                printf("\b\b\b\b%2d %%", (int)(((n_z_dot + 1.0) * j + k) / n_all * 100.0));
                switch(integrator)
                {
                    case 1:
                        runge_kutta(N, nt, m, h, sitnikov_eq, x0, res);
                        break;
                    case 2:
                        collo(N, nt, m, h, sitnikov_eq, x0, s, c, res);
                        break;
                }
                poincare_map_loc_loc = poincare_map_loc + nm * N * k;
                for(l = 1; l < nm; l++)
                {
                    for(ll = 0; ll < N; ll++)
                    {
                        poincare_map_loc_loc[(l - 1) * N + ll] = res[l * N + ll];
                    }
                }
            }
        }
        printf("\b\b\b\b100 %%\n");
        if(integrator == 2)
        {
            free(c);
        }
        free(res);
        write_poincare_map("poincare_map.dat", N, n_z, n_z_dot, nm - 1, poincare_map);
        free(poincare_map);
    }
    else
    {
        for(i = 0; i < N_PARAMS_PORTRAIT; i++)
        {
            portrait_values[i] = (char*)malloc(100 * sizeof(char));
        }
        read_config(CONFIG_NAME, N_PARAMS_PORTRAIT, portrait_parameters, portrait_values, portrait_is_req);
        sscanf(portrait_values[0], "%lf", &z0);
        sscanf(portrait_values[1], "%lf", &dz0);
        if(!strcmp(portrait_values[2], "#ND"))
        {
            m = 1;
        }
        else
        {
            sscanf(portrait_values[2], "%d", &m);
        }
        for(i = 0; i < N_PARAMS_PORTRAIT; i++)
        {
            free(portrait_values[i]);
        }
        nt = (int)(periods * 2.0 * M_PI / h);
        nm = nt / m;
        x0[0] = z0;
        x0[1] = dz0;
        res = (double*)malloc((nm + 1) * N * sizeof(double));
        switch(integrator)
        {
            case 2:
                c = (double*)malloc(s * sizeof(double));
                lobatto(s, c);
                collo(N, nt, m, h, sitnikov_eq, x0, s, c, res);
                free(c);
                break;
            case 1:
                runge_kutta(N, nt, m, h, sitnikov_eq, x0, res);
                break;
        }
        write_result("result.dat", N, nm, h * m, res);
        free(res);
    }
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

int write_poincare_map(const char* filename, int n, int n_x, int n_y, int n_periods, double* res)
{
    FILE* out;
    int i, j, k, l;
    double* map_loc;
    double* map_loc_loc;
    if(!(out = fopen(filename, "w")))
    {
        return 1;
    }
    for(i = 0; i <= n_x; i++)
    {
        map_loc = res + (n_y + 1) * n_periods * n * i;
        for(j = 0; j <= n_y; j++)
        {
            map_loc_loc = map_loc + n_periods * n * j;
            for(k = 0; k < n_periods; k++)
            {
                for(l = 0; l < n; l++)
                {
                    fprintf(out, "%.10f ", map_loc_loc[k * n + l]);
                }
                fprintf(out, "\n");
            }
        }
    }
    fclose(out);
}
