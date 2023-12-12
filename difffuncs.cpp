#include "difffuncs.h"
#include "linfuncs.h"

double* solve_bvp(int n, const Bound &left_bound, const Bound &right_bound, const Diffeq2 &eq_funcs, SOLVER_MODE mode)
{
    double h = (right_bound.x-left_bound.x)/n;
    short imzl = (left_bound.μ == 0), imzr = (right_bound.μ == 0);
    short n_eqs = n + 1 - (imzl + imzr);
    double *l = new double[n_eqs-1], *m = new double[n_eqs], *u = new double[n_eqs-1], *f = new double[n_eqs];

    double rh2 = 1/(h*h);
    double p_i;

    double* result = new double[n+1];
    if (imzl) {result[0] = left_bound.ψ/left_bound.λ;}
    if (imzr) {result[n] = right_bound.ψ/right_bound.λ;}

    //запись центральных уравнений
    double x = left_bound.x+(2-imzl)*h;
    for (int i = 1+imzl; i <= n-1-imzr; i++)
    {
        p_i = eq_funcs.p(x);
        l[i-1-imzl] = -(rh2 - p_i/2/h);
        m[i-imzl] = (2*rh2-eq_funcs.q(x));
        u[i-imzl] = -(rh2 + p_i/2/h);
        f[i-imzl] = -eq_funcs.g(x);
        x += h;
    }

    //запись первого уравнения
    if (imzl)
    {
        x = left_bound.x+h;
        p_i = eq_funcs.p(x);
        m[0] = (2*rh2-eq_funcs.q(x));
        u[0] = -(rh2 + p_i/2/h);
        f[0] = (rh2 - p_i/2/h)*left_bound.ψ/left_bound.λ - eq_funcs.g(x);
    }
    else
    {
        switch (mode)
        {
        case SOLVER_MODE::DSD:
            x = left_bound.x+h;
            p_i = eq_funcs.p(x);
            m[0] = left_bound.μ*rh2 + (left_bound.μ*p_i - left_bound.λ)/h - left_bound.λ*p_i/2;
            u[0] = -(left_bound.μ*rh2 + left_bound.μ*p_i/h + left_bound.μ*eq_funcs.q(x)/2);
            f[0] = -(1/h + p_i/2)*left_bound.ψ - left_bound.μ/2*eq_funcs.g(x);
            break;
        
        case SOLVER_MODE::FK:
            x = left_bound.x;
            p_i = eq_funcs.p(x);
            m[0] = left_bound.μ*rh2 - left_bound.λ/h + (left_bound.λ*p_i - left_bound.μ*eq_funcs.q(x))/2;
            u[0] = -left_bound.μ*rh2;
            f[0] = -(1/h - p_i/2)*left_bound.ψ - left_bound.μ/2*eq_funcs.g(x);
            break;
        }
    }

    //запись последнего уравнения
    if (imzr)
    {
        x = right_bound.x-h;
        p_i = eq_funcs.p(x);
        l[n_eqs-2] = -(rh2 - p_i/2/h);
        m[n_eqs-1] = (2*rh2-eq_funcs.q(x));
        f[n_eqs-1] = (rh2 + p_i/2/h)*right_bound.ψ/right_bound.λ - eq_funcs.g(x);
    }
    else
    {
        switch (mode)
        {
        case SOLVER_MODE::DSD:
            x = right_bound.x-h;
            p_i = eq_funcs.p(x);
            l[n_eqs-2] = -(right_bound.μ*rh2 - right_bound.μ*p_i/h + right_bound.μ*eq_funcs.q(x)/2);
            m[n_eqs-1] = right_bound.μ*rh2 - (right_bound.μ*p_i - right_bound.λ)/h - right_bound.λ*p_i/2;
            f[n_eqs-1] = (1/h - p_i/2)*right_bound.ψ - right_bound.μ/2*eq_funcs.g(x);
            break;
        
        case SOLVER_MODE::FK:
            x = right_bound.x;
            p_i = eq_funcs.p(x);
            l[n_eqs-2] = -right_bound.μ*rh2;
            m[n_eqs-1] = right_bound.μ*rh2 + right_bound.λ/h + (right_bound.λ*p_i - right_bound.μ*eq_funcs.q(x))/2;
            f[n_eqs-1] = (1/h + p_i/2)*right_bound.ψ - right_bound.μ/2*eq_funcs.g(x);
            break;
        }
    }

    double* sol = solve_tma(l, m, u, f, n_eqs);
    for (int i = imzl; i < n+1-imzr; i++)
    {
        result[i] = sol[i-imzl];
    }
    delete[] sol;
    delete[] l; delete[] m; delete[] u; delete[] f;

    return result;
}

double simpson(int n, double h, double *u)
{
    double s = 0;
    for (int i = 1; i < n; i+=2)
    {
        s += u[i-1] + 4*u[i] + u[i+1];
    }
    return s*h/3;
}