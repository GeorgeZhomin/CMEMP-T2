#pragma once

enum class SOLVER_MODE{DSD, FK}; //несимметричная производная, фиктивный узел

struct Bound
{
    double x, λ, μ, ψ;
};

struct Diffeq2
{
    double (*p)(double);
    double (*q)(double);
    double (*g)(double);
};

//n - количество сегментов
double* solve_bvp(int n, const Bound &left_bound, const Bound &right_bound, const Diffeq2 &eq_funcs, SOLVER_MODE mode = SOLVER_MODE::DSD);

//интегрирование по Симпсону с n сегментами на сетке с шагом h
double simpson(int n, double h, double *u);