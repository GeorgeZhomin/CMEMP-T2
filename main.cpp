#include <iostream>
#include <fstream>
#include <cmath>
#include "linfuncs.h"
#include "difffuncs.h"
#include <algorithm>

int main()
{
    //Вариант 8
    const Bound left_bound{0, 1, 0, 0}, right_bound{1, 1, 2, 0};
    const Diffeq2 eq_funcs{[](double x){return x;}, [](double x){return -std::sqrt(x);}, [](double x){return -3*std::exp(-x);}};
    const double ns[] = {25, 50, 100, 200, 500, 1000, 2000, 4000, 8000};

    std::ofstream file_sol("output/solution.txt");
    std::ofstream file_norm("output/norm.txt");

    for (int n : ns)
    {
        //решение ДУ
        double *u = solve_bvp(n, left_bound, right_bound, eq_funcs, SOLVER_MODE::DSD);
        double *v = solve_bvp(n, left_bound, right_bound, eq_funcs, SOLVER_MODE::FK);

        //запись решения в файл
        double x = left_bound.x, h = (right_bound.x-left_bound.x)/n;
        for (int i = 0; i < n+1; i++)
        {
            file_sol << n << '\t' << x << '\t' << u[i] << '\t' << v[i] << '\n';
            x += h;
        }

        //вычисление норм
        double *difference = new double[n+1];
        for (int i = 0; i < n+1; i++)
        {
            difference[i] = std::abs(u[i] - v[i]);
        }
        delete[] u; delete[] v;

        double *squared_difference = map(difference, n+1, [](double x){return x*x;});
        file_norm << n << '\t' << *std::max_element(difference, difference+(n+1)) << '\t' << simpson(n, h, difference) << '\t' << std::sqrt(simpson(n, h, squared_difference)) << '\n';

        delete[] difference; delete[] squared_difference;
    }
    return 0;
}