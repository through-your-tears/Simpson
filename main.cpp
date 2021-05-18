#include <iostream>
#include <omp.h>
#include <math.h>
#define N 10000

using namespace std;
double f(double x) {
    return sin(x) + cos(x);
}

double Simpson(double a, double b, int n) {
    const double width = (b-a)/n;
    double simpson_integral = 0;
#pragma omp parallel
    {
        double local_sum = 0;
#pragma omp for
        for (int step = 0; step < n; step++) {
            const double x1 = a + step * width;
            const double x2 = a + (step + 1) * width;

            local_sum += (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
        }
#pragma omp atomic
        simpson_integral += local_sum;
    };

    return simpson_integral;
}


int main()
{
    double a = 0;
    const double pi = M_PI;
    double b = pi;
    double s = 0;
    cout << Simpson(a, b, N);
    return 0;
}
