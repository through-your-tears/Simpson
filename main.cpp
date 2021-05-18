#include <iostream>
#include <omp.h>
#include <math.h>
#define N 10000

using namespace std;
double f(double x) {
    return sin(x) + cos(x);
}

double Simpson(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum1 = 0;
    double sum2 = 0;
    for (int k = 1; k <= n; k++)
    {
        double xk = a + k * h;
        if (k <= n - 1)
        {
            sum1 += f(xk);
        }

        double xk_1 = a + (k - 1) * h;
        sum2 += f((xk + xk_1) / 2);
    }

    double result = h / 3 * (1 / 2 * f(a) + sum1 + 2 * sum2 + 1 / 2 * f(b));
    return result;
}


int main()
{
    double a = 0;
    const double pi = M_PI;
    double b = pi;
    double s = 0;
    cout << Simpson(a, b, N);
    return 0;
/*#pragma omp parallel num_threads(4)
    {
        cout << "Parallel";
    }
    return 0;
*/
}
