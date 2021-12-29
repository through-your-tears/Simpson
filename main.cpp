#include <iostream>
#include <omp.h>
#include <math.h>
#include <time.h>
#define N 100000000

using namespace std;
double f(double x) {
    return 4*x;
}

double I(double a, int step, double width) {
    double x1 = a + step * width;
    double x2 = a + (step + 1) * width;
    return (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
}

double DI(int n, double a, int step, double width) {
    double s = 0.0;
    step *= n;
    width /= n;
    for(int i = 0; i < n; ++i) {
        double x1 = a + (step + i) * width;
        double x2 = a + (step + i + 1) * width;
        s += (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
    }
    return s;
}

double Simpson(double a, double b, int n) {
    double width = (b-a) / n;
    double simpson_integral = 0;
    double eps = 0.001;
    for (int step = 0; step < n; step++) {
        int counter = 1;
        while(I(a, step, width) - DI(counter * 2, a, step, width) > eps && counter <= 4) {
            counter++;
        }
        // cout << I(a, step, width) << " ------ " << DI(counter * 2, a, step, width) << "\n";
        simpson_integral += DI(counter * 2, a, step, width);
    }
    return simpson_integral;
}

//double ompSimpson(double a, double b, int n) {
//    const double width = (b-a)/n;
//    double simpson_integral = 0;
//#pragma omp parallel
//    {
//        double local_sum = 0;
//#pragma omp for
//        for (int step = 0; step < n; step++) {
//            const double x1 = a + step * width;
//            const double x2 = a + (step + 1) * width;
//
//            local_sum += (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
//        }
//#pragma omp atomic
//        simpson_integral += local_sum;
//    };
//
//    return simpson_integral;
//}


int main()
{
    double a;
    double b;
    int n;
    cin >> a >> b >> n;
    cout << Simpson(a, b, 1000);
    return 0;
}