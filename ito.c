#include "weierstrass.h"
#include <stdio.h>
#include <math.h>
#include <stdint.h>

double standard_normal() {
    static int have = 0;
    static double z1;
    if (have) { have = 0; return z1; }
    double u1, u2;
    do { u1 = (double)rand() / INT32_MAX; } while (u1 <= 1e-12);
    u2 = (double)rand() / INT32_MAX;
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * PI * u2;
    z1 = r * cos(theta);
    have = 1;
    return r * sin(theta);
}

// calculate it
double ito_integral(double (*f)(double, double), double a, double b, unsigned long n) {
    double dt = (b - a) / n;
    double Bt = 0.0; // B0 = 0
    double s = 0.0;
    for(unsigned long i = 0; i < n; i++) {
        double dB = sqrt(dt) * standard_normal(); 
        s += f(i * dt, Bt) * dB;
        Bt += dB;
    }
    return s;
}

// Calculates the variance of solution
double solution_var(double(*sigma)(double), double t, unsigned long n)
{
    double s = 0;
    double step = t / n;
    double x = 0;
    for(unsigned long i = 0;i < n;i++) {
        double v = sigma(x);
        s += step * (v * v);
        x += step;
    }
    return x;
}

// Calculates the expected value of solution
// E[X_t] = X_0 + integral{0, t} miu(s) ds
double solution_exp(double(*miu)(double), double t, unsigned long n, double x0)
{
    double s = x0;
    double step = t / n;
    double x = 0;
    for(unsigned long i = 0;i < n;i++) {
        s += step * miu(x);
        x += step;
    }
    return x;
}

double solve_gbm(double miu, double sigma, double t, double x0, double Bt)
{
    return x0 * exp((miu - 0.5 * sigma * sigma) * t + sigma * Bt);
}

double stochastic_weierstrass(double a, double b)
{
    return weierstrass(b, a, b, 100);
}

int main()
{
    printf("Ito Integral of Weierstrass: %lf", ito_integral(stochastic_weierstrass, 2, 4, 100));
    return 0;
}
