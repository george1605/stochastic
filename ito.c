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

double stratonovich_integral(double (*f)(double, double), double a, double b, unsigned long n) {
    double dt = (b - a) / n;
    double Bt = 0.0; // B0 = 0
    double s = 0.0;
    for(unsigned long i = 0; i < n; i++) {
        double dB = sqrt(dt) * standard_normal(); 
        s += f((i + 0.5) * dt, Bt) * dB;
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

double quadratic_variation(double *Y, int n) {
    double qv = 0.0;
    for (int i = 0; i < n-1; i++) {
        double diff = Y[i + 1] - Y[i];
        qv += diff * diff;
    }
    return qv;
}

double covariance(double* X, double *Y, int n) {
    double cov = 0.0;
    for (int i = 0; i < n-1; i++) {
        double diff1 = Y[i + 1] - Y[i];
        double diff2 = X[i + 1] - X[i];
        cov += diff1 * diff2;
    }
    return cov;
}

double stochastic_exp(double* Y, int n) {
    return exp(Y[n - 1] - 0.5 * quadratic_variation(Y, n));
}

// The Euler Maruyama SDE Solver
// Does the X[i] = X[i - 1] + dX_t
// Where dX_t = miu * dt + sigma * dBt
void euler_maruyama(double* X, int n, double dt, double(*miu)(double, double), double(*sigma)(double, double)) {
    double t = 0; // X[0] = initial condition, at t0 = 0
    for(int i = 1;i < n;i++)
    {
        double dBt = sqrt(dt) * standard_normal();
        X[i] = X[i-1] + miu(X[i-1], t) * dt + sigma(X[i-1], t) * dBt;
        t += dt;
    }
}

int main()
{
    double Y[] = { 0.2, 1.1, 0.9, 2.4, 3.3 };
    printf("Stochastic Exp: %lf", stochastic_exp(Y, 5));
    return 0;
}
