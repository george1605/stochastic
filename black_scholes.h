// The Black Scholes formula for pricing options
// One of the well known cases of SDEs in action
#include <math.h>

double black_scholes_solution(double S, double(*N)(double), double T, double K, double r, double sigma)
{
    double d1 = (log(S / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    return (N(d1) * S - N(d2) * K * exp(-r * T));
}