// The Weierstrass function
#include <assert.h>
#include <math.h>
#include <stdint.h>
#ifndef PI
#define PI 3.1415926
#endif

double weierstrass(double x, double a, double b, unsigned long n)
{
    // assert(a * b > (1 + 3 / 2 * PI) && "Invalid a and b");
    double res = 0.0;
    for(unsigned long i = 0;i < n;i++)
    {
        res += pow(a, n) * cos(pow(b, n) * PI * x);
    }
    return res;
}