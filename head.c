#include"head.h"

#define EPS 1e-14
#define MAX_ITER 10000

int rand_uniform(int n) 
{
    int max = RAND_MAX - RAND_MAX % n;
    int x;
    do
    {
        x = rand();
    }while ( x >= max );
    return x % n;
}

double M_Pi() 
{
    int inside_circle = 0;
    double x,y;
    for (int i = 0; i < SAMPLE_SIZE; i++) 
    {
        x = (rand() % 10000) / 10000.0;
        y = (rand() % 10000) / 10000.0;
        if (x * x + y * y <= 1.0) 
        {
            inside_circle++;
        }
    }
    return 4.0 * inside_circle / SAMPLE_SIZE;
}

double rand_normal(double mean, double stddev)
{
    static int hasSpare = 0;
    static double spare;
    double u1, u2, mag, z0;
    static double pi = 0.0;

    if (pi == 0.0)
        pi = M_Pi(SAMPLE_SIZE);

    if (hasSpare) 
    {
        hasSpare = 0;
        return spare * stddev + mean;
    }
    do 
    {
        u1 = (double)rand() / RAND_MAX;
    } while (u1 <= 1e-10);
    u2 = (double)rand() / RAND_MAX;

    mag = sqrt(-2.0 * log(u1));
    z0 = mag * cos(2 * pi * u2);
    spare = mag * sin(2 * pi * u2);
    hasSpare = 1;

    return z0 * stddev + mean;
}

double Chi_Square() 
{
    double expected = SAMPLE_SIZE / N;
    double chi_sq = 0.0;
    for (int i = 0; i < N; i++) 
    {
        double diff = histogram[i] - expected;
        chi_sq += diff * diff / expected;
    }
    return chi_sq;
}

double Entropy()
{
    double ent = 0.0;
    for (int i = 0; i < N; i++) 
    {
        if (histogram[i] > 0) 
        {
            double p = histogram[i] / (double)SAMPLE_SIZE;
            ent -= p * log2(p);
        }
    }
    return ent;
}

double Autocorrelation()
{
    double mean = 0.0;
    double numerator = 0.0;
    double denominator = 0.0;

    for (int i = 0; i < SAMPLE_SIZE; i++) 
    {
        mean += seq[i];
    }
    mean /= SAMPLE_SIZE;
    for (int i = 0; i < SAMPLE_SIZE - 1; i++) 
    {
        numerator += (seq[i] - mean) * (seq[i + 1] - mean);
        denominator += (seq[i] - mean) * (seq[i] - mean);
    }

    return numerator / denominator;
}

int Freedom_Degrees(int n) 
{
    if (n <= 0) return 0;
    return n - 1;
}

// 计算对数伽马函数 ln(Γ(x))
static double ln_gamma(double x) 
{
    static double coeffs[6] = 
    {
        76.18009172947146,   -86.50532032941677,
        24.01409824083091,   -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5
    };
    double y = x;
    double tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    for (int j=0; j<6; j++) 
    {
        y += 1.0;
        ser += coeffs[j]/y;
    }
    return -tmp + log(2.5066282746310005*ser/x);
}

// 级数展开计算P(a,x)
static double gamma_p_series(double a, double x) 
{
    double sum = 1.0 / a;
    double del = sum;
    double ap = a;

    for (int n = 1; n <= MAX_ITER; n++) 
    {
        ap += 1;
        del *= x / ap;
        sum += del;
        if (fabs(del) < fabs(sum) * EPS) break;
    }
    return sum * exp(-x + a * log(x) - ln_gamma(a));
}

// 连分式法计算Q(a,x)
static double gamma_q_contfrac(double a, double x) 
{
    double b = x + 1 - a;
    double c = 1.0 / 1.0e-30;
    double d = 1.0 / b;
    double h = d;

    for (int i = 1; i <= MAX_ITER; i++) 
    {
        double an = -i * (i - a);
        b += 2;
        d = an * d + b;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = b + an / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double delta = d * c;
        h *= delta;
        if (fabs(delta - 1.0) < EPS) break;
    }
    return h * exp(-x + a * log(x) - ln_gamma(a));
}

// 计算卡方检验的p值（右尾概率）
double Chi_Square_P_Value(double chi_square_stat, int df) 
{
    if (chi_square_stat < 0 || df < 1) return 0.0;
    double a = df / 2.0;
    double x = chi_square_stat / 2.0;
    if (x < a + 1) {
        return 1.0 - gamma_p_series(a, x);
    } else {
        return gamma_q_contfrac(a, x);
    }
}