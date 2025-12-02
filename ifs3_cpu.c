#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4
#define ITERATIONS 40000

#define USE_SIMD 0
#define USE_OMP 0

#if !USE_OMP
#define SEED 0
#else
#define SEED omp_get_thread_num()
#endif

// TYPE DEFINITIONS

typedef struct
{
    double a[N];
    double b[N];
    double c[N];
    double d[N];
    double e[N];
    double f[N];
    double p[N];       // Probability for each function.
    double prob_T[N];  // Probability table for Alias Method.
    double alias_T[N]; // Alias table for Alias Method.
} ifs;

typedef struct
{
    double x;
    double y;
} pos;

// RANDOM STUFF

static unsigned int g_seed;

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline __attribute__((__always_inline__)) int fast_rand(void)
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

inline __attribute__((__always_inline__)) double fast_rand_d(void)
{
    return (double)fast_rand() / 0x7FFF;
}

// DID NOT WRITE THIS!
void alias_preprocess(ifs *system)
{
    double scaled[N];
    int small[N];
    int large[N];
    int small_top = 0, large_top = 0;

    // Scale probabilities by n
    for (int i = 0; i < N; i++)
    {
        scaled[i] = system->p[i] * N;
        if (scaled[i] < 1.0)
            small[small_top++] = i;
        else
            large[large_top++] = i;
    }

    // Work through small and large lists
    while (small_top > 0 && large_top > 0)
    {
        int s = small[--small_top];
        int l = large[--large_top];

        system->prob_T[s] = scaled[s];
        system->alias_T[s] = l;

        scaled[l] = (scaled[l] + scaled[s]) - 1.0;

        if (scaled[l] < 1.0)
            small[small_top++] = l;
        else
            large[large_top++] = l;
    }

    // Remaining entries have probability 1
    while (large_top > 0)
    {
        int l = large[--large_top];
        system->prob_T[l] = 1.0;
        system->alias_T[l] = l;
    }
    while (small_top > 0)
    {
        int s = small[--small_top];
        system->prob_T[s] = 1.0;
        system->alias_T[s] = s;
    }
}

// NEITHER THIS
inline __attribute__((__always_inline__)) int alias_sample(ifs system)
{
    int i = fast_rand() % N;
    double u = fast_rand_d();
    return (u < system.prob_T[i]) ? i : system.alias_T[i];
}

// SETTERS

void set_function(ifs *system, int index, double a, double b, double c, double d, double e, double f, double p)
{
    system->a[index] = a;
    system->b[index] = b;
    system->c[index] = c;
    system->d[index] = d;
    system->e[index] = e;
    system->f[index] = f;
    system->p[index] = p;
}

void set_random_ifs(ifs *system)
{
    set_function(system, 0, fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d());

    for (int i = 1; i < N; i++)
        set_function(system, i, fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d(), fast_rand_d(), system->p[i - 1] + fast_rand_d());
}

void set_ifs_example_1(ifs *system)
{
    set_function(system, 0, 0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.01);
    set_function(system, 1, 0.85, 0.04, -0.04, 0.85, 0.0, 1.6, 0.85);
    set_function(system, 2, 0.2, -0.26, 0.23, 0.11, 0.0, 1.6, 0.07);
    set_function(system, 3, -0.15, 0.28, 0.26, 0.24, 0.0, 0.44, 0.07);
}

// HELPERS

void print_ifs(ifs system)
{
    for (int i = 0; i < N; i++)
        printf("f%d:\ta:%.2f\tb:%.2f\tc:%.2f\td:%.2f\te:%.2f\tf:%.2f\tp:%.2f\n", i + 1, system.a[i], system.b[i], system.c[i], system.d[i], system.e[i], system.f[i], system.p[i]);
}

// CRITICAL STUFF

inline __attribute__((__always_inline__)) pos func(pos p, double a, double b, double c, double d, double e, double f)
{
    return (pos){a * p.x + b * p.y + e, c * p.x + d * p.y + f};
}

inline __attribute__((__always_inline__)) pos next_iteration(pos p, ifs system)
{
    int f = alias_sample(system);
    return func(p, system.a[f], system.b[f], system.c[f], system.d[f], system.e[f], system.f[f]);
}

int main(int argc, char *argv[])
{
    // Set seed.
    g_seed = SEED;

    ifs *system = malloc(sizeof(ifs));
    set_ifs_example_1(system);
    // print_ifs(*system);

    // OPEN OUTPUT FILE
    FILE *output = fopen("ifs3.out", "w");

    // SETUP ALIAS TABLES
    alias_preprocess(system);

    pos p = (pos){0.0, 0.0};

    for (int i = 0; i < ITERATIONS; i++)
    {
        p = next_iteration(p, *system);
        fprintf(output, "%f %f\n", p.x, p.y);
    }

    fclose(output);

    return EXIT_SUCCESS;
}