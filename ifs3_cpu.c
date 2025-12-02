#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4
#define ITERATIONS 20000

#define USE_SIMD 0
#define USE_OMP 0

#if !USE_OMP
#define SEED 0
#else
#define SEED omp_get_thread_num()
#endif

typedef struct
{
    double a[N];
    double b[N];
    double c[N];
    double d[N];
    double e[N];
    double f[N];
    double p[N]; // Interval value not to be used for reading.
} ifs;

typedef struct
{
    double x;
    double y;
} pos;

inline __attribute__((__always_inline__)) double ran()
{
    return (double)rand() / RAND_MAX;
}

inline __attribute__((__always_inline__)) pos func(pos p, double a, double b, double c, double d, double e, double f)
{
    return (pos){a * p.x + b * p.y + e, c * p.x + d * p.y + f};
}

// There may be a faster algorithm for this.
inline __attribute__((__always_inline__)) int choose_function(ifs system)
{
    double r = ran();
}

inline __attribute__((__always_inline__)) pos next_iteration(pos p, ifs system)
{
    int f = choose_function(system);
    return func(p, system.a[f], system.b[f], system.c[f], system.d[f], system.e[f], system.f[f]);
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
    set_function(system, 0, ran(), ran(), ran(), ran(), ran(), ran(), ran());

    for (int i = 1; i < N; i++)
        set_function(system, i, ran(), ran(), ran(), ran(), ran(), ran(), system->p[i - 1] + ran());
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

int main(int argc, char *argv[])
{
    // Set seed.
    srand(SEED);

    ifs *system = malloc(sizeof(ifs));
    set_ifs_example_1(system);
    // print_ifs(*system);

    pos p = (pos){0.0, 0.0};

    for (int i = 0; i < ITERATIONS; i++)
    {
        p = next_iteration(p, *system);
    }

    return EXIT_SUCCESS;
}
