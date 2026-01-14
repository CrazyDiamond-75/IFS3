#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// CONSTANTS

#define N 4
#define ITERATIONS 100000000
#define GRADIENT 0.5

#define WIDTH 1000
#define HEIGHT 1000
#define CENTERX 500
#define CENTERY 0
#define SCALEX 100
#define SCALEY 100

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

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline __attribute__((__always_inline__)) unsigned int fast_rand(unsigned int *seed)
{
    *seed = (214013 * *seed + 2531011);
    return (*seed >> 16) & 0x7FFF;
}

inline __attribute__((__always_inline__)) double fast_rand_d(unsigned int *seed)
{
    return (double)fast_rand(seed) / 0x7FFF;
}

// ALIAS METHOD
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

inline __attribute__((__always_inline__)) int alias_sample(ifs system, unsigned int *seed)
{
    int i = fast_rand(seed) % N;
    double u = fast_rand_d(seed);
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

void set_random_ifs(ifs *system, unsigned int *seed)
{
    set_function(system, 0, fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed));

    for (int i = 1; i < N; i++)
        set_function(system, i, fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), fast_rand_d(seed), system->p[i - 1] + fast_rand_d(seed));
}

void set_ifs_example_1(ifs *system)
{
    assert(N == 4);

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

double normalize_min_max(double value, double min, double max)
{
    return (value - min) / (max - min);
}

double normalize_standard(double value, double mean, double stddev)
{
    return (value - mean) / stddev;
}

double mean(int *data, size_t length)
{
    double sum = 0.0;
    for (int i = 0; i < length; i++)
        sum += (double)data[i];
    return sum / length;
}

double stddev(int *data, size_t length, double mean)
{
    double sum = 0.0;
    for (int i = 0; i < length; i++)
        sum += (data[i] - mean) * (data[i] - mean);
    return sqrt(sum / length);
}

unsigned int clamp_255(double value)
{
    if (value < 0.0)
        return 0u;
    if (value > 255.0)
        return 255u;
    return (unsigned int)value;
}

// CRITICAL STUFF

inline __attribute__((__always_inline__)) pos func(pos p, double a, double b, double c, double d, double e, double f)
{
    return (pos){a * p.x + b * p.y + e, c * p.x + d * p.y + f};
}

inline __attribute__((__always_inline__)) pos next_iteration(pos p, ifs system, unsigned int *seed)
{
    int f = alias_sample(system, seed);
    return func(p, system.a[f], system.b[f], system.c[f], system.d[f], system.e[f], system.f[f]);
}

inline __attribute__((__always_inline__)) int inbounds(pos p)
{
    return (p.x >= 0.0 && p.x < WIDTH && p.y >= 0.0 && p.y < HEIGHT);
}

inline __attribute__((__always_inline__)) pos scale_and_center(pos p)
{
    return (pos){CENTERX + p.x * SCALEX, CENTERY + p.y * SCALEY};
}

// MAIN

int main(int argc, char *argv[])
{
    printf("Configuring IFS system...\n");

    // SETUP IFS SYSTEM
    ifs *system = malloc(sizeof(ifs));
    set_ifs_example_1(system);

    // SETUP ALIAS TABLES
    alias_preprocess(system);

    int screen[WIDTH * HEIGHT];

    // Thread count.
    omp_set_num_threads(32);

    printf("Computing %llu iterations...\n", (unsigned long long)ITERATIONS);

#pragma omp parallel
    {
        // Starting position.
        pos p = (pos){0.0, 0.0};

        // Set seed.
        unsigned int seed = omp_get_thread_num();

#pragma omp for
        for (unsigned long long i = 0; i < ITERATIONS; i++)
        {
            p = next_iteration(p, *system, &seed);
            pos p_s = scale_and_center(p);
            if (inbounds(p_s))
                screen[(int)p_s.x + (int)p_s.y * WIDTH]++;
        }
    }

    printf("Done computing iterations.\n");

    // // Find maximum value
    // int max = 0;
    // for (int i = 0; i < WIDTH * HEIGHT; i++)
    //     if (max < screen[i])
    //         max = screen[i];

    // double avg = mean(screen, WIDTH * HEIGHT);
    // double sd = stddev(screen, WIDTH * HEIGHT, avg);

    // unsigned int raw_value = screen[i + j * WIDTH];
    // unsigned int value = normalize_min_max((double) raw_value, 0, max) * 255;
    // unsigned int value = normalize_standard((double)raw_value, avg, sd) * 32 + 128;

    // Print the screen into a file.
    FILE *output = fopen("ifs3.pnm", "w");
    fprintf(output, "P2\n%d %d\n%d\n\n", WIDTH, HEIGHT, 255);
    for (int j = 0; j < HEIGHT; j++)
    {
        for (int i = 0; i < WIDTH; i++)
        {
            unsigned int raw_value = screen[i + j * WIDTH];

            unsigned int value = log2((double)raw_value + 1) / log2((double)ITERATIONS * GRADIENT) * 255;
            fprintf(output, "%4d", clamp_255(value));
        }
        fprintf(output, "\n");
    }

    fclose(output);

    return EXIT_SUCCESS;
}