#include <stdio.h>
#include <stdlib.h>

#define N 1000000
#define SEED 0
#define randD (float) rand() / RAND_MAX

int get_index(double v, double p[N])
{
    int L,R,M;
    L = 0;
    R = N - 1;

    while (L <= R) {
        M = L + (R - L) / 2;
        if (p[M] < v)
            L = M + 1;
        else if (p[M] > v)
            R = M - 1;
        else
            return M;
    }

    return -1;
}

int main()
{
    srand(SEED);

    double test[N] = {0};
    for (int i = 0; i < N - 1; i++)
        test[i + 1] = test[i] + randD;

    for (int i = 0; i < N; i++)
        get_index(test[i], test);
}
