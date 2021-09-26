#ifndef ALGEBRA_POLAR_H
#define ALGEBRA_POLAR_H

#endif //ALGEBRA_POLAR_H
#include <stdlib.h>

void Polar_Decomposition (int n, double (*Q)[n], double *singular, double (*P)[n],
                          double (*O)[n], double (*S)[n])
{
    double (*temp)[n] = calloc(n * n, n * n * sizeof(double));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) //
        temp[i][j] = singular[i] * P[j][i];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) //
        for (int k = 0; k < n; k++) S[i][j] += P[i][k] * temp[k][j];
    free(temp);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) O[i][j] += Q[i][k] * P[j][k];
}