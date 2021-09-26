#ifndef ALGEBRA_LINEAR_H
#define ALGEBRA_LINEAR_H

#endif //ALGEBRA_LINEAR_H
#include <stdlib.h>

void Triangle_Solve_Low (int n, double (*a)[n], double *b, double *x) { // O(n^2)
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) sum += x[j] * a[i][j];
        x[i] = (b[i] - sum) / a[i][i];
    }
}

void Triangle_Solve_Upper (int n, double (*a)[n], double *b, double *x) {
    for (int i = n - 1; -1 < i; i--) {
        double sum = 0;
        for (int j = n - 1; i < j; j--) sum += x[j] * a[i][j];
        x[i] = (b[i] - sum) / a[i][i];
    }
}