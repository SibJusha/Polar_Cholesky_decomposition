#ifndef ALGEBRA_CHOLESKY_H
#define ALGEBRA_CHOLESKY_H

#endif //ALGEBRA_CHOLESKY_H

#include <stdio.h>
#include <stdlib.h>

double Diag_find(int n, double (*S)[n], double (*L)[n], double *diag, int idx) {
    double sum = 0;
    for (int i = 0; i < idx; ++i) {
        sum += L[idx][i] * L[idx][i] * diag[i];
    }
    return S[idx][idx] - sum;
}

double Lel_find(int n, double (*S)[n], double (*L)[n], double *diag, int idx_i, int idx_j) {
    double sum = 0;
    for (int i = 0; i < idx_j; ++i) {
        sum += L[idx_i][i] * L[idx_j][i] * diag[i];  //находим элементы L по формуле
    }
    return (1 / diag[idx_j]) * (S[idx_i][idx_j] - sum);
}

void cholesky(int n, double (*S)[n], double *diag, double (*L)[n]) {
    for (int j = 0; j < n; ++j) {
        diag[j] = Diag_find(n, S, L, diag, j);    //диагональные элементы вычисляются
        for (int i = 0; i < n; ++i) {
            if (i < j) {                    //вычисляются элементы L
                L[i][j] = 0;
            } else if (i == j) {
                L[i][j] = 1;
            } else {
                L[i][j] = Lel_find(n, S, L, diag, i, j);  //до сюда
            }
            if (L[i][j] == 0) L[i][j] = 0;
        }
    }
}