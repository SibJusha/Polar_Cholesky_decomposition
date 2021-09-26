#include <stdarg.h>
#include <stdio.h>
#ifndef ALGEBRA_STUFF_H
#define ALGEBRA_STUFF_H

#endif //ALGEBRA_STUFF_H

typedef enum side {LEFT, RIGHT} side;
typedef enum where {TEXT, CONSOLE} where;

void Square_Matrix_Multiplication (int n, double (*resulted)[n], double (*a)[n], double (*b)[n]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) resulted[i][j] = 0;

    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++)
        resulted[i][j] += a[i][k] * b[k][j];
}

void Square_Matrix_Transpose (int n, double (*res)[n], double (*a)[n]) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        res[i][j] = a[j][i];
}

void Diagonal_Matrix_Multiplication (int n, side s, double (*resulted)[n],
                                     double (*a)[n], double *diag)
{
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) resulted[i][j] = 0;
    if (s == LEFT) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            resulted[i][j] = diag[i] * a[i][j];
    } else if (s == RIGHT) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            resulted[i][j] = diag[j] * a[i][j];
    } else printf("WASTED");
}

void On_Vector_Multiplication (int n, double *resulted, double (*a)[n], double *vector) {
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) sum += vector[j] * a[i][j];
        resulted[i] = sum;
    }
}

void Print_Matrix (int n, int m, double (*A)[n], where b, FILE *c) {
    if (b == CONSOLE)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) printf("%.20lf ", A[i][j]);
        printf("\n");
    }
    else if (b == TEXT)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) fprintf(c, "%lf ", A[i][j]);
            fprintf(c, "\n");
        }
    printf("\n");
}

int Symmetric_Check (int n, double (*a)[n]) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        if (a[i][j] != a[j][i]) return 0;
    return 1;
}
