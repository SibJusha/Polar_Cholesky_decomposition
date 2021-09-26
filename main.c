#include "svd.h"
#include "stuff.h"
#include "polar.h"
#include "cholesky.h"
#include "linear.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
    FILE *inp, *out;
    inp = fopen("input.txt", "r");
    out = fopen("output.txt", "w");
    int n;
    fscanf(inp, "%d", &n);

    double (*A)[n] = malloc(n * n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    double (*Q)[n] = malloc(n * n * sizeof(double));
    double (*S)[n] = calloc(n * n, n * n * sizeof(double));
    double *diag = (double *) malloc(n * sizeof(double));
    double (*L)[n] = malloc(n * n * sizeof(double));
    double *y = calloc(n, n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fscanf(inp, "%lf", &A[i][j]);
    for (int i = 0; i < n; i++) {
        fscanf(inp, "%lf", b + i);
    }

    if (A == NULL || Q == NULL
        || S == NULL || diag == NULL || L == NULL || y == NULL) {
        printf("No memory available\n");
        exit(-1);
    }
    fclose(inp);

    if (Symmetric_Check(n, A) == 1) {
        cholesky(n, A, diag, L);

        Triangle_Solve_Low(n, L, b, y); // L*y = b
        Square_Matrix_Transpose(n, A, L); // A = Lt
        Diagonal_Matrix_Multiplication(n, LEFT, Q, A, diag); // Q = diag * Lt = diag * A
        Triangle_Solve_Upper(n, Q, y, diag); // diag * Lt * x = y
        for (int i = 0; i < n; i++) fprintf(out, "%lf ", diag[i]);
        fprintf(out, "\n");

        fclose(out);
        free(A);
        free(Q);
        free(S);
        free(diag);
        free(L);
        free(b);
        return 0;
    }
    double (*P)[n] = malloc(n * n * sizeof(double));
    double *singular_values = (double *) malloc(n*sizeof(double));
    double *dummy_array = (double*) malloc(n * sizeof(double));
    double (*O)[n] = calloc(n * n, n * n * sizeof(double));
    if (P == NULL || singular_values == NULL
       || dummy_array == NULL || O == NULL ) {
        printf("No memory available\n");
        exit(-1);
    }

    int err = Singular_Value_Decomposition((double*) A, n, n, (double*) Q,
                                        singular_values, (double*) P, dummy_array);
    if (err < 0) printf("Failed to converge\n");

    Polar_Decomposition(n, Q, singular_values, P, O, S);

    cholesky(n, S, diag, L);

    Square_Matrix_Transpose(n, P, O); // P = Ot

    On_Vector_Multiplication(n, singular_values, P, b); // S = Ot * b

    Triangle_Solve_Low(n, L, singular_values, y); // L*y = Ot * b
    Square_Matrix_Transpose(n, A, L); // A = Lt
    Diagonal_Matrix_Multiplication(n, LEFT, Q, A, diag); // Q = diag * Lt = diag * A
    Triangle_Solve_Upper(n, Q, y, diag); // diag * Lt * x = L*y
    for (int i = 0; i < n; i++) fprintf(out, "%lf.8 ", diag[i]);
    fprintf(out, "\n");

    fclose(out);
    free(dummy_array);
    free(A);
    free(Q);
    free(P);
    free(O);
    free(S);
    free(diag);
    free(L);
    free(singular_values);
    free(b);
    free(y);
    return 0;
}
