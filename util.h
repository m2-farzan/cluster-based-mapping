#pragma once

#include <armadillo>
typedef arma::mat mat;

mat generate_samples() {
    const int N_DIM = 2;
    const int N_POINTS = 200;

    arma::arma_rng::set_seed(0);

    mat A1(N_DIM, N_POINTS / 2, arma::fill::randn);
    mat COV(N_DIM, N_DIM);
    COV = {{ 1, 0.8 }, { 0.8, 1 }};
    A1 = COV * A1;

    mat A2(N_DIM, N_POINTS / 2, arma::fill::randn);
    COV = {{ 0.5, -0.3 }, { -0.3, 1 }};
    A2 = COV * A2;
    for (int i = 0; i < N_POINTS/2; i++) {
        A2(0, i) += -4;
        A2(1, i) += 5;
    }

    mat A = join_horiz(A1, A2);
    return A;
}

mat distance_matrix(mat A) {
    const int N_DIM = A.n_rows;
    const int N_POINTS = A.n_cols;
    mat dist(N_POINTS, N_POINTS);
    for (int i = 0; i < N_POINTS; i++) {
        for (int j = 0; j < N_POINTS; j++) {
            dist(i, j) = hypot(A(0, i) - A(0, j), A(1, i) - A(1, j));
        }
    }
    return dist;
}
