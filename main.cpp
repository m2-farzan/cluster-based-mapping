#include "hdbscan.h"

void print_clusters(vector<Cluster>& clusters, int current_id, int indent=0) {
    for (int i = 0; i < indent; i++) {
        cout << " ";
    }
    cout << "Cluster(" << current_id << ")" << endl;
    for (int child_id : clusters[current_id].children_clusters) {
        print_clusters(clusters, child_id, indent+2);
    }
}

int main() {
    cout << "Program started" << endl;
    mat A {
        {4, -5},
        {0, 10},
        {1, 10},
        {2, 10},
        {1, 9},
        {0, 0}, // 5
        {-4, -4},
        {-5, -4},
        {-4, -5},
        {2, -4},
        {2, -5},
        {3, -5},
    };
    A = A.t();
    mat D(12, 12, arma::fill::ones);
    D *= 1000;
    D(2,4) = 2;
    D(1,4) = 3;
    D(3,4) = 4;
    D(6,7) = 5;
    D(9,10) = 6;
    D(9,11) = 7;
    D(9, 0) = 8;
    D(6,8) = 9;
    D(6,5) = 15;
    D(9,5) = 18;
    D(5,4) = 20;
    // reflect the matrix from the diagonal
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            D(i, j) = D(i, j) == 1000 ? D(j, i) : D(i, j);
        }
    }
    HDBSCAN model(3);

    // A = generate_samples();
    // D = distance_matrix(A);
    // HDBSCAN model(40);

    model.fit(A, D);
    auto E = model._get_mst();
    auto SLT = model._get_slt();
    auto condensed_clusters = model._get_condensed_clusters();
    auto stable_clusters = model.get_clusters();

    fstream A_csv;
    A_csv.open("a.points", ios::out);
    for (int i = 0; i < A.n_cols; i++) {
        A_csv << A(0, i) << " " << A(1, i) << "\n";
    }

    fstream MST;
    MST.open("mst.edges", ios::out);
    for (int i = 0; i < A.n_cols; i++) {
        MST << A(0, i) << " " << A(0, E[i]) << " " << A(1, i) << " " << A(1, E[i]) << "\n";
    }

    for (int i = 0; i < SLT.size(); i++) {
        cout << i << ": " << SLT[i].a << ", " << SLT[i].b << ", " << SLT[i].delta << ", " << SLT[i].size << endl;
    }

    print_clusters(condensed_clusters, 0);

    fstream M;
    M.open("membership.points", ios::out);
    for (int i = 0; i < stable_clusters.size(); i++) {
        for (int j = 0; j < stable_clusters[i].points.size(); j++) {
            int point_index = stable_clusters[i].points[j];
            M << A(0, point_index) << " " << A(1, point_index) << " " << i << "\n";
        }
    }
}
