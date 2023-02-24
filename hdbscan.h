#include <iostream>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>

using namespace std;
// using namespace arma;
typedef arma::mat mat;
typedef arma::umat umat;

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

vector<int> mst(mat D) {
    const int N_POINTS = D.n_cols;
    vector<int> E(N_POINTS, -1);
    mat C(1, N_POINTS);  // minimum cost of connecting the node
    umat C_other(2, N_POINTS);  // how is the minimum umat
    E[0] = 0;
    C = D.row(0);
    for (int i = 0; i < N_POINTS; i++) {
        C_other(0, i) = i;  // the one that isn't already in the tree
        C_other(1, i) = 0;
    }

    // for (int p = 0; p < 5; p++) {
    for (;;) {
        int next_edge_a = -1;
        int next_edge_b = -1;
        double next_edge_cost = 1e30;
        for (int i = 1; i < N_POINTS; i++) {
            if (E[i] == -1 && C[i] < next_edge_cost) {
                next_edge_cost = C[i];
                next_edge_a = C_other(0, i);
                next_edge_b = C_other(1, i);
            }
        }
        if (next_edge_a == -1) {
            break;
        }
        E[next_edge_a] = next_edge_b;
        for (int i = 0; i < N_POINTS; i++) {
            if (D(next_edge_a, i) < C(i)) {
                C(i) = D(next_edge_a, i);
                C_other(0, i) = i;
                C_other(1, i) = next_edge_a;
            }
        }
    }
    return E;
}

struct Joint {
    int a;
    int b;
    double delta;
    int size;
};

struct Edge {
    int a;
    int b;
    double distance;
};

bool edgeComparator(const Edge& a, const Edge& b) {
    return a.distance < b.distance;
}

struct UnionFind {
    vector<int> parent;
    vector<int> size;
    int next_label;

    UnionFind(int N) :
        parent(2*N, -1),
        size(2*N, 0),
        next_label(N)
    {
        for (int i = 0; i < N; i++) {
            size[i] = 1;
        }
    }

    void union_(int a, int b) {
        size[next_label] = size[a] + size[b];
        parent[a] = next_label;
        parent[b] = next_label;
        size[next_label] = size[a] + size[b];
        next_label++;
    }

    int fast_find(int n) {
        int p = n;
        while (parent[n] != -1) {
            n = parent[n];
        }
        // cache
        while (parent[p] != -1 && parent[p] != n) {  // but how does the pyx even work??
            int next_p = parent[p];
            parent[p] = n;
            p = next_p;
        }
        return n;
    }
};

vector<Joint> single_linkage_tree(vector<int> mst, mat D) {
    int N = mst.size();
    vector<Edge> L(N - 1);
    for (int i = 1; i < N; i++) {
        L[i-1] = { i, mst[i], D(i, mst[i]) };  // OPTIMIZE: get rid of reallocation
    }
    sort(L.begin(), L.end(), edgeComparator);
    UnionFind u(N);
    vector<Joint> result(N - 1);
    for (int i = 0; i < N-1; i++) {
        int root_a = u.fast_find(L[i].a);
        int root_b = u.fast_find(L[i].b);
        result[i] = {
            root_a,
            root_b,
            L[i].distance,
            u.size[root_a] + u.size[root_b],
        };
        u.union_(root_a, root_b);
    }
    return result;
}

struct Cluster {
    vector<int> points;
    vector<double> point_membership_strengths;
    double stability = 0;
    // int size = 0;
    double lambda_birth;
    vector<int> children_clusters;
    // bool selected = false;
};

void resolve_children(vector<Joint> slt, int current_joint_no, int current_base_cluster, vector<Cluster>& clusters, int MIN_CLUSTER_SIZE) {
    int N = slt.size() + 1;
    int a = slt[current_joint_no-N].a;
    int b = slt[current_joint_no-N].b;

    if (a < N) {
        clusters[current_base_cluster].points.push_back(a);
        double point_lambda = 1.0 / slt[current_joint_no-N].delta;
        clusters[current_base_cluster].point_membership_strengths.push_back(point_lambda);
        clusters[current_base_cluster].stability += point_lambda - clusters[current_base_cluster].lambda_birth; // fix
        // clusters[current_base_cluster].size++;
    }
    if (b < N) {
        clusters[current_base_cluster].points.push_back(b);
        double point_lambda = 1.0 / slt[current_joint_no-N].delta;
        clusters[current_base_cluster].point_membership_strengths.push_back(point_lambda);
        clusters[current_base_cluster].stability += point_lambda - clusters[current_base_cluster].lambda_birth; // fix
        // clusters[current_base_cluster].size++;
    }
    if (a >= N) {
        if (slt[a-N].size >= MIN_CLUSTER_SIZE && b >= N && slt[b-N].size >= MIN_CLUSTER_SIZE) {
            Cluster staged_cluster;
            staged_cluster.lambda_birth = 1.0 / slt[current_joint_no-N].delta;
            clusters.push_back(staged_cluster);
            clusters[current_base_cluster].children_clusters.push_back(clusters.size() - 1);
            resolve_children(slt, a, clusters.size()-1, clusters, MIN_CLUSTER_SIZE);
        } else {
            resolve_children(slt, a, current_base_cluster, clusters, MIN_CLUSTER_SIZE);
        }
    }
    if (b >= N) {
        if (slt[b-N].size >= MIN_CLUSTER_SIZE && a >= N && slt[a-N].size >= MIN_CLUSTER_SIZE) {
            Cluster staged_cluster;
            staged_cluster.lambda_birth = 1.0 / slt[current_joint_no-N].delta;
            clusters.push_back(staged_cluster);
            clusters[current_base_cluster].children_clusters.push_back(clusters.size() - 1);
            resolve_children(slt, b, clusters.size()-1, clusters, MIN_CLUSTER_SIZE);
        } else {
            resolve_children(slt, b, current_base_cluster, clusters, MIN_CLUSTER_SIZE);
        }
    }
}

vector<Cluster> get_condensed_clusters(vector<Joint> slt, int MIN_CLUSTER_SIZE) {
    vector<Cluster> result(1);
    result[0].lambda_birth = 0; // FIXME
    resolve_children(slt, 2 * slt.size(), 0, result, MIN_CLUSTER_SIZE);
    return result;
}

vector<Cluster> make_selection(vector<Cluster>& clusters, Cluster& current_cluster) {
    vector<Cluster> selection;  // Optimize: https://stackoverflow.com/questions/354442/looking-for-c-stl-like-vector-class-but-using-stack-storage
    if (current_cluster.children_clusters.size() == 0) {
        selection.push_back(current_cluster);
        return selection;
    }
    double children_total_stability = 0;
    for (int i = 0; i < current_cluster.children_clusters.size(); i++) {
        vector<Cluster> selected_children = make_selection(clusters, clusters[current_cluster.children_clusters[i]]);
        children_total_stability += clusters[current_cluster.children_clusters[i]].stability;
        selection.insert(selection.end(), selected_children.begin(), selected_children.end());
    }
    if (current_cluster.stability < children_total_stability) {
        current_cluster.stability = children_total_stability;
        return selection;
    } else {
        selection.clear();
        selection.push_back(current_cluster);
        return selection;
    }
}

vector<Cluster> get_stable_clusters(vector<Cluster> condensed_clusters) {
    return make_selection(condensed_clusters, condensed_clusters[0]);
}

// vector<int> cluster(vector<int> mst, mat D) {

// }
