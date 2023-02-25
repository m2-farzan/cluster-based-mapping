#pragma once

#include <iostream>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include "unionfind.h"
#include "util.h"

using namespace std;
typedef arma::mat mat;
typedef arma::umat umat;

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

struct Cluster {
    vector<int> points;
    vector<double> point_membership_strengths;
    double stability = 0;
    // int size = 0;
    double lambda_birth;
    vector<int> children_clusters;
    // bool selected = false;
};


class HDBSCAN {
    private:
    mat P;
    mat D;
    vector<int> mst;
    vector<Joint> slt;
    vector<Cluster> clusters; // = condensed clusters
    vector<Cluster> stable_clusters;
    int MIN_CLUSTER_SIZE;

    void calc_mst() {
        const int N_POINTS = D.n_cols;
        mst = vector<int>(N_POINTS, -1);
        mat C(1, N_POINTS);  // minimum cost of connecting the node
        umat C_other(2, N_POINTS);  // how is the minimum umat
        mst[0] = 0;
        C = D.row(0);
        for (int i = 0; i < N_POINTS; i++) {
            C_other(0, i) = i;  // the one that isn't already in the tree
            C_other(1, i) = 0;
        }

        for (;;) {
            int next_edge_a = -1;
            int next_edge_b = -1;
            double next_edge_cost = 1e30;
            for (int i = 1; i < N_POINTS; i++) {
                if (mst[i] == -1 && C[i] < next_edge_cost) {
                    next_edge_cost = C[i];
                    next_edge_a = C_other(0, i);
                    next_edge_b = C_other(1, i);
                }
            }
            if (next_edge_a == -1) {
                break;
            }
            mst[next_edge_a] = next_edge_b;
            for (int i = 0; i < N_POINTS; i++) {
                if (D(next_edge_a, i) < C(i)) {
                    C(i) = D(next_edge_a, i);
                    C_other(0, i) = i;
                    C_other(1, i) = next_edge_a;
                }
            }
        }
    }

    void calc_single_linkage_tree() {
        int N = mst.size();
        vector<Edge> L(N - 1);
        for (int i = 1; i < N; i++) {
            L[i-1] = { i, mst[i], D(i, mst[i]) };  // OPTIMIZE: get rid of reallocation
        }
        sort(L.begin(), L.end(), edgeComparator);
        UnionFind u(N);
        slt = vector<Joint>(N - 1);
        for (int i = 0; i < N-1; i++) {
            int root_a = u.fast_find(L[i].a);
            int root_b = u.fast_find(L[i].b);
            slt[i] = {
                root_a,
                root_b,
                L[i].distance,
                u.size[root_a] + u.size[root_b],
            };
            u.union_(root_a, root_b);
        }
    }

    void resolve_children(int current_joint_no, int current_base_cluster) {
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
                resolve_children(a, clusters.size()-1);
            } else {
                resolve_children(a, current_base_cluster);
            }
        }
        if (b >= N) {
            if (slt[b-N].size >= MIN_CLUSTER_SIZE && a >= N && slt[a-N].size >= MIN_CLUSTER_SIZE) {
                Cluster staged_cluster;
                staged_cluster.lambda_birth = 1.0 / slt[current_joint_no-N].delta;
                clusters.push_back(staged_cluster);
                clusters[current_base_cluster].children_clusters.push_back(clusters.size() - 1);
                resolve_children(b, clusters.size()-1);
            } else {
                resolve_children(b, current_base_cluster);
            }
        }
    }

    void calc_condensed_clusters() {
        clusters = vector<Cluster>(1);
        clusters[0].lambda_birth = 0; // FIXME
        resolve_children(2 * slt.size(), 0);
    }

    vector<Cluster> make_selection(Cluster& current_cluster) {
        vector<Cluster> selection;  // Optimize: https://stackoverflow.com/questions/354442/looking-for-c-stl-like-vector-class-but-using-stack-storage
        if (current_cluster.children_clusters.size() == 0) {
            selection.push_back(current_cluster);
            return selection;
        }
        double children_total_stability = 0;
        for (int i = 0; i < current_cluster.children_clusters.size(); i++) {
            vector<Cluster> selected_children = make_selection(clusters[current_cluster.children_clusters[i]]);
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

    void calc_stable_clusters() {
        stable_clusters = make_selection(clusters[0]);
    }

    public:
    HDBSCAN(int MIN_CLUSTER_SIZE) :
        MIN_CLUSTER_SIZE(MIN_CLUSTER_SIZE)
    {}

    void fit(mat P) {
        D = distance_matrix(P);
        fit(P, D);
    }

    void fit(mat P, mat D) {
        this->P = P;
        this->D = D;
        calc_mst();
        calc_single_linkage_tree();
        calc_condensed_clusters();
        calc_stable_clusters();
    }

    vector<Cluster> get_clusters() {
        return stable_clusters;
    }

    auto _get_mst() { return mst; }
    auto _get_slt() { return slt; }
    auto _get_condensed_clusters() { return clusters; }

};
