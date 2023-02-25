#pragma once
#include <vector>
using namespace std;

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
