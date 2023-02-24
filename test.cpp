#include <gtest/gtest.h>
#include "hdbscan.h"

TEST(Test, SmallTree) {
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
    vector<int> E = mst(D);
    vector<Joint> SLT = single_linkage_tree(E, D);
    vector<Cluster> condensed_clusters = get_condensed_clusters(SLT, 3);
    vector<Cluster> stable_clusters = get_stable_clusters(condensed_clusters);

    EXPECT_EQ(SLT[0].a, 2);
    EXPECT_EQ(SLT[0].b, 4);
    EXPECT_EQ(SLT[0].delta, 2);
    EXPECT_EQ(SLT[0].size, 2);

    EXPECT_EQ(SLT[1].a, 1);
    EXPECT_EQ(SLT[1].b, 12);
    EXPECT_EQ(SLT[1].delta, 3);
    EXPECT_EQ(SLT[1].size, 3);

    EXPECT_EQ(SLT[2].a, 3);
    EXPECT_EQ(SLT[2].b, 13);
    EXPECT_EQ(SLT[2].delta, 4);
    EXPECT_EQ(SLT[2].size, 4);

    EXPECT_EQ(SLT[3].a, 7);
    EXPECT_EQ(SLT[3].b, 6);
    EXPECT_EQ(SLT[3].delta, 5);
    EXPECT_EQ(SLT[3].size, 2);

    EXPECT_EQ(SLT[4].a, 10);
    EXPECT_EQ(SLT[4].b, 9);
    EXPECT_EQ(SLT[4].delta, 6);
    EXPECT_EQ(SLT[4].size, 2);

    EXPECT_EQ(SLT[5].a, 11);
    EXPECT_EQ(SLT[5].b, 16);
    EXPECT_EQ(SLT[5].delta, 7);
    EXPECT_EQ(SLT[5].size, 3);

    EXPECT_EQ(SLT[6].a, 17);
    EXPECT_EQ(SLT[6].b, 0);
    EXPECT_EQ(SLT[6].delta, 8);
    EXPECT_EQ(SLT[6].size, 4);

    EXPECT_EQ(SLT[7].a, 8);
    EXPECT_EQ(SLT[7].b, 15);
    EXPECT_EQ(SLT[7].delta, 9);
    EXPECT_EQ(SLT[7].size, 3);

    EXPECT_EQ(SLT[8].a, 19);
    EXPECT_EQ(SLT[8].b, 5);
    EXPECT_EQ(SLT[8].delta, 15);
    EXPECT_EQ(SLT[8].size, 4);

    EXPECT_EQ(SLT[9].a, 20);
    EXPECT_EQ(SLT[9].b, 18);
    EXPECT_EQ(SLT[9].delta, 18);
    EXPECT_EQ(SLT[9].size, 8);

    EXPECT_EQ(SLT[10].a, 14);
    EXPECT_EQ(SLT[10].b, 21);
    EXPECT_EQ(SLT[10].delta, 20);
    EXPECT_EQ(SLT[10].size, 12);

    EXPECT_EQ(stable_clusters.size(), 3);

    EXPECT_EQ(stable_clusters[0].points[0], 3);
    EXPECT_EQ(stable_clusters[0].points[1], 1);
    EXPECT_EQ(stable_clusters[0].points[2], 2);
    EXPECT_EQ(stable_clusters[0].points[3], 4);

    EXPECT_EQ(stable_clusters[1].points[0], 5);
    EXPECT_EQ(stable_clusters[1].points[1], 8);
    EXPECT_EQ(stable_clusters[1].points[2], 7);
    EXPECT_EQ(stable_clusters[1].points[3], 6);

    EXPECT_EQ(stable_clusters[2].points[0], 0);
    EXPECT_EQ(stable_clusters[2].points[1], 11);
    EXPECT_EQ(stable_clusters[2].points[2], 10);
    EXPECT_EQ(stable_clusters[2].points[3], 9);
}

