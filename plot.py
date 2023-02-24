#!/bin/env python3

import matplotlib.pyplot as plt


with open('a.points') as f:
    for line in f.readlines():
        x, y = [float(d) for d in line.split(' ')]
        plt.plot(x, y, 'k.')

# plt.show()

# with open('mst.edges') as f:
#     for line in f.readlines():
#         x0, x1, y0, y1 = [float(d) for d in line.split(' ')]
#         plt.plot([x0, x1], [y0, y1], 'r-o')

# plt.show()


colors = ['r', 'g', 'b', 'k', 'y']
with open('membership.points') as f:
    for line in f.readlines():
        x, y, cluster_no = [float(d) for d in line.split(' ')]
        print(int(cluster_no))
        plt.plot(x, y, colors[int(cluster_no) % 5] + 'o')

plt.show()
