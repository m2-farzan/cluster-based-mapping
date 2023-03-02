Test:
```
make test
./test
```

Visualize:
```
make
./main
./plot.py
```

# Description

This project is based on a raw idea and the main reason I did it is to refresh my memory of the C++ syntax. Anyway, the idea is to develop a navigation code that applies a clustering algorithm to the point cloud, so that we can reduce the dimensionality and extract useful information.

The point cloud is a set of points with known positions (relative to the robot). So each point (P) has two or three dimensions, Px, Py, Pz. We can extend these dimensions by appending extra information about the points. For example, we can do a Hough Line Transform or simple Kernel Convolution to the point cloud to assign each point a "straightness" value, Ps. Similarly, we can have roundness, clutterriness, angle, etc. If we add vision information, we can have "tree-ness", "rock-ness", "human-ness" and so on. We can also distinguish moving objects from stationary ones.

All these values can be assigned to every point in the cloud as a dimension. The idea is to feed these pre-processed points to a live (streaming) clustering code which reduces the dimensonality and extracts information about objects. This information can be used for behavioral decision-making, as well as localization. A robot may be able to infer its position when it sees a rock and a tree in special relative positions.

The roadmap:
- [x] Basic HDBSCAN
- [ ] Streaming HDBSCAN
  - [ ] Sparsify distance matrix
  - [ ] Optimize minimum spanning tree updates
- [ ] Use Bayesian inference to assign probabilities to clusters
- [ ] Use a reduction scheme to merge nearby points to a new point with doubled weight
- [ ] Add some demo pre-processing methods
  - [ ] LIDAR convolutional processing
- [ ] Make ROS package
- [ ] Testing & Improvement
