# Robust K-median clustering in C++

#### This is a C++ implementation of the robust K-median clustering algorithm.

K-median clustering is a type of unsupervised learning, which is used when you have unlabeled data (i.e., data without defined categories or groups). The goal of this algorithm is to find groups in the data, with the number of groups represented by the variable K. The algorithm works iteratively to assign each data point to one of K groups based on the features that are provided. Data points are clustered based on feature similarity.

    This implementation now contains multithreading support, which speeds up calculations for large vectors and thousands of points by parallelization.

## Instructions:
The input supports any number of points and any number of dimensions. Make the "input.csv" file accordingly.

* Run the robust_kmedian.cpp file or download the binary file "robust_kmedian.exe" from the repository.
* Run the robust_kmedian binary with the 2 input files name, number of clusters and output directory as command line arguments

Run this command :
`
./robust_kmedian CompleteData_miss=40_theta=0.05.csv Delta_miss=40_theta=0.05.csv 3 results
`

Output :
![Bash Output](image/cmd.PNG)
* The output will be 2 files in the output directory (eg: results) which contains the center point of each cluster, and cluster for each point.

![Clusters File Syntax](image/centroid.PNG)

-------------------------------
