#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <climits>

using namespace std;

// Structure of a point in high dimensional space
struct Point {
    vector<double> coordinates;
};

// Calculate the distance between two points
double distance(Point p1, Point p2) {
    double sum = 0;
    for (int i = 0; i < p1.coordinates.size(); i++) {
        sum += (p1.coordinates[i] - p2.coordinates[i]) * (p1.coordinates[i] - p2.coordinates[i]);
    }
    return sqrt(sum);
}

// Find the median of a set of points
Point findMedian(vector<Point> points) {
    int n = points.size();
    int dimensions = points[0].coordinates.size();
    Point median;
    median.coordinates.resize(dimensions);

    // Sort each dimension independently
    for (int d = 0; d < dimensions; d++) {
        sort(points.begin(), points.end(), [d](const Point& p1, const Point& p2) {
            return p1.coordinates[d] < p2.coordinates[d];
        });
        if (n % 2 == 0) {
            median.coordinates[d] = (points[n / 2 - 1].coordinates[d] + points[n / 2].coordinates[d]) / 2;
        }
        else {
            median.coordinates[d] = points[n / 2].coordinates[d];
        }
    }
    return median;
}

// Find the k-median of a set of points
Point findKMedian(vector<Point> points, int k) {
    // If there are fewer than k points, return the median of the set
    if (points.size() <= k) {
        return findMedian(points);
    }

    // Split the points into two sets
    vector<Point> set1, set2;
    Point median = findMedian(points);
    for (Point p : points) {
        if (distance(p, median) < distance(p, set1.empty() ? median : findMedian(set1))) {
            set1.push_back(p);
        }
        else {
            set2.push_back(p);
        }
    }

    // Recursively find the k-median of each set
    if (set1.size() >= k) {
        return findKMedian(set1, k);
    }
    else {
        return findKMedian(set2, k - set1.size());
    }
}

int main() {
    vector<Point> points = {{{1, 2, 3}}, {{3, 4, 5}}, {{5, 6, 7}}, {{7, 8, 9}}, {{9, 10, 11}}};
    int k = 3;
    Point kMedian = findKMedian(points, k);

    return 0;
}