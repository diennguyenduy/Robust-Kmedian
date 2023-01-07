#include <omp.h> // library for parallel programming in the SMP
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

// Convert a line that is read from input file to a vector of value
vector<double> lineToVec(string &line) {
    vector<double> values;
    string tmp = "";
    int length = (int)line.length();

    for (int i = 0; i < (int)line.length(); i++)
    {
        if ((48 <= int(line[i]) && int(line[i])  <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
        {
            tmp += line[i];
        }
        else if (tmp.length() > 0)
        {
            values.push_back(stod(tmp)); //convert string to float and then adds to the end of the vector.
            tmp = "";
        }
    }

    if (tmp.length() > 0) {
        values.push_back(stod(tmp));
        tmp = "";
    }

    return values;
}

class Point
{
private:
    int pointId, clusterId;
    int dimensions; // number of dimensions (features)
    vector<double> values; //vector of values/features of 1 data point

public:
    Point(int id, string line1)
    {
        pointId = id;
        values = lineToVec(line1);
        dimensions = values.size();
        clusterId = 0; // Initially not assigned to any cluster
    }

    int getDimensions() { return dimensions; }

    int getCluster() { return clusterId; }

    int getID() { return pointId; }

    void setCluster(int val) { clusterId = val; }

    double getVal(int pos) { return values[pos]; }
};

class Cluster
{
private:
    int clusterId;
    vector<double> centroid;
    vector<Point> points;

public:
    Cluster(int clusterId, Point centroid)
    {
        this->clusterId = clusterId;
        for (int i = 0; i < centroid.getDimensions(); i++)
        {
            this->centroid.push_back(centroid.getVal(i));
        }
        this->addPoint(centroid);
    }

    void addPoint(Point p)
    {
        p.setCluster(this->clusterId);
        points.push_back(p);
    }

    bool removePoint(int pointId)
    {
        int size = points.size();

        for (int i = 0; i < size; i++)
        {
            if (points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    void removeAllPoints() { points.clear(); }

    int getId() { return clusterId; }

    Point getPoint(int pos) { return points[pos]; }

    int getSize() { return points.size(); }

    double getCentroidByPos(int pos) { return centroid[pos]; }

    void setCentroidByPos(int pos, double val) { this->centroid[pos] = val; }
};

class KMedian
{
private:
    int K, iters, dimensions, total_points;
    vector<Cluster> clusters;
    string output_dir;

    void clearClusters()
    {
        for (int i = 0; i < K; i++)
        {
            clusters[i].removeAllPoints();
        }
    }

    int getNearestClusterId(Point point, vector<double> delta)
    {
        double sum = 0.0, min_dist;
        int NearestClusterId;

        if(dimensions==1) {
            min_dist = abs(clusters[0].getCentroidByPos(0) - point.getVal(0));
        }
        else
        {
          for (int i = 0; i < dimensions; i++)
          {
            sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i) + delta[i]);
          }
          min_dist = sum;
        }
        NearestClusterId = clusters[0].getId();

        for (int i = 1; i < K; i++)
        {
            double dist;
            sum = 0.0;

            if(dimensions==1) {
                  dist = abs(clusters[i].getCentroidByPos(0) - point.getVal(0));
            }
            else {
              for (int j = 0; j < dimensions; j++)
              {
                sum += (abs(clusters[i].getCentroidByPos(j) - point.getVal(j)) + delta[j]);
              }
                dist = sum;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                NearestClusterId = clusters[i].getId();
            }
        }

        return NearestClusterId;
    }

public:
    KMedian(int K, int iterations, string output_dir)
    {
        this->K = K;
        this->iters = iterations;
        this->output_dir = output_dir;
    }

    void run(vector<Point> &all_points, vector<vector<double>> &delta)
    {
        total_points = all_points.size();
        dimensions = all_points[0].getDimensions();

        // Initializing Clusters
        vector<int> used_pointIds;

        for (int i = 1; i <= K; i++)
        {
            while (true)
            {
                int index = rand() % total_points;

                if (find(used_pointIds.begin(), used_pointIds.end(), index) ==
                    used_pointIds.end())
                {
                    used_pointIds.push_back(index);
                    all_points[index].setCluster(i);
                    Cluster cluster(i, all_points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
        cout << "Clusters initialized = " << clusters.size() << endl
             << endl;

        cout << "Running K-Medians Clustering.." << endl;

        int iter = 1;
        while (true)
        {
            cout << "Iter - " << iter << "/" << iters << endl;
            bool done = true;

            // Add all points to their nearest cluster
            #pragma omp parallel for reduction(&&: done) num_threads(16)
            for (int i = 0; i < total_points; i++)
            {
                int currentClusterId = all_points[i].getCluster();
                int nearestClusterId = getNearestClusterId(all_points[i], delta[i]);

                if (currentClusterId != nearestClusterId)
                {
                    all_points[i].setCluster(nearestClusterId);
                    done = false;
                }
            }

            // clear all existing clusters
            clearClusters();

            // reassign points to their new clusters
            for (int i = 0; i < total_points; i++)
            {
                // cluster index is ID-1
                clusters[all_points[i].getCluster() - 1].addPoint(all_points[i]);
            }

            // Recalculating the center of each cluster
            for (int i = 0; i < K; i++)
            {
                int ClusterSize = clusters[i].getSize();

                for (int j = 0; j < dimensions; j++)
                {
                    double sum = 0.0;
                    if (ClusterSize > 0)
                    {
                        #pragma omp parallel for reduction(+: sum) num_threads(16)
                        for (int p = 0; p < ClusterSize; p++)
                        {
                            sum += clusters[i].getPoint(p).getVal(j);
                        }
                        clusters[i].setCentroidByPos(j, sum / ClusterSize);
                    }
                }
            }

            if (done || iter >= iters)
            {
                cout << "Clustering completed in iteration : " << iter << endl
                     << endl;
                break;
            }
            iter++;
        }

        //Write class of point to file
        ofstream pointsFile;
        pointsFile.open(output_dir + "/" + to_string(K) + "-points-KMedians.txt", ios::out);
        for (int i = 0; i < total_points; i++)
        {
            pointsFile << all_points[i].getCluster() << endl;
        }

        pointsFile.close();

        // Write clusters's centroid to file
        ofstream outfile;
        outfile.open(output_dir + "/" + to_string(K) + "-clusters-Kmedians.txt");
        if (outfile.is_open())
        {
            for (int i = 0; i < K; i++)
            {
                cout << "Cluster " << clusters[i].getId() << " centroid : ";
                for (int j = 0; j < dimensions; j++)
                {
                    cout << clusters[i].getCentroidByPos(j) << " ";    // Output to console
                    outfile << clusters[i].getCentroidByPos(j) << " "; // Output to file
                }
                cout << endl;
                outfile << endl;
            }
            outfile.close();
        }
        else
        {
            cout << "Error: Unable to write to clusters.txt";
        }
    }
};

int main(int argc, char **argv)
{
    // Need 4 arguments (except filename) to run, else exit
    if (argc != 5)
    {
        cout << "Error: command-line argument count mismatch. \n ./robust_kmedian <INPUT> <DELTA> <K> <OUT-DIR>" << endl;
        return 1;
    }

    string output_dir = argv[4];

    // Fetching number of clusters
    int K = atoi(argv[3]);

    // Open file for fetching points
    string filename1 = argv[1];
    ifstream infile1(filename1.c_str());

    if (!infile1.is_open())
    {
        cout << "Error: Failed to open file." << endl;
        return 1;
    }

    // Fetching points from file
    int pointId = 1;
    vector<Point> all_points;
    string line1;

    getline(infile1, line1); //skip the hearder when reading
    while (getline(infile1, line1))
    {
        Point point(pointId, line1);
        all_points.push_back(point);
        pointId++;
    }

    infile1.close();
    cout << "\nData fetched successfully!" << endl
         << endl;

    // Return if number of clusters > number of points
    if ((int)all_points.size() < K)
    {
        cout << "Error: Number of clusters greater than number of points." << endl;
        return 1;
    }

    // Open file for fetching delta values
    vector<vector<double>> delta;
    string filename2 = argv[2];
    ifstream infile2(filename2.c_str());
    string line2;

    getline(infile2, line2); //skip the hearder when reading
    while (getline(infile2, line2))
    {
        vector<double> values;
        values = lineToVec(line2);
        delta.push_back(values);
    }

    infile2.close();

    // Running K-Means Clustering
    int iters = 100;

    KMedian KMedian(K, iters, output_dir);
    KMedian.run(all_points, delta);

    return 0;
}
