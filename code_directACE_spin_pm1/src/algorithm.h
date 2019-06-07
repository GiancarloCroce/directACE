#ifndef ALGORITHM_H
#define ALGORITHM_H


// Data structures
#include "dataStructures.h"



// Cluster-related
void makeSingleCluster(Cluster &, int, double &);
void makeCluster(Cluster &, const Key &, int, bool, bool, double &);
void computeDFandDC(int[], int, std::vector<double> &, double &, bool, bool, double &);

// Main algorithm
void selectClusters_lax(std::set<Key> &, double,int, bool, std::set<Key> &, int);
void selectClusters_strict(std::set<Key> &, double,int, bool, std::set<Key> &, int);
void findNewClusters(const std::vector<Key> &);
void run(RunParameters &);




#endif
