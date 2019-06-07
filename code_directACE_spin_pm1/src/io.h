#ifndef IO_H
#define IO_H


#include <vector>
#include <math.h>
#include <stdio.h>
#include <string>
#include <cstring>

// Data structures
#include "dataStructures.h"
using namespace std;

void getCouplings(FILE *, std::vector<double> &);
void printCorrelation(std::string, const Vector &);
void printFree(FILE *, double,  const std::vector<double> &, Number, int, int,unsigned long, unsigned long, float, double );


// Given the size of a coupling or correlation vector, return the corresponding number of spins

inline int sizetolength(int size) {
    
    return ((sqrt(1 + 8 * size) - 1) / 2);
    
}


#endif
