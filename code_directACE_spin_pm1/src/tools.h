#ifndef TOOLS_H
#define TOOLS_H


#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>

// Data structures
#include "dataStructures.h"


// Typedefs
//typedef double Number;
//typedef std::vector<float> Vector;


//int pow(int b, int exp){ 

double innerProduct(double[], double[], int);
double L1(double[], int);
double L1(const std::vector<double> &);
double L2(double[], int);
double L2(const std::vector<double> &);
double LInfinity(double[], int);
double LInfinity(const std::vector<double> &);

void quickSort(double[][2], int, int);
void quickSort(std::vector<int> &, int, int);
bool checkSize(const Key &, const Key &);
void getSuperset(const Key &, const Key &, Key &);

void modifiedCholeskyInverse(double[], int);
double modifiedCholeskyInverseDet(double[], int);
double symmetricDet(double[], int);

void givens(double, double, double &, double &);
void householder(double[], int, double[], double &);
void tridiagonalize(double[], int, std::vector<std::vector<double> > &);
void QRStep(double[], int, std::vector<std::vector<double> > &);
void symmetricQR(double[], int, double[], std::vector<std::vector<double> > &);


// Checks to see whether the ith bit of an integer in a key is occupied

inline int checkBit(unsigned  long &key, int i) {
  
  return (key & (unsigned  long)1<<i);
    
}


// Flip the ith bit of an integer

inline void flipBit(unsigned  long &key, int i) {
   // cout << "FlipBit-key="<<key << " i=" <<i  << endl;
  //cout <<"*ok*" << endl;  
  key ^= (unsigned  long)1<<i;
    
}


// Count the number of bits set in a number (see http://www-graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan )

inline int countBits(unsigned  long v) {
   //cout << "CountBITS " << endl;
    unsigned long one=1;
    int c; // c accumulates the total bits set in v
    
    for (c=0; v; c++) v &= v - one; // clear the least significant bit set
   //cout <<"BITS " << c << endl;
    return c;
    
}

inline int countBitsCheap(unsigned long v) {
    
    //unsigned long one=1;
    int c=0; // c accumulates the total bits set in v

    if(v) { v &= v-1; c++; } else return c;
    if(v) { v &= v-1; c++; } else return c;
    if(v) c++;
    return c;

}

/*
// Count the number of bits set in a number (see http://www-graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan )

inline int countBits(char v) {
    
    //unsigned long one=1;
    int c; // c accumulates the total bits set in v
    
    for (c=0; v; c++) v &= v - 1; // clear the least significant bit set
    
    return c;
    
}
 */

// For j>i, the pair {i,j} in the list {{0,1}, {0,2}, ... } is located at offset(i,length)+j

inline int offset(int i, int length) {
    
    return (length + i * (length - 2) - (i * (i - 1)) / 2 - 1);
    
}


// For j>i, the pair {i,j} in the list {{0,0}, {0,1}, ... } is located at hindex(i,j,length)
// This routine is necessary when the diagonal term is also included

inline int hindex(int i, int j, int length) {
    
    return (i * length - (i * (i + 1)) / 2 + j);
    
}


// Returns as above, but with a check to make sure that j>i

inline int safeHindex(int i, int j, int length) {
    
    if (j>i) return (i * length - (i * (i + 1)) / 2 + j);
    else     return (j * length - (j * (j + 1)) / 2 + i);
    
}


// For j>i, the pair {i,j} in the list {{0},{1},...,{length-1},{0,1}, {0,2}, ... } is located at index(i,j,length)

inline int index(int i, int j, int length) {
    
    return (length + i * (length - 2) - (i * (i - 1)) / 2 - 1 + j);
    
}


// For k>j>i, return the index of {i,j,k} in the full list

inline int index(int i, int j, int k, int length) {
    
    return ( ( (i * (i + 1) * ( i + 2 )) / 3 
                  - (length-1) * (i * (i + 2) - 2 * j + 1) 
	       + (length-1)*(length-1) * (i - 1) 

                  - (j * (j + 1) - 2 * k)                ) / 2 );
    
}


// For l>k>j>i, return the index of {i,j,k,l} in the full list

inline int index(int i, int j, int k, int l, int length) {
    
    return ( (-i * (i + 1) * (i + 2) * (i + 3) 
              + 4 * (j * (1 + j) * (j + 2) - 3 * k * (k + 1) + 6 * l - 5 * (length-1)) 
              + 2 * (i * (11 + i * (9 + 2 * i)) - 6 * j * (j + 2) + 12 * k) * (length-1) 
              - 6 * (i * (i + 3) - 2 * j) * (length-1)* (length-1) 
              + 4 * (i - 1) * (length-1)* (length-1)* (length-1)                                           ) / 24 );
    
}


// Returns the sign of an integer x (and 0 if x==0)

inline int sign(int x) {
                
    static const int lookup_table[] = { -1, 1, 0, 0 };
    return lookup_table[(x==0) * 2 + (x>0)];

}


// Returns the sign of a double x (and 0 if x==0)

inline double sign(double x) {
    
    //double eps=pow(10.0,-20.0);
    
    static const double lookup_table[] = { -1, 1, 0, 0 };
    return lookup_table[(x==0) * 2 + (x>0)];
    
}


// Given components (a, b), this function computes parameters for a Givens rotation matrix G = ((c, s), (-s, c)),
// such that G^T (a, b) = (r, 0). See Matrix Computations (Golub) for details.

inline void givens(double a, double b, double &c, double &s) {
 
    if (b==0) { c=1; s=0; }
    
    else {
        
        if (fabs(b) > fabs(a)) { s = 1 / sqrt(1 + pow(a,2) / pow(b,2)); c = -a * s / b; }
        else                   { c = 1 / sqrt(1 + pow(b,2) / pow(a,2)); s = -b * c / a; }
        
    }
    
}


#endif
