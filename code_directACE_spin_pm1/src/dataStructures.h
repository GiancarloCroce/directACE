#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H


#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>
#include <iostream>
#include <bitset>
#include <stdio.h>

using namespace std;

// Typedefs
typedef double Number;
typedef std::vector<double> Vector;


// STRING MANIPULATION

// Converts generic to string

template <class T>

inline std::string tostring (const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}


// Converts a string to an integer

inline int strtoint(const std::string &s) {
    
    std::istringstream i(s);
    int x;
    
    if (!(i >> x)) printf("String to integer conversion failed!");
    return x;
    
}


// Converts a string to a double

inline double strtodouble(const std::string &s) {
    
    std::istringstream i(s);
    double x;
    
    if (!(i >> x)) printf("String to double conversion failed!");
    return x;

}


// CLUSTER CLASSES

// Key for referencing clusters in the clusterIndex map; this is made using an int array so that
// temporary keys can be built on the stack, rather than having to use std::vector<int> on the heap

typedef std::vector<unsigned long> Key;
// Cluster class, which holds the information for a given cluster

class Cluster {
    
public:
	
  Number dF;
  Vector dC;
  	
  Cluster() {}
  Cluster(int size) {

    dF = 0;
    dC.resize((size * (size + 1)) / 2);

  }
  ~Cluster() {}
};


// Functor which is used to determine whether or not a cluster should
// be counted among the significant clusters

class Significant {
    
public:
    
    double theta;
    double kmin;
    
    Significant(double t, int k){
        theta = t;
        kmin = k;
        };

    
    bool operator() (const Cluster &c) {
        
        int clust_size = 0;
        clust_size = (-1 + sqrt(1+8*c.dC.size())) / 2;
        bool sign_or_not = ((fabs(c.dF) >= theta) || clust_size <= kmin);
        //cout << "clust_size: " <<clust_size << " kmin:" <<kmin  << " "<< sign_or_not << endl;
        return sign_or_not;
        //return ((fabs(c.dF) > theta) || (c.dC.size()==1) || clust_size <= kmin);
        //return ( fabs(c.dF) >=theta);
    }
    
};


// PROGRAM SETTINGS

// This class holds the parameters needed for running the cluster algorithm

class RunParameters {
    
public:
    
    std::string directory;      // Path to the directory where the inut file is located
                                // Output will also be sent to this directory
    std::string infile;         // Input file from which correlations are to be read
    std::string outfile;        // Output file (prefix) where data is to be written
    
    int kmax;                   // Maximum cluster size before computation is truncated
    int kmin;                   // Minimal cluster size to be considered as significant
    double theta;               // The (starting) value of the cutoff
    double thetaMin;            // The minimum cutoff value (loop starts here)
    double beta;                // Inverse temperature 

    //TO BE IMPLEMENTED
    double thetaMax;            // The maximum cutoff value (loop ends here)  
    double thetaStep;           // Size of logarithmic steps when looping over theta

    bool lax;                // LAX selection rule
 
    bool useSupplementaryOutput;// If true, print supplementary information
    bool useVerbose;            // If true, print extra information while program is running

    
    RunParameters() {
        
        directory="";
        infile="in.corr";
        outfile="out";
        kmax=0;
        kmin=0;
        theta=pow(10.0,0);
        thetaMin=0.0;
        thetaMax=pow(10.0,0);
        thetaStep=1.05;
        beta=1;
        lax=false;
        useSupplementaryOutput=false;
        useVerbose=false;
        
    }
    std::string getInfile()                        { std::string s=(directory+"/"+infile+".hj");                         return s; }
    std::string getCouplingsOutfile()              { std::string s=(directory+"/"+outfile);                              return s; }
    std::string getFullOutfile()                   { std::string s=(directory+"/"+outfile+".full");                      return s; }
    std::string getCorrelationsOutfile()           { std::string s=(directory+"/"+outfile+".outCorr");                   return s; }
    std::string getSupplementaryOutfile()          { std::string s=(directory+"/"+outfile+".sce");                       return s; }
    std::string getCouplingsOutfile(int number)    { std::string s=(directory+"/"+outfile+"_"+tostring(number));         return s; }
    std::string getFullOutfile(int number)         { std::string s=(directory+"/"+outfile+"_"+tostring(number)+".full"); return s; }
    std::string getCorrelationsOutfile(int number) { std::string s=(directory+"/"+outfile+"_"+tostring(number)+".corr"); return s; }
    ~RunParameters() {}
    
};


#endif

