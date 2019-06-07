#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <fstream>

#include <fstream>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

// Data input and output
#include "io.h"

// Numerical tools
#include "tools.h"

using namespace std;

// Read the number of spins and the parameters of the model (couplings and fields) from the input file
void getCouplings(FILE *input, std::vector<double> &hj) {

  double n;
  fscanf(input,"%le",&n);
	
  hj.resize((n*(n+1))/2);
  
  int cont=0;	  
  for (int i=0;i<n;i++) {
    double p1; 		
    fscanf(input,"%le",&p1);
    hj[cont]=p1;
    cont++;	
  }

  
  int count=n;
  for(int i=0;i<n-1;i++){
    for(int j=i+1;j<n;j++){
      double  p2;
      fscanf(input,"%le",&p2);
      hj[count]=p2;
      
      count++;		
    }
  }

}

//print magentization (file .ms) and correlation (file .ccs) corresponding to the last theta
void printCorrelation(std::string file, const Vector &C) {

    int n = sizetolength(C.size());
    
    std::string mstring=file+".ms";
    std::string ccstring=file+".ccs";
    
    FILE *mout,*ccout;
    
    mout=fopen(mstring.c_str(),"w");
	for (int i=0;i<n;i++) fprintf(mout,"%s%lf %d\n",C[i]<0 ? " " : "  ",C[i],i+1);
    fflush(mout);
    
    int count=n;
    
    ccout=fopen(ccstring.c_str(),"w");
    for (int i=0;i<n;i++) {
        
        for (int j=i+1;j<n;j++) {
            
            fprintf(ccout,"%s%lf %d %d\n",C[count]<0 ? " " : "  ",C[count],i+1,j+1);
            count++;
            
        }
        
    }
    fflush(ccout);
	
}




// Print the free entropy (i.e. log(Z) as a function of theta)
// theta - error on magnetizations - error on correlations - error on free entropy - size of max constructed cluster - size of the max sign cluster -num of all constructed clusters - num of significant clusters - time


void printFree(FILE *output, double theta, const std::vector<double> &error, Number finalS, int maxClusterSize, int max_sign_size, unsigned long numClusters, unsigned long numSignificantClusters, float time, double num_op) {
    
    fprintf(output,"%le\t%le\t%le\t%le\t%le\t%d\t%d\t%lu\t%lu\t%le\t%le\n",theta,error[0],error[1], error[2],finalS,maxClusterSize,max_sign_size,numClusters,numSignificantClusters, time, num_op);
    fflush(output);
	
}



