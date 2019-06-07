#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>

// Algorithm declarations
#include "algorithm.h"

// Data input and output
#include "io.h"

// Numerical tools
#include "tools.h"

// Computations of cluster 
#include "direct.h"

using namespace std;

// GLOBAL VARIABLES
int N=0;      // total number of spins
int storageSize=8*sizeof(unsigned long); // number of spins which can be held in a single storage position of the key
int keySize=0;                           // number of entries necessary on the key to hold information on all spins


void (*computeF0andC0single_ptr)(double[], int,  std::vector<double> &, double &);
void (*computeF0andC0_ptr)(double[], int,  std::vector<double> &, double &);
void (*computeFandC_ptr)(double[], int,  std::vector<double> &,  double &, double &);

std::vector<Number> couplings;                                   // the set of all single and pair couplings
std::map<Key,Cluster> clusterIndex;                              // map from a key representing the spins in the cluster to the cluster
double oldF;                                          		       // previous free energy
std::vector<double> oldCorr;                                     // previous theta magnetization and correlations



// to be removed
void computeInfinityNorm(const std::vector<Number> &finalC,const std::vector<Number> &prevCorr,const double &finalF,const double &oldF, std::vector<double> &error,int N){

    error[0]=0.0;
    error[1]=0.0;
    error[2]=0.0;

    //prendo il max
    double delta=0;
    for(int n=0; n<N;n++){
        delta = fabs(finalC[n] - prevCorr[n]);
        if( delta > error[0]) error[0]=delta;
    }
    delta=0;
    for(int n=N; n<finalC.size();n++){
        delta = fabs(finalC[n] - prevCorr[n]);
        if( delta > error[1]) error[1]=delta;
    }
    delta=0;
    delta = fabs(finalF - oldF);
    error[2]=delta; 
}




// Calculate dF and dC for a given cluster

void computeDFandDC(int spins[], int clusterSize, std::vector<double> &dC, double &dF, bool useVerbose, bool lax, double &b) {

    int subset_mask[clusterSize];
    subset_mask[0]=1;
    for (int i=1;i<clusterSize;i++) subset_mask[i]=0;

    Key subsetKey(keySize,0);

    unsigned int twoToClusterSizeMinusTwo= (int(1) << clusterSize) - 2;


    for (unsigned int n=0;n<twoToClusterSizeMinusTwo;n++) {

        // Get the subset (if not found, make a new cluster)

        for (int i=0;i<keySize;i++) subsetKey[i]=0;

        int subsetSize=0;

        for (int i=0;i<clusterSize;i++) { if (subset_mask[i]) {

            subsetKey[spins[i]/storageSize] |= (unsigned long) 1 << (spins[i] % storageSize);
            subsetSize++;

        }} 



        // using the LAX rule some subclusters may have not been created: check it and create
        if(lax){ 
            if ((clusterIndex).count(subsetKey)==0) {
                Cluster c(subsetSize);
                makeCluster(c, subsetKey, subsetSize, false, lax,b);
                (clusterIndex)[subsetKey]=c;

            }
        }


        Cluster const &subCluster=clusterIndex[subsetKey];

        dF += subCluster.dF;

        // Map into the larger cluster and subtract the contribution to dC

        int dC_count=0;

        for (int i=0;i<clusterSize;i++) {
            if (subset_mask[i]){   		
                dC[i] += subCluster.dC[dC_count];  
                dC_count++;; 			
            }
        }


        for (int i=0;i<clusterSize;i++) {	
            if (subset_mask[i]){	  
                int off=offset(i,clusterSize);	  
                for (int j=i+1;j<clusterSize;j++) {	    
                    if (subset_mask[j]){ 
                        dC[off + j] += subCluster.dC[dC_count];
                        dC_count++;	      
                    }	    
                }	  
            }	
        }

        int  i;
        for (i=0;subset_mask[i];i++) { 
            subset_mask[i]=0;
            flipBit( subsetKey[spins[i]/storageSize] , spins[i]%storageSize ) ;
        }
        subset_mask[i]=1;
        flipBit( subsetKey[spins[i]/storageSize] , spins[i]%storageSize ) ;

    }

}

/*
 *** USELESS ****

// function to create cluster made by a single spin
void makeSingleCluster(Cluster &cluster, const Key &key, int clusterSize,bool useVerbose, bool lax, double &b) {
// Record individual spins in the key

int spins[clusterSize];
int n=0;

for (int i=0;i<keySize && n<clusterSize;i++) { for (int j=0;j<storageSize && n<clusterSize;j++) {

if (key[i] & (unsigned long) 1<<j) { spins[n] = storageSize * i + j; n++; }

} }

double hj[cluster.dC.size()];

// Get the set of couplings for this cluster

for (int i=0;i<clusterSize;i++) {
hj[i]=couplings[spins[i]];
//cout << hj[i] << endl;
}

std::vector<double> dC(cluster.dC.size());
double dF=0;

double Z=0;
Z=1+ exp(-b*hj[0]);

cluster.dF = log(Z);
cluster.dC = Vector(1, exp(-b*hj[0])/Z);
}
*/



// Makes a new cluster

void makeCluster(Cluster &cluster, const Key &key, int clusterSize, bool useVerbose, bool lax, double &b) {

    // Record the individual spins in the key

    int spins[clusterSize];
    int n=0;

    // read the key and check which are the corresponding spins
    for (int i=0;i<keySize && n<clusterSize;i++) { 
        for (int j=0;j<storageSize && n<clusterSize;j++) { 
            if (key[i] & (unsigned long) 1<<j) {
                spins[n] = storageSize * i + j;
                n++;
            }
        }
    }


    double hj[cluster.dC.size()];

    // Get the set of couplings for this cluster

    for (int i=0;i<clusterSize;i++) {
        hj[i]=couplings[spins[i]];
        //cout << spins[i] << "  " << couplings[spins[i]] << endl;
    }

    for (int i=0;i<clusterSize-1;i++) {         
        int off=offset(spins[i],N);
        for (int j=i+1;j<clusterSize;j++) {
            hj[n]=couplings[off + spins[j]];
            n++;
        }
    }

    // Compute values and assign them to the cluster

    std::vector<double> dC(cluster.dC.size()); double dF;

    (*computeF0andC0_ptr)(hj,clusterSize,dC, dF);

    // compute dF and dC of the subclusters
    if (clusterSize>1) computeDFandDC(spins,clusterSize,dC,dF,useVerbose, lax,b); 

    std::vector<double> C(dC.size());
    double F;

    (*computeFandC_ptr)(hj,clusterSize,C,F,b); 
    for (int i=0;i<cluster.dC.size();i++) {cluster.dC[i]=C[i]-dC[i];}

    cluster.dF=F-dF;

    //if dF is zero
    if(fabs(cluster.dF)<0.000000001)  cluster.dF=0; 

    // Print the key and the spin of the new clusters

    if(useVerbose){
        if(fabs(cluster.dF)>0){
            cout << "Key:( ";
            for(int m=0;m<keySize; m++) cout  << key[m]  <<"-";
            cout << ")   spins:";
            for( int k=0; k<clusterSize; k++) cout << spins[k]<< "-";
            cout << "   " << clusterSize;
            // cout << " ->"; 
            cout << "   F:" << F << " dF:" << dF <<" cluster.dF " << cluster.dF << endl;

        }
    }
}





//************** STRICT selection rule ********************* 
// Runs the cluster selection and creation algorithm

void selectClusters_strict(std::set<Key> &clusters, double theta, int clusterSize, bool useVerbose, std::set<Key> &Non_sign, int kmin) {

    std::vector<Key> significantClusters;
    significantClusters.reserve(clusters.size());

    Significant isSignificant(theta, kmin);
    if(useVerbose) cout << " num_non_significant: "; 

    int cont_non_sign=0;
    for (std::set<Key>::iterator i=clusters.begin();i!=clusters.end();++i) {

        if ( isSignificant((clusterIndex)[*i])) significantClusters.push_back(*i);
        else {  
            Non_sign.insert(Non_sign.end(),(*i));
            int nonSignSize=(sqrt(1 + 8 * ((clusterIndex)[*i]).dC.size()) - 1) / 2;
            if (clusterSize==nonSignSize) cont_non_sign+=1;
        }
    }
    if(useVerbose){  cout << cont_non_sign<< endl;}    

    clusters.clear();



    // Build list of potential new cluster elements

    std::map<int,int> spinMap;
    std::set<int> spinSet;

    std::map<std::vector<int>,int> pairMap;
    std::set<std::vector<int> > pairSet;
    std::vector<int> pair(2,0);

    for (int i=0;i<significantClusters.size();i++) {

        // Extract list of spins

        int spins[clusterSize];
        int n=0;

        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {

            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }

        } }
        /*	
            if(useVerbose){
            cout << "  Key: ";
            for(int h=0; h<keySize; h++) cout << significantClusters[i][h] << " ";
            cout << " spins ";
            for(int i=0; i<clusterSize; i++) cout << spins[i] << " ";
            cout << endl;
            }
            */

        // Count singles and pairs

        for (int j=0;j<clusterSize;j++) { 

            if (spinMap.count(spins[j])==0) spinMap[spins[j]]=1;
            else                            spinMap[spins[j]]+=1;

            pair[0]=spins[j];

            for (int k=j+1;k<clusterSize;k++) { 
                pair[1]=spins[k];
                //cout << "  pair " << pair[0] << " " << pair[1] ;
                if (pairMap.count(pair)==0) pairMap[pair]=1;
                else                        pairMap[pair]+=1;

                //cout << " ->" << pairMap[pair] << endl; 

            }


        }
    }      

    // Select potential singles

    int spinCutSize=clusterSize-1;
    // if (clusterSize==2) spinCutSize=0;

    for (std::map<int,int>::iterator i=spinMap.begin();i!=spinMap.end();++i) {

        if ((*i).second>spinCutSize) spinSet.insert(spinSet.end(),(*i).first); 

    }

    //cout << "           SpinSet ";
    // for(std::set<int>::iterator i=spinSet.begin(); i!= spinSet.end(); ++i) cout << *i << " ";
    // cout << endl;

    spinMap.clear();

    // Select potential pairs

    int pairCutSize=clusterSize-2;

    for (std::map<std::vector<int>,int>::iterator i=pairMap.begin();i!=pairMap.end();++i) {

        if ((*i).second>pairCutSize) pairSet.insert(pairSet.end(),(*i).first);

    }


    //cout << "pairSet.size() " << pairSet.size() << endl;
    pairMap.clear();


    // Find supersets of significant clusters and reassign them to clusters

    Key supersetKey(keySize,0);

    std::map<Key,int> buildNum;

    for (int i=0;i<significantClusters.size();i++) {

        // Extract list of spins, skip if single or pair is absent

        int spins[clusterSize];
        int n=0;
        bool halt=false;

        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {

            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }

        } }

        for (int j=0;j<clusterSize && !halt;j++) {              

            if (spinSet.count(spins[j])==0) { halt=true; break; }

            pair[0]=spins[j];

            for (int k=j+1;k<clusterSize && !halt;k++) {

                pair[1]=spins[k];
                if (pairSet.count(pair)==0) { halt=true; break; }

            }

        }

        if (halt) continue;      

        // Add potential new clusters      
        for (int j=0;j<keySize;j++) supersetKey[j]=significantClusters[i][j];
        for (std::set<int>::iterator j=spinSet.begin();j!=spinSet.end();++j) { 

            // If cluster is unchanged, skip. Else the spin considered is not in the cluster 
            if (supersetKey[(*j)/storageSize] & (unsigned long) 1 << ((*j) % storageSize) )  continue;  

            // Test for potential cluster

            bool isPossible=true;
            for (int k=0;k<clusterSize;k++) { 
                if (spins[k]<(*j)) { pair[0]=spins[k]; pair[1]=(*j);     }
                else               { pair[0]=(*j);     pair[1]=spins[k]; }
                // cout << pair[0] << " " << pair[1] << endl;
                if (pairSet.count(pair)==0) { isPossible=false; break; }

            }

            // Add to count

            if (isPossible){

                supersetKey[(*j)/storageSize] |= (unsigned long) 1 << ((*j) % storageSize);

                if (buildNum.count(supersetKey)==0) buildNum[supersetKey]=1;
                else                                buildNum[supersetKey]+=1;


                supersetKey[(*j)/storageSize]  = significantClusters[i][(*j)/storageSize];
            }
        }

    }


    int cutSize=clusterSize;
    //if (clusterSize==2) cutSize=1;
    for (std::map<Key,int>::iterator i=buildNum.begin();i!=buildNum.end();++i) {
        if ((*i).second>cutSize) clusters.insert(clusters.end(),(*i).first);
    }
}


//*************** Regola di selezione LAX ***********

void selectClusters_lax(std::set<Key> &clusters,  double theta,int clusterSize, bool useVerbose,  std::set<Key> &Non_sign , int kmin) {


    std::vector<Key> significantClusters;
    significantClusters.reserve(clusters.size());

    Significant isSignificant(theta, kmin );
    if(useVerbose) cout << " Non_significant: "; 
    for (std::set<Key>::iterator i=clusters.begin();i!=clusters.end();++i) {

        if ( isSignificant((clusterIndex)[*i])) significantClusters.push_back(*i);
        else { /*  if(useVerbose) { for(int k=0; k<keySize; k++) cout << " " << (*i)[k] << " "<<"dF:"<<clusterIndex[*i].dF;
                   cout << " -";} */ 
            Non_sign.insert(Non_sign.end(),(*i));
        }
    }
    if(useVerbose){  cout << Non_sign.size()  << endl;}    

    clusters.clear(); 

    // Build list of potential new cluster elements

    std::map<int,int> spinMap;
    std::set<int> spinSet;

    for (int i=0;i<significantClusters.size();i++) {

        // Extract list of spins

        int spins[clusterSize];
        int n=0;

        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {

            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }

        } }



        //if(useVerbose){
        //cout << "  Key: ";
        //for(int h=0; h<keySize; h++) cout << significantClusters[i][h] << " ";
        //cout << " spins ";
        //for(int i=0; i<clusterSize; i++) cout << spins[i] << " ";
        //cout << endl;
        //}

        // Select potential singles, count singles

        for (int j=0;j<clusterSize;j++) {

            if (spinMap.count(spins[j])==0) spinMap[spins[j]]=1;
            else                            spinMap[spins[j]]+=1;
        }
    }

    int spinCutSize=0;  //2;
    //if (clusterSize==2) spinCutSize=1;

    for (std::map<int,int>::iterator i=spinMap.begin();i!=spinMap.end();++i) {

        if ((*i).second>spinCutSize) spinSet.insert(spinSet.end(),(*i).first);

    }
    /*
       cout << "   spinset: ";
       for(std::set<int>::iterator i= spinSet.begin(); i!=spinSet.end(); i++) {
       cout << *i <<" ";}
       cout << endl;
       */

    spinMap.clear();

    for (int i=0;i<significantClusters.size();i++) {

        // Extract list of spins, skip if single or pair is absent

        int spins[clusterSize];
        int n=0;
        bool halt=false;

        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {

            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }

        } }

        for (int j=0;j<clusterSize && !halt;j++) {

            if (spinSet.count(spins[j])==0) { halt=true; break; }
        }
    }

    // Find supersets of significant clusters and reassign them to clusters

    Key supersetKey(keySize,0);

    std::map<Key,int> buildNum;

    for (int i=0;i<significantClusters.size();i++) {

        // Extract list of spins, skip if single or pair is absent

        int spins[clusterSize];
        int n=0;
        bool halt=false;

        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {

            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }

        } }

        for (int j=0;j<clusterSize && !halt;j++) {// if spin not in spinSet, skip 

            if (spinSet.count(spins[j])==0) { halt=true; break;}
        }

        if (halt) continue;

        // Add potential new clusters

        for (int j=0;j<keySize;j++) supersetKey[j]=significantClusters[i][j];

        for (std::set<int>::iterator j=spinSet.begin();j!=spinSet.end();++j) {

            // If cluster is unchanged, skip

            if (supersetKey[(*j)/storageSize] & (unsigned long) 1 << ((*j) % storageSize)) continue;

            // Add to count

            supersetKey[(*j)/storageSize] |= (unsigned long) 1 << ((*j) % storageSize);

            if (buildNum.count(supersetKey)==0) buildNum[supersetKey]=1;
            else                                buildNum[supersetKey]+=1;

            supersetKey[(*j)/storageSize]  = significantClusters[i][(*j)/storageSize];

        }

    }

    int cutSize=1;

    for (std::map<Key,int>::iterator i=buildNum.begin();i!=buildNum.end();++i) {

        if ((*i).second>cutSize) clusters.insert(clusters.end(),(*i).first);

    }

}




// Run the program

void run(RunParameters &r) {
    printf("\n");
    // Retrieve couplings from file and set system, key sizes

    FILE *datain=fopen(r.getInfile().c_str(),"r");

    if (datain!=NULL) {

        getCouplings(datain, couplings);
        //for(int n=0;n<couplings.size();n++) cout << couplings[n] << "  ";

        N = (sqrt(1 + 8 * couplings.size()) - 1) / 2;
        keySize = (N + storageSize - 1) / storageSize;
    }

    else printf("Problem retrieving data from file: %s",r.getInfile().c_str());

    cout << "Number of spin: "<< N << endl;
    cout << "Inverse temperature: " << r.beta << endl;
    //random starting point for free energy and correlation
    oldF=double(rand())/RAND_MAX; 
    oldCorr.resize(couplings.size(),0.0);
    for(int n=0;n<couplings.size();n++) oldCorr[n]=double(rand())/RAND_MAX; 

    // Open supplementary output file, if used
    FILE *supout=fopen(r.getSupplementaryOutfile().c_str(),"w");

    // Set kmax to the system size, if no specific maximum size is given 

    if (r.kmax<=0) r.kmax=N;

    cout << "Max cluster length: " << r.kmax << endl;

    // See kmin to the system size, if no specific minimum size is given 

    if (r.kmin<=0) r.kmin=1;

    cout << "Minimum cluster size: " << r.kmin << endl;

    // Sanity checks on theta and gamma input

    if (r.thetaMin<0) { RunParameters rDefault; r.thetaMin  = rDefault.thetaMin;  }

    computeFandC_ptr=&computeFandC; 
    computeF0andC0single_ptr=&computeF0andC0_Empty;
    computeF0andC0_ptr=&computeF0andC0_Empty;

    printf("Computing Ising model correlations using theta = %.8e \n",r.theta);
    if(r.lax) cout <<"LAX selection rule"<< endl;
    else cout <<"STRICT selection rule"<< endl;
    printf("Storage size = %d, key size = %d\n",storageSize,keySize);


    // MAIN ALGORITHM

    // list of non significant clusters
    std::set<Key> Non_sign;

    // Clusters of size one 
    int maxClusterSize=1; 

    for (int i=0;i<N;i++) {

        Key clusterKey(keySize,0);
        clusterKey[i/storageSize] |= (unsigned long) 1 << (i % storageSize);

        Cluster c(1);
        makeCluster(c,clusterKey,maxClusterSize,r.useVerbose, r.lax,r.beta);
        clusterIndex[clusterKey]=c;

    }

    if (r.useVerbose) printf("Computing all clusters of size 1. Found %d cluster(s).\n", N);

    std::set<Key> clusters;

    for (int i=0;i<N;i++) {
        for (int j=i+1;j<N;j++) {	
            Key clusterKey(keySize,0);
            clusterKey[i/storageSize] |= (unsigned long) 1 << (i % storageSize);
            clusterKey[j/storageSize] |= (unsigned long) 1 << (j % storageSize);        
            clusters.insert(clusterKey);
        }
    }


    // Continue to general algorithm
    while ((clusters.size()>0 && maxClusterSize!=r.kmax)  ) {

        maxClusterSize++;

        for (std::set<Key>::iterator i=clusters.begin();i!=clusters.end();++i) {

            Cluster c(maxClusterSize);
            makeCluster(c,*i,maxClusterSize, r.useVerbose, r.lax,r.beta);
            clusterIndex[*i]=c;

            //if shorter then kmin then use them for construction
            //if(maxClusterSize < r.kmin) {
            //   significantClusters.push_back(*i);
            //}
        }

        if (r.lax) selectClusters_lax(clusters,r.theta,maxClusterSize, r.useVerbose, Non_sign, r.kmin );
        else     selectClusters_strict(clusters,r.theta,maxClusterSize, r.useVerbose, Non_sign, r.kmin );

    }


    // GET FINAL F, C

    double *hj=new double[couplings.size()];
    for (int i=0;i<couplings.size();i++) hj[i]=couplings[i];

    std::vector<Number> finalC0(couplings.size(),0.0);
    Number finalF0=0.0;

    (*computeF0andC0_ptr)(hj,N,finalC0,  finalF0);

    delete[] hj;

    Vector finalC(finalC0.size());
    Number finalF=finalF0;


    unsigned long numSignificantClusters=0;
    unsigned long numClusters=clusterIndex.size();


    for (std::map<Key,Cluster>::iterator i=clusterIndex.begin();i!=clusterIndex.end();++i) {

        Significant isSignificant(r.theta, r.kmin );

        if ( isSignificant((*i).second) ) {

            // Count the number of significant clusters

            numSignificantClusters++;

            // Add the cluster's contribution to dF

            finalF += (*i).second.dF;

            // Map the cluster to the whole system and add the contribution to dC

            int clusterSize=(sqrt(1 + 8 * (*i).second.dC.size()) - 1) / 2;
            int spins[clusterSize];

            int n=0;


            for (int j=0;j<keySize && n<clusterSize;j++) { 
                for (int k=0;k<storageSize && n<clusterSize;k++) { 
                    if ((*i).first[j] & (unsigned long) 1<<k) {
                        spins[n] = storageSize * j + k;
                        n++;                        
                    }
                }
            }


            for (int j=0;j<clusterSize;j++) { //prova
                //cout << "prova " <<(*i).second.dC[spins[j]] << endl; 
                finalC[spins[j]] += (*i).second.dC[j];
                //cout << "finalC["<<spins[j]<<"] " << finalC[spins[j]] << endl;
            }

            for (int j=0;j<clusterSize-1;j++) {
                int off=offset(spins[j], N);
                for (int k=j+1;k<clusterSize;k++) { 
                    finalC[off + spins[k]] += (*i).second.dC[n];
                    n++;
                }
            }

        }    
    }


    std::vector<double> error(3,0.0);
    computeInfinityNorm(finalC,oldCorr,finalF,oldF,error,N);   //  compiute the distance (infinity norm) between old and new F
    oldF=finalF;
    for(int n=0;n<oldCorr.size();n++) oldCorr[n]=finalC[n];


    unsigned long lastNumSignificantClusters=numSignificantClusters;

    //max dF of non significant clusters
    double max=0;  
    double temp=0;
    for (std::set<Key>::iterator i=Non_sign.begin();i!=Non_sign.end();++i){
        temp=fabs(clusterIndex[*i].dF);
        if(temp > max) max=temp;
    }
    if(r.useVerbose)     {cout <<"Number of non significant: " << Non_sign.size() << " |max_dF| " << max << endl;}



    int max_sign_size=1;
    double num_op=0.0;

    //MAIN WHILE CYCLE
    while ( ( (numSignificantClusters <= (((unsigned long)(1) << N )-1)) &&  maxClusterSize<=r.kmax   && r.theta> r.thetaMin)  ) {

        Non_sign.clear();

        //compute the num_of_operation
    for (std::map<Key,Cluster>::iterator i=clusterIndex.begin();i!=clusterIndex.end();++i) {
            int clusterSize=(sqrt(1 + 8 * (*i).second.dC.size()) - 1) / 2;
            num_op += clusterSize<<2 ;
        }
   


        // First record data
        printFree(supout,r.theta,error,finalF,maxClusterSize,max_sign_size,numClusters,numSignificantClusters, float(clock())/CLOCKS_PER_SEC, num_op);
        printCorrelation(r.getCouplingsOutfile().c_str(),finalC);

        // Then lower threshold in such a way at each step new clusters are included
        r.theta=0.99*max;  

        //CHECK HERE!!!
        //if(r.theta==0) break; //se il |max| è 0 -> esci dal while




        if(r.useVerbose){ cout << endl;
            cout << "***** Theta=" << r.theta << endl;}
        // Rerun steps above but do not compute cluster if it already exists

        maxClusterSize=1; 

        clusters.clear();

        for (int i=0;i<N;i++) {
            for (int j=i+1;j<N;j++) {
                Key clusterKey(keySize,0);
                clusterKey[i/storageSize] |= (unsigned long) 1 << (i % storageSize);
                clusterKey[j/storageSize] |= (unsigned long) 1 << (j % storageSize);
                clusters.insert(clusterKey);
            }        
        }

        // Continue to general algorithm    

        while ((clusters.size()>0 && maxClusterSize!=r.kmax)  ) {
            maxClusterSize++;

            if(r.useVerbose) { cout << "Cluster of lenght " << maxClusterSize << ", num_created: " << clusters.size();
            }
            for (std::set<Key>::iterator i=clusters.begin();i!=clusters.end();++i) { if (clusterIndex.count(*i)==0) {
                Cluster c(maxClusterSize);
                makeCluster(c,*i,maxClusterSize, r.useVerbose, r.lax, r.beta);
                clusterIndex[*i]=c;
            }
            }
            //only if i have already reached kmin I can select, otherwise construct all clusters
            if (r.lax) selectClusters_lax(clusters,r.theta,maxClusterSize, r.useVerbose, Non_sign, r.kmin );
            else     selectClusters_strict(clusters,r.theta,maxClusterSize, r.useVerbose, Non_sign, r.kmin);            
        }

        // Get final Free energy and correlations 

        finalF=finalF0;
        for (int i=0;i<finalC.size();i++) finalC[i]=finalC0[i];

        numSignificantClusters=0;
        numClusters=clusterIndex.size();


        // Length max of constructed clusters
        max_sign_size=0;

        for (std::map<Key,Cluster>::iterator i=clusterIndex.begin();i!=clusterIndex.end();++i) {
            Significant isSignificant(r.theta, r.kmin);
            if ( isSignificant((*i).second)){

                numSignificantClusters++;  // Count the number of significant clusters


                // Add the cluster's contribution to dF
                finalF += (*i).second.dF;

                // Map the cluster to the whole system and add the contribution to dC

                int clusterSize=(sqrt(1 + 8 * (*i).second.dC.size()) - 1) / 2;
                int spins[clusterSize];
                int n=0;


                for (int j=0;j<keySize && n<clusterSize;j++) { 
                    for (int k=0;k<storageSize && n<clusterSize;k++) { 
                        if ((*i).first[j] & (unsigned long) 1<<k) {
                            spins[n] = storageSize * j + k;
                            n++;                
                        }
                    }
                }

                for (int j=0;j<clusterSize;j++) {finalC[spins[j]] += (*i).second.dC[j];}// cout << (*i).first[0] << "->"<< (*i).second.dC[j] << " ";}
            //	cout << endl;

            for (int j=0;j<clusterSize-1;j++) {     
                int off=offset(spins[j], N);
                for (int k=j+1;k<clusterSize;k++) {      
                    finalC[off + spins[k]] += (*i).second.dC[n];
                    n++;               
                }            
            }    

            //max size of significat clusters
            if(clusterSize>max_sign_size) max_sign_size=clusterSize;	
        }        
    }



    // compute max dF of NON significant clusters
    max=0;
    double temp=0;
    for (std::set<Key>::iterator i=Non_sign.begin();i!=Non_sign.end();++i){
        temp=fabs(clusterIndex[*i].dF);
        if(temp > max) max=temp;
    }
    if(r.useVerbose)     {cout <<"Tot_num_of_non_significant: " << Non_sign.size() << " |max_dF| " << max << endl;}


    if (numSignificantClusters!=lastNumSignificantClusters){   
        computeInfinityNorm(finalC,oldCorr,finalF,oldF,error,N);
        oldF=finalF;
        for(int n=0;n<oldCorr.size();n++) oldCorr[n]=finalC[n];
    }

    lastNumSignificantClusters=numSignificantClusters;

}

// Record final data and exit
printCorrelation(r.getCouplingsOutfile().c_str(),finalC);
printFree(supout,r.theta,error,finalF,maxClusterSize,max_sign_size,numClusters,numSignificantClusters, float(clock())/CLOCKS_PER_SEC, num_op);

}



/*********************************************************************

  COMMAND LINE INPUT FORMAT

  Command line instructions tell the program where to look for input
  files and where to send output, as well as the setting of various
  parameters (theta, etc) and flags (verbose, etc).

  -d, -dir, -directory: string
  Path to the directory where the data file is located, and where
  output will be written.

  -i, -in, -input: string
  The location of the file containing a set of couplings

  -o, -out, -output: string
  The location of the file where output is to be sent.

  -lax: none
  Use LAX selection rule

  -k, -kmax: integer
  Maximum cluster size.

  -t, -theta: real number
  Stariting point for the cutoff theta .

  -tmin, -thetamin: real number
  The minimum value of the cutoff, when the algorithm is to loop 
  over a range of possible cutoff values.

  -v, -verbose: 

 *********************************************************************/

// MAIN PROGRAM

int main(int argc, char *argv[]) {

    RunParameters r;

    // Process command line input

    for (int i=1;i<argc;i++) {

        if (strcmp(argv[i],"-d")==0 || strcmp(argv[i],"-dir")==0 || strcmp(argv[i],"-directory")==0)   { if (++i==argc) break; else r.directory=argv[i]; }
        else if (strcmp(argv[i],"-i")==0 || strcmp(argv[i],"-in")==0 || strcmp(argv[i],"-input")==0)   { if (++i==argc) break; else r.infile=argv[i];    }
        else if (strcmp(argv[i],"-o")==0 || strcmp(argv[i],"-out")==0 || strcmp(argv[i],"-output")==0) { if (++i==argc) break; else r.outfile=argv[i];   }

        else if (strcmp(argv[i],"-t")==0 || strcmp(argv[i],"-theta")==0)        { if (++i==argc) break; else r.theta=strtodouble(argv[i]);      }
        else if (strcmp(argv[i],"-tmin")==0 || strcmp(argv[i],"-thetamin")==0)  { if (++i==argc) break; else r.thetaMin=strtodouble(argv[i]);   }

        else if (strcmp(argv[i],"-b")==0  || strcmp(argv[i],"-beta")==0)       {if (++i==argc) break; else r.beta=strtodouble(argv[i]);         }  

        else if (strcmp(argv[i],"-kmax")==0)   {if (++i==argc) break; else  r.kmax=strtodouble(argv[i]);      }
        else if (strcmp(argv[i],"-kmin")==0)   {if (++i==argc) break; else  r.kmin=strtodouble(argv[i]);      }

        else if (strcmp(argv[i],"-lax")==0)   {  r.lax=true; }
        else if (strcmp(argv[i],"-v")==0  || strcmp(argv[i],"-verbose")==0) {  r.useVerbose=true;             }  

        else printf("Unrecognized command! '%s'\n",argv[i]);

    }

    run(r);

    return 0;

}

