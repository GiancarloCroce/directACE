#include <vector>
// Algorithm for binary search - can do own implementation later
#include <algorithm>

// Numerical tools
#include "tools.h"

// Computations of cluster entropies
#include "direct.h"


using namespace std;


//compute log(Z) and correlation

void computeFandC(double hj[], int length, std::vector<double> &C,  double &F, double &b) {

    // Import couplings and fields
    std::vector<double> coupling(C.size());
    for (int i=0;i<coupling.size();i++) {
        coupling[i]=hj[i];
    }

    // print them
    //cout << "coupling: ";
    //for(int i=0; i<coupling.size(); i++) cout << coupling[i] << " ";
    //cout << endl; 

    std::vector<double> grad(C.size(),0);

    //compute log(Z), magnetization and correlation  
    optimizeF(coupling, F, grad, length, b);

    for (int i=0;i<C.size();i++) {
        C[i]=grad[i];
        //cout << "C[i] " << grad[i] << " "  << C[i] << endl;
    }

}

// modified for calculation using spin \pm1

void optimizeF(const std::vector<double> &J, double &F, std::vector<double> &grad, int length, double &b) {


    int spins[length];
    long double Z=0;

    // Calculate Z

    int i;
    for (i=0;i<length;i++) { spins[i]=0; grad[i]=0; }
    for (i=length;i<J.size();i++)     {  grad[i]=0;}

    int foundGrad[J.size()];
    double sign[J.size()];  //keep the sign for \pm 1 spin

    unsigned long twoToLength = (unsigned long) 1 << length;	
    for (i=0;i<twoToLength;i++) { 

        int lastGrad=0;
        long double temp=1;

        int j;
        for (j=0;j<length;j++) {
            //cout << "length " << length << " j: " << j<< endl;

            if (spins[j]) { //lo spin j è down
                foundGrad[lastGrad]=j;
                sign[lastGrad]=1;
                lastGrad++;

                temp *= exp(-b*J[j]); 
                int offj=offset(j,length);

                for (int k=j+1;k<length;k++) {

                    if (spins[k]) { //spin k is down

                        temp *= exp(-b*J[offj+k]);  
                        foundGrad[lastGrad]=offj + k;
                        sign[lastGrad]=1;
                        lastGrad++;

                    }
                    if (!spins[k]) { //spin k is up


                        temp *= exp(b*J[offj+k]);  ;  

                        foundGrad[lastGrad]=offj + k;
                        sign[lastGrad]=-1;
                        lastGrad++;

                    }          
                }

            }

            if (!spins[j]) { // spin j is up

                foundGrad[lastGrad]=j;
                sign[lastGrad]=-1;
                lastGrad++;
                temp *= exp(b*J[j]); 

                int offj=offset(j,length);

                for (int k=j+1;k<length;k++) {

                    if (spins[k]) { //spin k is down

                        temp *= exp(b*J[offj+k]);  
                        foundGrad[lastGrad]=offj + k;
                        sign[lastGrad]=-1;
                        lastGrad++;

                    }
                    if (!spins[k]) { //lo spin k is up

                        temp *= exp(-b*J[offj+k]); 
                        foundGrad[lastGrad]=offj + k;
                        sign[lastGrad]=1;
                        lastGrad++; 

                    }          
                }

            }

        }	
        //cout << "temp " << temp << endl;
        Z += temp; 

        /*
           cout << "spins: ";
           for(int k=0; k< length; k++) cout << spins[k] << " ";
           cout << endl;

           cout << "found: ";
           for (j=0;j<lastGrad;j++) cout <<  foundGrad[j]  << " ";
           cout << endl;
           cout  << "last grad " << lastGrad << endl;
           */

        for (j=0;j<lastGrad;j++){
            // grad[foundGrad[j]] += (2*spins[foundGrad[j]]-1)* temp;
            grad[foundGrad[j]] += sign[j]* temp;
        }

        for (j=0;j<length && spins[j];j++) spins[j]=0;
        if (j<length) spins[j]=1;

    }

    for (i=0;i<grad.size();i++) {grad[i] /= Z;}


    //compute connected correlations <sigma_i * sigma_j> - <sigma_i><sigma_j>

    for (i=0;i<length;i++) {
        int off=offset(i,length);
        for (int j=i+1;j<length;j++) {
	    grad[off+j] -= grad[i] * grad[j];
        }
    }

    //cout << "lungh " << length << endl;
    //for (int i=0;i<grad.size();i++) {    cout << "grad[i] " << grad[i] << endl;} 

    //cout << " Z " << Z << endl;  
    F = log(Z) ;	

}


// A reference entropy call that simply returns 0 for all couplings, fields, and the entropy

void computeF0andC0_Empty(double hj[], int length, std::vector<double> &C0, double &F0) {

    for (int i=0;i<C0.size();i++) C0[i]=0;
    F0=0;

}



//OLD STUFF FOR INVERSE PROBLEM
/*
   void optimizeFastF(const std::vector<double> &J, double &F, double grad[], int length) {

   int spins[length];
   int offJvect[length]; //compute it once for all
   double expJ[J.size()];

   double Z=1;  //start with the all spins=0 already counted
   int Zterms= (1 << length); //I don't like pow function :)
   int foundGrad[J.size()];

// Calculate Z
// During the computation keep track of elements to calculate the gradient and hessian


spins[0]=1;
expJ[0] = exp(J[0]); 
grad[0]=0;
offJvect[0]=length-1;


for (int i=1;i<length;i++) { spins[i]=0; expJ[i] = exp(J[i]); grad[i]=0; offJvect[i]=offset(i,length);} //the for starts from 1!!
for (int i=length;i<J.size();i++) {      expJ[i] = exp(J[i]); grad[i]=0; }

int firstSpin=0;  //the first non zero spin
int lastSpin =1;  //the last up spin is lastSpin-1

for (int i=1;i<Zterms;i++) {//the for starts from 1!!

// Compute contribution to the partition function

double temp=1;
int lastGrad=0;


for (int j=firstSpin;j<lastSpin;j++) {// the first spin=1 is in position firstSpin
if (spins[j]) {
temp *= expJ[j];
int offj = offJvect[j];
foundGrad[lastGrad]=j;
lastGrad++;

for (int k=j+1;k<lastSpin;k++) {
if (spins[k]){
temp *= expJ[offj + k];
foundGrad[lastGrad]=offj + k;
lastGrad++;

}
}
}
}
Z += temp;
for (int j=0;j<lastGrad;j++) grad[foundGrad[j]] += temp;


// Increment the counter

for(firstSpin=0; spins[ firstSpin ]; firstSpin++) spins[ firstSpin ]=0;
spins[ firstSpin ]=1;
if(firstSpin==lastSpin) lastSpin++;
}

// Normalize

for (int i=0;i<J.size();i++) { grad[i] /= Z; }

//connessi

for (int i=0;i<length;i++) {
int off=offset(i,length);
for (int j=i+1;j<length;j++) {
    grad[off+j] -= grad[i] * grad[j];
}
}


// Return the result

F = log(Z) ;  
//cout << "Z " << Z << " F " << F << endl; 
}


//Extra reference

void computeF0andC0_SingleExtra(double hj[], int length, double unused, std::vector<double> &C0, double mZero[], double &F0){
    C0[0] = mZero[0]; 
    F0 = - C0[0] * log(C0[0]) - (1.0-C0[0]) * log(1.0 - C0[0]) + C0[0] * hj[0];
}

void computeF0andC0_Extra(double hj[], int length, double unused, std::vector<double> &C0, double mZero[], double &F0){

    computeF0andC0loop(hj, length, false, C0, mZero, F0);

}




//Intra reference


void computeF0andC0_SingleIntra(double hj[], int length, double unused, std::vector<double> &C0, double mZero[], double &F0){\

    // C0[0] = 1.0/(1.0 + exp(- hj[0]) );
    F0= -log(1. - C0[0] );

    //C0[0] = mZero[0]; 
    //F0 = - C0[0] * log(C0[0]) - (1.0-C0[0]) * log(1.0 - C0[0]) + C0[0] * hj[0];
}

void computeF0andC0_Intra(double hj[], int length, double unused, std::vector<double> &C0, double useless[], double &F0){

    double mZero[length];

    computeM0(mZero,length,hj);

    computeF0andC0loop(hj, length, true, C0, mZero, F0);

}


//Reference Functions

void computeM0( double mZero[], int n, double finalJ[]){

    vector <int> off(n);
    vector <double> expMinHloc(n);
    vector < vector <double> > J(n);
    vector < vector <int> > neighbor(n);
    vector <double> mZeroOld(n);

    for(int i=0;i<n;i++) {
        mZero[i] = 1.0/(1.0 + exp(- finalJ[i]) );
        //mZero[i] = 0.5;
        off[i]=offset(i,n);
        expMinHloc[i] = exp(-finalJ[i]); 
    }

    for (int i=0;i<n;i++) {
        for (int j=i+1;j<n;j++) {
            if (fabs(finalJ[off[i] + j])>0) {
                neighbor[i].push_back(j);
                neighbor[j].push_back(i);

                J[i].push_back(  finalJ[off[i] + j]  );
                J[j].push_back(  finalJ[off[i] + j]  );                
            }            
        }
    }


    for(int i=0;i<n;i++){
        for(int j=0; j < neighbor[i].size(); j++){
            expMinHloc[i] *= exp( - J[i][j]*  mZero[ neighbor[i][j] ] );  			
        }
    }

    double check=10.0;
    while(check > 0.00001){
        check=0.0;
        for(int i=0;i<n;i++){
            mZeroOld[i] = mZero[i];
            mZero[i] = 1.0/(1.0 + expMinHloc[i]);
            check += fabs(mZero[i] -mZeroOld[i]);
        }
        for(int i=0;i<n;i++){
            for(int j=0; j < neighbor[i].size(); j++){
                expMinHloc[i] *= exp( - J[i][j] * ( mZero[ neighbor[i][j] ] -  mZeroOld[ neighbor[i][j] ] ) );  
            }
        }	
        //cout << check << "  " << mZero[0] << endl;

    }
}


void computeF0andC0loop(double hj[], int length, bool deltaP, std::vector<double> &C0, double mZero[], double &F0){

    vector < vector < double> > b(length);
    vector < vector < double> > bEigVect(length);	
    vector < vector < double> > bINV(length);
    vector <double> ll(length);
    vector <int> off(length);
    vector <double> bEigVal(length);
    int nRot;


    F0=1.0;

    for(int i=0; i< length; i++){
        off[i]=offset(i,length);
        b[i].resize(length);
        bINV[i].resize(length,0.0);
        bEigVect[i].resize(length,0.0);
        b[i][i]= 1.0/ mZero[i]/(1-mZero[i]);
        ll[i]= sqrt( mZero[i]*(1.0- mZero[i]) );
        F0 *=  mZero[i]*(1-mZero[i]);
        C0[i] = mZero[i];
    }

    F0 = - 0.5*log(F0);

    for(int i=0; i< length; i++){
        F0 += hj[i] * mZero[i];
        for(int j=i+1; j<length; j++){
            F0 +=  hj[ off[i] + j ]* mZero[i]*mZero[j];	

            b[j][i] = - hj[off[i]+j];
            b[i][j] = b[j][i];
        }
    }

    for(int i=0; i< length; i++) F0 += - mZero[i]* log(mZero[i]) - (1.0-mZero[i])* log(1.0-mZero[i]);

    jacobi(b,length,bEigVal,bEigVect,&nRot);  //si puo' ottimizzare???

    double detB(1.0);

    for(int i=0; i< length; i++){
        detB *= bEigVal[i];
        for(int j=i+1; j<length; j++){
            for(int k=0;k<length;k++){
                bINV[i][j] += bEigVect[i][k]*bEigVect[j][k]/ bEigVal[k];
            }
            //bINV[j][i] = bINV[i][j]; 
        }
    }

    F0 -= 0.5 * log(detB);

    for(int i=0; i< length; i++){
        for(int j=i+1; j<length; j++){
            C0[off[i]+j] = bINV[i][j]  ;
        }
    }

    if(deltaP){
        for(int i=0; i< length; i++){
            for(int k=0;k<length;k++){
                C0[i] += 0.5 * bINV[k][i] *(1- 2.0* mZero[k]) / mZero[k]/(1-mZero[k]) * ( bINV[k][k]/mZero[k]/(1.0-mZero[k]) -1.0   ) ;
            }
        }
    }



}
*/





/*
   void computeM0( vector <double> &mZero, int n, vector <double> &finalJ){

   vector <int> off(n);
   vector <double> expMinHloc(n);
   vector < vector <double> > J(n);
   vector < vector <int> > neighbor(n);
   vector <double> mZeroOld(n);

   for(int i=0;i<n;i++) {
   mZero[i] = 1.0/(1.0 + exp(- finalJ[i]) );
//mZero[i] = 0.5;
off[i]=offset(i,n);
expMinHloc[i] = exp(-finalJ[i]); 
}

for (int i=0;i<n;i++) {
for (int j=i+1;j<n;j++) {
if (fabs(finalJ[off[i] + j])>0) {
neighbor[i].push_back(j);
neighbor[j].push_back(i);

J[i].push_back(  finalJ[off[i] + j]  );
J[j].push_back(  finalJ[off[i] + j]  );                
}            
}
}


for(int i=0;i<n;i++){
for(int j=0; j < neighbor[i].size(); j++){
expMinHloc[i] *= exp( - J[i][j]*  mZero[ neighbor[i][j] ] );  			
}
}

double check=0.0;
while(check > 0.00001){
check=0.0;
for(int i=0;i<n;i++){
mZeroOld[i] = mZero[i];
mZero[i] = 1.0/(1.0 + expMinHloc[i]);
check += fabs(mZero[i] -mZeroOld[i]);
}
for(int i=0;i<n;i++){
for(int j=0; j < neighbor[i].size(); j++){
expMinHloc[i] *= exp( - J[i][j] * ( mZero[ neighbor[i][j] ] -  mZeroOld[ neighbor[i][j] ] ) );  
}
}	
}

}

*/




