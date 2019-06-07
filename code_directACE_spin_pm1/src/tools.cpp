#include <vector>
#include <math.h>

// Numerical tools
#include "tools.h"


// Returns the inner product of two arrays of doubles

double innerProduct(double list1[], double list2[], int length) {

    double product=0;
    
    for (int i=0;i<length;i++) product += list1[i] * list2[i];
    
    return product;
    
}


// Computes the L1 norm of an array of doubles

double L1(double list[], int length) {
    
    double norm=0;
    
    for (int i=0;i<length;i++) norm += fabs(list[i]);
    
    return norm;
    
}


// Computes the L1 norm of a vector of doubles

double L1(const std::vector<double> &list) {
    
    double norm=0;
    
    for (int i=0;i<list.size();i++) norm += fabs(list[i]);
    
    return norm;
    
}


// Computes the L2 norm of an array of doubles

double L2(double list[], int length) {
    
    double norm=0;
    
    //for (int i=0;i<length;i++) norm += pow(list[i],2.0);
    for (int i=0;i<length;i++) norm += list[i]*list[i];

    
    return sqrt(norm);
    
}




// Computes the L2 norm of a vector of doubles

double L2(const std::vector<double> &list) {
    
    double norm=0;
    
    //for (int i=0;i<list.size();i++) norm += pow(list[i],2.0);
    for (int i=0;i<list.size();i++) norm += list[i]*list[i];
    
    return sqrt(norm);
    
}


// Computes the LInfinity norm of an array of doubles

double LInfinity(double list[], int length) {
    
    double norm=0;
    
    for (int i=0;i<length;i++) { if (fabs(list[i])>norm) norm = fabs(list[i]); }
    
    return norm;
    
}


// Computes the L1 norm of a vector of doubles

double LInfinity(const std::vector<double> &list) {
    
    double norm=0;
    
    for (int i=0;i<list.size();i++) { if (fabs(list[i])>norm) norm = fabs(list[i]); }
    
    return norm;
    
}


// Sorts the list of doubles (with corresponding index also arranged properly) using the quickSort algorithm
// Sort order is FROM HIGH TO LOW

void quickSort(double a[][2], int left, int right) {
    
    int i=left;
    int j=right;
    
    double temp0, temp1;
    double pivot = a[(left + right) / 2][0];
    
    // Partition
    
    while (i <= j) {
        
        while (a[i][0] > pivot) i++;
        while (a[j][0] < pivot) j--;
        
        if (i <= j) {
            
            temp0 = a[i][0];
            temp1 = a[i][1];
            
            a[i][0] = a[j][0];
            a[i][1] = a[j][1];
            
            a[j][0] = temp0;
            a[j][1] = temp1;
            
            i++;
            j--;

        }

    }

    // Recursion
    
    if (left < j)  quickSort(a, left, j);
    if (i < right) quickSort(a, i, right);

}


// Sorts the vector of integers (with corresponding index also arranged properly) using the quickSort algorithm
// Sort order is FROM LOW TO HIGH

void quickSort(std::vector<int> &a, int left, int right) {
    
    int i=left;
    int j=right;
    
    int temp;
    int pivot = a[(left + right) / 2];
    
    // Partition
    
    while (i <= j) {
        
        while (a[i] < pivot) i++;
        while (a[j] > pivot) j--;
        
        if (i <= j) {
            
            temp = a[i];
            a[i] = a[j];
            a[j] = temp;
            
            i++;
            j--;
            
        }
        
    }
    
    // Recursion
    
    if (left < j)  quickSort(a, left, j);
    if (i < right) quickSort(a, i, right);
    
}


// Check to see if two sets of spins of size k have a superset of size k+1

bool checkSize(const Key &a, const Key &b) {
    
    int i=0,j=0,error=0;
    
    while (error<3) {
        
        if      (i==a.size()) { error += b.size()-j; break; }
		else if (j==b.size()) { error += a.size()-i; break; }
        
		else if (a[i]<b[j])   { i++; error++; continue; }
		else if (a[i]>b[j])   { j++; error++; continue; }
        else if (a[i]==b[j])  { i++; j++;     continue; }
		else break;
        
	}
    
    if (error!=2) return false;
    else          return true;
    
}


// Returns the superset of two sets of spins

void getSuperset(const Key &a, const Key &b, Key &superset) {
	
	int i=0,j=0,index=0;
	
	while (true) {
		
		if      (i==a.size()) { for (int k=j;k<b.size();k++) { superset[index]=b[k]; index++; } break; }
		else if (j==b.size()) { for (int k=i;k<a.size();k++) { superset[index]=a[k]; index++; } break; }
        
		else if (a[i]==b[j]) { superset[index]=a[i]; index++; i++; j++; continue; }
		else if (a[i]<b[j])  { superset[index]=a[i]; index++; i++;      continue; }
        else if (a[i]>b[j])  { superset[index]=b[j]; index++;      j++; continue; }
		else break;
        
	}
	
}


// This function performs the modified Cholesky factorization on a positive-definite, but
// possibly near singular matrix and uses this to invert the matrix. The inversion is regularized
// to insure that the condition number of the matrix does not exceed a certain maximum value.

void modifiedCholeskyInverse(double m[], int length) {
    
    double c[length][length];
    double d[length];
    double delta=pow(10.0,-10.0);
    double beta =pow(10.0,  2.0);
    // Note: comp should be = beta * sqrt(delta)
    double comp =pow(10.0, -3.0);
    
    // First compute modified LDLT Cholesky factorization
    
    for (int i=0;i<length;i++) {
        
        c[i][i] = m[hindex(i,i,length)];
        
        //for (int j=0;j<i;j++) c[i][i] -= d[j] * pow(c[j][i],2);
        for (int j=0;j<i;j++) c[i][i] -= d[j] * c[j][i]*c[j][i];
        
        double max=comp;
        
        for (int j=i+1;j<length;j++) {
            
            c[i][j] = m[hindex(i,j,length)];
            
            for (int k=0;k<i;k++) c[i][j] -= d[k] * c[k][i] * c[k][j];
            
            if (fabs(c[i][j])>max) max=fabs(c[i][j]);
            
        }
        
        //max = pow(max/beta,2);
        max = max/beta*max/beta;
        
        if (fabs(c[i][i])>max) d[i] = fabs(c[i][i]);
        else                   d[i] = max;
        
        for (int j=i;j<length;j++) c[i][j] /= d[i];
            
    }
    
    // Now take the inverse
    
    double n[length][length];
    
    // First get D^-1 and L^-1
    
    for (int i=0;i<length;i++) {
        
        n[i][i] = 1;
        
        for (int j=i+1;j<length;j++) {
            
            n[i][j] = -c[i][j];
            
            for (int k=i+1;k<j;k++) n[i][j] -= n[i][k] * c[k][j];
            
        }
        
    }
    
    // Multiply out (L^-1)T D^-1 L^-1 to get the inverse
    
    for (int i=0;i<length;i++) {
        
        for (int j=i;j<length;j++) {
            
            m[hindex(i,j,length)]=0;
            
            for (int k=j;k<length;k++) m[hindex(i,j,length)] += n[i][k] * n[j][k] / d[k];
            
        }
        
    }
    
}


// Same as above, but also returns the determinant

double modifiedCholeskyInverseDet(double m[], int length) {
    
    double c[length][length];
    double d[length];
    double delta=pow(10.0,-10.0);
    double beta =pow(10.0,  2.0);
    // Note: comp should be = beta * sqrt(delta)
    double comp =pow(10.0, -3.0);
    
    // First compute modified LDLT Cholesky factorization
    
    for (int i=0;i<length;i++) {
        
        c[i][i] = m[hindex(i,i,length)];
        
        //for (int j=0;j<i;j++) c[i][i] -= d[j] * pow(c[j][i],2);
        for (int j=0;j<i;j++) c[i][i] -= d[j] * c[j][i]*c[j][i];
        
        double max=comp;
        
        for (int j=i+1;j<length;j++) {
            
            c[i][j] = m[hindex(i,j,length)];
            
            for (int k=0;k<i;k++) c[i][j] -= d[k] * c[k][i] * c[k][j];
            
            if (fabs(c[i][j])>max) max=fabs(c[i][j]);
            
        }
        
        //max = pow(max/beta,2);
        max = max/beta*max/beta;        
        if (fabs(c[i][i])>max) d[i] = fabs(c[i][i]);
        else                   d[i] = max;
        
        for (int j=i;j<length;j++) c[i][j] /= d[i];
        
    }
    
    // Compute the determinant
    
    double det=1;
    
    for (int i=0;i<length;i++) det *= d[i];
    
    // Now take the inverse
    
    double n[length][length];
    
    // First get D^-1 and L^-1
    
    for (int i=0;i<length;i++) {
        
        n[i][i] = 1;
        
        for (int j=i+1;j<length;j++) {
            
            n[i][j] = -c[i][j];
            
            for (int k=i+1;k<j;k++) n[i][j] -= n[i][k] * c[k][j];
            
        }
        
    }
    
    // Multiply out (L^-1)T D^-1 L^-1 to get the inverse
    
    for (int i=0;i<length;i++) {
        
        for (int j=i;j<length;j++) {
            
            m[hindex(i,j,length)]=0;
            
            for (int k=j;k<length;k++) m[hindex(i,j,length)] += n[i][k] * n[j][k] / d[k];
            
        }
        
    }
    
    return det;
    
}


// Returns the determinant of a symmetric matrix

double symmetricDet(double m[], int length) {
    
    double c[length][length];
    double d[length];
    double delta=pow(10.0,-10.0);
    double beta =pow(10.0,  2.0);
    // Note: comp should be = beta * sqrt(delta)
    double comp =pow(10.0, -3.0);
    
    // First compute modified LDLT Cholesky factorization
    
    for (int i=0;i<length;i++) {
        
        c[i][i] = m[hindex(i,i,length)];
        
        //for (int j=0;j<i;j++) c[i][i] -= d[j] * pow(c[j][i],2);
        for (int j=0;j<i;j++) c[i][i] -= d[j] * c[j][i]*c[j][i];        
        double max=comp;
        
        for (int j=i+1;j<length;j++) {
            
            c[i][j] = m[hindex(i,j,length)];
            
            for (int k=0;k<i;k++) c[i][j] -= d[k] * c[k][i] * c[k][j];
            
            if (fabs(c[i][j])>max) max=fabs(c[i][j]);
            
        }
        
        //max = pow(max/beta,2);
        max = max/beta*max/beta;
        
        if (fabs(c[i][i])>max) d[i] = fabs(c[i][i]);
        else                   d[i] = max;
        
        for (int j=i;j<length;j++) c[i][j] /= d[i];
        
    }
    
    // Compute the determinant
    
    double det=1;
    
    for (int i=0;i<length;i++) det *= d[i];
    
    return det;
    
}


// Given a vector x of length n, this function computes parameters for a Householder matrix H,
// such that H = 1 - beta v v^T is orthogonal and H x = L2(x) e_1. See Matrix Computations (Golub) for details.

void householder(double x[], int length, double v[], double &beta) {
    
    double sigma = 0;
    
    v[0] = 1;
    
    for (int i=1;i<length;i++) { sigma += pow(x[i],2); v[i] = x[i]; }
    
    if (sigma==0) beta=0;
    
    else {
        
        double mu = sqrt(pow(x[0],2) + sigma);
        
        if (x[0]<=0) v[0] = x[0] - mu;
        else         v[0] = -sigma / (x[0] + mu);
        
        beta = 2 * pow(v[0],2) / (sigma + pow(v[0],2));
        
        for (int i=length-1;i>=0;i--)  v[i] /= v[0];
        
    }
    
}


// Given a symmetric matrix A, this algorithm overwrites A with T = Q^T A Q, with T tridiagonal,
// and Q orthogonal. See Matrix Computations (Golub) for details. This algorithm assumes Q
// has already been initialized.

void tridiagonalize(double A[], int length, std::vector<std::vector<double> > &Q) {
    
    //for (int i=0;i<(length*(length+1))/2;i++) Q[i]=0;
    //    for (int i=0;i<length;i++)                Q[index(i,i,length)]=1;
    
    for (int k=0;k<length-2;k++) {
        
        // Compute the Householder matrix
        
        double v[length-(k+1)], a[length-(k+1)];
        double beta;
        
        for (int i=k+1;i<length;i++) a[i-(k+1)] = A[hindex(k,i,length)];
        
        householder(a,length-(k+1),v,beta);
        
        // Update A
        
        double p[length-(k+1)], w[length-(k+1)];
        double pv=0;
        
        for (int i=k+1;i<length;i++) {
            
            p[i-(k+1)]                             = 0;
            for (int j=k+1;j<i;j++)    p[i-(k+1)] += A[hindex(j,i,length)] * v[j-(k+1)];
            for (int j=i;j<length;j++) p[i-(k+1)] += A[hindex(i,j,length)] * v[j-(k+1)];
            p[i-(k+1)]                            *= beta;
            
            w[i-(k+1)]=p[i-(k+1)];
            
            pv += p[i-(k+1)] * v[i-(k+1)];
            
        }
        
        for (int i=0;i<length-(k+1);i++) w[i] -= beta * pv * v[i] / 2;
        
        A[hindex(k,k+1,length)]=L2(a,length-(k+1));
        
        for (int i=k+1;i<length;i++) {
            
            for (int j=i;j<length;j++)  A[hindex(i,j,length)] -= v[i-(k+1)] * w[j-(k+1)] + v[j-(k+1)] * w[i-(k+1)];
            
        }
        
        // Update Q
        
        double ww[length];
        
        for (int i=0;i<length;i++) {
            
            //            w[i-(k+1)]                               = Q[hindex(i,i,length)] * v[i-(k+1)];
            //            for (int j=i+1;j<length;j++) w[i-(k+1)] += Q[hindex(i,j,length)] * v[j-(k+1)];
            //            w[i-(k+1)]                              *= beta;
            
            //w[i-(k+1)]                             = Q[hindex(i,i,length)] * v[i-(k+1)];
            //            for (int j=k+1;j<i;j++)    w[i-(k+1)] += Q[hindex(j,i,length)] * v[j-(k+1)];
            //            for (int j=i;j<length;j++) w[i-(k+1)] += Q[hindex(i,j,length)] * v[j-(k+1)];
            //            w[i-(k+1)]                            *= beta;
            
            ww[i]                               = 0;
            for (int j=k+1;j<length;j++) ww[i] += Q[i][j] * v[j-(k+1)];
            ww[i]                              *= beta;
            
            
        }
        
        for (int i=0;i<length;i++) {
            
            for (int j=k+1;j<length;j++)  Q[i][j] -= ww[i] * v[j-(k+1)];
            
        }
        
        // Manually zero out entries on the kth row
        //DEBUG
        //for (int i=k+2;i<length;i++) A[hindex(k,i,length)]=0;
        
    }
    
    //DEBUG
    for (int i=0;i<length;i++) {
        
        for (int j=i+2;j<length;j++) A[hindex(i,j,length)]=0;
        
    }
    
}


// Given an unreduced symmetric tridiagonal matrix T, this algorithm overwrites T with Z^T T Z,
// where Z is a product of Givens rotations. See Matrix Computations (Golub) for details.

void QRStep(double T[], int length, std::vector<std::vector<double> > &Z) {
    
    double d = (T[hindex(length-2,length-2,length)] - T[hindex(length-1,length-1,length)]) / 2;
    
    double mu = T[hindex(length-1,length-1,length)] - pow(T[hindex(length-2,length-1,length)],2)
                     / (d * (1 + sqrt(pow(d,2) + pow(T[hindex(length-2,length-1,length)],2)) / fabs(d) ) );
    
    double x = T[0] - mu;
    double z = T[hindex(0,1,length)];
    
    for (int k=0;k<length-1;k++) {
        
        double c, s;
        givens(x, z, c, s);
        
        double G[4][4] = {{1, 0, 0, 0}, {0, c, s, 0}, {0, -s, c, 0}, {0, 0, 0, 1}};
        double newT[(length*(length+1))/2];
        
        // Multiply T by G on the right and G^T on the left
        
        // Start with i=k-1 --> m=k-1 
        
        if (k>0) {
        
            for (int j=k-1;j<k+2;j++) {
            
                                          newT[hindex(k-1,j,length)]  = 0;
                for (int n=k-1;n<k+2;n++) newT[hindex(k-1,j,length)] += T[hindex(k-1,n,length)] * G[n-(k-1)][j-(k-1)];
                
            }
                
        }
        
        // i=k --> m=k, k+1
        
        for (int j=k;j<k+3 && j<length;j++) {
            
                                                  newT[hindex(k,j,length)]  = 0;
            //if (k>0)                              newT[hindex(k,j,length)] +=  c * T[hindex(k-1,k,length)] * G[0][j-(k-1)]; 
            for (int n=k;n<k+3 && n<length;n++)   newT[hindex(k,j,length)] +=  c * T[hindex(k,n,length)]   * G[n-(k-1)][j-(k-1)];
                                                  newT[hindex(k,j,length)] += -s * T[hindex(k,k+1,length)] * G[1][j-(k-1)];
            for (int n=k+1;n<k+3 && n<length;n++) newT[hindex(k,j,length)] += -s * T[hindex(k+1,n,length)] * G[n-(k-1)][j-(k-1)];
            
        }
        
        // i=k+1 --> m=k, k+1
        
        for (int j=k+1;j<k+3 && j<length;j++) {
            
                                                  newT[hindex(k+1,j,length)]  = 0;
            //if (k>0)                              newT[hindex(k+1,j,length)] +=  s * T[hindex(k-1,k,length)] * G[0][j-(k-1)];
            for (int n=k;n<k+3 && n<length;n++)   newT[hindex(k+1,j,length)] += s * T[hindex(k,n,length)]   * G[n-(k-1)][j-(k-1)];
                                                  newT[hindex(k+1,j,length)] += c * T[hindex(k,k+1,length)] * G[1][j-(k-1)];
            for (int n=k+1;n<k+3 && n<length;n++) newT[hindex(k+1,j,length)] += c * T[hindex(k+1,n,length)] * G[n-(k-1)][j-(k-1)];
            
        }
        
        // Copy values to T
        
        if (k>0) { for (int j=k-1;j<k+2;j++)  T[hindex(k-1,j,length)] = newT[hindex(k-1,j,length)]; }
        for (int j=k;  j<k+3 && j<length;j++) T[hindex(k,j,length)]   = newT[hindex(k,j,length)];
        for (int j=k+1;j<k+3 && j<length;j++) T[hindex(k+1,j,length)] = newT[hindex(k+1,j,length)];
        
        // Update Z
        
        for (int i=0;i<length;i++) {
            
            double z1 = Z[i][k];
            double z2 = Z[i][k+1];
            
            Z[i][k]   = c * z1 - s * z2;
            Z[i][k+1] = s * z1 + c * z2;
            
        }
        
        // Adjust x and z
        
        if (k<length-2) {
            
            x = T[hindex(k,k+1,length)];
            z = T[hindex(k,k+2,length)];
            
        }
        
    }
    
}


// Given a symmetric matrix A, this algorithm computes its eigenvalues and eigenvectors
// via the (symmetric) QR algorithm. See Matrix Computations (Golub) for details. This function
// assumes eigenvectors is preset to the identity matrix before beginning the algorithm.

void symmetricQR(double A[], int length, double eigenvalues[], std::vector<std::vector<double> > &eigenvectors) {
    
    double eps = pow(10,-8);    // error tolerance
    //TEST
    //double min = pow(10,-10);    // minimum allowed d
    
    //DEBUG
    //printf("\nstart A = {");
    //for (int i=0;i<(length*(length+1))/2;i++) printf(" %lf",A[i]);
    //printf(" }\n");
    
    // Tridiagonalize the matrix
    
    tridiagonalize(A, length, eigenvectors);
    
    //DEBUG
    //printf("T = {");
    //for (int i=0;i<(length*(length+1))/2;i++) printf(" %lf",A[i]);
    //printf(" }\n");
    
    // Determine the unreduced subset
    
    int q=0;
    
    for (int i=length-2;i>=0 && fabs(A[hindex(i,i+1,length)])==0;i--) q++;
    
    int p=length-q-2;
    
    for (int i=p-1;i>=0 && fabs(A[hindex(i,i+1,length)])!=0;i--) p--;
    
    //DEBUG
    /*
    if (length==6) {
        printf("\nstart A = {");
        for (int i=0;i<(length*(length+1))/2;i++) printf(" %lf",A[i]);
        printf(" }\n");
        printf("p = %d, q= %d\n",p,q);
    }
    //DEBUG */
    
    // Loop until convergence
    
    //DBUG
    int loop=0;
    
    while (q<length-1 /*&& loop<10000*/) {
        
        //DEBUG
        loop++;
        
        // Zero small off-diagonal terms
        
        //for (int i=0;i<length-1;i++) { if (fabs(A[hindex(i,i+1,length)]) <= eps * (fabs(A[hindex(i,i,length)]) + fabs(A[hindex(i,i+1,length)]))) A[hindex(i,i+1,length)]=0; }
        
        // Apply the QR step to the unreduced part of the matrix between p and length-q-1 (inclusive)
        
        int n=length-q-p;
        
        double d = (A[hindex(p+n-2,p+n-2,length)] - A[hindex(p+n-1,p+n-1,length)]) / 2;
        double mu;
        
        if (d==0) mu = A[hindex(p+n-1,p+n-1,length)] - fabs(A[hindex(p+n-2,p+n-1,length)]);
        else      mu = A[hindex(p+n-1,p+n-1,length)] - pow(A[hindex(p+n-2,p+n-1,length)],2)
                       / (d + sign(d) * sqrt(pow(d,2) + pow(A[hindex(p+n-2,p+n-1,length)],2)));
        
        double x = A[hindex(p,p,length)] - mu;
        double z = A[hindex(p,p+1,length)];
        
        for (int k=0;k<n-1;k++) {
            
            double c, s;
            givens(x, z, c, s);
            
            //if (c==1) break;
            
            double G[4][4] = {{1, 0, 0, 0}, {0, c, s, 0}, {0, -s, c, 0}, {0, 0, 0, 1}};
            double newA[(n*(n+1))/2];
            
            // Multiply A by G on the right and G^T on the left
            
            // Start with i=k-1
            
            if (k>0) {
                
                for (int j=k-1;j<k+2;j++) {
                    
                    newA[hindex(k-1,j,n)]                            = 0;
                    for (int m=k-1;m<k+2;m++) newA[hindex(k-1,j,n)] += A[hindex(p+k-1,p+m,length)] * G[m-(k-1)][j-(k-1)];
                    
                }
                
            }
            
            // i=k
            
            for (int j=k;j<k+3 && j<n;j++) {
                
                newA[hindex(k,j,n)]                                   = 0;
                for (int m=k;m<k+3 && m<n;m++)   newA[hindex(k,j,n)] +=  c * A[hindex(p+k,p+m,length)]   * G[m-(k-1)][j-(k-1)];
                newA[hindex(k,j,n)]                                  += -s * A[hindex(p+k,p+k+1,length)] * G[1][j-(k-1)];
                for (int m=k+1;m<k+3 && m<n;m++) newA[hindex(k,j,n)] += -s * A[hindex(p+k+1,p+m,length)] * G[m-(k-1)][j-(k-1)];
                
            }
            
            // i=k+1
            
            for (int j=k+1;j<k+3 && j<n;j++) {
                
                newA[hindex(k+1,j,n)]                                   = 0;
                for (int m=k;m<k+3 && m<n;m++)   newA[hindex(k+1,j,n)] += s * A[hindex(p+k,p+m,length)]   * G[m-(k-1)][j-(k-1)];
                newA[hindex(k+1,j,n)]                                  += c * A[hindex(p+k,p+k+1,length)] * G[1][j-(k-1)];
                for (int m=k+1;m<k+3 && m<n;m++) newA[hindex(k+1,j,n)] += c * A[hindex(p+k+1,p+m,length)] * G[m-(k-1)][j-(k-1)];
                
            }
            
            // Copy values to A
            
            if (k>0) { for (int j=k-1;j<k+2;j++) A[hindex(p+k-1,p+j,length)] = newA[hindex(k-1,j,n)]; }
            for (int j=k;  j<k+3 && j<n;j++)     A[hindex(p+k,p+j,length)]   = newA[hindex(k,j,n)];
            for (int j=k+1;j<k+3 && j<n;j++)     A[hindex(p+k+1,p+j,length)] = newA[hindex(k+1,j,n)];
            
            // Update matrix of eigenvectors
            
            for (int i=0;i<length;i++) {
                
                double z1 = eigenvectors[i][p+k];
                double z2 = eigenvectors[i][p+k+1];
                
                eigenvectors[i][p+k]   = c * z1 - s * z2;
                eigenvectors[i][p+k+1] = s * z1 + c * z2;
                
            }
            
            // Adjust x and z
            //TEST
            if (k<(n-1)/* && (p+k+2)<length*/) {
                
                x = A[hindex(p+k,p+k+1,length)];
                z = A[hindex(p+k,p+k+2,length)];
                
            }
            
            //DEBUG
            /*
             if (length==6){
                 printf("loop: ");
                printf("\nd = %.16e, mu = %.16e\n",d,mu);
                printf("x = %.16e, z = %.16e\n",x,z);
                printf("c = %.16e, s = %.16e\n",c,s);
                printf("A = {\n");
                for (int i=0;i<(length*(length+1))/2;i++) printf(" %.8e",A[i]);
                printf(" }\n");
                printf("p = %d, q= %d\n",p,q);
                loop=0;
            }//*/
            
        }
        
        // Determine the remaining unreduced subset
        
        q=0;
        
        //for (int i=length-2;i>=0 && fabs(A[hindex(i,i+1,length)])==0;i--) q++;
        for (int i=length-2;i>=0 && fabs(A[hindex(i,i+1,length)]) <= eps * (fabs(A[hindex(i,i,length)]) + fabs(A[hindex(i,i+1,length)])); i--) {
        
            q++;
            A[hindex(i,i+1,length)]=0;
            
        }
        
        p=length-q-2;
        
        for (int i=p-1;i>=0 && fabs(A[hindex(i,i+1,length)])!=0;i--) {
            
            p--;
            
            if (fabs(A[hindex(i,i+1,length)]) <= eps * (fabs(A[hindex(i,i,length)]) + fabs(A[hindex(i,i+1,length)]))) {
                
                A[hindex(i,i+1,length)]=0;
                p++;
                break;
                
            }
            
        }
        
    }
    
    for (int i=0;i<length;i++) eigenvalues[i]=A[hindex(i,i,length)];
    
    //DEBUG
    //if (loop==10000) { for (int i=0;i<length;i++) eigenvalues[i]=0; }
    
}
