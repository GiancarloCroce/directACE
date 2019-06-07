#ifndef DIRECT_H
#define DIRECT_H


// Typedefs
//typedef double Number;
//typedef std::vector<float> Vector;



void computeFandC(double[], int,  std::vector<double> &, double &, double &);
void optimizeF(const std::vector<double> &, double &, std::vector<double> &, int, double &);
void optimizeFastF(const std::vector<double> &, double &, double[], int);

/*
//Mean Field
void computeM0( double[], int, double[]);
void computeF0andC0loop(double[], int, bool, std::vector<double> &, double[], double &);
//void computeP0loop( double[], int, double[]){}
*/


//No Reference
void computeF0andC0_Empty(double[], int, std::vector<double> &, double &);

/*
//Extra Reference
void computeF0andC0_Extra(double[], int, double, std::vector<double> &,double[], double &);
void computeF0andC0_SingleExtra(double[], int, double, std::vector<double> &,double[], double &);


//Intra Reference
void computeF0andC0_Intra(double[], int, double, std::vector<double> &,double[], double &);
void computeF0andC0_SingleIntra(double[], int, double, std::vector<double> &,double[], double &);
*/

/*
// Auxiliary
void computeF0andC0_Empty(double[], int, double, std::vector<double> &,double[], double &);
void computeF(const std::vector<double> &, double &, int, double[]);
void computeFandC_S(double[], int, double, std::vector<double> &, const double[], double &);
void optimizeF_S(const std::vector<double> &, double &, double[], double[], int, double[]);

// No regularization
void computeF0andC0single(double[], int, double, std::vector<double> &, double[],double &);
void computeF0andC0(double[], int, double, std::vector<double> &, double[],double &);
void computeFandC(double[], int, double, std::vector<double> &, const double[], double &);
void optimizeF(const std::vector<double> &, double &, double[], double[], int, double[]);

*/



#endif
