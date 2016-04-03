#include <iostream>
#include <stdlib.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

#define THRESHOLD 0.00001

using namespace Eigen;


class Interpolation{
private:
    MatrixXd inter;
    MatrixXd kernelv;
    MatrixXd order;
    MatrixXd func;
    int deg;
    int k;
    int rank;
    int maxa,maxb,genmax;

public:
    Interpolation(MatrixXd,int);
    ~Interpolation();
    bool solver();
    MatrixXd generatePoly(int);

    void calcmax();

    MatrixXd get_inter(){ return inter; }
    MatrixXd get_kernel() { return kernelv; }
    int get_deg(){ return deg; }
    int get_k(){ return k; }
    int get_rank(){ return rank; }

    int get_genmax(){ return genmax; }
    int get_maxa(){ return maxa; }
    int get_maxb(){ return maxb; }


};
