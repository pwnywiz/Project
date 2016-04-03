#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "sy.h"

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <lapacke.h>

using namespace std;
using namespace Eigen;

class Generalized {
private:
    int m,d;
    MatrixXd L0;
    MatrixXd L1;

    MatrixXd multiys;
    int multiy;

public:
    Generalized(PolyMatrixSy&);
    ~Generalized();
    MatrixXd solver();
    int get_multiy() { return multiy; }
	MatrixXd get_multiys() { return multiys; }
};
