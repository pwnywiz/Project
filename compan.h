#include <iostream>
#include <stdlib.h>
#include "sy.h"

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;

class Companion{
private:
	MatrixXd comparray;
	MatrixXd multiys;
	int m,d;
	int multiy;

public:
	Companion(PolyMatrixSy&);
	Companion(UniPoly*,int);
	~Companion();
	MatrixXd solver (void);
	MatrixXd unisolver(void);
	int get_multiy() { return multiy; }
	MatrixXd get_multiys() { return multiys; }


};
