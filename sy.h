#ifndef SY_H_
#define SY_H_

#include "sylvester.h"

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;

class MatrixCoeff{
private:
    double **array;
    int size;

    double *eigenarray;
public:
    MatrixCoeff(double**&,int&);
    ~MatrixCoeff();

    int get_size(){return size;}
    void printme();
    void multiplyme(double**&,double);
    void getmymd(double[]);
};

/*****************************************************************************
******************************************************************************/

class PolyMatrixSy{
private:
    MatrixCoeff** Mcoeffs;
    MatrixXd invertedmd;

    int **triangle;
    int mcoeffsize;
    int power;
    int mccounter;
    double kappa;
    double t1,t2,t3,t4;

public:
    PolyMatrixSy(Sylv&);
    PolyMatrixSy(int x, int y);
    ~PolyMatrixSy();

    void printall();
    void printindi(int);
    int printdimen();
    int numofmatr(){return mccounter;}
    bool computeDet(int);
    double findkappa();
    void invert();
    void fillMd(double **array, int md);
    void pascal(int x);
    void varchange(PolyMatrixSy&);
    //void multiplymd(int**&,int,int);
    MatrixXd multiply(int md);
    MatrixXd MdToMatrix(int x);
    double multi_tafs(double);

    void lapGEP();
};

#endif /* MC_H_ */
