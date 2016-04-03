#include "interpolation.h"

using namespace std;

/******************************************************************************
******                     Interpolation Functions                       ******
******************************************************************************/

Interpolation::Interpolation(MatrixXd points,int d){
//Interpolation::Interpolation(){

    int a;
    int b;

    rank = 0;
    deg = d;
    k = ( ((deg+2)*(deg+1)) / 2 ) - 1;
    inter.resize(k,k+1);

    order.resize(k+1,2);


  	bool flag;

    // cout << "The (" << k << ") Selected points are : " << endl << points << endl;

    //Interpolation matrix initialization
    for(int i = 0; i < k; i++){
        a = 0;
        b = 0;
        flag = false;

        for(int j = 0; j <= k; j++){
            if( b <= deg && (a+b) <= deg){
                if(flag == true){b++; flag=false;}
                if(i == 0){ order(j,0) = a; order(j,1) = b; }
                inter(i,j) = pow(points(i,0),a)*pow(points(i,1),b);
                b++;
            }
            else{
                a++;
                b=0;
                flag = true;
                if(i == 0){ order(j,0) = a; order(j,1) = b; }
                inter(i,j) = pow(points(i,0),a)*pow(points(i,1),b);
            }
        }
    }

    // cout << "Interpolation = \n" << inter << endl;

}

/*****************************************************************************
******************************************************************************/

Interpolation::~Interpolation(){}

/*****************************************************************************
******************************************************************************/

bool Interpolation::solver() {

    FullPivLU<MatrixXd> luA(inter);
    rank = luA.rank();

    // cout << endl << "RANK = " << rank << endl ;

    if (rank != k) {
        // cout << "The problem is infeasible or the solution is numerically unstable" << endl;
        return false;
    }

    kernelv = luA.kernel();
    // cout << "KERNEL = " << endl << kernelv << endl;

    return true;
}

/*****************************************************************************
******************************************************************************/

MatrixXd Interpolation::generatePoly(int main) {

  int mlim,hlim;

  if (main == 1) {
  	mlim = maxa;
    hlim = maxb;
  }
  else {
  	mlim = maxb;
    hlim = maxa;
  }

  func.resize(mlim+1,hlim+1);
  func.setZero(mlim+1,hlim+1);

  for( int i = 0; i <= k; i++){
    if (order(i,0) <= maxa && order(i,1) <= maxb) {
      if (main == 1) {
      	func(order(i,0),order(i,1)) = kernelv(i);
      }
      else func(order(i,1),order(i,0)) = kernelv(i);
    }
  }

  return func;
}
/*****************************************************************************
******************************************************************************/

void Interpolation::calcmax(){

    maxa = 0;
    maxb = 0;

    for(int i = 0; i<=k; i++){
        if(order(i,0) > maxa && kernelv(i) != 0){
            maxa = order(i,0);
        }

        if(order(i,1) > maxb && kernelv(i) != 0){
            maxb = order(i,1);
        }
    }

    if(maxb > maxa) genmax = 2;
    else            genmax = 1;
}
