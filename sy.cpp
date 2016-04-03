#include "sy.h"
#include "initfuncs.h"

#include <lapacke.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace Eigen;

/******************************************************************************
******                      PolyMatrixSy Functions                       ******
******************************************************************************/

//Constructor which takes a sylvester matrix as refference
PolyMatrixSy::PolyMatrixSy(Sylv& Slv){

    //We have to create so many matrix-coeffs as is the max exponent of our
    //hidden variable (y) + 1 for the constant term
    mccounter=Slv.get_hidden() + 1;

    //We need <mccounter> class instances of Mcoeffs class type
    Mcoeffs = new MatrixCoeff*[mccounter];

    //Each matrix-coeff internal array will have the same size
    //with Sylvester matrix
    mcoeffsize = Slv.get_total();//returns size of sylvester array

    double** buffer; //buffer array for obtaining elements from Sylvester matrix

    buffer = new double*[mcoeffsize];

    for(int i = 0;i < mcoeffsize;i++){
        buffer[i] = new double[mcoeffsize];
    }

    for(int i = 0; i<mccounter; i++){

        for(int k = 0 ; k<mcoeffsize; k++){
            for(int n=0; n<mcoeffsize; n++){
                buffer[k][n] = Slv.specific(k,n,i); //row collumn part
                //We obatin from Sylvester matrix the element which is in
                //row k,collumn n and in position i in the array
                //of the UniPoly class
            }
        }

        //new MatrixCoeff object which is constructed
        //with the buffer array as initialization
        Mcoeffs[i] = new MatrixCoeff(buffer,mcoeffsize);
    }

    //We dont need anymore the buffer array
    for (int i = 0; i < mcoeffsize; i++) {
        delete [] buffer[i];
    }
    delete [] buffer;
}

/*****************************************************************************
******************************************************************************/

//New Constructor which is used in variable change
PolyMatrixSy::PolyMatrixSy(int mccount, int matrixsize) {
    Mcoeffs = new MatrixCoeff*[mccount];
    mccounter = mccount;
    mcoeffsize = matrixsize;
}

/*****************************************************************************
******************************************************************************/

PolyMatrixSy::~PolyMatrixSy(){
    for (int i = 0; i < mccounter; i++) {
        delete Mcoeffs[i];
    }
    delete [] Mcoeffs;

    //Each PolyMatrixSy Holds its Euclidean triangle
    for (int i = 0; i <= power; i++) {
      delete [] triangle[i];
    }
    delete [] triangle;
}

/*****************************************************************************
******************************************************************************/

// Used to fill the Matrix Coeffs that are required for the change of variable
void PolyMatrixSy::fillMd(double **array, int md) {
    Mcoeffs[md] = new MatrixCoeff(array,mcoeffsize);
}

/*****************************************************************************
******************************************************************************/

//It prints out every Matrix-Coeff object of the Sy
void PolyMatrixSy::printall(){
    for(int i = 0; i<mccounter; i++){
        if(i!=0){
            cout<<"Printing Coefficient matrix of y^"<<i<<endl;
        }
        else{
            cout<<"Printing constant Coefficient matrix"<<endl;
        }
        Mcoeffs[i]->printme();
    }
}

/*****************************************************************************
******************************************************************************/

//It prints out the individual Matrix-Coeff object
//that we give as parameter (which)
void PolyMatrixSy::printindi(int which){
    /*
    if(which!=0){
        cout<<"Printing Coefficient matrix of y^"<<which<<endl;
    }
    else{
        cout<<"Printing constant Coefficient matrix"<<endl;
    }
    */
    Mcoeffs[which]->printme();
}

/*****************************************************************************
******************************************************************************/

int PolyMatrixSy::printdimen() {
    return Mcoeffs[0]->get_size();
}

/*****************************************************************************
******************************************************************************/

//Function to compute Kappa for Md matrix of M(y)
double PolyMatrixSy::findkappa(void){

    int size = Mcoeffs[mccounter-1]->get_size();

    double eigenarrayTmp[size*size];
    Mcoeffs[mccounter-1]->getmymd(eigenarrayTmp);

    // for(int i = 0; i < size*size; i++){
    //     cout<<eigenarrayTmp[i]<<" , ";
    // }

    Map<MatrixXd> m(eigenarrayTmp,size,size);
    JacobiSVD<MatrixXd> svd(m,Eigen::ComputeFullU);
    JacobiSVD<MatrixXd>::SingularValuesType singular = svd.singularValues();

    //cout << "The singular values are : "<<svd.singularValues()<<endl;

    // for(int i = 0; i<singular.rows(); i++){
    //     cout << "singular Values" << i <<":"<<singular(i)<<endl;
    // }

    return singular(0)/singular(singular.rows()-1);

}

/*****************************************************************************
******************************************************************************/
//Function to invert the Md Matrix-Coefficient
void PolyMatrixSy::invert(){

    int size = Mcoeffs[mccounter-1]->get_size();
    double eigenarrayTmp[size*size];
    Mcoeffs[mccounter-1]->getmymd(eigenarrayTmp);

    Map<MatrixXd> m(eigenarrayTmp,size,size);

    invertedmd = m.inverse();

    // Map<MatrixXd>( eigenarrayTmp, invertedmd.rows(), invertedmd.cols() ) = invertedmd;
}

/*****************************************************************************
******************************************************************************/

//Fucntion which implements the multiplication of Md  with another Matrix
MatrixXd PolyMatrixSy::multiply(int md){

    MatrixXd multipliedmd;
    int size = Mcoeffs[mccounter-1]->get_size();
    double eigenarrayTmp[size*size];

    Mcoeffs[md]->getmymd(eigenarrayTmp);

    Map<MatrixXd> m(eigenarrayTmp,size,size);

    multipliedmd = -1*invertedmd*m;

    return multipliedmd;
}

/*****************************************************************************
******************************************************************************/
//Function which convers the array of doubles to MatrixXd Array
MatrixXd PolyMatrixSy::MdToMatrix(int x) {
    int size = Mcoeffs[mccounter-1]->get_size();
    double temp[size*size];
    Mcoeffs[x]->getmymd(temp);
    Map<MatrixXd> m(temp,size,size);

    return m;
}

/*****************************************************************************
******************************************************************************/

//BONUS - 1 - Computation of det of the given Yo Matrix-Coefficient
bool PolyMatrixSy::computeDet(int yo){

    int size = Mcoeffs[mccounter-1]->get_size();
    double mtemp[size*size];

    Mcoeffs[yo]->getmymd(mtemp);

    Map<MatrixXd> compdet(mtemp,size,size);

    double det = compdet.determinant();

    cout <<"\nDeterminant of : M"<<yo<<" = "<<det<<endl;

    if(det == 0) return true;
    else return false;

}

/*****************************************************************************
******************************************************************************/

//For every PolyMatrixSy we have a Pascal triangle implementation
//Ideal for use in variable change
void PolyMatrixSy::pascal(int x) {
  int j;

  power = x;
  triangle = new int*[power+1];
  for(int i = 0; i <= power; ++i)
      triangle[i] = new int[power+1];

  for (int i = 0; i <= power; i++) {
    for (int k = 0; k <= power; k++) {
      triangle[i][k] = 0;
    }
  }

  for (int i = 0; i <= power; i++) {
    j = 0;
    while (j <= i) {
      if (j == 0 || j == i) {
        triangle[i][j] = 1;
      }
      else {
        triangle[i][j] = triangle[i-1][j-1] + triangle[i-1][j];
      }

      j++;
    }
  }
}

/*****************************************************************************
******************************************************************************/

//Task4 - Variable change with random numbers
void PolyMatrixSy::varchange(PolyMatrixSy &New) {
  t1 =((double)rand()/(double)RAND_MAX) * ((rand() % 20) -10);
  t2 =((double)rand()/(double)RAND_MAX) * ((rand() % 20) -10);
  t3 =((double)rand()/(double)RAND_MAX) * ((rand() % 20) -10);
  t4 =((double)rand()/(double)RAND_MAX) * ((rand() % 20) -10);
  int d = mccounter - 1;
  int k;
  double num[2];
  double den[2];

  num[0] = t2;
  num[1] = t1;
  den[0] = t4;
  den[1] = t3;

  double *devel1 = new double[d+1];//Holds the numerator coeff analysis
  double *devel2 = new double[d+1];//Holds the denominator coeff analysis

  double **final;

  final = new double*[d+1];

  for (int i = 0; i <= d; i++) {
      final[i] = new double[d+1];
  }

  for (int i = 0; i <= d; i++) {
      for (int j = 0; j <= d; j++) {
          final[i][j] = 0;
      }
  }

  //Analyse every coefficient of M MatrixCoeff
  for (int m = 0; m <= d; m++) {
      for (k = 0; k <= d-m; k++) {
          devel1[k] = pow(num[0],(d-m)-k)*pow(num[1],k)*triangle[d-m][k];
          // devel2[k] = pow(den[0],k)*pow(den[1],d-k)*triangle[m][k];
      }

      //Fill whatever is left with zeros
      while (k <= d) {
          devel1[k] = 0;
          k++;
      }

      for (k = 0; k <= m; k++) {
          devel2[k] = pow(den[0],m-k)*pow(den[1],k)*triangle[m][k];
      }
      while (k <= d) {
          devel2[k] = 0;
          k++;
      }

      for (int i = 0; i <= d; i++) {
          for (int j = 0; j <= d; j++) {
              if ((i+j) <= d) {
                  final[m][i+j] += devel1[i]*devel2[j];
              }
          }
      }
  }

  //We have created the final matrix which holds every coefficient for
  //every matrix we have.Now we multiply every coeff of final with the
  //correct Matrix and then we add them all together to produce the new M'
  double **temp;
  temp = new double*[mcoeffsize];

  for (int i = 0; i < mcoeffsize; i++) {
      temp[i] = new double[mcoeffsize];
  }
  for (int i = 0; i < mcoeffsize; i++) {
      for (int j = 0; j < mcoeffsize; j++) {
          temp[i][j] = 0;
      }
  }

  double **mega;
  mega = new double*[mcoeffsize];

  for (int i = 0; i < mcoeffsize; i++) {
      mega[i] = new double[mcoeffsize];
  }
  for (int i = 0; i < mcoeffsize; i++) {
      for (int j = 0; j < mcoeffsize; j++) {
          mega[i][j] = 0;
      }
  }

  for (int m = 0; m <= d; m++) {
      for(int j = 0; j <= d; j++){
          Mcoeffs[d-j]->multiplyme(temp,final[j][m]);

          for(int l = 0; l < mcoeffsize; l++){
              for(int i = 0; i < mcoeffsize; i++){
                  mega[l][i] += temp[l][i];
              }
          }
      }

      New.fillMd(mega,m);

      for (int i = 0; i < mcoeffsize; i++) {
          for (int j = 0; j < mcoeffsize; j++) {
              mega[i][j] = 0;
          }
      }
  }

  delete [] devel1;
  delete [] devel2;

  for (int i = 0; i <= d; i++) {
    delete [] final[i];
  }
  delete [] final;

  for (int i = 0; i < mcoeffsize; i++) {
    delete [] temp[i];
  }
  delete [] temp;

  for (int i = 0; i < mcoeffsize; i++) {
    delete [] mega[i];
  }
  delete [] mega;
}

/*****************************************************************************
******************************************************************************/

double PolyMatrixSy::multi_tafs(double psi){

    double temp;
    temp = (t1*psi + t2)/(t3*psi +t4);
    return temp;
}

/******************************************************************************
******                      MatrixCoeff Functions                       *******
******************************************************************************/

//Constructor which gets as parameter a refference to
//a 2d array (buffer) and its size in order to initialize
//its own privare array
MatrixCoeff::MatrixCoeff(double**& buffer,int& matrixsize){

    size=matrixsize;
    int i,j;

    array = new double*[size];
    eigenarray = new double[size*size];

    for(i = 0;i < size;i++){
        array[i] = new double[(size)];
    }

    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            array[i][j] = buffer [i][j];
            eigenarray[i*size + j] = buffer[j][i];
        }
    }
}

/*****************************************************************************
******************************************************************************/

MatrixCoeff::~MatrixCoeff(){
    for (int i = 0; i < size; i++) {
        delete [] array[i];
    }
    delete [] array;
    delete [] eigenarray;
}

/*****************************************************************************
******************************************************************************/

//It prints out itself
void MatrixCoeff::printme(){
    for(int i = 0; i<size; i++){
        for(int j = 0; j<size; j++){
            cout<<array[i][j]<<" ";
        }
        cout<<endl;
    }
}

/*****************************************************************************
******************************************************************************/

void MatrixCoeff::getmymd(double temparray[]){

    for(int i = 0; i<size*size; i++){
        temparray[i] = eigenarray[i];
    }
}

/*****************************************************************************
******************************************************************************/

void MatrixCoeff::multiplyme(double**& temp,double multi){

    for(int i = 0; i<size; i++){
        for(int j = 0; j<size; j++){
            temp[i][j] = array[i][j] * multi;
        }
    }

}
