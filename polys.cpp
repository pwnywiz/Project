#include "polys.h"
#include "sylvester.h"
#include <iostream>

using namespace std;

/******************************************************************************
******                      Polynomials Functions                        ******
******************************************************************************/

BivPoly::BivPoly(int x)
{
  int a,b;

  // Generate the degree of x and y
  a = rand() % (x+1);
  b = x - a + (rand() % (a+1)); // The total degree must be at least x (plus a random acceptable number)

  // Decide which will be the main degree
  if (a > b || a == b) {
    main = a;
    hidden = b;
  } else {
    main = b;
    hidden = a;
  }

  // 2D Array initialization
  func = new double*[main+1];
  for (int i = 0; i <= main; i++) {
    func[i] = new double[hidden+1];
  }

  // Randomly fill in the array
  for (int i = 0; i <= main; i++) {
    for (int j = 0; j <= hidden; j++) {
      if ((i + j) <= x) {
        func[i][j] = ((double)rand()/(double)RAND_MAX) * ((rand() % 100) -50 ); // Random values from -50 to 50 as requested
      } else {
        func[i][j] = 0;
      }
      //cout << func[i][j] << " ";
    }
    cout << endl;
  }
  cout<<endl;
}

/*****************************************************************************
******************************************************************************/

// Pretty much the same with above, accepting user input instead
BivPoly::BivPoly(int*& expmain,int*& exphid,double*& coeff,int& maxm,int& maxh,int& monoms){
    int i,j;

    main = maxm;    //main max degree
    hidden = maxh; //hidden max degree

    // 2D Array initialization
    func = new double*[(maxm+1)];
    for(i = 0;i <= maxm;i++){
        func[i] = new double[(maxh+1)];
    }

    for(i = 0; i <=maxm; i++){
        for(j = 0; j <=maxh; j++){
            func[i][j]=0;
        }
    }

    for(i = 0; i<monoms; i++){
        func[expmain[i]][exphid[i]] = coeff[i];
    }

    // cout<<endl;
}

/*****************************************************************************
******************************************************************************/

BivPoly::BivPoly(MatrixXd function){

    main = function.rows() - 1;
    hidden = function.cols() - 1;

    func = new double*[function.rows()];

    for(int i = 0; i < function.rows(); i++){
        func[i] = new double[function.cols()];
    }

    for(int i = 0; i < function.rows(); i++){
        for(int j = 0; j < function.cols(); j++){
            /* THRESHOLD ENABLE FOR INTERPOLATION USAGE */

            // if(abs(function(i,j)) < THRESHOLD) func[i][j] = 0;
            //  else func[i][j] = function(i,j);

            func[i][j] = function(i,j);
        }
    }
}

/*****************************************************************************
******************************************************************************/

void BivPoly::fill_line(double *line, int k)
{
  for (int i = 0; i <= hidden; i++) {
    line[i] = func[k][i];
  }
}

/*****************************************************************************
******************************************************************************/

bool BivPoly::solveReal(double x,double y){

    double result = 0;

    for(int i = 0; i<=main; i++){
        for(int j = 0; j<=hidden; j++){
            if(func[i][j] != 0){
                result+= func[i][j] * pow(x,(double)i) * pow(y,(double)j);
            }
        }
    }

    // cout << "apotelesma epilisis me x = " << x << "y  = " << y << "   === "<< result <<endl << endl;

    if(abs(result) < pow(10,-5) ) return true;
    else return false;

}

/*****************************************************************************
******************************************************************************/

MatrixXd BivPoly::solvey(double y) {

  double temp1 = 0;
  double maxcoef;

  polyx.resize(main+1, 1);
  for (int i = 0; i <= main; i++) {
    for (int j = 0; j <= hidden; j++) {
      temp1 += func[i][j]*pow(y,j);
    }
    polyx(i) = temp1;
    temp1 = 0;
  }


  for(int i = main; i >= 0; i--){
      if( polyx(i) != 0){
          maxcoef = polyx(i);
          maxuni = i;
          break;
      }
  }
  for (int i = 0; i <= main; i++) {
    // cout << "THA DIAIRESW TO " << polyx(i) << " ME TO " << maxcoef << endl;
    polyx(i) /= maxcoef;
    // cout << "DIVISION = " << polyx(i) << endl;
  }

  return polyx;
}

/*****************************************************************************
******************************************************************************/

QString BivPoly::printfunc(void){
    QString str;
    /*
    for(int i = 0; i <= main; i++){
        for(int j = 0; j <= hidden; j++){
            if(func[i][j] != 0){
                if(func[i][j] > 0 ) cout <<"+ ";
                if(i == 0 && j == 0) cout << func[i][j] << " ";
                else if(i == 0) cout << func[i][j] <<"*y^" << j <<" ";
                else if(j == 0) cout << func[i][j] <<"*x^" << i <<" ";
                else cout << func[i][j] <<"*x^" << i <<"*y^" << j <<" ";
            }
        }
    }
    */
    for(int i = 0; i <= main; i++){
        for(int j = 0; j <= hidden; j++){
            if(func[i][j] != 0){
                if(func[i][j] > 0 ) str.append("+ ");
                if(i == 0 && j == 0) {
                    str.append(QString::number(func[i][i]));
                    str.append(" ");
                }
                else if(i == 0) {
                    str.append(QString::number(func[i][j]));
                    str.append("*y^");
                    str.append(QString::number(j));
                    str.append(" ");
                }
                else if(j == 0) {
                    str.append(QString::number(func[i][j]));
                    str.append("*x^");
                    str.append(QString::number(i));
                    str.append(" ");
                }
                else {
                    str.append(QString::number(func[i][j]));
                    str.append("*x^");
                    str.append(QString::number(i));
                    str.append("*y^");
                    str.append(QString::number(j));
                    str.append(" ");
                }
            }
        }
    }

    // cout << endl << endl;
    return str;
}
/*****************************************************************************
******************************************************************************/

BivPoly::~BivPoly()
{
  for (int i = 0; i <= main; i++) {
    delete [] func[i];
  }
  delete [] func;
}

/*****************************************************************************
*******        Unipoly/Univariate polynomial functions                ********
******************************************************************************/

UniPoly::UniPoly(int x)
{
    // 1D Array allocation
    line = new double[x+1];
    for (int i = 0; i <= x; i++) {
      line[i] = 0; // Fill in the array with zeros
    }
}

/*****************************************************************************
******************************************************************************/

UniPoly::UniPoly(BivPoly *cls, int k, int max)
{
  line = new double[max+1];
  for (int i = 0; i <= max; i++) {
    line[i] = 0;
  }
  cls->fill_line(line, k); // Send the array with a line request to the 2D struct for filling
}

/*****************************************************************************
******************************************************************************/

UniPoly::UniPoly(MatrixXd uniarray){
    int size = uniarray.rows();

    line = new double[size];
    for(int i = 0; i < size; i++){
        line[i] = uniarray(i);
    }
}

/*****************************************************************************
******************************************************************************/

void UniPoly::print(int y)
{
  cout << "[";
  for (int i = 0; i <= y; i++) {
    cout << line[i];
    if (i < y) {
      cout << ",";
    }
  }
  cout << "]";
}

/*****************************************************************************
******************************************************************************/

UniPoly::~UniPoly() {
  delete [] line;
}
