#include "sylvester.h"
#include "polys.h"
#include <iostream>

using namespace std;

/******************************************************************************
******                      Sylvester Functions                          ******
******************************************************************************/

Sylv::Sylv(BivPoly *clsA, BivPoly *clsB){

  // Store the individual function degrees
  main1 = clsA->get_main();
  hidden1 = clsA->get_hidden();
  main2 = clsB->get_main();
  hidden2 = clsB->get_hidden();

  // Store the max main and hidden degree
  if (main1 > main2) {
    maxMain = main1;
  } else {
    maxMain = main2;
  }
  if (hidden1 > hidden2) {
    maxHidden = hidden1;
  } else {
    maxHidden = hidden2;
  }

  // 1D Array of pointers to UniPoly object for each function
  // in order to seperate the hidden variable
  func1 = new UniPoly*[main1+1];
  for (int i = 0; i <= main1; i++) {
    func1[main1-i] = new UniPoly(clsA, i, maxHidden);
  }

  func2 = new UniPoly*[main2+1];
  for (int i = 0; i <= main2; i++) {
    func2[main2-i] = new UniPoly(clsB, i, maxHidden);
  }

  // 1D array of zeros for the matrix filling
  zeros = new UniPoly(maxHidden);

  // 2D array of pointers to UniPoly object since the values
  // of each line are repeated instead of doing multiple saves (memory sufficient)
  sylvester = new UniPoly**[main1+main2];
  for (int i = 0; i < main1+main2; ++i) {
    sylvester[i] = new UniPoly*[main1+main2];
  }

  // Each element array points to a UniPoly object
  for (int i = 0; i < main1+main2; i++) {
    if (i < main2) {
      for (int j = 0; j < main1+main2; j++) {
        if (j < i) {
          sylvester[i][j] = zeros;
        } else if (j-i < (main1+1)) {
          sylvester[i][j] = func1[j-i];
        } else {
          sylvester[i][j] = zeros;
        }
      }
    } else {
      for (int j = 0; j < main1+main2; j++) {
        if (j < (i - main2)) {
          sylvester[i][j] = zeros;
        } else if (j-(i-main2) < (main2+1)) {
          sylvester[i][j] = func2[j-(i-main2)];
        } else {
          sylvester[i][j] = zeros;
        }
      }
    }
  }
}

/*****************************************************************************
******************************************************************************/

void Sylv::printSylv()
{
  for (int i = 0; i < main1+main2; i++) {
    for (int j = 0; j < main1+main2; j++) {
      sylvester[i][j]->print(maxHidden);
    }
    cout << endl;
  }
  cout<<endl;
}

/*****************************************************************************
******************************************************************************/

double Sylv::specific(int row, int col, int part)
{
  return sylvester[row][col]->get_val(part);
}

/*****************************************************************************
******************************************************************************/

Sylv::~Sylv()
{
  for (int i = 0; i <= main1; i++) {
    delete func1[i];
  }
  delete [] func1;

  for (int i = 0; i <= main2; i++) {
    delete func2[i];
  }
  delete [] func2;

  delete zeros;

  for (int i = 0 ; i < main1+main2 ; i++) {
    delete[] sylvester[i];
  }
  delete[] sylvester;
}
