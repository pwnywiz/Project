#ifndef POLYS_H_
#define POLYS_H_

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>

#include <QString>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

#define THRESHOLD 0.00001

using namespace Eigen;


// Creates a 2D array with number of rows = main degree
// and columns = hidden degree
class BivPoly
{
	private:
		double **func;
		int main;
		int hidden;
		MatrixXd polyx;
		int maxuni;

	public:
		BivPoly(int x); // Constructor for generated function
		BivPoly(int*&,int*&,double*&,int&,int&,int&); // Constructor for user input function
		BivPoly(MatrixXd);

		~BivPoly();

		void fill_line(double *line, int k); // Fills a 1D array given a specific line
		int get_main() { return main; } // Returns the main degree of the function
		int get_hidden() { return hidden; } // Returns the hidden degree of the function
		bool solveReal(double,double);
		MatrixXd solvey(double);
		int get_maxuni(){ return maxuni; }
        QString printfunc(void);
};

// Creates a 1D array for each function member
// Instances of it are created inside sylvester
class UniPoly
{
	private:
		double *line;

	public:
		UniPoly(int x); // Constructor used for zeros
		UniPoly(BivPoly *cls, int k, int max);
		UniPoly(MatrixXd uniarray);

		~UniPoly();

		void print(int y);
		double get_val(int x) { return line[x]; } // Returns a specific variable
};

#endif /* POLYS_H_ */
