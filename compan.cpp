#include "sy.h"
#include "compan.h"

#include <stdio.h>

using namespace std;

/******************************************************************************
******                      Companion Functions                       	 ******
******************************************************************************/

Companion::Companion(PolyMatrixSy& sy){

	d = sy.numofmatr()-1;
	m = sy.printdimen();

	MatrixXd temp;

	comparray.resize(m*d,m*d);
	sy.invert();

	int diag,flag;
	diag = m;

	if(d > 1){
		for (int i=0; i<((m*d)-m); i++){
			flag=1;
			for(int j=0; j<m*d;	j++){
				if(diag==j && flag ){
					comparray(i,j) = 1;
					diag ++;
					flag=0;
				}else{
					comparray(i,j) = 0;
				}
			}
		}
	}

	for (int i = 0; i < d; i++) {
		temp = sy.multiply(i);
		comparray.block((m*d)-m,i*m,m,m) = temp;
	}
}

/*****************************************************************************
******************************************************************************/

Companion::Companion(UniPoly* polyx,int size){

  comparray.resize(size,size);

  int diag,flag;
	diag = 1;

	if(size  > 1){
		for (int i = 0; i < (size-1); i++){
			flag=1;
			for(int j = 0; j < size; j++){
				if(diag == j && flag ){
					comparray(i,j) = 1;
					diag ++;
					flag=0;
				}else{
					comparray(i,j) = 0;
				}
			}
		}
	}

	for (int i = 0; i < size; i++) {
    	comparray(size-1, i) = -1*polyx->get_val(i);
	}
}

/*****************************************************************************
******************************************************************************/

Companion::~Companion(){}

/*****************************************************************************
******************************************************************************/

MatrixXd Companion::solver(){

	int i,j;
    // double x;

	EigenSolver<MatrixXd> es(comparray);

	int solutions = es.eigenvalues().rows();

	//cout << "Eigenvalues = \n" << es.eigenvalues() << endl;
	//cout << "Eigenvectors = \n" << es.eigenvectors() << endl;

	MatrixXd realSolVal = es.eigenvalues().real();  //holds only real values
	MatrixXd imgSolVal = es.eigenvalues().imag();
	MatrixXd solVec = es.eigenvectors().real();

	int multicheck[solutions];   //checks y multiplicity

    // int realcounter = 0;
	int counter = 0;

	multiy = 0;
	bool multiyflg = false;

	//MatrixXd multiys(solutions,1);
	multiys.resize(solutions,1);
	MatrixXd realSols(solutions,2);

	for( i = 0; i < solutions; i++){
		multicheck[i] = 1;
	}

	for(i = 0; i < solutions; i++){
		for(j = 0; j < solutions; j++){
			//Check the values that appear more than once
			if(fabs(realSolVal(i)-realSolVal(j))<(pow(10,-5)) && i!=j){
				multicheck[i]++;
			}
		}
	}

	for(i = 0; i < solutions; i++){
		if(multicheck[i] > 1 && (fabs(imgSolVal(i)) < pow(10,-5))){
			//printf("y = %lf %lfi with Multiplicity = %d\n",realSolVal(i),imgSolVal(i),multicheck[i]);
			//printf("y = %lf with Multiplicity = %d\n",realSolVal(i),multicheck[i]);

			for(int v = 0; v < multiy; v++){
				multiyflg = false;
				if((fabs(multiys(v) - realSolVal(i))) < pow(10,-5)){
					multiyflg = true;
					break;
				}
			}
			if(!multiyflg){
                //printf("y = %lf with Multiplicity = %d\n",realSolVal(i),multicheck[i]);
				multiys(multiy) = realSolVal(i);
				multiy++;
				multiyflg = false;
			}
		}
		else if(imgSolVal(i) == 0 && realSolVal(i) != INFINITY && realSolVal(i) != -INFINITY && !(realSolVal(i)!=realSolVal(i))){

			//cout << "y = " << realSolVal(i) << " is a real solution." <<endl;

			if(solVec(solutions-1,i) != 0 && solVec(solutions-1,i) != INFINITY && !(solVec(solutions-1,i)!=solVec(solutions-1,i))){
				realSols(counter,0) = realSolVal(i);
				realSols(counter,1) = solVec(solutions-2,i) / solVec(solutions-1,i);
				counter++;
			}
		}
		else{
			//MANAGE COMPLEX ROOTS HERE - [BONUS 2]
		}
	}

	// realSols.resize(counter,2);
	multiys.conservativeResize(multiy,1);
	realSols.conservativeResize(counter,2);

    // DEBUG DEBUG DEBUG DEBUG

     // cout << "Real solutions given from Eigen are : ( y - x )"<<endl;
     // cout << endl << realSols << endl;

     // cout << "Y with multiplicity > 1 are :"<<endl;
     // cout << multiys << endl;

	return realSols;
}

/*****************************************************************************
******************************************************************************/

MatrixXd Companion::unisolver() {

    //cout << "COMPANION ARRAY = " << endl << comparray << endl;

	EigenSolver<MatrixXd> es(comparray);

	int solutions = es.eigenvalues().rows();

    //cout << "Eigenvalues = \n" << es.eigenvalues() << endl;
    //cout << "Eigenvectors = \n" << es.eigenvectors() << endl;

	MatrixXd realSolVal = es.eigenvalues().real();  //holds only real values
	MatrixXd imgSolVal = es.eigenvalues().imag();

	int counter = 0;

	MatrixXd realSols(solutions,1);

	for(int i = 0; i < solutions; i++){
        if((fabs(imgSolVal(i)) < pow(10,-5)) && realSolVal(i) != INFINITY && realSolVal(i) != -INFINITY && !(realSolVal(i)!=realSolVal(i))){

            realSols(counter) = realSolVal(i);
            counter++;
		}
		else{
			//MANAGE COMPLEX ROOTS HERE - [BONUS 2]
		}
	}

	realSols.conservativeResize(counter,1);

     // DEBUG DEBUG DEBUG DEBUG

     // cout << "Real solutions given from Eigen are : ( x )"<<endl;
     // cout << endl << realSols << endl;

     // cout << "Y with multiplicity > 1 are :"<<endl;
     // cout << multiys << endl;

	return realSols;
}
