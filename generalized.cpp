#include "generalized.h"

/******************************************************************************
******                      Generalized Functions                       ******
******************************************************************************/

Generalized::Generalized(PolyMatrixSy& sy) {
    d = sy.numofmatr()-1;
	m = sy.printdimen();

    if (d == 1) {   //In this case L0=M0,L1=-M1
        L0 = sy.MdToMatrix(0);
        L1 = -1*sy.MdToMatrix(1);

        // cout << "L0 Matrix is\n" << L0 << endl;
        // cout << "\n-L1 Matrix is\n" << L1 << endl;
    }
    else {
        int diag,flag;
    	diag = m;

        L0.resize(m*d,m*d);
        L1.resize(m*d,m*d);

        for (int i=0; i<((m*d)-m); i++){
            flag=1;
            for(int j=0; j<m*d; j++){
                if(diag==j && flag ){
                    L0(i,j) = -1;
                    diag ++;
                    flag=0;
                }else{
                    L0(i,j) = 0;
                }
            }
        }

        diag = 0;
        for (int i=0; i<((m*d)-m); i++){
            flag=1;
            for(int j=0; j<m*d; j++){
                if(diag==j && flag ){
                    L1(i,j) = -1;
                    diag ++;
                    flag=0;
                }else{
                    L1(i,j) = 0;
                }
            }
        }

        for (int i = 0; i < d; i++) {
    		L0.block((m*d)-m,i*m,m,m) = sy.MdToMatrix(i);
            if (i < (d-1))
                L1.block((m*d)-m,i*m,m,m) = MatrixXd::Zero(m, m);
            else
                L1.block((m*d)-m,i*m,m,m) = sy.MdToMatrix(i+1);
            L1 = -1*L1;
    	}

        // cout << "L0 Matrix is\n" << L0 << endl;
        // cout << "\n-L1 Matrix is\n" << L1 << endl;

    }
}

/*****************************************************************************
******************************************************************************/

Generalized::~Generalized() {}

/*****************************************************************************
******************************************************************************/

MatrixXd Generalized::solver(){

    int N = L0.cols();

    MatrixXd v(N,N);
    MatrixXd lambda(N,3);

    int LDA = L0.outerStride();
    int LDB = L1.outerStride();
    int LDV = v.outerStride();
    int INFO = 0;

    double* alphar = lambda.col(0).data();
    double* alphai = lambda.col(1).data();
    double* beta = lambda.col(2).data();

    INFO = LAPACKE_dggev(LAPACK_COL_MAJOR,'N','V',N,L0.data(),LDA,L1.data(),LDB,alphar,alphai,beta,0,LDV,v.data(),LDV);

    // cout << "LAMBDA IS :\n"<<lambda << endl<< endl<< endl;
    // cout << "VECTOR IS :\n"<< v << endl;

    int solutions ;

    solutions = lambda.rows();

    // double x,y;
    int i,j;
    int multicheck[solutions];   //checks y multiplicity

    MatrixXd yvalues(solutions,1);
    MatrixXd yvaluesImg(solutions,1);
    MatrixXd realSols(solutions,2);

    multiy = 0;
    bool multiyflg = false;
    //MatrixXd multiys(solutions,1);

    multiys.resize(solutions,1);


    for( i = 0; i < solutions; i++){
        multicheck[i] = 1;//Init to one
        yvalues(i,0) = lambda(i,0)/lambda(i,2); //Holds real y values
        yvaluesImg(i,0) = lambda(i,1);  //Holds Imaginary y values
    }

    //cout << endl << yvalues << endl;

    for(i = 0; i < solutions; i++){
        for(j = 0; j < solutions; j++){
            if(abs(yvalues(i,0)-yvalues(j,0))<(pow(10,-5)) && i!=j){
                multicheck[i]++;
            }
        }
    }

    int counter = 0;

    for(i = 0; i < solutions; i++){
        if(multicheck[i] > 1 && (fabs(yvaluesImg(i)) < pow(10,-5))){
            //printf("y = %lf %lfi with Multiplicity = %d\n",realSolVal(i),imgSolVal(i),multicheck[i]);
            //printf("y = %lf with Multiplicity = %d\n",realSolVal(i),multicheck[i]);

            for(int v = 0; v < multiy; v++){
                multiyflg = false;
                if((fabs(multiys(v) - yvalues(i))) < pow(10,-5)){
                    multiyflg = true;
                    break;
                }
            }
            if(!multiyflg){
                // printf("y = %lf with Multiplicity = %d\n",yvalues(i),multicheck[i]);
                multiys(multiy) = yvalues(i);
                multiy++;
                multiyflg = false;
            }
        }
        else if(yvaluesImg(i) == 0 && yvalues(i) != INFINITY && yvalues(i) != -INFINITY && !(yvalues(i)!=yvalues(i))){

            //cout << "y = " << realSolVal(i) << " is a real solution." <<endl;

            if(v(solutions-1,i) != 0 && v(solutions-1,i) != INFINITY && !(v(solutions-1,i)!=v(solutions-1,i))){
                realSols(counter,0) = yvalues(i);
                realSols(counter,1) = v(solutions-2,i) / v(solutions-1,i);
                counter++;
            }
        }
        else{
            //MANAGE COMPLEX ROOTS HERE - [BONUS 2]
        }
    }

    multiys.conservativeResize(multiy,1);
    realSols.conservativeResize(counter,2);

    // DEBUG DEBUG DEBUG DEBUG

    // cout << "Real solutions given from Eigen are : ( y - x )"<<endl;
    // cout << endl << realSols << endl;

    // cout << "Y with multiplicity > 1 are :"<<endl;
    // cout << multiys << endl;

    return realSols;
}
