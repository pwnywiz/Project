#include "initfuncs.h"
#include "polyfuncs.h"

/******************************************************************************
******           Various functions for polynomials                       ******
******************************************************************************/

//Reads from terminal the bivariate polynomial
//and its defgree
string readfterm(int& deg){

    string function;

    cout<<"Give the polynomial function : \n";
    getline(cin,function);
    cout<<"and now give its' degree : \n";
    cin>>deg;
    cin.ignore();
    return function;
}

/*****************************************************************************
******************************************************************************/

//Reads a pair of polynomials from file,line by line
void readffile(char*& infile,string& f1,string& f2){

    ifstream inputfile;
    inputfile.open(infile);

    getline(inputfile,f1);
    getline(inputfile,f2);

    //inputfile.close();
}

/*****************************************************************************
******************************************************************************/

//It gets a bivariate polynomial in string
//and counts its mononyms,action which is used
//for better analysis of the polynomial
int countmononyms(string& f){
    int i,monoms=1;
    int fsize=f.size();

    for(i=1;i<fsize;i++){
        if(f[i]=='+' || f[i]=='-') monoms++;
    }
    return monoms;
}

/*****************************************************************************
******************************************************************************/

//It gets an array and returns its
//max element
int findmax(int*& array,int size){
    int max=array[0];

    for(int i=1;i<size;i++){
        if(array[i]>max) max=array[i];
    }
    return max;
}

/*****************************************************************************
******************************************************************************/

//Analyses bivariate polynomial from string into an array of integers
bool polanalysis(string& f1,int& d1,int& monoms1,int*& exponx,int*& expony,double*& coeff){

    setlocale (LC_NUMERIC, "en_US.utf8");
    int k,i,j;

    int checkd=-1;

    int symbolpos,digitcount,digitpos,expo,monomsize;
    string symbol,previous;

    int f1size=f1.size();

    j=0;k=0;
    string monomStr[monoms1]; //holds each mononym

    for(i=1;i<f1size;++i){
        if(f1[i] == '+' || f1[i] == '-')
		{
			monomStr[k++] = f1.substr(j, i - j);
			j = i;
		}
    }

    monomStr[k]=f1.substr(j,f1size-j+1); // adds the last mononym to array
    monomStr[k].push_back(' '); // adds a space char to the end of the last

    // cout<<f1<<endl;

    //prints out the mononyms of the polynoial function

    // for(i=0;i<monoms1;i++) cout<<monomStr[i]<<"|";
    // cout<<endl;

    //sets default values to the arrays
    for(i=0;i<monoms1;i++){
        coeff[i]=1; exponx[i]=0; expony[i]=0;
    }

    //for each mononym of the polynomial function
    //runs the following alanysis
    for(i=0;i<monoms1;i++){

        monomsize=monomStr[i].size();

        symbolpos=-1;
        symbol=" ";
        previous=" ";
        digitcount=0;
        expo=0;
        digitpos=-1;

        for(j=0;j<monomsize;j++){
            if(symbolpos==-1){ //define sign
                if((monomStr[i][j]!='+' && monomStr[i][j]!='-') || monomStr[i][j]=='+') symbol="+";
                else symbol="-";
                symbolpos=1; //flag for symbol when found
            }

            //if we have digit,continue to read
            if(isdigit(monomStr[i][j]) || monomStr[i][j]=='.'){
                digitcount++;
                if(digitcount==1) digitpos=j;
                if(isdigit(monomStr[i][j+1]) || monomStr[i][j]=='.') continue;
            }

            //if the exponent is flagged..
            if(expo==1){
                if(previous=="x"){
                    exponx[i]=atoi((monomStr[i].substr(digitpos,digitcount)).c_str());
                }
                else if(previous=="y"){
                    expony[i]=atoi((monomStr[i].substr(digitpos,digitcount)).c_str());
                }
                expo=0;
                digitcount=0;
                digitpos=-1;
            }

            if(monomStr[i][j]=='*'){
                if(digitcount!=0){  //coeff found
                    if(symbol=="-"){
                        coeff[i]=(-1)*atof((monomStr[i].substr(digitpos,j)).c_str());
                    }
                    else{
                        coeff[i]=atof((monomStr[i].substr(digitpos,j)).c_str());
                    }
                    digitcount=0;
                    digitpos=-1;
                }
                else if(monomStr[i][j-1]=='x') exponx[i]=1; //x without ^
                else if(monomStr[i][j-1]=='y') expony[i]=1; //y without ^
            }

            else if(monomStr[i][j]=='x' || monomStr[i][j]=='y'){
                previous=monomStr[i][j];
                if(coeff[i]==1 && symbol=="-"){
                    coeff[i]=-1; //case with no given coeff but with - in front
                }
            }
            else if(monomStr[i][j]=='^') expo=1; //flag the exponent


            if(monomStr[i][j]==' '){
                if(digitcount!=0){
                    if(symbol=="-"){
                        coeff[i]=(-1)*atof((monomStr[i].substr(digitpos,j)).c_str());
                    }
                    else{
                        coeff[i]=atof((monomStr[i].substr(digitpos,j)).c_str());
                    }
                    exponx[i]=0;
                    expony[i]=0;
                    digitpos=-1;
                    digitcount=0;
                }
                else if(monomStr[i][j-1]=='x') exponx[i]=1; //x without ^
                else if(monomStr[i][j-1]=='y') expony[i]=1; //y without ^
            }

        }

        if(checkd < (exponx[i] + expony[i])){
            checkd = exponx[i] + expony[i];
        }

    }

    //it prints out the analysis of the polynomial
    //  for(i=0;i<monoms1;i++){
    //        cout<<i+1<<") coef = "<<coeff[i]<<" || exponx = "<<exponx[i]<<" || expony = "<<expony[i]<<endl;
    //    }
    
    if(checkd!=d1){
        d1=checkd;
        return false;
    }
    else return true;
}
