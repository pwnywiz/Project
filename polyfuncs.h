#ifndef POLYFUNCS_H_
#define POLYFUNCS_H_

#include <string>
#include <locale.h>

using namespace std;

string readfterm(int&);
void readffile(char*& ,string& ,string& );

int countmononyms(string&);
int findmax(int*&,int);

bool polanalysis(string& ,int& ,int&,int*&,int*&,double*&);

#endif /* POLYFUNCS_H_ */
