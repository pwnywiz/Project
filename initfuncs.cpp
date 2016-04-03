#include "initfuncs.h"

/******************************************************************************
******      Initialize functions for use with command line execusion     ******
******************************************************************************/

//Function which prints out the usage of our programm
//it's called when -man is given as cl argument
void use_info(char* name){
    cout <<endl<< "##### Command line parameters usage methods are #####" <<endl<<endl;
    cout <<"1) "<<name<<" -generate -d1 <D1> -d2 <D2> : Random generation of two polynomials with <D1> and <D2> degree each."<<endl<<endl;
    cout <<"2) "<<name<<" -read -console : User inserts two polynomials directly on console."<<endl<<endl;
    cout <<"3) "<<name<<" -read -i <file_path> -d1 <D1> -d2 <D2> : User inserts two polynomials to be read from file and also gives their degrees."<<endl<<endl;
    cout <<"By adding -solve <B> in command line arguments,the system's solutions are calculated."<<endl<<endl;
}

/*****************************************************************************
******************************************************************************/

//analyses the command line arguments for each case
void initfunc(int& argc,char** argv,int& method,int& d1,int& d2,int& B,char*& file){

    //in this case we need to print out the manual
    if((string)argv[1]=="-man"){
        method=-1;
        use_info(argv[0]);
    }

    //in this case we have random generation of polynomials
    else if((string)argv[1]=="-generate") {
        method=1;
        if((string)argv[2]=="-d1"){
            d1=atoi(argv[3]);
            d2=atoi(argv[5]);
        }
        else{
            d1=atoi(argv[5]);
            d2=atoi(argv[3]);
        }
        if(argc >= 7){
            if((string)argv[6]=="-solve"){
                if(argc == 8)   B=atoi(argv[7]);
                else            B=7;

                cout<<"\nGeneration method with d1 = "<<d1<<" and d2 = "<<d2<<" and solving with B = "<<B<<endl;
            }
        }
        else{
            cout<<"\nGeneration method with d1 = "<<d1<<" and d2 = "<<d2<<endl;
        }

    }
    else if((string)argv[1]=="-read"){

        if((string)argv[2]=="-console"){

            method=2;

            if(argc >= 4){
                if((string)argv[3]=="-solve"){
                    if(argc == 5)   B=atoi(argv[4]);
                    else            B=7;

                    cout<<"\nReading from command line method and solving with B = "<<B<<endl;
                }
            }
            else{
                cout<<"\nReading from command line method . . ."<<endl<<endl;
            }
        }
        else if((string)argv[2]=="-i"){
            method=3;
            file=argv[3];
            if((string)argv[4]=="-d1"){
                d1=atoi(argv[5]);
                d2=atoi(argv[7]);
            }
            else{
                d1=atoi(argv[7]);
                d2=atoi(argv[5]);
            }

            if(argc >= 9){
                if((string)argv[8]=="-solve"){
                    if(argc == 10)   B=atoi(argv[9]);
                    else            B=7;

                    cout<<"\nReading from file: "<<file<<" with d1 = "<<d1<<" and d2 = "<<d2<<" and solving with B = "<<B<<endl<<endl;
                }
            }
            else{
                cout<<"\nReading from file: "<<file<<" with d1 = "<<d1<<" and d2 = "<<d2<<endl<<endl;
            }
        }
    }
}
