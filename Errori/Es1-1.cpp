//Gruppo: 
//Dellepiane Emanuele, 4876072
//Manini Filippo, 4798004
//Miggiano Davide, 4840761 
#include <iostream>
#include <cmath>
using namespace std;


double findA(int d0, int i)
{
    double a = (double)(d0 + 1) * pow(10,i);
    return a;
}

double findB(int d1)
{
    double b = (double)(d1 + 1) * pow(10,20);
    return b;
}

double findC(double b)
{
    return b*(-1);
}

double f1(double a, double b, double c)
{
    double ris = (a + b) + c;
    return ris;
}

double f2(double a, double b, double c)
{
    double ris = a + (b + c);
    return ris;
}

int main()
{
    int d0 = 7;
    int d1 = 2;
    double a = 0.0 ,b = 0.0 ,c = 0.0, risf1 = 0.0, risf2 = 0.0;

    cout<< "d0: "<< d0 << "\nd1: "<< d1 <<endl;

    //Calcolo b e c
    b = findB(d1);
	cout<< "b: "<< b << endl;

    c = findC(b);
    cout<< "c: "<< c << endl;

    for(int i=0; i < 7; i++)
    {
        a = findA(d0,i);
        cout<< "a["<< i << "]" << ": "<< a <<endl;
        risf1 = f1(a,b,c);
        risf2 = f2(a,b,c);
        cout<< "(a + b) + c: "<<  risf1 <<endl;
        cout<< "a + (b + c): "<<  risf2 <<endl;

    }
    return 0;
}

