//Gruppo: 
//Dellepiane Emanuele, 4876072
//Manini Filippo, 4798004
//Miggiano Davide, 4840761 
#include <iostream>
#include <cmath>
#include <limits>
using namespace std;

//eps = 2^-d
//d : il più grande intero positivo t.c 1+2^(-d) > 1

//const int MAX=numeric_limits<int>::max();

//powf : pow per i float
float singleEps()
{
    int d = 0;
    float eps = powf(2,-d);
    while(1+eps > 1.0)
    {
        d++;
        eps = powf(2,-d);
    }
    return d;
}

double doubleEps()
{
    int d = 0;
    double eps = powf(2,-d);
    while(1+eps > 1.0)
    {
        d++;
        eps = powf(2,-d);
    }
    return d;
}

int main()
{
    cout<<"Singola precisione: "<< singleEps() << "\n"<<endl;
    cout<<"Doppia precisione: "<< doubleEps() << "\n" <<endl;
    return 0;
}
